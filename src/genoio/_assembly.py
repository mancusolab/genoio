# pattern: Functional Core

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy import sparse as scipy_sparse

SparseMatrixResult = scipy_sparse.csc_matrix | scipy_sparse.csr_matrix
MatrixResult = NDArray[Any] | SparseMatrixResult
ReadResult = MatrixResult | tuple[MatrixResult, pl.DataFrame] | tuple[MatrixResult, pl.DataFrame, pl.DataFrame]

_SAMPLE_COLUMNS = ["fid", "iid", "father", "mother", "sex", "phenotype"]
_VARIANT_COLUMNS = [
    "chrom",
    "pos",
    "id",
    "a0",
    "a1",
]
_DENSE_LAYOUT_SAMPLE_MAJOR = "sample_major"
_DENSE_LAYOUT_VARIANT_MAJOR = "variant_major"


MetadataColumns = dict[str, Sequence[Any]]


def samples_frame(columns: MetadataColumns, *, include_haplotype_columns: bool = False) -> pl.DataFrame:
    r"""Build the public sample metadata frame from Rust sample columns.

    **Arguments:**

    - `columns`: column-oriented sample payload returned by the Rust extension.
    - `include_haplotype_columns`: preserve haplotype row mapping columns even
      when the retained sample set is empty.

    **Returns:**

    Polars DataFrame in source sample order.
    """
    frame = pl.DataFrame(
        {column: columns[column] for column in _SAMPLE_COLUMNS},
        schema=_SAMPLE_COLUMNS,
    )
    haplotype_indices = columns.get("haplotype_index")
    if include_haplotype_columns or (haplotype_indices and all(index is not None for index in haplotype_indices)):
        source_sample_indices = columns.get("source_sample_index", [None] * frame.height)
        haplotype_indices = haplotype_indices or [None] * frame.height
        return frame.with_columns(
            pl.Series("source_sample_index", source_sample_indices),
            pl.Series("haplotype_index", haplotype_indices),
        )
    return frame


def variants_frame(columns: MetadataColumns) -> pl.DataFrame:
    r"""Build the public variant metadata frame from Rust variant columns.

    **Arguments:**

    - `columns`: column-oriented variant payload returned by the Rust extension.

    **Returns:**

    Polars DataFrame in source variant order. The public schema is deliberately
    limited to the columns needed to interpret matrix columns and counted
    allele orientation.
    """
    # Rust returns in-memory columns across the PyO3 boundary. Eager DataFrame
    # assembly is the adapter boundary here because there is no file scan or
    # deferred query plan for Polars to optimize, and source order is already
    # the downstream contract.
    return pl.DataFrame(
        {column: columns[column] for column in _VARIANT_COLUMNS},
        schema=_VARIANT_COLUMNS,
    )


def dense_array_from_rust(
    *,
    values: Sequence[float] | NDArray[Any],
    shape: tuple[int, int],
    dtype: np.dtype[Any],
    values_layout: str = "sample_major",
) -> NDArray[Any]:
    r"""Convert a flat Rust dense matrix payload into a NumPy array.

    Rust returns Python-owned NumPy buffers through the extension boundary with
    dense missing-data policy already applied. This function handles dtype
    conversion and layout views.

    **Arguments:**

    - `values`: flat matrix values in `values_layout` order.
    - `shape`: `(n_samples, n_variants)`.
    - `dtype`: output NumPy dtype.
    - `values_layout`: flat buffer layout, either `"sample_major"` or
      `"variant_major"`.

    **Returns:**

    Dense NumPy matrix with shape `shape`.
    """
    values_layout = _validate_dense_layout(values_layout)
    values_array = np.asarray(values, dtype=dtype)
    return _reshape_dense_payload(values_array, shape, values_layout)


def _reshape_dense_payload(array: NDArray[Any], shape: tuple[int, int], values_layout: str) -> NDArray[Any]:
    """Return the public sample-by-variant view for a Rust dense payload."""
    expected_size = shape[0] * shape[1]
    if array.size != expected_size:
        raise AssertionError(f"dense values length {array.size} does not match shape {shape}")
    if values_layout == _DENSE_LAYOUT_SAMPLE_MAJOR:
        return array.reshape(shape)
    if values_layout == _DENSE_LAYOUT_VARIANT_MAJOR:
        n_samples, n_variants = shape
        # Variant-major buffers are already Python-owned NumPy arrays here; the
        # transpose is a strided view, not another Rust-side materialization.
        return array.reshape((n_variants, n_samples)).T
    raise AssertionError(f"unvalidated dense value layout: {values_layout}")


def _validate_dense_layout(values_layout: str) -> str:
    """Return a known dense layout tag or fail at the private bridge boundary."""
    if values_layout in {_DENSE_LAYOUT_SAMPLE_MAJOR, _DENSE_LAYOUT_VARIANT_MAJOR}:
        return values_layout
    raise AssertionError(f"unvalidated dense value layout: {values_layout}")


def sparse_matrix_from_rust(
    *,
    indptr: list[int],
    indices: list[int],
    data: list[float],
    shape: tuple[int, int],
    dtype: np.dtype[Any],
    sparse_format: str,
) -> SparseMatrixResult:
    r"""Convert Rust CSC arrays into the requested SciPy sparse format.

    Rust always emits CSC because variants are accumulated column-wise. CSR is
    a Python-side view conversion after the validated CSC arrays cross the
    extension boundary.

    **Arguments:**

    - `indptr`: CSC column pointer array.
    - `indices`: CSC row index array.
    - `data`: nonzero genotype values.
    - `shape`: `(n_samples, n_variants)`.
    - `dtype`: output value dtype.
    - `sparse_format`: `"csc"` or `"csr"`.

    **Returns:**

    SciPy sparse genotype matrix.
    """
    matrix = scipy_sparse.csc_matrix(
        (
            np.asarray(data, dtype=dtype),
            np.asarray(indices, dtype=np.int64),
            np.asarray(indptr, dtype=np.int64),
        ),
        shape=shape,
    )
    if sparse_format == "csc":
        return matrix
    if sparse_format == "csr":
        return matrix.tocsr()
    raise AssertionError(f"unvalidated sparse format: {sparse_format}")


def read_result_tuple(
    genotype_matrix: MatrixResult,
    samples: pl.DataFrame | None,
    variants: pl.DataFrame | None,
    *,
    return_samples: bool,
    return_variants: bool,
) -> ReadResult:
    r"""Attach optional metadata frames to a matrix result.

    `samples` and `variants` are optional at this assembly boundary because
    Rust metadata columns are only converted when the corresponding return flag
    is set.

    **Returns:**

    Matrix alone or tuple with requested metadata frames.
    """
    if return_samples and return_variants:
        # These asserts document the flag/data invariant for static checkers;
        # callers build metadata conditionally from the same flags.
        assert samples is not None
        assert variants is not None
        return genotype_matrix, samples, variants
    if return_samples:
        assert samples is not None
        return genotype_matrix, samples
    if return_variants:
        assert variants is not None
        return genotype_matrix, variants
    return genotype_matrix
