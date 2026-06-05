# pattern: Functional Core

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy import sparse as scipy_sparse

from ._errors import MissingDataError

SparseMatrixResult = scipy_sparse.spmatrix | scipy_sparse.sparray
MatrixResult = NDArray[Any] | SparseMatrixResult

_SAMPLE_COLUMNS = ["fid", "iid", "father", "mother", "sex", "phenotype"]
_VARIANT_COLUMNS = [
    "chrom",
    "pos",
    "id",
    "a0",
    "a1",
]


MetadataColumns = dict[str, Sequence[Any]]


def samples_frame(columns: MetadataColumns) -> pl.DataFrame:
    r"""Build the public sample metadata frame from Rust sample columns.

    **Arguments:**

    - `columns`: column-oriented sample payload returned by the Rust extension.

    **Returns:**

    Polars DataFrame in source sample order.
    """
    frame = pl.DataFrame(
        {column: columns[column] for column in _SAMPLE_COLUMNS},
        schema=_SAMPLE_COLUMNS,
    )
    haplotype_indices = columns["haplotype_index"]
    if haplotype_indices and all(index is not None for index in haplotype_indices):
        return frame.with_columns(
            pl.Series("source_sample_index", columns["source_sample_index"]),
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
    values: list[float],
    shape: tuple[int, int],
    missing_mask: list[bool],
    missing: str,
    dtype: np.dtype[Any],
) -> NDArray[Any]:
    r"""Convert a flat Rust dense matrix payload into a NumPy array.

    Rust returns Python-owned NumPy buffers through the extension boundary; this
    function applies the public dtype and missing-data policy.

    **Arguments:**

    - `values`: flat sample-major matrix values.
    - `shape`: `(n_samples, n_variants)`.
    - `missing_mask`: flat boolean mask aligned with `values`.
    - `missing`: validated missing-data policy.
    - `dtype`: output NumPy dtype.

    **Returns:**

    Dense NumPy matrix with shape `shape`.
    """
    array = np.asarray(values, dtype=dtype).reshape(shape)
    mask = np.asarray(missing_mask, dtype=bool).reshape(shape)
    if missing == "nan":
        array[mask] = np.nan
        return array
    if missing == "raise":
        if mask.any():
            raise MissingDataError("missing genotype calls are present in retained data")
        return array
    if missing == "impute":
        return _impute_missing_by_variant(array, mask)
    raise AssertionError(f"unvalidated missing-data policy: {missing}")


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
) -> MatrixResult | tuple[MatrixResult, pl.DataFrame] | tuple[MatrixResult, pl.DataFrame, pl.DataFrame]:
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


def _impute_missing_by_variant(array: NDArray[Any], mask: NDArray[np.bool_]) -> NDArray[Any]:
    if not mask.any():
        return array

    called_mask = ~mask
    called_counts = called_mask.sum(axis=0)
    all_missing_columns = np.flatnonzero(mask.any(axis=0) & (called_counts == 0))
    if all_missing_columns.size:
        raise MissingDataError(f"cannot impute all-missing variant at column {int(all_missing_columns[0])}")

    called_sums = np.where(mask, 0, array).sum(axis=0)
    means = np.divide(
        called_sums,
        called_counts,
        out=np.zeros_like(called_sums),
        where=called_counts != 0,
    )

    # Each missing entry receives the mean of its own variant column. The
    # indexed assignment avoids per-variant Python loops on large matrices.
    missing_rows, missing_columns = np.nonzero(mask)
    imputed = array.copy()
    imputed[missing_rows, missing_columns] = means[missing_columns]
    return imputed
