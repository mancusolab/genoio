# pattern: Functional Core

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy import sparse as scipy_sparse

from ._errors import MissingDataError

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
    missing: str,
    dtype: np.dtype[Any],
    values_layout: str = "sample_major",
    missing_mask: Sequence[bool] | NDArray[Any] | None = None,
    missing_indices: Sequence[int] | NDArray[Any] | None = None,
) -> NDArray[Any]:
    r"""Convert a flat Rust dense matrix payload into a NumPy array.

    Rust returns Python-owned NumPy buffers through the extension boundary; this
    function applies the public dtype and missing-data policy.

    **Arguments:**

    - `values`: flat matrix values in `values_layout` order.
    - `shape`: `(n_samples, n_variants)`.
    - `missing`: validated missing-data policy.
    - `dtype`: output NumPy dtype.
    - `values_layout`: flat buffer layout, either `"sample_major"` or
      `"variant_major"`.
    - `missing_mask`: optional flat boolean mask aligned with `values`.
    - `missing_indices`: optional flat indices aligned with `values`.

    **Returns:**

    Dense NumPy matrix with shape `shape`.
    """
    values_array = np.asarray(values, dtype=dtype)
    array = _reshape_dense_payload(values_array, shape, values_layout)
    if missing_mask is not None and missing_indices is not None:
        raise AssertionError("dense missing payload must use mask or indices, not both")
    if missing_indices is not None:
        return _apply_dense_missing_indices(
            values_array,
            array,
            shape,
            values_layout,
            missing,
            missing_indices,
        )
    if missing_mask is None:
        return array
    mask = _dense_missing_mask_from_flat(
        np.asarray(missing_mask, dtype=bool),
        values_array.size,
        shape,
        values_layout,
    )
    return _apply_dense_missing_mask(array, mask, missing)


def _apply_dense_missing_mask(
    array: NDArray[Any],
    mask: NDArray[np.bool_],
    missing: str,
) -> NDArray[Any]:
    """Apply a dense missing mask after both payloads have public matrix shape."""
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


def _dense_missing_mask_from_flat(
    mask_values: NDArray[np.bool_],
    values_size: int,
    shape: tuple[int, int],
    values_layout: str,
) -> NDArray[np.bool_]:
    """Return a public-shaped missing mask aligned with a flat Rust values buffer."""
    if mask_values.size != values_size:
        raise AssertionError("dense missing mask length does not match values length")
    return _reshape_dense_payload(mask_values, shape, values_layout)


def _apply_dense_missing_indices(
    flat_values: NDArray[Any],
    array: NDArray[Any],
    shape: tuple[int, int],
    values_layout: str,
    missing: str,
    missing_indices: Sequence[int] | NDArray[Any],
) -> NDArray[Any]:
    """Apply sparse missing-value positions aligned with the flat Rust buffer."""
    indices = _dense_missing_indices_from_flat(missing_indices)
    if np.any(indices < 0) or np.any(indices >= flat_values.size):
        raise AssertionError("dense missing index is outside values buffer")
    if missing == "nan":
        flat_values[indices] = np.nan
        return array
    if missing == "raise":
        if indices.size:
            raise MissingDataError("missing genotype calls are present in retained data")
        return array
    if missing == "impute":
        return _impute_missing_by_variant_indices(array, indices, shape, values_layout)
    raise AssertionError(f"unvalidated missing-data policy: {missing}")


def _dense_missing_indices_from_flat(
    missing_indices: Sequence[int] | NDArray[Any],
) -> NDArray[np.int64]:
    """Return validated int64 indices into the flat Rust values buffer."""
    indices = np.asarray(missing_indices)
    if indices.ndim != 1:
        raise AssertionError("dense missing indices must be one-dimensional")
    if not np.issubdtype(indices.dtype, np.integer):
        raise AssertionError("dense missing indices must be integers")
    return indices.astype(np.int64, copy=False)


def _impute_missing_by_variant_indices(
    array: NDArray[Any],
    flat_indices: NDArray[np.int64],
    shape: tuple[int, int],
    values_layout: str,
) -> NDArray[Any]:
    """Impute sparse missing positions without expanding them to a dense mask."""
    if not flat_indices.size:
        return array
    if np.unique(flat_indices).size != flat_indices.size:
        raise AssertionError("dense missing indices must be unique")

    missing_rows, missing_columns = _flat_indices_to_matrix_coordinates(
        flat_indices,
        shape,
        values_layout,
    )
    n_samples, n_variants = shape

    missing_counts = np.bincount(missing_columns, minlength=n_variants)
    called_counts = n_samples - missing_counts
    all_missing_columns = np.flatnonzero((missing_counts != 0) & (called_counts == 0))
    if all_missing_columns.size:
        first_all_missing_column = int(all_missing_columns[0])
        raise MissingDataError(f"cannot impute all-missing variant at column {first_all_missing_column}")

    # Missing positions can contain backend-specific placeholder values. Column
    # means must be based only on called entries, so subtract sparse missing
    # contributions instead of first materializing a dense missing mask.
    column_sums = array.sum(axis=0)
    missing_sums = np.bincount(
        missing_columns,
        weights=array[missing_rows, missing_columns],
        minlength=n_variants,
    )
    called_sums = column_sums - missing_sums
    means = np.divide(
        called_sums,
        called_counts,
        out=np.zeros_like(called_sums),
        where=called_counts != 0,
    )

    array[missing_rows, missing_columns] = means[missing_columns]
    return array


def _flat_indices_to_matrix_coordinates(
    indices: NDArray[np.int64],
    shape: tuple[int, int],
    values_layout: str,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Map Rust flat-buffer indices to public sample-by-variant coordinates."""
    n_samples, n_variants = shape
    if values_layout == "sample_major":
        return indices // n_variants, indices % n_variants
    if values_layout == "variant_major":
        return indices % n_samples, indices // n_samples
    raise AssertionError(f"unvalidated dense value layout: {values_layout}")


def _reshape_dense_payload(array: NDArray[Any], shape: tuple[int, int], values_layout: str) -> NDArray[Any]:
    """Return the public sample-by-variant view for a Rust dense payload."""
    if values_layout == "sample_major":
        return array.reshape(shape)
    if values_layout == "variant_major":
        n_samples, n_variants = shape
        # Variant-major buffers are already Python-owned NumPy arrays here; the
        # transpose is a strided view, not another Rust-side materialization.
        return array.reshape((n_variants, n_samples)).T
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
