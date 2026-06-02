# pattern: Functional Core

from __future__ import annotations

from typing import Any

import numpy as np
import polars as pl
from scipy import sparse as scipy_sparse

from ._errors import MissingDataError

_SAMPLE_COLUMNS = ["fid", "iid", "father", "mother", "sex", "phenotype"]
_VARIANT_COLUMNS = [
    "chrom",
    "pos",
    "id",
    "a0",
    "a1",
    "ref",
    "alt",
    "source_a0",
    "source_a1",
    "flipped",
    "af",
    "maf",
    "mac",
    "missing_rate",
    "n_called",
]


def samples_frame(records: list[dict[str, Any]]) -> pl.DataFrame:
    frame = pl.DataFrame(
        {column: [record.get(column) for record in records] for column in _SAMPLE_COLUMNS},
        schema=_SAMPLE_COLUMNS,
    )
    if records and all(record.get("haplotype_index") is not None for record in records):
        return frame.with_columns(
            pl.Series("source_sample_index", [record["source_sample_index"] for record in records]),
            pl.Series("haplotype_index", [record["haplotype_index"] for record in records]),
        )
    return frame


def variants_frame(records: list[dict[str, Any]]) -> pl.DataFrame:
    # Rust returns in-memory metadata records across the PyO3 boundary. Eager
    # DataFrame assembly is the adapter boundary here because there is no file
    # scan or deferred query plan for Polars to optimize, and source order is
    # already the downstream contract.
    return pl.DataFrame(
        {
            "chrom": [record.get("chrom") for record in records],
            "pos": [record.get("pos") for record in records],
            "id": [record.get("id") for record in records],
            "a0": [record.get("a0") for record in records],
            "a1": [record.get("a1") for record in records],
            "ref": [record.get("ref_allele") for record in records],
            "alt": [record.get("alt_allele") for record in records],
            "source_a0": [record.get("source_a0") for record in records],
            "source_a1": [record.get("source_a1") for record in records],
            "flipped": [record.get("flipped", False) for record in records],
            "af": [record.get("af") for record in records],
            "maf": [record.get("maf") for record in records],
            "mac": [record.get("mac") for record in records],
            "missing_rate": [record.get("missing_rate") for record in records],
            "n_called": [record.get("n_called") for record in records],
        },
        schema=_VARIANT_COLUMNS,
    )


def dense_array_from_rust(
    *,
    values: list[float],
    shape: tuple[int, int],
    missing_mask: list[bool],
    missing: str,
    dtype: np.dtype[Any],
) -> np.ndarray:
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
) -> Any:
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
    genotype_matrix: Any,
    samples: pl.DataFrame,
    variants: pl.DataFrame,
    *,
    return_samples: bool,
    return_variants: bool,
) -> Any:
    if return_samples and return_variants:
        return genotype_matrix, samples, variants
    if return_samples:
        return genotype_matrix, samples
    if return_variants:
        return genotype_matrix, variants
    return genotype_matrix


def _impute_missing_by_variant(array: np.ndarray, mask: np.ndarray) -> np.ndarray:
    imputed = array.copy()
    for variant_index in range(imputed.shape[1]):
        missing_rows = mask[:, variant_index]
        if not missing_rows.any():
            continue
        called_rows = ~missing_rows
        if not called_rows.any():
            raise MissingDataError(f"cannot impute all-missing variant at column {variant_index}")
        imputed[missing_rows, variant_index] = imputed[called_rows, variant_index].mean()
    return imputed
