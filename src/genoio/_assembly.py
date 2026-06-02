# pattern: Functional Core

from __future__ import annotations

from typing import Any

import polars as pl

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
    return pl.DataFrame(
        {column: [record.get(column) for record in records] for column in _SAMPLE_COLUMNS},
        schema=_SAMPLE_COLUMNS,
    )


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
            "af": [None for _ in records],
            "maf": [None for _ in records],
            "mac": [None for _ in records],
            "missing_rate": [None for _ in records],
            "n_called": [None for _ in records],
        },
        schema=_VARIANT_COLUMNS,
    )
