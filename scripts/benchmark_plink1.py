#!/usr/bin/env python
# pattern: Mixed

from __future__ import annotations

import argparse
import importlib
from pathlib import Path
from typing import Any, cast

import numpy as np
from bench_common import benchmark, compare_summaries, positive_int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare genoio PLINK1 reads against pandas-plink.")
    parser.add_argument("--prefix", type=Path, default=Path("data/chr22_hg38"))
    parser.add_argument("--max-variants", type=positive_int, default=1_000)
    parser.add_argument("--repeats", type=positive_int, default=3)
    parser.add_argument("--backend", choices=["both", "genoio", "pandas_plink"], default="both")
    parser.add_argument("--pandas-ref", choices=["a0", "a1"], default="a0")
    parser.add_argument("--no-compare", action="store_true")
    return parser.parse_args()


def read_genoio(args: argparse.Namespace) -> np.ndarray:
    import genoio

    return next(
        genoio.bfile(args.prefix).iter_blocks(
            args.max_variants,
            missing="nan",
            dtype=np.float32,
        )
    )


def _patch_dask_memmap_tokenization() -> None:
    """Avoid hashing the whole BED mmap while pandas-plink builds its Dask graph."""
    dask_tokenize = cast(Any, importlib.import_module("dask.tokenize"))

    # Force Dask's lazy NumPy token handlers to register before overriding memmap.
    _ = dask_tokenize.normalize_token.dispatch(np.memmap)
    dask_tokenize.normalize_token.register(
        np.memmap,
        lambda mmap: (
            "memmap",
            str(mmap.filename),
            mmap.dtype.str,
            mmap.mode,
            mmap.offset,
            mmap.shape,
            mmap.strides,
        ),
    )


def read_pandas_plink(args: argparse.Namespace) -> np.ndarray:
    from pandas_plink import read_plink1_bin  # type: ignore[import-not-found]

    _patch_dask_memmap_tokenization()
    data = read_plink1_bin(
        str(args.prefix.with_suffix(".bed")),
        str(args.prefix.with_suffix(".bim")),
        str(args.prefix.with_suffix(".fam")),
        verbose=False,
        ref=args.pandas_ref,
    )
    return np.asarray(data.isel(variant=slice(0, args.max_variants)).values, dtype=np.float32)


def main() -> None:
    args = parse_args()
    genoio_matrix = None
    pandas_plink_matrix = None
    if args.backend in {"both", "genoio"}:
        genoio_matrix = benchmark("genoio_plink1", lambda: read_genoio(args), args.repeats)
    if args.backend in {"both", "pandas_plink"}:
        pandas_plink_matrix = benchmark("pandas_plink", lambda: read_pandas_plink(args), args.repeats)
    if not args.no_compare and genoio_matrix is not None and pandas_plink_matrix is not None:
        compare_summaries("genoio_plink1", genoio_matrix, "pandas_plink", pandas_plink_matrix)


if __name__ == "__main__":
    main()
