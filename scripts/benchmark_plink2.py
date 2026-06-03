#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from bench_common import benchmark, compare_summaries, plink2_prefix_with_uncompressed_pvar, positive_int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare genoio PLINK2 reads against pgenlib.PgenReader.")
    parser.add_argument("--prefix", type=Path, default=Path("data/chr22_hg38"))
    parser.add_argument("--max-variants", type=positive_int, default=1_000)
    parser.add_argument("--repeats", type=positive_int, default=3)
    parser.add_argument("--backend", choices=["both", "genoio", "pgenlib"], default="both")
    parser.add_argument(
        "--pgenlib-path",
        type=Path,
        default=None,
        help="Optional path to plink-ng/2.0/Python when pgenlib is built in-place but not installed.",
    )
    parser.add_argument("--no-compare", action="store_true")
    return parser.parse_args()


def read_genoio(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    return next(
        genoio.open(prefix, format="plink2").blocks(
            max_variants,
            missing="nan",
            dtype=np.float32,
        )
    )


def import_pgenlib(pgenlib_path: Path | None):
    if pgenlib_path is not None:
        sys.path.insert(0, str(pgenlib_path))
        sys.path.insert(0, str(pgenlib_path / "src"))
    try:
        import pgenlib
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "pgenlib is not importable; install it with `pip install Pgenlib` or pass "
            "`--pgenlib-path plink-ng/2.0/Python` after building it in-place"
        ) from error

    return pgenlib


def read_pgenlib(args: argparse.Namespace) -> np.ndarray:
    pgenlib = import_pgenlib(args.pgenlib_path)
    reader = pgenlib.PgenReader(bytes(args.prefix.with_suffix(".pgen")))
    try:
        sample_ct = reader.get_raw_sample_ct()
        variant_ct = min(args.max_variants, reader.get_variant_ct())
        values = np.empty((sample_ct, variant_ct), dtype=np.int8)
        reader.read_range(0, variant_ct, values, sample_maj=True)
    finally:
        reader.close()
    matrix = values.astype(np.float32, copy=False)
    matrix[matrix == -9] = np.nan
    return matrix


def main() -> None:
    args = parse_args()
    genoio_matrix = None
    pgenlib_matrix = None
    if args.backend in {"both", "genoio"}:
        with plink2_prefix_with_uncompressed_pvar(args.prefix) as genoio_prefix:
            genoio_matrix = benchmark(
                "genoio_plink2",
                lambda: read_genoio(genoio_prefix, args.max_variants),
                args.repeats,
            )
    if args.backend in {"both", "pgenlib"}:
        pgenlib_matrix = benchmark("pgenlib_pgenreader", lambda: read_pgenlib(args), args.repeats)
    if not args.no_compare and genoio_matrix is not None and pgenlib_matrix is not None:
        compare_summaries("genoio_plink2", genoio_matrix, "pgenlib_pgenreader", pgenlib_matrix)


if __name__ == "__main__":
    main()
