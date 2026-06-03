#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from bench_common import benchmark, compare_summaries, nonnegative_float, positive_int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare genoio VCF reads against cyvcf2.")
    parser.add_argument("--vcf", type=Path, default=Path("data/chr22_hg38.vcf.gz"))
    parser.add_argument("--max-variants", type=positive_int, default=1_000)
    parser.add_argument("--repeats", type=positive_int, default=3)
    parser.add_argument("--backend", choices=["both", "genoio", "cyvcf2"], default="both")
    parser.add_argument("--qual-min", type=nonnegative_float, default=None)
    parser.add_argument("--maf-max", type=nonnegative_float, default=None)
    parser.add_argument("--no-compare", action="store_true")
    return parser.parse_args()


def genoio_filter(args: argparse.Namespace):
    import genoio

    expr = genoio.biallelic()
    if args.qual_min is not None:
        expr = expr & genoio.qual(min=args.qual_min)
    if args.maf_max is not None:
        expr = expr & genoio.maf(max=args.maf_max)
    return expr


def read_genoio(args: argparse.Namespace) -> np.ndarray:
    import genoio

    return next(
        genoio.open(args.vcf).blocks(
            args.max_variants,
            variants=genoio_filter(args),
            missing="nan",
            dtype=np.float32,
        )
    )


def read_cyvcf2(args: argparse.Namespace) -> np.ndarray:
    import cyvcf2

    columns = []
    for variant in cyvcf2.VCF(str(args.vcf), gts012=True, strict_gt=True):
        if len(variant.ALT or []) != 1:
            continue
        if args.qual_min is not None and (variant.QUAL is None or variant.QUAL < args.qual_min):
            continue
        if args.maf_max is not None:
            alt_freq = float(np.atleast_1d(variant.aaf)[0])
            if min(alt_freq, 1.0 - alt_freq) > args.maf_max:
                continue
        values = variant.gt_types.astype(np.float32, copy=True)
        values[values == 3] = np.nan
        columns.append(values)
        if len(columns) == args.max_variants:
            break
    if not columns:
        return np.empty((0, 0), dtype=np.float32)
    return np.stack(columns, axis=1)


def main() -> None:
    args = parse_args()
    genoio_matrix = None
    cyvcf2_matrix = None
    if args.backend in {"both", "genoio"}:
        genoio_matrix = benchmark("genoio_vcf", lambda: read_genoio(args), args.repeats)
    if args.backend in {"both", "cyvcf2"}:
        cyvcf2_matrix = benchmark("cyvcf2_vcf", lambda: read_cyvcf2(args), args.repeats)
    if not args.no_compare and genoio_matrix is not None and cyvcf2_matrix is not None:
        compare_summaries("genoio_vcf", genoio_matrix, "cyvcf2_vcf", cyvcf2_matrix)


if __name__ == "__main__":
    main()
