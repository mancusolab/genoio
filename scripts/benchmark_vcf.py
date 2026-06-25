#!/usr/bin/env python
# pattern: Mixed

from __future__ import annotations

import argparse
import gzip
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
from bench_common import benchmark, compare_summaries, matrix_summary, nonnegative_float, positive_int, print_result

SCENARIOS = (
    "metadata",
    "matrix-only",
    "with-variants",
    "sample-filtered",
    "genotype-filtered",
    "indexed-region",
    "indexed-region-sample-filtered",
)
KINDS = ("geno", "dosage", "haplo")
_last_variant_metadata_length: int | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare genoio VCF reads against cyvcf2.")
    parser.add_argument("--vcf", type=Path, default=Path("data/chr22_hg38.vcf.gz"))
    parser.add_argument("--max-variants", type=positive_int, default=1_000)
    parser.add_argument("--repeats", type=positive_int, default=3)
    parser.add_argument("--backend", choices=["both", "genoio", "cyvcf2"], default="both")
    parser.add_argument("--scenario", choices=[*SCENARIOS, "all"], default="matrix-only")
    parser.add_argument(
        "--kind",
        choices=KINDS,
        default="geno",
        help=(
            "Matrix kind to time. Defaults to genotype hardcalls; "
            '"dosage" reads FORMAT/DS; "haplo" times phased haplotype hardcalls.'
        ),
    )
    parser.add_argument("--sparse", action="store_true", help="Time genoio sparse CSC output.")
    parser.add_argument("--region", default="22:20000000-21000000")
    parser.add_argument(
        "--samples",
        default=None,
        help="Optional comma-separated sample IDs for sample-filtered reads.",
    )
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


def _read_options(kind: str, sparse: bool) -> dict[str, object]:
    options: dict[str, object] = {
        "missing": "raise" if sparse else "nan",
        "dtype": np.float32,
        "sparse": "csc" if sparse else False,
    }
    if kind == "haplo":
        options["kind"] = "haplo"
    if kind == "dosage":
        options["dosage"] = "dosage"
    return options


def _genoio_label(kind: str, scenario: str, sparse: bool) -> str:
    suffix = scenario.replace("-", "_")
    kind_part = kind if kind != "geno" else ""
    sparse_part = "sparse" if sparse else ""
    parts = ["genoio_vcf", kind_part, sparse_part, suffix]
    return "_".join(part for part in parts if part)


def _requested_sample_ids(args: argparse.Namespace) -> list[str]:
    if args.samples:
        sample_ids = [sample for sample in args.samples.split(",") if sample]
        if not sample_ids:
            raise RuntimeError("--samples must contain at least one sample ID")
        return sample_ids
    sample_ids = _read_vcf_sample_ids(args.vcf)
    keep_count = max(1, len(sample_ids) // 2)
    return sample_ids[:keep_count]


def read_genoio_matrix_only(args: argparse.Namespace) -> Any:
    import genoio

    matrix = next(
        genoio.vcf(args.vcf).iter_blocks(
            args.max_variants,
            variants=genoio_filter(args),
            **_read_options(args.kind, args.sparse),
        )
    )
    return matrix


def read_genoio_metadata(args: argparse.Namespace) -> np.ndarray:
    import genoio

    dataset = genoio.vcf(args.vcf)
    samples = dataset.samples()
    variants = dataset.variants()
    return np.array([samples.height, variants.height], dtype=np.int64)


def benchmark_metadata(name: str, fn: Callable[[], Any], repeats: int) -> Any:
    start = time.perf_counter()
    first = fn()
    cold_time = time.perf_counter() - start
    summary = matrix_summary(first)
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        observed = fn()
        times.append(time.perf_counter() - start)
        observed_summary = matrix_summary(observed)
        if observed_summary["shape"] != summary["shape"]:
            raise RuntimeError(f"{name} shape changed between repeats")
    print_result(name, summary, times)
    print(f"  cold={cold_time:.4f}s")
    return first


def read_genoio_with_variants(args: argparse.Namespace) -> Any:
    import genoio

    global _last_variant_metadata_length
    matrix, variants = next(
        genoio.vcf(args.vcf).iter_blocks(
            args.max_variants,
            variants=genoio_filter(args),
            return_variants=True,
            **_read_options(args.kind, args.sparse),
        )
    )
    _last_variant_metadata_length = variants.height
    return matrix


def read_genoio_sample_filtered(args: argparse.Namespace) -> Any:
    import genoio

    matrix = next(
        genoio.vcf(args.vcf).iter_blocks(
            args.max_variants,
            variants=genoio_filter(args),
            samples=_requested_sample_ids(args),
            **_read_options(args.kind, args.sparse),
        )
    )
    return matrix


def read_genoio_genotype_filtered(args: argparse.Namespace) -> Any:
    import genoio

    matrix = next(
        genoio.vcf(args.vcf).iter_blocks(
            args.max_variants,
            variants=genoio_filter(args) & genoio.maf(min=0.01),
            **_read_options(args.kind, args.sparse),
        )
    )
    return matrix


def read_genoio_indexed_region(args: argparse.Namespace) -> Any:
    import genoio

    matrix, variants = next(
        genoio.vcf(args.vcf).iter_blocks(
            args.max_variants,
            variants=genoio.region(args.region) & genoio_filter(args),
            return_variants=True,
            **_read_options(args.kind, args.sparse),
        )
    )
    _validate_region_variants(variants, args.region)
    return matrix


def read_genoio_indexed_region_sample_filtered(args: argparse.Namespace) -> Any:
    import genoio

    matrix, variants = next(
        genoio.vcf(args.vcf).iter_blocks(
            args.max_variants,
            variants=genoio.region(args.region) & genoio_filter(args),
            samples=_requested_sample_ids(args),
            return_variants=True,
            **_read_options(args.kind, args.sparse),
        )
    )
    _validate_region_variants(variants, args.region)
    return matrix


def read_cyvcf2(args: argparse.Namespace) -> np.ndarray:
    import cyvcf2  # type: ignore[import-not-found]

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


def _read_vcf_sample_ids(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                sample_ids = line.rstrip("\n").split("\t")[9:]
                if sample_ids:
                    return sample_ids
                break
    raise RuntimeError(f"{path} does not contain VCF sample IDs")


def _validate_region_variants(variants: Any, region: str) -> None:
    if not variants.height:
        return
    chrom, coords = region.split(":", 1)
    start, end = (int(value) for value in coords.split("-", 1))
    observed_min = int(variants["pos"].min())
    observed_max = int(variants["pos"].max())
    if variants["chrom"].unique().to_list() != [chrom] or observed_min < start or observed_max > end:
        raise RuntimeError(f"indexed region result escaped requested region: {region}")


def selected_scenarios(scenario: str) -> tuple[str, ...]:
    if scenario == "all":
        return SCENARIOS
    return (scenario,)


def benchmark_genoio_scenario(scenario: str, args: argparse.Namespace) -> Any:
    if scenario == "metadata":
        return benchmark_metadata(
            _genoio_label(args.kind, scenario, args.sparse),
            lambda: read_genoio_metadata(args),
            args.repeats,
        )
    if scenario == "matrix-only":
        return benchmark(
            _genoio_label(args.kind, scenario, args.sparse),
            lambda: read_genoio_matrix_only(args),
            args.repeats,
        )
    if scenario == "with-variants":
        variant_metadata_length = None
        global _last_variant_metadata_length
        _last_variant_metadata_length = None

        def read_matrix() -> Any:
            nonlocal variant_metadata_length
            result = read_genoio_with_variants(args)
            variant_metadata_length = _last_variant_metadata_length
            return result

        matrix = benchmark(_genoio_label(args.kind, scenario, args.sparse), read_matrix, args.repeats)
        print(f"  variant_metadata length={variant_metadata_length}")
        return matrix
    if scenario == "sample-filtered":
        return benchmark(
            _genoio_label(args.kind, scenario, args.sparse),
            lambda: read_genoio_sample_filtered(args),
            args.repeats,
        )
    if scenario == "genotype-filtered":
        return benchmark(
            _genoio_label(args.kind, scenario, args.sparse),
            lambda: read_genoio_genotype_filtered(args),
            args.repeats,
        )
    if scenario == "indexed-region":
        return benchmark(
            _genoio_label(args.kind, scenario, args.sparse),
            lambda: read_genoio_indexed_region(args),
            args.repeats,
        )
    if scenario == "indexed-region-sample-filtered":
        return benchmark(
            _genoio_label(args.kind, scenario, args.sparse),
            lambda: read_genoio_indexed_region_sample_filtered(args),
            args.repeats,
        )
    raise ValueError(f"unknown scenario: {scenario}")


def print_cyvcf2_skip(scenario: str, args: argparse.Namespace) -> None:
    if scenario == "metadata":
        message = "skipped cyvcf2 comparison for metadata: benchmark only compares matrix-only genotype reads"
    elif args.kind != "geno":
        message = f"skipped cyvcf2 comparison for {args.kind} {scenario}: benchmark only compares genotype hardcalls"
    elif args.sparse:
        message = f"skipped cyvcf2 comparison for sparse {scenario}: comparison backend returns dense genotypes"
    else:
        message = f"skipped cyvcf2 comparison for {scenario}: benchmark only implements matrix-only genotype reads"
    print(message)


def main() -> None:
    args = parse_args()
    for scenario in selected_scenarios(args.scenario):
        genoio_matrix = None
        cyvcf2_matrix = None
        if args.backend in {"both", "genoio"}:
            genoio_matrix = benchmark_genoio_scenario(scenario, args)
        if scenario == "matrix-only" and args.kind == "geno" and not args.sparse and args.backend in {"both", "cyvcf2"}:
            cyvcf2_matrix = benchmark("cyvcf2_vcf", lambda: read_cyvcf2(args), args.repeats)
        elif args.backend in {"both", "cyvcf2"}:
            print_cyvcf2_skip(scenario, args)
        if (
            scenario == "matrix-only"
            and args.kind == "geno"
            and not args.sparse
            and not args.no_compare
            and genoio_matrix is not None
            and cyvcf2_matrix is not None
        ):
            compare_summaries("genoio_vcf_matrix_only", genoio_matrix, "cyvcf2_vcf", cyvcf2_matrix)


if __name__ == "__main__":
    main()
