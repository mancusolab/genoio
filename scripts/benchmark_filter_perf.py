#!/usr/bin/env python
# pattern: Mixed

"""Benchmark Rust-side genotype filters against NumPy post-filtering.

The script is intentionally mixed: it owns CLI parsing, dataset I/O, and the
small NumPy mirror used to compare filter execution strategies.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
from bench_common import benchmark, nonnegative_float, positive_int, read_first_block

SOURCE_FORMATS = ("vcf", "bfile", "pfile", "bgen")
SCENARIOS = ("rust", "numpy")
WINDOW_MODES = ("retained", "source")
PREDICATES = ("maf", "mac", "missing_rate", "polymorphic")
FILTER_SHAPES = (
    "polymorphic",
    "mac",
    "mac_min",
    "mac_max",
    "maf",
    "maf_range",
    "missing_rate",
    "mac_missing_rate",
    "maf_missing_rate_polymorphic",
    "generic_fallback",
)
BASE_FILTERS = ("auto", "none", "biallelic", "not_multiallelic", "non_multiallelic", "snp")


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be nonnegative")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark Rust-side variant filtering against NumPy post-filtering.")
    parser.add_argument(
        "--source-format",
        choices=SOURCE_FORMATS,
        required=True,
        help="Input backend to benchmark.",
    )
    parser.add_argument(
        "--path",
        type=Path,
        required=True,
        help="Input path or prefix accepted by the selected genoio source constructor.",
    )
    parser.add_argument("--max-variants", type=positive_int, default=1_000)
    parser.add_argument("--repeats", type=positive_int, default=3)
    parser.add_argument("--scenario", choices=[*SCENARIOS, "both"], default="both")
    parser.add_argument(
        "--window-mode",
        choices=WINDOW_MODES,
        default="retained",
        help=(
            "retained reads until each scenario returns up to max passing variants; "
            "source filters the same first max source/base-filter variants"
        ),
    )
    parser.add_argument(
        "--predicate",
        choices=PREDICATES,
        default="maf",
        help="Legacy genotype-stat predicate to benchmark when --filter-shape is omitted.",
    )
    parser.add_argument(
        "--filter-shape",
        choices=FILTER_SHAPES,
        default=None,
        help="Named genotype-filter shape to benchmark.",
    )
    parser.add_argument(
        "--dosage",
        choices=["auto", "hardcall", "dosage"],
        default="auto",
        help="Genotype value source. BGEN uses dosage when auto is selected.",
    )
    parser.add_argument(
        "--base-filter",
        choices=BASE_FILTERS,
        default="auto",
        help="Metadata filter applied before genotype-stat filtering. Auto uses biallelic for VCF.",
    )
    parser.add_argument("--maf-min", type=nonnegative_float, default=0.01)
    parser.add_argument("--maf-max", type=nonnegative_float, default=None)
    parser.add_argument("--mac-min", type=nonnegative_int, default=1)
    parser.add_argument("--mac-max", type=nonnegative_int, default=None)
    parser.add_argument("--missing-rate-max", type=nonnegative_float, default=None)
    return parser.parse_args()


def selected_scenarios(scenario: str) -> tuple[str, ...]:
    if scenario == "both":
        return SCENARIOS
    return (scenario,)


def active_filter_shape(args: argparse.Namespace) -> str:
    return args.filter_shape or args.predicate


def read_options(args: argparse.Namespace) -> dict[str, object]:
    dosage = "dosage" if args.dosage == "auto" and args.source_format == "bgen" else args.dosage
    options: dict[str, object] = {
        "dtype": np.float32,
        "missing": "nan",
    }
    if dosage != "auto":
        options["dosage"] = dosage
    return options


def dataset_for_args(args: argparse.Namespace) -> Any:
    import genoio

    constructors = {
        "vcf": genoio.vcf,
        "bfile": genoio.bfile,
        "pfile": genoio.pfile,
        "bgen": genoio.bgen,
    }
    return constructors[args.source_format](args.path)


def base_filter(args: argparse.Namespace) -> Any:
    import genoio

    filter_name = args.base_filter
    if filter_name == "auto":
        filter_name = "biallelic" if args.source_format == "vcf" else "none"
    if filter_name == "none":
        return None
    if filter_name in {"biallelic", "not_multiallelic", "non_multiallelic"}:
        return genoio.biallelic()
    if filter_name == "snp":
        return genoio.snp()
    raise ValueError(f"unknown base filter: {filter_name}")


def stat_filter(args: argparse.Namespace) -> Any:
    import genoio

    shape = active_filter_shape(args)
    if shape in {"maf", "maf_min"}:
        expr = genoio.maf(min=args.maf_min, max=args.maf_max)
    elif shape == "maf_range":
        if args.maf_max is None:
            raise ValueError("--maf-max is required when --filter-shape maf_range")
        expr = genoio.maf(min=args.maf_min, max=args.maf_max)
    elif shape in {"mac", "mac_min"}:
        expr = genoio.mac(min=args.mac_min, max=args.mac_max)
    elif shape == "mac_max":
        if args.mac_max is None:
            raise ValueError("--mac-max is required when --filter-shape mac_max")
        expr = genoio.mac(max=args.mac_max)
    elif shape == "missing_rate":
        if args.missing_rate_max is None:
            raise ValueError("--missing-rate-max is required when benchmarking missing_rate")
        expr = genoio.missing_rate(max=args.missing_rate_max)
    elif shape == "polymorphic":
        expr = genoio.polymorphic()
    elif shape == "mac_missing_rate":
        if args.missing_rate_max is None:
            raise ValueError("--missing-rate-max is required when --filter-shape mac_missing_rate")
        expr = genoio.mac(min=args.mac_min, max=args.mac_max) & genoio.missing_rate(max=args.missing_rate_max)
    elif shape == "maf_missing_rate_polymorphic":
        if args.maf_max is None:
            raise ValueError("--maf-max is required when --filter-shape maf_missing_rate_polymorphic")
        if args.missing_rate_max is None:
            raise ValueError("--missing-rate-max is required when --filter-shape maf_missing_rate_polymorphic")
        expr = (
            genoio.maf(min=args.maf_min, max=args.maf_max)
            & genoio.missing_rate(max=args.missing_rate_max)
            & genoio.polymorphic()
        )
    elif shape == "generic_fallback":
        if args.missing_rate_max is None:
            raise ValueError("--missing-rate-max is required when --filter-shape generic_fallback")
        expr = genoio.mac(min=args.mac_min, max=args.mac_max) | genoio.missing_rate(max=args.missing_rate_max)
    else:
        raise ValueError(f"unknown filter shape: {shape}")
    if args.filter_shape is None and args.predicate != "missing_rate" and args.missing_rate_max is not None:
        expr = expr & genoio.missing_rate(max=args.missing_rate_max)
    return expr


def combined_filter(args: argparse.Namespace) -> Any:
    base = base_filter(args)
    stats = stat_filter(args)
    return stats if base is None else base & stats


def source_window_filter(args: argparse.Namespace, variant_ids: list[str]) -> Any:
    import genoio

    validate_source_window_variant_ids(variant_ids)
    base = base_filter(args)
    stats = stat_filter(args)
    source_ids = genoio.id_in(variant_ids)
    expr = source_ids & stats
    return expr if base is None else base & expr


def validate_source_window_variant_ids(variant_ids: list[str]) -> None:
    if len(variant_ids) != len(set(variant_ids)):
        raise ValueError(
            "source-window Rust filtering cannot use duplicate variant IDs; "
            "use --window-mode source --scenario numpy for duplicate-ID sources"
        )


def read_rust_filtered(args: argparse.Namespace) -> np.ndarray:
    dataset = dataset_for_args(args)
    return read_first_block(
        dataset,
        args.max_variants,
        variants=combined_filter(args),
        **read_options(args),
    )


def read_rust_source_window_filtered(args: argparse.Namespace, variant_ids: list[str]) -> np.ndarray:
    dataset = dataset_for_args(args)
    return read_first_block(
        dataset,
        args.max_variants,
        variants=source_window_filter(args, variant_ids),
        **read_options(args),
    )


def read_base_block(args: argparse.Namespace) -> np.ndarray:
    dataset = dataset_for_args(args)
    return read_first_block(
        dataset,
        args.max_variants,
        variants=base_filter(args),
        **read_options(args),
    )


def read_base_block_with_variants(args: argparse.Namespace) -> tuple[np.ndarray, Any]:
    dataset = dataset_for_args(args)
    return read_first_block(
        dataset,
        args.max_variants,
        variants=base_filter(args),
        return_variants=True,
        **read_options(args),
    )


def source_window_variant_ids(args: argparse.Namespace) -> list[str]:
    """Return IDs from the first source/base-filter block.

    Rust source-window filtering uses these IDs to force the public retained
    window API to operate on the same source variants as NumPy. This only works
    when source IDs are unique.
    """
    try:
        _, variants = read_base_block_with_variants(args)
    except StopIteration:
        return []
    return variant_ids_from_frame(variants)


def variant_ids_from_frame(variants: Any) -> list[str]:
    if hasattr(variants, "get_column"):
        return [str(value) for value in variants.get_column("id").to_list()]
    return [str(value) for value in variants["id"]]


def read_numpy_source_window_postfiltered(args: argparse.Namespace) -> np.ndarray:
    """Filter the first source/base-filter block in NumPy."""
    matrix = np.asarray(read_base_block(args))
    mask = numpy_variant_mask_for_args(matrix, args)
    return matrix[:, mask]


def read_numpy_retained_postfiltered(args: argparse.Namespace) -> np.ndarray:
    """Read source blocks until NumPy has up to ``max_variants`` passing columns."""
    dataset = dataset_for_args(args)
    pieces: list[np.ndarray] = []
    empty: np.ndarray | None = None
    retained = 0
    iterator = dataset.iter_blocks(
        args.max_variants,
        variants=base_filter(args),
        **read_options(args),
    )
    with iterator as blocks:
        for block in blocks:
            matrix = np.asarray(block)
            if empty is None:
                empty = matrix[:, :0]
            mask = numpy_variant_mask_for_args(matrix, args)
            filtered = matrix[:, mask]
            if filtered.shape[1] == 0:
                continue
            remaining = args.max_variants - retained
            pieces.append(filtered[:, :remaining])
            retained += min(filtered.shape[1], remaining)
            if retained >= args.max_variants:
                break
    if pieces:
        return np.concatenate(pieces, axis=1)
    if empty is not None:
        return empty
    return np.empty((0, 0), dtype=np.float32)


def numpy_variant_mask_for_args(matrix: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    return numpy_variant_mask(
        matrix,
        filter_shape=active_filter_shape(args),
        maf_min=args.maf_min,
        maf_max=args.maf_max,
        mac_min=args.mac_min,
        mac_max=args.mac_max,
        missing_rate_max=args.missing_rate_max,
    )


def numpy_variant_mask(
    matrix: np.ndarray,
    *,
    filter_shape: str | None = None,
    predicate: str | None = None,
    maf_min: float,
    maf_max: float | None,
    mac_min: int,
    mac_max: int | None,
    missing_rate_max: float | None,
) -> np.ndarray:
    """Return the NumPy mask for the requested benchmark filter shape."""
    shape = filter_shape or predicate
    if shape is None:
        raise ValueError("filter_shape is required")
    missing = np.isnan(matrix)
    called = np.count_nonzero(~missing, axis=0)
    total_samples = matrix.shape[0]
    allele_sum = np.nansum(matrix, axis=0, dtype=np.float64)

    with np.errstate(invalid="ignore", divide="ignore"):
        af = allele_sum / (2.0 * called)
    maf = np.minimum(af, 1.0 - af)
    mac = np.minimum(allele_sum, (2.0 * called) - allele_sum)
    missing_rate = np.zeros_like(maf, dtype=np.float64)
    if total_samples:
        missing_rate = np.count_nonzero(missing, axis=0) / total_samples

    if shape in {"maf", "maf_min"}:
        keep = called > 0
        keep &= maf >= maf_min
        if maf_max is not None:
            keep &= maf <= maf_max
    elif shape == "maf_range":
        if maf_max is None:
            raise ValueError("maf_range filter shape requires maf_max")
        keep = called > 0
        keep &= maf >= maf_min
        keep &= maf <= maf_max
    elif shape in {"mac", "mac_min"}:
        keep = called > 0
        keep &= mac >= mac_min
        if mac_max is not None:
            keep &= mac <= mac_max
    elif shape == "mac_max":
        if mac_max is None:
            raise ValueError("mac_max filter shape requires mac_max")
        keep = called > 0
        keep &= mac <= mac_max
    elif shape == "missing_rate":
        if missing_rate_max is None:
            raise ValueError("missing_rate predicate requires missing_rate_max")
        keep = missing_rate <= missing_rate_max
    elif shape == "polymorphic":
        called_alleles = 2.0 * called
        keep = called > 0
        keep &= allele_sum > 0.0
        keep &= allele_sum < called_alleles
    elif shape == "mac_missing_rate":
        if missing_rate_max is None:
            raise ValueError("mac_missing_rate filter shape requires missing_rate_max")
        keep = called > 0
        keep &= mac >= mac_min
        if mac_max is not None:
            keep &= mac <= mac_max
        keep &= missing_rate <= missing_rate_max
    elif shape == "maf_missing_rate_polymorphic":
        if maf_max is None:
            raise ValueError("maf_missing_rate_polymorphic filter shape requires maf_max")
        if missing_rate_max is None:
            raise ValueError("maf_missing_rate_polymorphic filter shape requires missing_rate_max")
        called_alleles = 2.0 * called
        keep = called > 0
        keep &= allele_sum > 0.0
        keep &= allele_sum < called_alleles
        keep &= maf >= maf_min
        keep &= maf <= maf_max
        keep &= missing_rate <= missing_rate_max
    elif shape == "generic_fallback":
        if missing_rate_max is None:
            raise ValueError("generic_fallback filter shape requires missing_rate_max")
        keep = called > 0
        keep &= mac >= mac_min
        if mac_max is not None:
            keep &= mac <= mac_max
        keep |= missing_rate <= missing_rate_max
    else:
        raise ValueError(f"unknown filter shape: {shape}")
    if filter_shape is None and predicate != "missing_rate" and missing_rate_max is not None:
        keep &= missing_rate <= missing_rate_max
    return keep


def benchmark_scenario(scenario: str, args: argparse.Namespace) -> np.ndarray:
    label = f"genoio_{args.source_format}_{active_filter_shape(args)}_{args.window_mode}_window_filter_{scenario}"
    if scenario == "rust":
        if args.window_mode == "source":
            variant_ids = source_window_variant_ids(args)
            return benchmark(
                label,
                lambda: read_rust_source_window_filtered(args, variant_ids),
                args.repeats,
            )
        return benchmark(label, lambda: read_rust_filtered(args), args.repeats)
    if scenario == "numpy":
        if args.window_mode == "source":
            return benchmark(label, lambda: read_numpy_source_window_postfiltered(args), args.repeats)
        return benchmark(label, lambda: read_numpy_retained_postfiltered(args), args.repeats)
    raise ValueError(f"unknown scenario: {scenario}")


def main() -> None:
    args = parse_args()
    if args.window_mode == "source":
        print("window_semantics", "source=max source/base-filter variants before genotype filtering")
    else:
        print("window_semantics", "retained=max retained variants after genotype filtering")
    for scenario in selected_scenarios(args.scenario):
        benchmark_scenario(scenario, args)


if __name__ == "__main__":
    main()
