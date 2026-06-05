#!/usr/bin/env python
# pattern: Mixed

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import numpy as np
from bench_common import benchmark, compare_summaries, positive_int

SCENARIOS = ("matrix-only", "with-variants", "sample-filtered", "genotype-filtered", "indexed-region")
_last_variant_metadata_length: int | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark genoio BGEN dosage reads.")
    parser.add_argument("--prefix", type=Path, default=Path("data/chr22_hg38"))
    parser.add_argument("--max-variants", type=positive_int, default=1_000)
    parser.add_argument("--repeats", type=positive_int, default=3)
    parser.add_argument("--backend", choices=["both", "genoio", "bgen_reader"], default="both")
    parser.add_argument("--scenario", choices=[*SCENARIOS, "all"], default="matrix-only")
    parser.add_argument("--region", default="22:20000000-21000000")
    parser.add_argument("--no-compare", action="store_true")
    return parser.parse_args()


def read_genoio_matrix_only(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    return next(
        genoio.bgen(prefix).blocks(
            max_variants,
            dosage="dosage",
            missing="nan",
            dtype=np.float32,
        )
    )


def read_genoio_with_variants(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    global _last_variant_metadata_length
    matrix, variants = next(
        genoio.bgen(prefix).blocks(
            max_variants,
            dosage="dosage",
            missing="nan",
            dtype=np.float32,
            return_variants=True,
        )
    )
    _last_variant_metadata_length = variants.height
    return matrix


def read_genoio_sample_filtered(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    sample_ids = _read_bgen_sample_ids(prefix.with_suffix(".sample"))
    keep_count = max(1, len(sample_ids) // 2)
    return next(
        genoio.bgen(prefix).blocks(
            max_variants,
            dosage="dosage",
            missing="nan",
            dtype=np.float32,
            samples=sample_ids[:keep_count],
        )
    )


def read_genoio_genotype_filtered(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    return next(
        genoio.bgen(prefix).blocks(
            max_variants,
            dosage="dosage",
            missing="nan",
            dtype=np.float32,
            variants=genoio.maf(min=0.01),
        )
    )


def read_genoio_indexed_region(prefix: Path, max_variants: int, region: str) -> np.ndarray:
    import genoio

    matrix, variants = next(
        genoio.bgen(prefix).blocks(
            max_variants,
            dosage="dosage",
            missing="nan",
            dtype=np.float32,
            variants=genoio.region(region),
            return_variants=True,
        )
    )
    if variants.height:
        chrom, coords = region.split(":", 1)
        start, end = (int(value) for value in coords.split("-", 1))
        observed_min = int(variants["pos"].min())
        observed_max = int(variants["pos"].max())
        if variants["chrom"].unique().to_list() != [chrom] or observed_min < start or observed_max > end:
            raise RuntimeError(f"indexed region result escaped requested region: {region}")
    return matrix


def read_bgen_reader_expected_dosage(prefix: Path, max_variants: int) -> np.ndarray:
    _configure_bgen_reader_cache()
    from bgen_reader import open_bgen  # type: ignore[import-not-found]
    from cbgen import bgen_file  # type: ignore[import-not-found]

    bgen_path = prefix.with_suffix(".bgen")
    bgen = open_bgen(bgen_path, verbose=False)
    columns = []
    with bgen_file(str(bgen_path)) as cbgen:
        for variant_index in range(max_variants):
            probabilities = np.asarray(cbgen.read_probability(int(bgen._vaddr[variant_index]), 64))
            weights = _expected_dosage_weights(probabilities.shape[1])
            columns.append(probabilities @ weights)
    if not columns:
        return np.empty((len(bgen.samples), 0), dtype=np.float32)
    return np.stack(columns, axis=1).astype(np.float32, copy=False)


def _configure_bgen_reader_cache() -> None:
    cache = Path(tempfile.gettempdir()) / "genoio-bgen-reader-cache"
    os.environ.setdefault("BGEN_CACHE_HOME", str(cache))
    os.environ.setdefault("BGEN_READER_CACHE_HOME", str(cache))


def _expected_dosage_weights(probability_count: int) -> np.ndarray:
    if probability_count == 3:
        return np.array([0.0, 1.0, 2.0], dtype=np.float64)
    if probability_count == 4:
        # Phased biallelic diploid BGEN stores per-haplotype allele
        # probabilities: hap0 A0, hap0 A1, hap1 A0, hap1 A1.
        return np.array([0.0, 1.0, 0.0, 1.0], dtype=np.float64)
    raise RuntimeError(f"unsupported bgen_reader probability count: {probability_count}")


def _read_bgen_sample_ids(path: Path) -> list[str]:
    with path.open() as handle:
        header = handle.readline().split()
        if not header:
            raise RuntimeError(f"{path} does not contain a sample header")
        sample_index = header.index("ID_2") if "ID_2" in header else 0
        handle.readline()
        sample_ids = [fields[sample_index] for line in handle if (fields := line.split())]
    if sample_ids:
        return sample_ids
    raise RuntimeError(f"{path} does not contain sample IDs")


def selected_scenarios(scenario: str) -> tuple[str, ...]:
    if scenario == "all":
        return SCENARIOS
    return (scenario,)


def benchmark_genoio_scenario(
    scenario: str,
    prefix: Path,
    max_variants: int,
    repeats: int,
    region: str,
) -> np.ndarray:
    if scenario == "matrix-only":
        return benchmark(
            "genoio_bgen_matrix_only",
            lambda: read_genoio_matrix_only(prefix, max_variants),
            repeats,
        )
    if scenario == "with-variants":
        variant_metadata_length = None
        global _last_variant_metadata_length
        _last_variant_metadata_length = None

        def read_matrix() -> np.ndarray:
            nonlocal variant_metadata_length
            result = read_genoio_with_variants(prefix, max_variants)
            variant_metadata_length = _last_variant_metadata_length
            return result

        matrix = benchmark("genoio_bgen_with_variants", read_matrix, repeats)
        print(f"  variant_metadata length={variant_metadata_length}")
        return matrix
    if scenario == "sample-filtered":
        return benchmark(
            "genoio_bgen_sample_filtered",
            lambda: read_genoio_sample_filtered(prefix, max_variants),
            repeats,
        )
    if scenario == "genotype-filtered":
        return benchmark(
            "genoio_bgen_genotype_filtered",
            lambda: read_genoio_genotype_filtered(prefix, max_variants),
            repeats,
        )
    if scenario == "indexed-region":
        return benchmark(
            "genoio_bgen_indexed_region",
            lambda: read_genoio_indexed_region(prefix, max_variants, region),
            repeats,
        )
    raise ValueError(f"unknown scenario: {scenario}")


def benchmark_bgen_reader_scenario(scenario: str, prefix: Path, max_variants: int, repeats: int) -> np.ndarray | None:
    if scenario == "matrix-only":
        return benchmark(
            "bgen_reader_expected_dosage",
            lambda: read_bgen_reader_expected_dosage(prefix, max_variants),
            repeats,
        )
    print(f"skipped bgen_reader comparison for {scenario}: benchmark only implements matrix-only expected dosage")
    return None


def main() -> None:
    args = parse_args()
    for scenario in selected_scenarios(args.scenario):
        genoio_matrix = None
        bgen_reader_matrix = None
        if args.backend in {"both", "genoio"}:
            genoio_matrix = benchmark_genoio_scenario(
                scenario, args.prefix, args.max_variants, args.repeats, args.region
            )
        if args.backend in {"both", "bgen_reader"}:
            bgen_reader_matrix = benchmark_bgen_reader_scenario(
                scenario,
                args.prefix,
                args.max_variants,
                args.repeats,
            )
        if (
            scenario == "matrix-only"
            and not args.no_compare
            and genoio_matrix is not None
            and bgen_reader_matrix is not None
        ):
            compare_summaries(
                "genoio_bgen_matrix_only", genoio_matrix, "bgen_reader_expected_dosage", bgen_reader_matrix
            )


if __name__ == "__main__":
    main()
