#!/usr/bin/env python
# pattern: Mixed

from __future__ import annotations

import argparse
import os
import tempfile
from itertools import islice
from pathlib import Path

import numpy as np
from bench_common import benchmark, compare_summaries, positive_int, read_first_block

SCENARIOS = ("matrix-only", "with-variants", "sample-filtered", "genotype-filtered", "indexed-region")
KINDS = ("geno", "haplo")
BACKENDS = ("both", "all", "genoio", "bgen_reader", "bgen")
_last_variant_metadata_length: int | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark genoio BGEN dosage reads.")
    parser.add_argument("--prefix", type=Path, default=Path("data/chr22_hg38"))
    parser.add_argument("--max-variants", type=positive_int, default=1_000)
    parser.add_argument("--repeats", type=positive_int, default=3)
    parser.add_argument("--backend", choices=BACKENDS, default="both")
    parser.add_argument("--scenario", choices=[*SCENARIOS, "all"], default="matrix-only")
    parser.add_argument(
        "--kind",
        choices=KINDS,
        default="geno",
        help='Matrix kind to time. Defaults to genotype dosage; "haplo" times dense BGEN haplotype dosage.',
    )
    parser.add_argument("--region", default="22:20000000-21000000")
    parser.add_argument("--no-compare", action="store_true")
    return parser.parse_args()


def _read_options(kind: str) -> dict[str, object]:
    options: dict[str, object] = {
        "dosage": "dosage",
        "missing": "nan",
        "dtype": np.float32,
    }
    if kind == "haplo":
        options["kind"] = "haplo"
    return options


def _genoio_label(kind: str, scenario: str) -> str:
    suffix = scenario.replace("-", "_")
    kind_part = "haplo" if kind == "haplo" else ""
    parts = ["genoio_bgen", kind_part, suffix]
    return "_".join(part for part in parts if part)


def read_genoio_matrix_only(prefix: Path, max_variants: int, kind: str = "geno") -> np.ndarray:
    import genoio

    return read_first_block(
        genoio.bgen(prefix),
        max_variants,
        **_read_options(kind),
    )


def read_genoio_with_variants(prefix: Path, max_variants: int, kind: str = "geno") -> np.ndarray:
    import genoio

    global _last_variant_metadata_length
    matrix, variants = read_first_block(
        genoio.bgen(prefix),
        max_variants,
        return_variants=True,
        **_read_options(kind),
    )
    _last_variant_metadata_length = variants.height
    return matrix


def read_genoio_sample_filtered(prefix: Path, max_variants: int, kind: str = "geno") -> np.ndarray:
    import genoio

    sample_ids = _read_bgen_sample_ids(prefix.with_suffix(".sample"))
    keep_count = max(1, len(sample_ids) // 2)
    return read_first_block(
        genoio.bgen(prefix),
        max_variants,
        samples=sample_ids[:keep_count],
        **_read_options(kind),
    )


def read_genoio_genotype_filtered(prefix: Path, max_variants: int, kind: str = "geno") -> np.ndarray:
    import genoio

    return read_first_block(
        genoio.bgen(prefix),
        max_variants,
        variants=genoio.maf(min=0.01),
        **_read_options(kind),
    )


def read_genoio_indexed_region(prefix: Path, max_variants: int, region: str, kind: str = "geno") -> np.ndarray:
    import genoio

    matrix, variants = read_first_block(
        genoio.bgen(prefix),
        max_variants,
        variants=genoio.region(region),
        return_variants=True,
        **_read_options(kind),
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


def read_bgen_package_matrix_only(prefix: Path, max_variants: int) -> np.ndarray:
    matrix, _ = _read_bgen_package_alt_dosage_matrix(prefix, max_variants)
    return matrix


def read_bgen_package_with_variants(prefix: Path, max_variants: int) -> np.ndarray:
    matrix, variant_count = _read_bgen_package_alt_dosage_matrix(prefix, max_variants, read_variant_metadata=True)
    global _last_variant_metadata_length
    _last_variant_metadata_length = variant_count
    return matrix


def read_bgen_package_sample_filtered(prefix: Path, max_variants: int) -> np.ndarray:
    from bgen import BgenReader  # type: ignore[import-not-found]

    sample_path = _bgen_package_sample_path(prefix)
    with BgenReader(prefix.with_suffix(".bgen"), sample_path=sample_path, delay_parsing=True) as bgen:
        keep_count = max(1, len(bgen.samples) // 2)
        columns = [
            np.asarray(variant.alt_dosage[:keep_count], dtype=np.float32) for variant in islice(bgen, max_variants)
        ]
    return _stack_bgen_package_columns(columns, keep_count)


def read_bgen_package_genotype_filtered(prefix: Path, max_variants: int) -> np.ndarray:
    from bgen import BgenReader  # type: ignore[import-not-found]

    sample_path = _bgen_package_sample_path(prefix)
    with BgenReader(prefix.with_suffix(".bgen"), sample_path=sample_path, delay_parsing=True) as bgen:
        sample_count = len(bgen.samples)
        columns = []
        for variant in bgen:
            dosage = np.asarray(variant.alt_dosage, dtype=np.float32)
            if _dosage_maf(dosage) >= 0.01:
                columns.append(dosage)
                if len(columns) == max_variants:
                    break
    return _stack_bgen_package_columns(columns, sample_count)


def read_bgen_package_indexed_region(prefix: Path, max_variants: int, region: str) -> np.ndarray:
    from bgen import BgenReader  # type: ignore[import-not-found]

    sample_path = _bgen_package_sample_path(prefix)
    chrom, coords = region.split(":", 1)
    start, end = (int(value) for value in coords.split("-", 1))
    with BgenReader(prefix.with_suffix(".bgen"), sample_path=sample_path, delay_parsing=True) as bgen:
        sample_count = len(bgen.samples)
        columns = []
        for variant in islice(bgen.fetch(chrom, start, end), max_variants):
            columns.append(np.asarray(variant.alt_dosage, dtype=np.float32))
    return _stack_bgen_package_columns(columns, sample_count)


def _read_bgen_package_alt_dosage_matrix(
    prefix: Path,
    max_variants: int,
    *,
    read_variant_metadata: bool = False,
) -> tuple[np.ndarray, int]:
    from bgen import BgenReader  # type: ignore[import-not-found]

    sample_path = _bgen_package_sample_path(prefix)
    with BgenReader(prefix.with_suffix(".bgen"), sample_path=sample_path, delay_parsing=True) as bgen:
        sample_count = len(bgen.samples)
        columns = []
        variant_count = 0
        for variant in islice(bgen, max_variants):
            if read_variant_metadata:
                # Force the same metadata strings/coordinates that genoio returns
                # with `return_variants=True`; the values are not materialized here.
                _ = (variant.varid, variant.rsid, variant.chrom, variant.pos, variant.alleles)
            columns.append(np.asarray(variant.alt_dosage, dtype=np.float32))
            variant_count += 1
    matrix = _stack_bgen_package_columns(columns, sample_count)
    return matrix, variant_count


def _bgen_package_sample_path(prefix: Path) -> str:
    sample_path = prefix.with_suffix(".sample")
    return str(sample_path) if sample_path.exists() else ""


def _stack_bgen_package_columns(columns: list[np.ndarray], sample_count: int) -> np.ndarray:
    if not columns:
        return np.empty((sample_count, 0), dtype=np.float32)
    return np.stack(columns, axis=1).astype(np.float32, copy=False)


def _dosage_maf(dosage: np.ndarray) -> float:
    called = dosage[~np.isnan(dosage)]
    if called.size == 0:
        return 0.0
    allele_frequency = float(called.sum()) / (2.0 * float(called.size))
    return min(allele_frequency, 1.0 - allele_frequency)


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
    kind: str,
    prefix: Path,
    max_variants: int,
    repeats: int,
    region: str,
) -> np.ndarray:
    if scenario == "matrix-only":
        return benchmark(
            _genoio_label(kind, scenario),
            lambda: read_genoio_matrix_only(prefix, max_variants, kind),
            repeats,
        )
    if scenario == "with-variants":
        variant_metadata_length = None
        global _last_variant_metadata_length
        _last_variant_metadata_length = None

        def read_matrix() -> np.ndarray:
            nonlocal variant_metadata_length
            result = read_genoio_with_variants(prefix, max_variants, kind)
            variant_metadata_length = _last_variant_metadata_length
            return result

        matrix = benchmark(_genoio_label(kind, scenario), read_matrix, repeats)
        print(f"  variant_metadata length={variant_metadata_length}")
        return matrix
    if scenario == "sample-filtered":
        return benchmark(
            _genoio_label(kind, scenario),
            lambda: read_genoio_sample_filtered(prefix, max_variants, kind),
            repeats,
        )
    if scenario == "genotype-filtered":
        return benchmark(
            _genoio_label(kind, scenario),
            lambda: read_genoio_genotype_filtered(prefix, max_variants, kind),
            repeats,
        )
    if scenario == "indexed-region":
        return benchmark(
            _genoio_label(kind, scenario),
            lambda: read_genoio_indexed_region(prefix, max_variants, region, kind),
            repeats,
        )
    raise ValueError(f"unknown scenario: {scenario}")


def benchmark_bgen_reader_scenario(
    scenario: str,
    kind: str,
    prefix: Path,
    max_variants: int,
    repeats: int,
) -> np.ndarray | None:
    if kind == "haplo":
        print(
            f"skipped bgen_reader comparison for haplo {scenario}: "
            "comparison backend only computes diploid expected dosage"
        )
        return None
    if scenario == "matrix-only":
        return benchmark(
            "bgen_reader_expected_dosage",
            lambda: read_bgen_reader_expected_dosage(prefix, max_variants),
            repeats,
        )
    print(f"skipped bgen_reader comparison for {scenario}: benchmark only implements matrix-only expected dosage")
    return None


def benchmark_bgen_package_scenario(
    scenario: str,
    kind: str,
    prefix: Path,
    max_variants: int,
    repeats: int,
    region: str,
) -> np.ndarray | None:
    if kind == "haplo":
        print(f"skipped bgen package comparison for haplo {scenario}: backend returns diploid dosage")
        return None
    if scenario == "matrix-only":
        return benchmark(
            "bgen_package_alt_dosage",
            lambda: read_bgen_package_matrix_only(prefix, max_variants),
            repeats,
        )
    if scenario == "with-variants":
        variant_metadata_length = None
        global _last_variant_metadata_length
        _last_variant_metadata_length = None

        def read_matrix() -> np.ndarray:
            nonlocal variant_metadata_length
            result = read_bgen_package_with_variants(prefix, max_variants)
            variant_metadata_length = _last_variant_metadata_length
            return result

        matrix = benchmark("bgen_package_with_variants", read_matrix, repeats)
        print(f"  variant_metadata length={variant_metadata_length}")
        return matrix
    if scenario == "sample-filtered":
        return benchmark(
            "bgen_package_sample_filtered",
            lambda: read_bgen_package_sample_filtered(prefix, max_variants),
            repeats,
        )
    if scenario == "genotype-filtered":
        return benchmark(
            "bgen_package_genotype_filtered",
            lambda: read_bgen_package_genotype_filtered(prefix, max_variants),
            repeats,
        )
    if scenario == "indexed-region":
        return benchmark(
            "bgen_package_indexed_region",
            lambda: read_bgen_package_indexed_region(prefix, max_variants, region),
            repeats,
        )
    raise ValueError(f"unknown scenario: {scenario}")


def main() -> None:
    args = parse_args()
    for scenario in selected_scenarios(args.scenario):
        genoio_matrix = None
        bgen_reader_matrix = None
        bgen_package_matrix = None
        if args.backend in {"both", "all", "genoio"}:
            genoio_matrix = benchmark_genoio_scenario(
                scenario, args.kind, args.prefix, args.max_variants, args.repeats, args.region
            )
        if args.backend in {"both", "all", "bgen_reader"}:
            bgen_reader_matrix = benchmark_bgen_reader_scenario(
                scenario,
                args.kind,
                args.prefix,
                args.max_variants,
                args.repeats,
            )
        if args.backend in {"all", "bgen"}:
            bgen_package_matrix = benchmark_bgen_package_scenario(
                scenario,
                args.kind,
                args.prefix,
                args.max_variants,
                args.repeats,
                args.region,
            )
        if scenario == "matrix-only" and args.kind == "geno" and not args.no_compare and genoio_matrix is not None:
            if bgen_reader_matrix is not None:
                compare_summaries(
                    "genoio_bgen_matrix_only", genoio_matrix, "bgen_reader_expected_dosage", bgen_reader_matrix
                )
            if bgen_package_matrix is not None:
                compare_summaries(
                    "genoio_bgen_matrix_only", genoio_matrix, "bgen_package_alt_dosage", bgen_package_matrix
                )


if __name__ == "__main__":
    main()
