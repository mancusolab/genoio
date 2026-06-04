#!/usr/bin/env python
# pattern: Mixed

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from bench_common import benchmark, compare_summaries, plink2_prefix_with_uncompressed_pvar, positive_int

SCENARIOS = ("matrix-only", "with-variants", "sample-filtered", "genotype-filtered")
_last_variant_metadata_length: int | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare genoio PLINK2 reads against pgenlib.PgenReader.")
    parser.add_argument("--prefix", type=Path, default=Path("data/chr22_hg38"))
    parser.add_argument("--max-variants", type=positive_int, default=1_000)
    parser.add_argument("--repeats", type=positive_int, default=3)
    parser.add_argument("--backend", choices=["both", "genoio", "pgenlib"], default="both")
    parser.add_argument("--scenario", choices=[*SCENARIOS, "all"], default="matrix-only")
    parser.add_argument(
        "--pgenlib-path",
        type=Path,
        default=None,
        help="Optional path to plink-ng/2.0/Python when pgenlib is built in-place but not installed.",
    )
    parser.add_argument("--no-compare", action="store_true")
    return parser.parse_args()


def read_genoio_matrix_only(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    return next(
        genoio.pfile(prefix).blocks(
            max_variants,
            missing="nan",
            dtype=np.float32,
        )
    )


def read_genoio_with_variants(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    global _last_variant_metadata_length
    matrix, variants = next(
        genoio.pfile(prefix).blocks(
            max_variants,
            missing="nan",
            dtype=np.float32,
            return_variants=True,
        )
    )
    _last_variant_metadata_length = variants.height
    return matrix


def read_genoio_sample_filtered(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    sample_ids = _read_psam_sample_ids(prefix.with_suffix(".psam"))
    keep_count = max(1, len(sample_ids) // 2)
    return next(
        genoio.pfile(prefix).blocks(
            max_variants,
            missing="nan",
            dtype=np.float32,
            samples=sample_ids[:keep_count],
        )
    )


def read_genoio_genotype_filtered(prefix: Path, max_variants: int) -> np.ndarray:
    import genoio

    return next(
        genoio.pfile(prefix).blocks(
            max_variants,
            missing="nan",
            dtype=np.float32,
            variants=genoio.maf(min=0.01),
        )
    )


def _read_psam_sample_ids(path: Path) -> list[str]:
    with path.open() as handle:
        header: list[str] | None = None
        sample_index = 0
        sample_ids: list[str] = []
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if stripped.startswith("#"):
                header = [field.removeprefix("#") for field in fields]
                sample_column = "IID" if "IID" in header else header[0]
                sample_index = header.index(sample_column)
                continue
            sample_ids.append(fields[sample_index])
    if sample_ids:
        return sample_ids
    raise RuntimeError(f"{path} does not contain sample IDs")


def import_pgenlib(pgenlib_path: Path | None):
    if pgenlib_path is not None:
        sys.path.insert(0, str(pgenlib_path))
        sys.path.insert(0, str(pgenlib_path / "src"))
    try:
        import pgenlib  # type: ignore[import-not-found]
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


def selected_scenarios(scenario: str) -> tuple[str, ...]:
    if scenario == "all":
        return SCENARIOS
    return (scenario,)


def benchmark_genoio_scenario(scenario: str, prefix: Path, max_variants: int, repeats: int) -> np.ndarray:
    if scenario == "matrix-only":
        return benchmark(
            "genoio_plink2_matrix_only",
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
            if isinstance(result, tuple):
                matrix, variant_metadata_length = result
                return matrix
            variant_metadata_length = _last_variant_metadata_length
            return result

        matrix = benchmark("genoio_plink2_with_variants", read_matrix, repeats)
        print(f"  variant_metadata length={variant_metadata_length}")
        return matrix
    if scenario == "sample-filtered":
        return benchmark(
            "genoio_plink2_sample_filtered",
            lambda: read_genoio_sample_filtered(prefix, max_variants),
            repeats,
        )
    if scenario == "genotype-filtered":
        return benchmark(
            "genoio_plink2_genotype_filtered",
            lambda: read_genoio_genotype_filtered(prefix, max_variants),
            repeats,
        )
    raise ValueError(f"unknown scenario: {scenario}")


def print_pgenlib_skip(scenario: str) -> None:
    message = f"skipped pgenlib comparison for {scenario}: pgenlib does not provide the same metadata/filter contract"
    print(message)


def main() -> None:
    args = parse_args()
    for scenario in selected_scenarios(args.scenario):
        genoio_matrix = None
        pgenlib_matrix = None
        if args.backend in {"both", "genoio"}:
            with plink2_prefix_with_uncompressed_pvar(args.prefix) as genoio_prefix:
                genoio_matrix = benchmark_genoio_scenario(scenario, genoio_prefix, args.max_variants, args.repeats)
        if scenario == "matrix-only" and args.backend in {"both", "pgenlib"}:
            pgenlib_matrix = benchmark("pgenlib_pgenreader", lambda: read_pgenlib(args), args.repeats)
        elif scenario != "matrix-only" and args.backend in {"both", "pgenlib"}:
            print_pgenlib_skip(scenario)
        if (
            scenario == "matrix-only"
            and not args.no_compare
            and genoio_matrix is not None
            and pgenlib_matrix is not None
        ):
            compare_summaries("genoio_plink2_matrix_only", genoio_matrix, "pgenlib_pgenreader", pgenlib_matrix)


if __name__ == "__main__":
    main()
