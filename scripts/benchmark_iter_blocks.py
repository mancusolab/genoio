#!/usr/bin/env python
# pattern: Mixed (unavoidable)
# Reason: Timed dataset I/O and low-overhead shape accounting must share one benchmark boundary.

"""Benchmark sustained ``Dataset.iter_blocks`` throughput.

The module is intentionally mixed because the timed operation must keep the
dataset I/O loop and its low-overhead correctness accounting at one benchmark
boundary. It never concatenates yielded matrices or scans their values.
"""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from bench_common import positive_int

SOURCE_FORMATS = ("vcf", "bfile", "pfile", "bgen")
SCENARIOS = ("matrix-only", "with-variants", "genotype-filtered")
KINDS = ("hardcall", "dosage", "haplo-hardcall", "haplo-dosage")
DEFAULT_BLOCK_SIZES = (128, 512, 2048)


@dataclass(frozen=True)
class ScanSummary:
    """Shape accounting for one fixed-prefix streaming scan."""

    blocks: int
    variants: int
    rows: int


@dataclass(frozen=True)
class BenchmarkResult:
    """Timing samples and validated output for one benchmark case."""

    scenario: str
    block_size: int
    summary: ScanSummary
    seconds: tuple[float, ...]

    @property
    def median_seconds(self) -> float:
        return statistics.median(self.seconds)

    @property
    def min_seconds(self) -> float:
        return min(self.seconds)

    @property
    def variants_per_second(self) -> float:
        return self.summary.variants / self.median_seconds


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be nonnegative")
    return parsed


def block_sizes(value: str) -> tuple[int, ...]:
    try:
        parsed = tuple(int(part) for part in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError("block sizes must be comma-separated integers") from error
    if not parsed or any(size < 1 for size in parsed):
        raise argparse.ArgumentTypeError("block sizes must be positive integers")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark sustained genoio iter_blocks throughput across fixed-size "
            "blocks while holding the retained-variant workload constant."
        )
    )
    parser.add_argument("--source-format", choices=SOURCE_FORMATS, required=True)
    parser.add_argument(
        "--path",
        type=Path,
        required=True,
        help="Input path or prefix accepted by the selected genoio constructor.",
    )
    parser.add_argument(
        "--label",
        required=True,
        help="Explicit revision/run label used in human and JSON output.",
    )
    parser.add_argument(
        "--block-sizes",
        type=block_sizes,
        default=DEFAULT_BLOCK_SIZES,
        help="Comma-separated block sizes. Each must divide --max-variants.",
    )
    parser.add_argument(
        "--max-variants",
        type=positive_int,
        default=16_384,
        help="Exact number of retained variants consumed by every timed scan.",
    )
    parser.add_argument("--repeats", type=positive_int, default=5)
    parser.add_argument("--warmups", type=nonnegative_int, default=1)
    parser.add_argument("--scenario", choices=[*SCENARIOS, "all"], default="matrix-only")
    parser.add_argument(
        "--kind",
        choices=KINDS,
        default="hardcall",
        help="Matrix representation. BGEN dosage benchmarks should use dosage.",
    )
    parser.add_argument("--sparse", action="store_true", help="Request sparse CSC matrices.")
    parser.add_argument("--output-json", type=Path, default=None)
    return parser.parse_args()


def validate_workload(requested_block_sizes: tuple[int, ...], max_variants: int) -> None:
    """Reject sweeps that do not perform directly comparable work."""
    if len(requested_block_sizes) != len(set(requested_block_sizes)):
        raise ValueError("block sizes must be distinct")
    incompatible = [size for size in requested_block_sizes if max_variants % size]
    if incompatible:
        rendered = ", ".join(str(size) for size in incompatible)
        raise ValueError(f"--max-variants must be divisible by every block size; incompatible: {rendered}")


def dataset_for_source(source_format: str, path: Path) -> Any:
    import genoio

    constructors = {
        "vcf": genoio.vcf,
        "bfile": genoio.bfile,
        "pfile": genoio.pfile,
        "bgen": genoio.bgen,
    }
    return constructors[source_format](path)


def selected_scenarios(scenario: str) -> tuple[str, ...]:
    if scenario == "all":
        return SCENARIOS
    return (scenario,)


def read_options_for_case(
    source_format: str,
    scenario: str,
    kind: str,
    sparse: bool,
) -> dict[str, object]:
    """Build public read options for one benchmark scenario."""
    import genoio

    options: dict[str, object] = {
        "dtype": np.float32,
        "missing": "raise" if sparse else "nan",
    }
    if sparse:
        options["sparse"] = "csc"
    if kind == "dosage":
        options["dosage"] = "dosage"
    elif kind == "haplo-hardcall":
        options["kind"] = "haplo"
        options["dosage"] = "hardcall"
    elif kind == "haplo-dosage":
        options["kind"] = "haplo"
        options["dosage"] = "dosage"

    variant_filter = genoio.biallelic() if source_format == "vcf" else None
    if scenario == "with-variants":
        options["return_variants"] = True
    if scenario == "genotype-filtered":
        genotype_filter = genoio.maf(min=0.01)
        variant_filter = genotype_filter if variant_filter is None else variant_filter & genotype_filter
    if variant_filter is not None:
        options["variants"] = variant_filter
    return options


def consume_blocks(
    dataset: Any,
    *,
    block_size: int,
    max_variants: int,
    read_options: dict[str, object],
) -> ScanSummary:
    """Consume and validate an exact retained-variant prefix without copying it."""
    return_variants = read_options.get("return_variants") is True
    block_count = 0
    variant_count = 0
    row_count: int | None = None

    with dataset.iter_blocks(block_size, **read_options) as blocks:
        for yielded in blocks:
            if return_variants:
                matrix, variants = yielded
            else:
                matrix = yielded
                variants = None

            rows, width = matrix.shape
            if width < 1 or width > block_size:
                raise RuntimeError(f"iter_blocks yielded width {width}; expected between 1 and {block_size}")
            if row_count is None:
                row_count = rows
            elif rows != row_count:
                raise RuntimeError(f"matrix rows changed between blocks: {row_count} then {rows}")

            if variants is not None and variants.height != width:
                raise RuntimeError(f"variant metadata rows ({variants.height}) do not match matrix columns ({width})")

            block_count += 1
            variant_count += width
            if variant_count > max_variants:
                raise RuntimeError(
                    "iter_blocks exceeded the exact workload; use block sizes that divide --max-variants"
                )
            if variant_count == max_variants:
                assert row_count is not None
                return ScanSummary(
                    blocks=block_count,
                    variants=variant_count,
                    rows=row_count,
                )

    raise ValueError(f"source ended after {variant_count} retained variants; benchmark requires {max_variants}")


def benchmark_case(
    dataset: Any,
    *,
    scenario: str,
    block_size: int,
    max_variants: int,
    read_options: dict[str, object],
    repeats: int,
    warmups: int,
) -> BenchmarkResult:
    """Run warm-cache timing samples for one scenario and block size."""
    expected: ScanSummary | None = None
    for _ in range(warmups):
        observed = consume_blocks(
            dataset,
            block_size=block_size,
            max_variants=max_variants,
            read_options=read_options,
        )
        expected = _validate_summary(expected, observed)

    samples: list[float] = []
    for _ in range(repeats):
        start = time.perf_counter()
        observed = consume_blocks(
            dataset,
            block_size=block_size,
            max_variants=max_variants,
            read_options=read_options,
        )
        samples.append(time.perf_counter() - start)
        expected = _validate_summary(expected, observed)

    assert expected is not None
    return BenchmarkResult(
        scenario=scenario,
        block_size=block_size,
        summary=expected,
        seconds=tuple(samples),
    )


def _validate_summary(expected: ScanSummary | None, observed: ScanSummary) -> ScanSummary:
    if expected is not None and observed != expected:
        raise RuntimeError(f"stream summary changed between runs: {expected} then {observed}")
    return observed


def print_result(result: BenchmarkResult) -> None:
    print(f"scenario={result.scenario} block_size={result.block_size}")
    print(
        "  stream",
        f"blocks={result.summary.blocks}",
        f"variants={result.summary.variants}",
        f"rows={result.summary.rows}",
    )
    print(
        "  time",
        f"median={result.median_seconds:.4f}s",
        f"min={result.min_seconds:.4f}s",
        f"variants/s={result.variants_per_second:,.0f}",
        "runs=" + " ".join(f"{value:.4f}" for value in result.seconds),
    )


def result_record(result: BenchmarkResult) -> dict[str, object]:
    return {
        "scenario": result.scenario,
        "block_size": result.block_size,
        **asdict(result.summary),
        "median_seconds": result.median_seconds,
        "min_seconds": result.min_seconds,
        "variants_per_second": result.variants_per_second,
        "seconds": list(result.seconds),
    }


def write_json_report(
    output_path: Path,
    args: argparse.Namespace,
    results: list[BenchmarkResult],
) -> None:
    import genoio

    report = {
        "label": args.label,
        "genoio_version": genoio.__version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "source_format": args.source_format,
        "path": str(args.path.resolve()),
        "max_variants": args.max_variants,
        "block_sizes": list(args.block_sizes),
        "repeats": args.repeats,
        "warmups": args.warmups,
        "scenario": args.scenario,
        "kind": args.kind,
        "sparse": args.sparse,
        "results": [result_record(result) for result in results],
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n")


def prepare_json_destination(output_path: Path) -> None:
    """Verify that a JSON destination is writable before timed work begins."""
    with output_path.open("a", encoding="utf-8"):
        pass


def print_run_context(args: argparse.Namespace) -> None:
    """Print the run identity needed to distinguish benchmark revisions."""
    print(
        f"label={args.label}",
        f"source_format={args.source_format}",
        f"path={args.path.resolve()}",
    )


def main() -> None:
    import genoio

    args = parse_args()
    try:
        validate_workload(args.block_sizes, args.max_variants)
        if args.output_json is not None:
            prepare_json_destination(args.output_json)

        print_run_context(args)
        dataset = dataset_for_source(args.source_format, args.path)
        results = []
        for scenario in selected_scenarios(args.scenario):
            read_options = read_options_for_case(
                args.source_format,
                scenario,
                args.kind,
                args.sparse,
            )
            for block_size in args.block_sizes:
                result = benchmark_case(
                    dataset,
                    scenario=scenario,
                    block_size=block_size,
                    max_variants=args.max_variants,
                    read_options=read_options,
                    repeats=args.repeats,
                    warmups=args.warmups,
                )
                print_result(result)
                results.append(result)

        if args.output_json is not None:
            write_json_report(args.output_json, args, results)
    except genoio.InternalError:
        raise
    except (genoio.GenoioError, OSError, ValueError) as error:
        raise SystemExit(f"error: {error}") from None


if __name__ == "__main__":
    main()
