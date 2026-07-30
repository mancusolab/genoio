# pattern: Mixed (unavoidable)
# Reason: CLI shell tests share real lightweight integration fixtures with pure helper tests.

from __future__ import annotations

import argparse
import io
import json
import platform
import sys
from collections.abc import Iterable
from contextlib import redirect_stdout
from pathlib import Path
from typing import Any, cast

import numpy as np
import polars as pl
import pytest
from script_loader import load_benchmark_script

import genoio


def _benchmark_iter_blocks() -> Any:
    return cast(Any, load_benchmark_script("benchmark_iter_blocks"))


class _ContextBlocks:
    def __init__(self, blocks: Iterable[Any]) -> None:
        self._blocks = iter(blocks)
        self.entered = False
        self.closed = False

    def __iter__(self) -> _ContextBlocks:
        return self

    def __next__(self) -> Any:
        return next(self._blocks)

    def __enter__(self) -> _ContextBlocks:
        self.entered = True
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def close(self) -> None:
        self.closed = True


def test_read_first_block_uses_iterator_context_manager() -> None:
    bench_common = cast(Any, load_benchmark_script("bench_common"))
    expected = np.zeros((3, 2), dtype=np.float32)
    blocks = _ContextBlocks([expected])

    class Dataset:
        def iter_blocks(self, size: int, **options: object) -> _ContextBlocks:
            assert size == 2
            assert options == {"missing": "nan"}
            return blocks

    observed = bench_common.read_first_block(Dataset(), 2, missing="nan")

    assert observed is expected
    assert blocks.entered is True
    assert blocks.closed is True


def test_parse_args_accepts_streaming_workload_options(monkeypatch) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_iter_blocks.py",
            "--source-format",
            "vcf",
            "--path",
            "cohort.vcf.gz",
            "--label",
            "candidate",
            "--block-sizes",
            "128,512,2048",
            "--max-variants",
            "4096",
            "--repeats",
            "7",
            "--warmups",
            "2",
            "--scenario",
            "all",
            "--output-json",
            "results.json",
        ],
    )

    args = benchmark_iter_blocks.parse_args()

    assert args.source_format == "vcf"
    assert args.path == Path("cohort.vcf.gz")
    assert args.label == "candidate"
    assert args.block_sizes == (128, 512, 2048)
    assert args.max_variants == 4096
    assert args.repeats == 7
    assert args.warmups == 2
    assert args.scenario == "all"
    assert args.output_json == Path("results.json")


def test_json_report_records_reproducible_run_context(tmp_path: Path) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    source_path = tmp_path / "cohort.vcf.gz"
    output_path = tmp_path / "candidate.json"
    args = argparse.Namespace(
        source_format="vcf",
        path=source_path,
        label="candidate-05bc6c1",
        max_variants=4,
        block_sizes=(2,),
        repeats=2,
        warmups=1,
        scenario="matrix-only",
        kind="hardcall",
        sparse=False,
    )
    result = benchmark_iter_blocks.BenchmarkResult(
        scenario="matrix-only",
        block_size=2,
        summary=benchmark_iter_blocks.ScanSummary(blocks=2, variants=4, rows=3),
        seconds=(0.5, 0.25),
    )

    benchmark_iter_blocks.write_json_report(output_path, args, [result])

    report = json.loads(output_path.read_text())
    assert report["label"] == "candidate-05bc6c1"
    assert report["genoio_version"] == genoio.__version__
    assert report["python_version"] == platform.python_version()
    assert report["platform"] == platform.platform()
    assert report["machine"] == platform.machine()
    assert report["path"] == str(source_path.resolve())


@pytest.mark.parametrize("source_format", ["vcf", "bfile", "pfile", "bgen"])
def test_dataset_for_source_dispatches_to_selected_constructor(
    monkeypatch,
    source_format: str,
) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    sentinel = object()
    calls: list[Path] = []

    def constructor(path: Path) -> object:
        calls.append(path)
        return sentinel

    monkeypatch.setattr(genoio, source_format, constructor)

    observed = benchmark_iter_blocks.dataset_for_source(source_format, Path("cohort"))

    assert observed is sentinel
    assert calls == [Path("cohort")]


@pytest.mark.parametrize(
    ("block_sizes", "max_variants", "message"),
    [
        ((128, 512), 1000, "divisible"),
        ((128, 128), 1024, "distinct"),
    ],
)
def test_validate_workload_rejects_biased_sweeps(
    block_sizes: tuple[int, ...],
    max_variants: int,
    message: str,
) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()

    with pytest.raises(ValueError, match=message):
        benchmark_iter_blocks.validate_workload(block_sizes, max_variants)


def test_sparse_case_uses_supported_missing_value_policy() -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()

    options = benchmark_iter_blocks.read_options_for_case(
        "bfile",
        "matrix-only",
        "hardcall",
        True,
    )

    assert options["sparse"] == "csc"
    assert options["missing"] == "raise"


def test_consume_blocks_counts_exact_prefix_and_closes_early() -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    blocks = _ContextBlocks(
        [
            np.zeros((3, 2), dtype=np.float32),
            np.ones((3, 2), dtype=np.float32),
            np.full((3, 2), 2.0, dtype=np.float32),
        ]
    )

    class Dataset:
        def iter_blocks(self, size: int, **options: object) -> _ContextBlocks:
            assert size == 2
            assert options == {"missing": "nan", "dtype": np.float32}
            return blocks

    summary = benchmark_iter_blocks.consume_blocks(
        Dataset(),
        block_size=2,
        max_variants=4,
        read_options={"missing": "nan", "dtype": np.float32},
    )

    assert summary == benchmark_iter_blocks.ScanSummary(blocks=2, variants=4, rows=3)
    assert blocks.entered is True
    assert blocks.closed is True


def test_consume_blocks_validates_variant_metadata_alignment() -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()

    class Dataset:
        def iter_blocks(self, size: int, **options: object) -> _ContextBlocks:
            assert size == 2
            assert options["return_variants"] is True
            return _ContextBlocks(
                [
                    (
                        np.zeros((3, 2), dtype=np.float32),
                        pl.DataFrame({"id": ["rs1"]}),
                    )
                ]
            )

    with pytest.raises(RuntimeError, match="metadata rows"):
        benchmark_iter_blocks.consume_blocks(
            Dataset(),
            block_size=2,
            max_variants=2,
            read_options={"return_variants": True},
        )


def test_consume_blocks_reports_short_source_as_workload_error() -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()

    class Dataset:
        def iter_blocks(self, size: int, **options: object) -> _ContextBlocks:
            return _ContextBlocks([np.zeros((3, 1), dtype=np.float32)])

    with pytest.raises(ValueError, match="source ended after 1 retained variants"):
        benchmark_iter_blocks.consume_blocks(
            Dataset(),
            block_size=1,
            max_variants=2,
            read_options={},
        )


def test_consume_blocks_reads_multiple_real_vcf_blocks() -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    fixture = Path(__file__).parent / "fixtures" / "vcf" / "tiny.vcf"
    dataset = benchmark_iter_blocks.dataset_for_source("vcf", fixture)
    read_options = benchmark_iter_blocks.read_options_for_case(
        "vcf",
        "with-variants",
        "hardcall",
        False,
    )

    summary = benchmark_iter_blocks.consume_blocks(
        dataset,
        block_size=1,
        max_variants=2,
        read_options=read_options,
    )

    assert summary == benchmark_iter_blocks.ScanSummary(blocks=2, variants=2, rows=3)


def _main_argv(
    tmp_path: Path,
    *,
    output_json: Path | None = None,
) -> list[str]:
    argv = [
        "benchmark_iter_blocks.py",
        "--source-format",
        "vcf",
        "--path",
        str(tmp_path / "cohort.vcf.gz"),
        "--label",
        "candidate",
        "--block-sizes",
        "2",
        "--max-variants",
        "2",
        "--repeats",
        "1",
        "--warmups",
        "0",
    ]
    if output_json is not None:
        argv.extend(["--output-json", str(output_json)])
    return argv


def _fixed_result(benchmark_iter_blocks: Any) -> Any:
    return benchmark_iter_blocks.BenchmarkResult(
        scenario="matrix-only",
        block_size=2,
        summary=benchmark_iter_blocks.ScanSummary(blocks=1, variants=2, rows=3),
        seconds=(0.25,),
    )


def test_main_prints_run_identity_and_result(monkeypatch, tmp_path: Path) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    monkeypatch.setattr(sys, "argv", _main_argv(tmp_path))
    monkeypatch.setattr(benchmark_iter_blocks, "dataset_for_source", lambda source_format, path: object())
    monkeypatch.setattr(
        benchmark_iter_blocks, "benchmark_case", lambda *args, **kwargs: _fixed_result(benchmark_iter_blocks)
    )

    stream = io.StringIO()
    with redirect_stdout(stream):
        benchmark_iter_blocks.main()

    output = stream.getvalue()
    assert "label=candidate" in output
    assert f"path={(tmp_path / 'cohort.vcf.gz').resolve()}" in output
    assert "source_format=vcf" in output
    assert "scenario=matrix-only block_size=2" in output


def test_main_writes_json_report(monkeypatch, tmp_path: Path) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    output_path = tmp_path / "candidate.json"
    monkeypatch.setattr(sys, "argv", _main_argv(tmp_path, output_json=output_path))
    monkeypatch.setattr(benchmark_iter_blocks, "dataset_for_source", lambda source_format, path: object())
    monkeypatch.setattr(
        benchmark_iter_blocks, "benchmark_case", lambda *args, **kwargs: _fixed_result(benchmark_iter_blocks)
    )

    benchmark_iter_blocks.main()

    report = json.loads(output_path.read_text())
    assert report["label"] == "candidate"
    assert report["results"] == [
        {
            "scenario": "matrix-only",
            "block_size": 2,
            "blocks": 1,
            "variants": 2,
            "rows": 3,
            "median_seconds": 0.25,
            "min_seconds": 0.25,
            "variants_per_second": 8.0,
            "seconds": [0.25],
        }
    ]


def test_main_validates_json_destination_before_benchmarking(
    monkeypatch,
    tmp_path: Path,
) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    calls: list[str] = []
    monkeypatch.setattr(sys, "argv", _main_argv(tmp_path, output_json=tmp_path))
    monkeypatch.setattr(
        benchmark_iter_blocks,
        "dataset_for_source",
        lambda source_format, path: calls.append("dataset"),
    )

    with pytest.raises(SystemExit, match=r"^error: .*") as raised:
        benchmark_iter_blocks.main()

    assert raised.value.code != 0
    assert calls == []


@pytest.mark.parametrize(
    "error",
    [
        genoio.InvalidSourceError("source unavailable"),
        genoio.UnsupportedRepresentation("representation unavailable"),
        OSError("input read failed"),
        ValueError("invalid workload"),
    ],
)
def test_main_reports_expected_benchmark_errors_without_traceback(
    monkeypatch,
    tmp_path: Path,
    error: Exception,
) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    monkeypatch.setattr(sys, "argv", _main_argv(tmp_path))
    monkeypatch.setattr(
        benchmark_iter_blocks,
        "dataset_for_source",
        lambda source_format, path: (_ for _ in ()).throw(error),
    )

    with pytest.raises(SystemExit, match=rf"^error: {error}$") as raised:
        benchmark_iter_blocks.main()

    assert raised.value.__cause__ is None


def test_main_does_not_hide_internal_errors(monkeypatch, tmp_path: Path) -> None:
    benchmark_iter_blocks = _benchmark_iter_blocks()
    monkeypatch.setattr(sys, "argv", _main_argv(tmp_path))
    error = genoio.InternalError("backend invariant")
    monkeypatch.setattr(
        benchmark_iter_blocks,
        "dataset_for_source",
        lambda source_format, path: (_ for _ in ()).throw(error),
    )

    with pytest.raises(genoio.InternalError, match="backend invariant"):
        benchmark_iter_blocks.main()
