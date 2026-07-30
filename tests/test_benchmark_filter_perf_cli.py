# pattern: Mixed

from __future__ import annotations

import sys
from collections.abc import Iterable
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import polars as pl
from script_loader import load_benchmark_script

benchmark_filter_perf = cast(Any, load_benchmark_script("benchmark_filter_perf"))


class _ContextBlocks:
    def __init__(self, blocks: Iterable[Any]) -> None:
        self._blocks = iter(blocks)
        self.closed = False

    def __iter__(self) -> _ContextBlocks:
        return self

    def __next__(self) -> Any:
        return next(self._blocks)

    def __enter__(self) -> _ContextBlocks:
        return self

    def __exit__(self, *args: object) -> None:
        self.closed = True


def test_parse_args_accepts_polymorphic_predicate(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_filter_perf.py",
            "--source-format",
            "vcf",
            "--path",
            "cohort.vcf",
            "--predicate",
            "polymorphic",
            "--scenario",
            "numpy",
        ],
    )

    args = benchmark_filter_perf.parse_args()

    assert args.predicate == "polymorphic"
    assert args.source_format == "vcf"
    assert args.scenario == "numpy"


def test_parse_args_accepts_mac_and_non_multiallelic_alias(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_filter_perf.py",
            "--source-format",
            "vcf",
            "--path",
            "cohort.vcf",
            "--predicate",
            "mac",
            "--mac-min",
            "2",
            "--mac-max",
            "5",
            "--base-filter",
            "not_multiallelic",
        ],
    )

    args = benchmark_filter_perf.parse_args()

    assert args.predicate == "mac"
    assert args.mac_min == 2
    assert args.mac_max == 5
    assert args.base_filter == "not_multiallelic"


def test_parse_args_accepts_missing_rate_predicate(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_filter_perf.py",
            "--source-format",
            "bfile",
            "--path",
            "cohort",
            "--predicate",
            "missing_rate",
            "--missing-rate-max",
            "0.25",
        ],
    )

    args = benchmark_filter_perf.parse_args()

    assert args.predicate == "missing_rate"
    assert args.missing_rate_max == 0.25


def test_parse_args_accepts_composite_filter_shape(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_filter_perf.py",
            "--source-format",
            "pfile",
            "--path",
            "cohort",
            "--filter-shape",
            "mac_missing_rate",
            "--mac-min",
            "2",
            "--missing-rate-max",
            "0.05",
        ],
    )

    args = benchmark_filter_perf.parse_args()

    assert args.filter_shape == "mac_missing_rate"
    assert benchmark_filter_perf.active_filter_shape(args) == "mac_missing_rate"


def test_parse_args_accepts_source_window_mode(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_filter_perf.py",
            "--source-format",
            "pfile",
            "--path",
            "cohort",
            "--filter-shape",
            "mac_max",
            "--mac-max",
            "2",
            "--window-mode",
            "source",
        ],
    )

    args = benchmark_filter_perf.parse_args()

    assert args.window_mode == "source"


def test_numpy_retained_window_reads_until_enough_passing_variants(monkeypatch) -> None:
    chunks = [
        np.array([[0.0, 1.0, 1.0], [0.0, 1.0, 1.0]], dtype=np.float32),
        np.array([[2.0, 1.0], [2.0, 1.0]], dtype=np.float32),
    ]
    blocks = _ContextBlocks(chunks)

    class Dataset:
        def iter_blocks(self, size: int, **options: object) -> _ContextBlocks:
            assert size == 2
            return blocks

    args = SimpleNamespace(
        max_variants=2,
        filter_shape="mac_max",
        predicate="mac",
        maf_min=0.01,
        maf_max=None,
        mac_min=1,
        mac_max=1,
        missing_rate_max=None,
    )
    monkeypatch.setattr(benchmark_filter_perf, "dataset_for_args", lambda args: Dataset())
    monkeypatch.setattr(benchmark_filter_perf, "base_filter", lambda args: None)
    monkeypatch.setattr(benchmark_filter_perf, "read_options", lambda args: {})

    matrix = benchmark_filter_perf.read_numpy_retained_postfiltered(args)

    np.testing.assert_array_equal(matrix, np.array([[0.0, 2.0], [0.0, 2.0]], dtype=np.float32))
    assert blocks.closed is True


def test_source_window_variant_ids_come_from_first_base_block(monkeypatch) -> None:
    variants = pl.DataFrame({"id": ["rs1", "rs2", "rs3"]})
    blocks = _ContextBlocks([(np.zeros((2, 3), dtype=np.float32), variants)])

    class Dataset:
        def iter_blocks(self, size: int, **options: object) -> _ContextBlocks:
            assert size == 3
            assert options["return_variants"] is True
            return blocks

    args = SimpleNamespace(max_variants=3)
    monkeypatch.setattr(benchmark_filter_perf, "dataset_for_args", lambda args: Dataset())
    monkeypatch.setattr(benchmark_filter_perf, "base_filter", lambda args: None)
    monkeypatch.setattr(benchmark_filter_perf, "read_options", lambda args: {})

    assert benchmark_filter_perf.source_window_variant_ids(args) == ["rs1", "rs2", "rs3"]
    assert blocks.closed is True


def test_source_window_filter_rejects_duplicate_variant_ids() -> None:
    try:
        benchmark_filter_perf.validate_source_window_variant_ids(["rs1", "rs1"])
    except ValueError as error:
        assert "duplicate variant IDs" in str(error)
    else:
        raise AssertionError("expected duplicate IDs to be rejected")


def test_numpy_variant_mask_supports_polymorphic_predicate() -> None:
    matrix = np.array(
        [
            [0.0, 0.0, 2.0, np.nan],
            [0.0, 1.0, 2.0, np.nan],
            [0.0, 0.0, 2.0, np.nan],
        ],
        dtype=np.float32,
    )

    mask = benchmark_filter_perf.numpy_variant_mask(
        matrix,
        predicate="polymorphic",
        maf_min=0.01,
        maf_max=None,
        mac_min=1,
        mac_max=None,
        missing_rate_max=None,
    )

    np.testing.assert_array_equal(mask, np.array([False, True, False, False]))


def test_numpy_variant_mask_supports_mac_max_shape() -> None:
    matrix = np.array(
        [
            [0.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float32,
    )

    mask = benchmark_filter_perf.numpy_variant_mask(
        matrix,
        filter_shape="mac_max",
        maf_min=0.01,
        maf_max=None,
        mac_min=1,
        mac_max=2,
        missing_rate_max=None,
    )

    np.testing.assert_array_equal(mask, np.array([True, True, False]))


def test_numpy_variant_mask_supports_mac_predicate() -> None:
    matrix = np.array(
        [
            [0.0, 0.0, 2.0, np.nan],
            [0.0, 1.0, 2.0, 2.0],
            [1.0, 0.0, 2.0, 2.0],
        ],
        dtype=np.float32,
    )

    mask = benchmark_filter_perf.numpy_variant_mask(
        matrix,
        filter_shape="mac",
        maf_min=0.01,
        maf_max=None,
        mac_min=1,
        mac_max=2,
        missing_rate_max=None,
    )

    np.testing.assert_array_equal(mask, np.array([True, True, False, False]))


def test_numpy_variant_mask_supports_composite_shapes() -> None:
    matrix = np.array(
        [
            [0.0, 0.0, 1.0, np.nan],
            [1.0, 0.0, 1.0, np.nan],
            [1.0, 1.0, 1.0, 2.0],
            [2.0, 1.0, 1.0, 2.0],
        ],
        dtype=np.float32,
    )

    mac_missing = benchmark_filter_perf.numpy_variant_mask(
        matrix,
        filter_shape="mac_missing_rate",
        maf_min=0.2,
        maf_max=0.4,
        mac_min=2,
        mac_max=None,
        missing_rate_max=0.25,
    )
    maf_missing_poly = benchmark_filter_perf.numpy_variant_mask(
        matrix,
        filter_shape="maf_missing_rate_polymorphic",
        maf_min=0.2,
        maf_max=0.4,
        mac_min=1,
        mac_max=None,
        missing_rate_max=0.25,
    )

    np.testing.assert_array_equal(mac_missing, np.array([True, True, True, False]))
    np.testing.assert_array_equal(maf_missing_poly, np.array([False, True, False, False]))


def test_numpy_variant_mask_supports_missing_rate_predicate() -> None:
    matrix = np.array(
        [
            [0.0, np.nan, 2.0],
            [0.0, 1.0, np.nan],
            [0.0, 0.0, np.nan],
            [0.0, 0.0, 2.0],
        ],
        dtype=np.float32,
    )

    mask = benchmark_filter_perf.numpy_variant_mask(
        matrix,
        filter_shape="missing_rate",
        maf_min=0.01,
        maf_max=None,
        mac_min=1,
        mac_max=None,
        missing_rate_max=0.25,
    )

    np.testing.assert_array_equal(mask, np.array([True, True, False]))
