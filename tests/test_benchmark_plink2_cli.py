# pattern: Mixed

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import benchmark_plink2  # noqa: E402


def _matrix(value: float = 1.0) -> np.ndarray:
    return np.array([[value, value + 1.0]], dtype=np.float32)


def test_parse_args_accepts_scenario_and_preserves_existing_options(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_plink2.py",
            "--scenario",
            "all",
            "--backend",
            "genoio",
            "--max-variants",
            "7",
            "--repeats",
            "2",
            "--pgenlib-path",
            "plink-ng/2.0/Python",
            "--no-compare",
        ],
    )

    args = benchmark_plink2.parse_args()

    assert args.scenario == "all"
    assert args.backend == "genoio"
    assert args.max_variants == 7
    assert args.repeats == 2
    assert args.pgenlib_path == Path("plink-ng/2.0/Python")
    assert args.no_compare is True


def test_all_scenario_names_each_genoio_reader_and_skips_pgenlib(monkeypatch, capsys) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_plink2.py",
            "--scenario",
            "all",
            "--backend",
            "both",
            "--max-variants",
            "3",
            "--repeats",
            "1",
        ],
    )
    monkeypatch.setattr(benchmark_plink2, "read_genoio_matrix_only", lambda prefix, max_variants: _matrix(1.0))
    monkeypatch.setattr(benchmark_plink2, "read_genoio_with_variants", lambda prefix, max_variants: (_matrix(2.0), 3))
    monkeypatch.setattr(benchmark_plink2, "read_genoio_sample_filtered", lambda prefix, max_variants: _matrix(3.0))
    monkeypatch.setattr(benchmark_plink2, "read_genoio_genotype_filtered", lambda prefix, max_variants: _matrix(4.0))
    monkeypatch.setattr(benchmark_plink2, "read_pgenlib", lambda args: _matrix(1.0))

    benchmark_plink2.main()

    output = capsys.readouterr().out
    assert "genoio_plink2_matrix_only" in output
    assert "pgenlib_pgenreader" in output
    assert "comparison" in output
    assert "genoio_plink2_with_variants" in output
    assert "variant_metadata length=3" in output
    assert "genoio_plink2_sample_filtered" in output
    assert "genoio_plink2_genotype_filtered" in output
    assert output.count("skipped pgenlib comparison") == 3


def test_matrix_only_scenario_includes_pgenlib_comparison(monkeypatch, capsys) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_plink2.py",
            "--scenario",
            "matrix-only",
            "--backend",
            "both",
            "--max-variants",
            "3",
            "--repeats",
            "1",
        ],
    )
    monkeypatch.setattr(benchmark_plink2, "read_genoio_matrix_only", lambda prefix, max_variants: _matrix(1.0))
    monkeypatch.setattr(benchmark_plink2, "read_pgenlib", lambda args: _matrix(1.0))

    benchmark_plink2.main()

    output = capsys.readouterr().out
    assert "genoio_plink2_matrix_only" in output
    assert "pgenlib_pgenreader" in output
    assert "  time median=" in output
    assert "comparison" in output
