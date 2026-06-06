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
            "--kind",
            "haplo-dosage",
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
    assert args.kind == "haplo-dosage"
    assert args.pgenlib_path == Path("plink-ng/2.0/Python")
    assert args.no_compare is True


def test_haplotype_kinds_dispatch_genoio_and_skip_pgenlib(monkeypatch, capsys) -> None:
    calls: list[str] = []

    def read_matrix(prefix, max_variants, kind) -> np.ndarray:
        calls.append(kind)
        return _matrix(1.0)

    monkeypatch.setattr(benchmark_plink2, "read_genoio_matrix_only", read_matrix)
    monkeypatch.setattr(benchmark_plink2, "read_pgenlib", lambda args: _matrix(1.0))

    for kind in ("haplo-hardcall", "haplo-dosage"):
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "benchmark_plink2.py",
                "--scenario",
                "matrix-only",
                "--kind",
                kind,
                "--backend",
                "both",
                "--max-variants",
                "3",
                "--repeats",
                "1",
            ],
        )

        benchmark_plink2.main()

    output = capsys.readouterr().out
    assert calls == ["haplo-hardcall", "haplo-hardcall", "haplo-dosage", "haplo-dosage"]
    assert "genoio_plink2_haplo_hardcall_matrix_only" in output
    assert "genoio_plink2_haplo_dosage_matrix_only" in output
    assert output.count("skipped pgenlib comparison") == 2
    assert "pgenlib_pgenreader" not in output


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
    monkeypatch.setattr(benchmark_plink2, "read_genoio_matrix_only", lambda prefix, max_variants, kind: _matrix(1.0))
    monkeypatch.setattr(
        benchmark_plink2,
        "read_genoio_with_variants",
        lambda prefix, max_variants, kind: (_matrix(2.0), 3),
    )
    monkeypatch.setattr(
        benchmark_plink2, "read_genoio_sample_filtered", lambda prefix, max_variants, kind: _matrix(3.0)
    )
    monkeypatch.setattr(
        benchmark_plink2, "read_genoio_genotype_filtered", lambda prefix, max_variants, kind: _matrix(4.0)
    )
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
    monkeypatch.setattr(benchmark_plink2, "read_genoio_matrix_only", lambda prefix, max_variants, kind: _matrix(1.0))
    monkeypatch.setattr(benchmark_plink2, "read_pgenlib", lambda args: _matrix(1.0))

    benchmark_plink2.main()

    output = capsys.readouterr().out
    assert "genoio_plink2_matrix_only" in output
    assert "pgenlib_pgenreader" in output
    assert "  time median=" in output
    assert "comparison" in output
