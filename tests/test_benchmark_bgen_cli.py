# pattern: Mixed

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import benchmark_bgen  # noqa: E402


def _matrix(value: float = 1.0) -> np.ndarray:
    return np.array([[value, value + 1.0]], dtype=np.float32)


def test_parse_args_accepts_scenario_and_preserves_existing_options(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_bgen.py",
            "--scenario",
            "all",
            "--backend",
            "bgen_reader",
            "--max-variants",
            "7",
            "--repeats",
            "2",
            "--region",
            "22:1-100",
            "--no-compare",
        ],
    )

    args = benchmark_bgen.parse_args()

    assert args.scenario == "all"
    assert args.backend == "bgen_reader"
    assert args.max_variants == 7
    assert args.repeats == 2
    assert args.region == "22:1-100"
    assert args.no_compare is True


def test_read_bgen_sample_ids_uses_id_2_and_skips_type_row(tmp_path) -> None:
    path = tmp_path / "cohort.sample"
    path.write_text(
        """\
ID_1 ID_2 missing sex
0 0 0 D
0 S1 0 1
0 S2 0 2
"""
    )

    assert benchmark_bgen._read_bgen_sample_ids(path) == ["S1", "S2"]


def test_all_scenario_names_each_genoio_reader(monkeypatch, capsys) -> None:
    def read_with_variants(prefix, max_variants) -> np.ndarray:
        benchmark_bgen._last_variant_metadata_length = 3
        return _matrix(2.0)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_bgen.py",
            "--scenario",
            "all",
            "--backend",
            "genoio",
            "--max-variants",
            "3",
            "--repeats",
            "1",
        ],
    )
    monkeypatch.setattr(benchmark_bgen, "read_genoio_matrix_only", lambda prefix, max_variants: _matrix(1.0))
    monkeypatch.setattr(benchmark_bgen, "read_genoio_with_variants", read_with_variants)
    monkeypatch.setattr(benchmark_bgen, "read_genoio_sample_filtered", lambda prefix, max_variants: _matrix(3.0))
    monkeypatch.setattr(benchmark_bgen, "read_genoio_genotype_filtered", lambda prefix, max_variants: _matrix(4.0))
    monkeypatch.setattr(benchmark_bgen, "read_genoio_indexed_region", lambda prefix, max_variants, region: _matrix(5.0))

    benchmark_bgen.main()

    output = capsys.readouterr().out
    assert "genoio_bgen_matrix_only" in output
    assert "genoio_bgen_with_variants" in output
    assert "variant_metadata length=3" in output
    assert "genoio_bgen_sample_filtered" in output
    assert "genoio_bgen_genotype_filtered" in output
    assert "genoio_bgen_indexed_region" in output


def test_matrix_only_scenario_compares_bgen_reader(monkeypatch, capsys) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_bgen.py",
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
    monkeypatch.setattr(benchmark_bgen, "read_genoio_matrix_only", lambda prefix, max_variants: _matrix(1.0))
    monkeypatch.setattr(benchmark_bgen, "read_bgen_reader_expected_dosage", lambda prefix, max_variants: _matrix(1.0))

    benchmark_bgen.main()

    output = capsys.readouterr().out
    assert "genoio_bgen_matrix_only" in output
    assert "bgen_reader_expected_dosage" in output
    assert "comparison" in output


def test_expected_dosage_weights_support_unphased_and_phased_probability_shapes() -> None:
    np.testing.assert_array_equal(benchmark_bgen._expected_dosage_weights(3), np.array([0.0, 1.0, 2.0]))
    np.testing.assert_array_equal(benchmark_bgen._expected_dosage_weights(4), np.array([0.0, 1.0, 0.0, 1.0]))
