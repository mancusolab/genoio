# pattern: Mixed

from __future__ import annotations

import contextlib
import io
import sys
from typing import Any, cast

import numpy as np
from scipy import sparse as scipy_sparse
from script_loader import load_benchmark_script

bench_common = cast(Any, load_benchmark_script("bench_common"))
benchmark_vcf = cast(Any, load_benchmark_script("benchmark_vcf"))


def _matrix(value: float = 1.0) -> np.ndarray:
    return np.array([[value, value + 1.0]], dtype=np.float32)


def _capture_stdout(func) -> str:
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        func()
    return output.getvalue()


def test_parse_args_accepts_scenario_kind_sparse_and_samples(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
            "--scenario",
            "all",
            "--backend",
            "genoio",
            "--max-variants",
            "7",
            "--repeats",
            "2",
            "--region",
            "22:1-100",
            "--kind",
            "haplo",
            "--sparse",
            "--samples",
            "S2,S1",
            "--no-compare",
        ],
    )

    args = benchmark_vcf.parse_args()

    assert args.scenario == "all"
    assert args.backend == "genoio"
    assert args.max_variants == 7
    assert args.repeats == 2
    assert args.region == "22:1-100"
    assert args.kind == "haplo"
    assert args.sparse is True
    assert args.samples == "S2,S1"
    assert args.no_compare is True


def test_read_vcf_sample_ids_supports_plain_vcf(tmp_path) -> None:
    path = tmp_path / "cohort.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
22\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1
"""
    )

    assert benchmark_vcf._read_vcf_sample_ids(path) == ["S1", "S2"]


def test_haplotype_kind_dispatches_genoio_without_cyvcf2_comparison(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
            "--scenario",
            "matrix-only",
            "--kind",
            "haplo",
            "--backend",
            "both",
            "--max-variants",
            "3",
            "--repeats",
            "1",
        ],
    )
    monkeypatch.setattr(benchmark_vcf, "read_genoio_matrix_only", lambda args: _matrix(1.0))
    monkeypatch.setattr(benchmark_vcf, "read_cyvcf2", lambda args: _matrix(1.0))

    output = _capture_stdout(benchmark_vcf.main)
    assert "genoio_vcf_haplo_matrix_only" in output
    assert "skipped cyvcf2 comparison for haplo matrix-only" in output
    assert "cyvcf2_vcf" not in output


def test_dosage_kind_dispatches_genoio_without_cyvcf2_comparison(monkeypatch) -> None:
    observed_options = []

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
            "--scenario",
            "matrix-only",
            "--kind",
            "dosage",
            "--backend",
            "both",
            "--max-variants",
            "3",
            "--repeats",
            "1",
        ],
    )

    def read_matrix(args) -> np.ndarray:
        observed_options.append(benchmark_vcf._read_options(args.kind, args.sparse))
        return _matrix(1.0)

    monkeypatch.setattr(benchmark_vcf, "read_genoio_matrix_only", read_matrix)
    monkeypatch.setattr(benchmark_vcf, "read_cyvcf2", lambda args: _matrix(1.0))

    output = _capture_stdout(benchmark_vcf.main)
    assert observed_options == [
        {"missing": "nan", "dtype": np.float32, "sparse": False, "dosage": "dosage"},
        {"missing": "nan", "dtype": np.float32, "sparse": False, "dosage": "dosage"},
    ]
    assert "genoio_vcf_dosage_matrix_only" in output
    assert "skipped cyvcf2 comparison for dosage matrix-only" in output
    assert "cyvcf2_vcf" not in output


def test_all_scenario_names_each_genoio_reader(monkeypatch) -> None:
    def read_with_variants(args) -> np.ndarray:
        benchmark_vcf._last_variant_metadata_length = 3
        return _matrix(2.0)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
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
    monkeypatch.setattr(benchmark_vcf, "read_genoio_matrix_only", lambda args: _matrix(1.0))
    monkeypatch.setattr(
        benchmark_vcf,
        "read_genoio_metadata",
        lambda args: np.array([2, 3], dtype=np.int64),
    )
    monkeypatch.setattr(benchmark_vcf, "read_genoio_with_variants", read_with_variants)
    monkeypatch.setattr(benchmark_vcf, "read_genoio_sample_filtered", lambda args: _matrix(3.0))
    monkeypatch.setattr(benchmark_vcf, "read_genoio_genotype_filtered", lambda args: _matrix(4.0))
    monkeypatch.setattr(benchmark_vcf, "read_genoio_indexed_region", lambda args: _matrix(5.0))
    monkeypatch.setattr(benchmark_vcf, "read_genoio_indexed_region_sample_filtered", lambda args: _matrix(6.0))

    output = _capture_stdout(benchmark_vcf.main)
    assert "genoio_vcf_metadata" in output
    assert "genoio_vcf_matrix_only" in output
    assert "genoio_vcf_with_variants" in output
    assert "variant_metadata length=3" in output
    assert "genoio_vcf_sample_filtered" in output
    assert "genoio_vcf_genotype_filtered" in output
    assert "genoio_vcf_indexed_region" in output
    assert "genoio_vcf_indexed_region_sample_filtered" in output


def test_metadata_scenario_reports_cold_time(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
            "--scenario",
            "metadata",
            "--backend",
            "genoio",
            "--repeats",
            "1",
        ],
    )
    monkeypatch.setattr(
        benchmark_vcf,
        "read_genoio_metadata",
        lambda args: np.array([2, 3], dtype=np.int64),
    )

    output = _capture_stdout(benchmark_vcf.main)
    assert "genoio_vcf_metadata" in output
    assert "cold=" in output


def test_haplotype_sparse_indexed_region_scenario_names_reader(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
            "--scenario",
            "indexed-region",
            "--backend",
            "genoio",
            "--kind",
            "haplo",
            "--sparse",
            "--max-variants",
            "3",
            "--repeats",
            "1",
        ],
    )
    monkeypatch.setattr(benchmark_vcf, "read_genoio_indexed_region", lambda args: _matrix(1.0))

    output = _capture_stdout(benchmark_vcf.main)
    assert "genoio_vcf_haplo_sparse_indexed_region" in output


def test_indexed_region_sample_filtered_scenario_names_reader(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
            "--scenario",
            "indexed-region-sample-filtered",
            "--backend",
            "genoio",
            "--max-variants",
            "3",
            "--repeats",
            "1",
        ],
    )
    monkeypatch.setattr(benchmark_vcf, "read_genoio_indexed_region_sample_filtered", lambda args: _matrix(1.0))

    output = _capture_stdout(benchmark_vcf.main)
    assert "genoio_vcf_indexed_region_sample_filtered" in output


def test_matrix_only_scenario_compares_cyvcf2(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
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
    monkeypatch.setattr(benchmark_vcf, "read_genoio_matrix_only", lambda args: _matrix(1.0))
    monkeypatch.setattr(benchmark_vcf, "read_cyvcf2", lambda args: _matrix(1.0))

    output = _capture_stdout(benchmark_vcf.main)
    assert "genoio_vcf_matrix_only" in output
    assert "cyvcf2_vcf" in output
    assert "comparison" in output


def test_sparse_dispatches_genoio_without_cyvcf2_comparison(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "benchmark_vcf.py",
            "--scenario",
            "matrix-only",
            "--backend",
            "both",
            "--sparse",
            "--max-variants",
            "3",
            "--repeats",
            "1",
        ],
    )
    monkeypatch.setattr(benchmark_vcf, "read_genoio_matrix_only", lambda args: _matrix(1.0))
    monkeypatch.setattr(benchmark_vcf, "read_cyvcf2", lambda args: _matrix(1.0))

    output = _capture_stdout(benchmark_vcf.main)
    assert "genoio_vcf_sparse_matrix_only" in output
    assert "skipped cyvcf2 comparison for sparse matrix-only" in output
    assert "cyvcf2_vcf" not in output


def test_matrix_summary_uses_sparse_storage_without_densifying(monkeypatch) -> None:
    matrix = scipy_sparse.csc_matrix(
        (
            np.array([1.0, 2.0, 3.0], dtype=np.float32),
            np.array([0, 2, 1], dtype=np.int64),
            np.array([0, 2, 3], dtype=np.int64),
        ),
        shape=(3, 2),
    )
    monkeypatch.setattr(
        matrix,
        "toarray",
        lambda: (_ for _ in ()).throw(AssertionError("sparse benchmark output should not be densified")),
    )

    summary = bench_common.matrix_summary(matrix)

    assert summary == {
        "shape": (3, 2),
        "dtype": "float32",
        "sum": 6.0,
        "missing": 0,
    }
