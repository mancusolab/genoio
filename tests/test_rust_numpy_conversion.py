# pattern: Imperative Shell

from pathlib import Path

import numpy as np

from genoio import _rust


def _write_tiny_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "tiny.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
"""
    )
    return path


def _write_missing_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "missing.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t./.
"""
    )
    return path


def _write_tiny_phased_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "tiny_phased.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|1
1\t20\trs2\tA\tG\t.\tPASS\t.\tGT\t1|0\t0|0
"""
    )
    return path


def _read_private_dense_values(path: Path, *, missing: str = "raise") -> dict[str, object]:
    return _rust.read_dense(
        "vcf",
        {"vcf": str(path)},
        {
            "samples": None,
            "variants": None,
            "variant_window": None,
            "dosage": "hardcall",
            "missing": missing,
            "return_samples": False,
            "return_variants": False,
            "matrix_only": True,
        },
    )


def _read_private_haplotype_dense_values(path: Path) -> dict[str, object]:
    return _rust.read_haplotypes_dense(
        "vcf",
        {"vcf": str(path)},
        {
            "samples": None,
            "variants": None,
            "variant_window": None,
            "dosage": "hardcall",
            "missing": "raise",
            "return_samples": False,
            "return_variants": False,
            "matrix_only": True,
        },
    )


def _read_private_sparse_values(path: Path) -> dict[str, object]:
    return _rust.read_sparse(
        "vcf",
        {"vcf": str(path)},
        {
            "samples": None,
            "variants": None,
            "variant_window": None,
            "dosage": "hardcall",
            "missing": "raise",
            "return_samples": False,
            "return_variants": False,
            "matrix_only": True,
        },
    )


def _is_bytearray_backed(array: np.ndarray) -> bool:
    return isinstance(getattr(array.base, "obj", None), bytearray)


def test_rust_f32_values_transfer_ownership_to_numpy_without_bytearray_base(tmp_path):
    result = _read_private_dense_values(_write_tiny_vcf(tmp_path))

    values = result["values"]

    assert isinstance(values, np.ndarray)
    assert values.dtype == np.dtype("float32")
    assert values.flags.writeable
    assert not _is_bytearray_backed(values)


def test_rust_dense_raise_policy_omits_missing_payload_when_no_calls_are_missing(tmp_path):
    result = _read_private_dense_values(_write_tiny_vcf(tmp_path))

    assert "missing_mask" not in result
    assert "missing_indices" not in result


def test_rust_dense_nan_policy_writes_nan_without_missing_payload(tmp_path):
    result = _read_private_dense_values(_write_missing_vcf(tmp_path), missing="nan")

    values = result["values"]

    assert "missing_mask" not in result
    assert "missing_indices" not in result
    np.testing.assert_array_equal(values, np.array([0.0, np.nan], dtype=np.float32))


def test_rust_dense_impute_policy_returns_missing_indices(tmp_path):
    result = _read_private_dense_values(_write_missing_vcf(tmp_path), missing="impute")

    indices = result["missing_indices"]

    assert "missing_mask" not in result
    assert isinstance(indices, np.ndarray)
    assert indices.dtype == np.dtype("int64")
    assert not _is_bytearray_backed(indices)
    np.testing.assert_array_equal(indices, np.array([1], dtype=np.int64))


def test_rust_dense_impute_policy_omits_missing_indices_when_no_calls_are_missing(tmp_path):
    result = _read_private_dense_values(_write_tiny_vcf(tmp_path), missing="impute")

    assert "missing_mask" not in result
    assert "missing_indices" not in result


def test_rust_sparse_indices_transfer_ownership_to_numpy_without_bytearray_base(tmp_path):
    result = _read_private_sparse_values(_write_tiny_vcf(tmp_path))

    indptr = result["indptr"]
    indices = result["indices"]

    assert isinstance(indptr, np.ndarray)
    assert isinstance(indices, np.ndarray)
    assert indptr.dtype == np.dtype("int64")
    assert indices.dtype == np.dtype("int64")
    assert not _is_bytearray_backed(indptr)
    assert not _is_bytearray_backed(indices)


def test_rust_variant_major_haplotype_values_report_layout_without_transpose(tmp_path):
    result = _read_private_haplotype_dense_values(_write_tiny_phased_vcf(tmp_path))

    values = result["values"]

    assert result["values_layout"] == "variant_major"
    np.testing.assert_array_equal(
        values,
        np.array([0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=np.float32),
    )
