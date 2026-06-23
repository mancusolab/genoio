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


def _read_private_dense_values(path: Path) -> dict[str, object]:
    return _rust.read_dense(
        "vcf",
        {"vcf": str(path)},
        {
            "samples": None,
            "variants": None,
            "variant_window": None,
            "dosage": "hardcall",
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


def test_rust_dense_missing_mask_transfers_ownership_to_numpy_without_bytearray_base(tmp_path):
    result = _read_private_dense_values(_write_tiny_vcf(tmp_path))

    missing_mask = result["missing_mask"]

    assert isinstance(missing_mask, np.ndarray)
    assert missing_mask.dtype == np.dtype("bool")
    assert not _is_bytearray_backed(missing_mask)


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
