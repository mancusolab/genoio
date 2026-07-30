# pattern: Imperative Shell

from pathlib import Path

import numpy as np
import pytest
from fixture_writers import (
    write_bgen_dosage,
    write_canonical_plink1,
    write_fixed_width_phased_dosage_plink2,
    write_fixed_width_plink2,
    write_fixed_width_plink2_dosage,
    write_phased_hardcall_plink2,
)

from genoio import _rust


def _write_block_vcf(tmp_path: Path, *, missing: bool = False) -> Path:
    second_call = "./." if missing else "1|1"
    path = tmp_path / ("missing.vcf" if missing else "blocks.vcf")
    path.write_text(
        f"""\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DS\t0|0:0.1\t0|1:0.9
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT:DS\t0|1:1.2\t{second_call}:1.8
"""
    )
    return path


def _options(
    *,
    samples: list[str] | None = None,
    dosage: str = "hardcall",
    missing: str = "raise",
) -> dict[str, object]:
    return {
        "samples": samples,
        "variants": None,
        "dosage": dosage,
        "missing": missing,
        "return_samples": False,
        "return_variants": False,
        "matrix_only": True,
    }


def _reader(
    path: Path,
    *,
    kind: str = "geno",
    sparse: bool = False,
    dosage: str = "hardcall",
    missing: str = "raise",
    block_size: int = 1,
):
    return _rust._BlockReader(
        "vcf",
        {"vcf": str(path)},
        kind,
        sparse,
        _options(dosage=dosage, missing=missing),
        block_size,
    )


@pytest.mark.parametrize(
    ("kind", "sparse", "expected_shape"),
    [
        pytest.param("geno", False, (2, 1), id="dense-genotype"),
        pytest.param("geno", True, (2, 1), id="sparse-genotype"),
        pytest.param("haplo", False, (4, 1), id="dense-haplotype"),
        pytest.param("haplo", True, (4, 1), id="sparse-haplotype"),
    ],
)
def test_pbr_py_native_001_private_reader_advances_supported_vcf_modes(
    tmp_path,
    kind,
    sparse,
    expected_shape,
):
    reader = _reader(_write_block_vcf(tmp_path), kind=kind, sparse=sparse)

    first = reader.next_block()
    second = reader.next_block()

    assert first["shape"] == expected_shape
    assert second["shape"] == expected_shape
    assert reader.next_block() is None
    assert reader.next_block() is None


def test_pbr_py_native_001_private_reader_supports_dense_dosage(tmp_path):
    reader = _reader(_write_block_vcf(tmp_path), dosage="dosage", block_size=2)

    block = reader.next_block()

    np.testing.assert_array_equal(
        block["values"],
        np.array([0.1, 0.9, 1.2, 1.8], dtype=np.float32),
    )
    assert reader.next_block() is None


@pytest.mark.parametrize(
    ("format", "kind", "sparse", "dosage", "samples", "expected_rows"),
    [
        pytest.param("bgen", "geno", False, "dosage", None, 2, id="bgen-genotype-dosage"),
        pytest.param("bgen", "haplo", False, "dosage", None, 4, id="bgen-haplotype-dosage"),
        pytest.param("plink1", "geno", False, "hardcall", ["S1", "S3", "S4"], 3, id="plink1-dense"),
        pytest.param("plink1", "geno", True, "hardcall", ["S1", "S3", "S4"], 3, id="plink1-sparse"),
        pytest.param("plink2", "geno", False, "hardcall", ["S1", "S3"], 2, id="plink2-dense"),
        pytest.param("plink2", "geno", True, "hardcall", ["S1", "S3"], 2, id="plink2-sparse"),
        pytest.param(
            "plink2",
            "haplo",
            False,
            "hardcall",
            ["S1", "S2"],
            4,
            id="plink2-dense-haplotype",
        ),
        pytest.param(
            "plink2",
            "haplo",
            True,
            "hardcall",
            ["S1", "S2"],
            4,
            id="plink2-sparse-haplotype",
        ),
        pytest.param("plink2", "geno", False, "dosage", None, 3, id="plink2-genotype-dosage"),
        pytest.param("plink2", "haplo", False, "dosage", None, 6, id="plink2-haplotype-dosage"),
    ],
)
def test_pbr_py_native_001_private_reader_advances_non_vcf_backend_modes(
    tmp_path,
    format,
    kind,
    sparse,
    dosage,
    samples,
    expected_rows,
):
    if format == "bgen":
        path = write_bgen_dosage(tmp_path, phased=kind == "haplo")
        members = {"bgen": str(path)}
    elif format == "plink1":
        prefix = write_canonical_plink1(tmp_path)
        members = {
            "bed": str(prefix.with_suffix(".bed")),
            "bim": str(prefix.with_suffix(".bim")),
            "fam": str(prefix.with_suffix(".fam")),
        }
    else:
        if kind == "haplo" and dosage == "hardcall":
            prefix = write_phased_hardcall_plink2(tmp_path)
        elif kind == "haplo":
            prefix = write_fixed_width_phased_dosage_plink2(tmp_path)
        elif dosage == "dosage":
            prefix = write_fixed_width_plink2_dosage(tmp_path)
        else:
            prefix = write_fixed_width_plink2(tmp_path)
        members = {
            "pgen": str(prefix.with_suffix(".pgen")),
            "pvar": str(prefix.with_suffix(".pvar")),
            "psam": str(prefix.with_suffix(".psam")),
        }

    reader = _rust._BlockReader(
        format,
        members,
        kind,
        sparse,
        _options(samples=samples, dosage=dosage),
        1,
    )

    block = reader.next_block()

    assert block is not None
    assert block["shape"] == (expected_rows, 1)


def test_pbr_py_lifecycle_001_close_is_idempotent_and_terminal(tmp_path):
    reader = _reader(_write_block_vcf(tmp_path))
    assert reader.next_block() is not None

    reader.close()
    reader.close()

    assert reader.next_block() is None


def test_pbr_py_error_001_unsupported_representation_keeps_native_error_type(tmp_path):
    with pytest.raises(
        _rust.RustUnsupportedRepresentationError,
        match="sparse dosage-backed genotype",
    ):
        _reader(_write_block_vcf(tmp_path), sparse=True, dosage="dosage")


def test_pbr_py_error_001_source_error_keeps_native_error_type(tmp_path):
    with pytest.raises(_rust.RustInvalidSourceError):
        _reader(tmp_path / "absent.vcf")


def test_pbr_py_error_001_missing_data_error_is_raised_by_affected_block(tmp_path):
    reader = _reader(_write_block_vcf(tmp_path, missing=True))
    assert reader.next_block() is not None

    with pytest.raises(_rust.RustMissingDataError, match="missing"):
        reader.next_block()

    assert reader.next_block() is None
