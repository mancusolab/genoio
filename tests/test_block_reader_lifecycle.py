# pattern: Imperative Shell

import gc
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
from fixture_writers import (
    write_bgen_dosage,
    write_canonical_plink1,
    write_fixed_width_phased_dosage_plink2,
    write_fixed_width_plink2,
    write_fixed_width_plink2_dosage,
    write_phased_dosage_plink2,
    write_phased_hardcall_plink2,
)

import genoio
from genoio import _rust


def _write_block_vcf(tmp_path: Path, *, missing: bool = False) -> Path:
    second_call = ".|." if missing else "1|1"
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


def _dense_native_block() -> dict[str, object]:
    return {
        "values": [0.0],
        "shape": (1, 1),
        "samples": {},
        "variants": {},
        "diagnostics": {},
    }


class _ControlledReader:
    def __init__(
        self,
        blocks: list[dict[str, object]] | None = None,
        *,
        next_error: BaseException | None = None,
        close_error: BaseException | None = None,
    ) -> None:
        self.blocks = [] if blocks is None else list(blocks)
        self.next_error = next_error
        self.close_error = close_error
        self.next_calls = 0
        self.close_calls = 0

    def next_block(self) -> dict[str, object] | None:
        self.next_calls += 1
        if self.next_error is not None:
            raise self.next_error
        if self.blocks:
            return self.blocks.pop(0)
        return None

    def close(self) -> None:
        self.close_calls += 1
        if self.close_error is not None:
            raise self.close_error


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


def test_pbr_py_error_001_public_header_error_is_deferred_to_first_next(tmp_path):
    path = tmp_path / "invalid_header.vcf"
    path.write_text("not a VCF header\n")
    dataset = genoio.vcf(path)

    iterator = dataset.iter_blocks(size=1)

    with pytest.raises(genoio.InvalidSourceError, match="header|VCF"):
        next(iterator)


@pytest.mark.parametrize(
    ("kind", "sparse"),
    [
        pytest.param("geno", False, id="dense-genotype"),
        pytest.param("geno", "csc", id="sparse-genotype"),
        pytest.param("haplo", False, id="dense-haplotype"),
        pytest.param("haplo", "csr", id="sparse-haplotype"),
    ],
)
def test_pbr_py_error_001_missing_error_matches_full_read_on_affected_block(
    tmp_path,
    kind,
    sparse,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path, missing=True))
    read_options = {
        "kind": kind,
        "sparse": sparse,
        "missing": "raise",
    }

    with pytest.raises(genoio.MissingDataError, match="missing"):
        cast(Any, dataset.read)(**read_options)
    iterator = dataset.iter_blocks(size=1, **read_options)
    assert next(iterator).shape[1] == 1

    with pytest.raises(genoio.MissingDataError, match="missing"):
        next(iterator)


@pytest.mark.parametrize(
    ("dataset_factory", "read_options", "match"),
    [
        pytest.param(
            lambda tmp_path: genoio.vcf(_write_block_vcf(tmp_path)),
            {"dosage": "dosage", "sparse": True},
            "sparse dosage",
            id="vcf-sparse-dosage",
        ),
        pytest.param(
            lambda tmp_path: genoio.bgen(write_bgen_dosage(tmp_path)),
            {"dosage": "dosage", "sparse": True},
            "sparse genotype",
            id="bgen-sparse-dosage",
        ),
        pytest.param(
            lambda tmp_path: genoio.bfile(write_canonical_plink1(tmp_path)),
            {"kind": "haplo"},
            "haplo",
            id="plink1-haplotype",
        ),
        pytest.param(
            lambda tmp_path: genoio.pfile(write_phased_dosage_plink2(tmp_path)),
            {"kind": "haplo", "dosage": "dosage", "sparse": True},
            "sparse haplotype",
            id="plink2-sparse-haplotype-dosage",
        ),
    ],
)
def test_pbr_py_error_001_unsupported_error_matches_full_read(
    tmp_path,
    dataset_factory,
    read_options,
    match,
):
    dataset = dataset_factory(tmp_path)

    with pytest.raises(genoio.UnsupportedRepresentation, match=match):
        cast(Any, dataset.read)(**read_options)
    with pytest.raises(genoio.UnsupportedRepresentation, match=match):
        dataset.iter_blocks(size=1, **read_options)


def test_pbr_py_error_001_constructor_error_maps_to_public_exception(
    tmp_path,
    monkeypatch,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path))

    def fail_constructor(*args, **kwargs):
        raise _rust.RustInvalidSourceError("constructor failed")

    monkeypatch.setattr(_rust, "_BlockReader", fail_constructor)
    iterator = dataset.iter_blocks(size=1)

    with pytest.raises(genoio.InvalidSourceError, match="constructor failed"):
        next(iterator)


def test_pbr_py_error_001_advance_error_maps_and_closes_reader(
    tmp_path,
    monkeypatch,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path))
    reader = _ControlledReader(
        next_error=_rust.RustMissingDataError("advance failed"),
    )
    monkeypatch.setattr(_rust, "_BlockReader", lambda *args, **kwargs: reader)
    iterator = dataset.iter_blocks(size=1)

    with pytest.raises(genoio.MissingDataError, match="advance failed"):
        next(iterator)

    assert reader.close_calls == 1


def test_pbr_py_error_001_already_public_error_passes_through_unchanged(
    tmp_path,
    monkeypatch,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path))
    public_error = genoio.InvalidSourceError("public failure")
    reader = _ControlledReader(next_error=public_error)
    monkeypatch.setattr(_rust, "_BlockReader", lambda *args, **kwargs: reader)

    with pytest.raises(genoio.InvalidSourceError, match="public failure") as raised:
        next(dataset.iter_blocks(size=1))

    assert raised.value is public_error
    assert reader.close_calls == 1


def test_pbr_py_error_001_explicit_close_failure_maps_and_propagates(
    tmp_path,
    monkeypatch,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path))
    reader = _ControlledReader(
        [_dense_native_block()],
        close_error=_rust.RustInvalidSourceError("close failed"),
    )
    monkeypatch.setattr(_rust, "_BlockReader", lambda *args, **kwargs: reader)
    iterator = dataset.iter_blocks(size=1)
    next(iterator)

    with pytest.raises(genoio.InvalidSourceError, match="close failed"):
        cast(Any, iterator).close()

    assert reader.close_calls == 1


def test_pbr_py_error_001_primary_error_is_preserved_with_close_failure_note(
    tmp_path,
    monkeypatch,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path))
    reader = _ControlledReader(
        next_error=_rust.RustMissingDataError("primary failed"),
        close_error=_rust.RustInvalidSourceError("close also failed"),
    )
    monkeypatch.setattr(_rust, "_BlockReader", lambda *args, **kwargs: reader)

    with pytest.raises(genoio.MissingDataError, match="primary failed") as raised:
        next(dataset.iter_blocks(size=1))

    assert reader.close_calls == 1
    assert any("close also failed" in note for note in raised.value.__notes__)


def test_pbr_py_error_001_normal_exhaustion_exposes_close_failure(
    tmp_path,
    monkeypatch,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path))
    reader = _ControlledReader(
        close_error=_rust.RustInvalidSourceError("exhaustion close failed"),
    )
    monkeypatch.setattr(_rust, "_BlockReader", lambda *args, **kwargs: reader)

    with pytest.raises(genoio.InvalidSourceError, match="exhaustion close failed"):
        next(dataset.iter_blocks(size=1))

    assert reader.close_calls == 1


@pytest.mark.parametrize(
    "exit_path",
    [
        pytest.param("exhaustion", id="exhaustion"),
        pytest.param("explicit-close", id="explicit-close"),
        pytest.param("read-error", id="read-error"),
        pytest.param("conversion-error", id="conversion-error"),
        pytest.param("finalization", id="finalization"),
    ],
)
def test_pbr_py_lifecycle_001_public_iterator_closes_reader_once(
    tmp_path,
    monkeypatch,
    exit_path,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path))
    if exit_path == "exhaustion":
        reader = _ControlledReader()
    elif exit_path == "read-error":
        reader = _ControlledReader(
            next_error=_rust.RustMissingDataError("read failed"),
        )
    elif exit_path == "conversion-error":
        reader = _ControlledReader([{"shape": (1, 1)}])
    else:
        reader = _ControlledReader([_dense_native_block()])
    monkeypatch.setattr(_rust, "_BlockReader", lambda *args, **kwargs: reader)
    iterator = dataset.iter_blocks(size=1)

    if exit_path == "exhaustion":
        assert list(iterator) == []
    elif exit_path == "explicit-close":
        next(iterator)
        cast(Any, iterator).close()
    elif exit_path == "read-error":
        with pytest.raises(genoio.MissingDataError, match="read failed"):
            next(iterator)
    elif exit_path == "conversion-error":
        with pytest.raises(KeyError, match="values"):
            next(iterator)
    else:
        next(iterator)
        del iterator
        gc.collect()

    assert reader.close_calls == 1


def test_pbr_py_iterator_001_pbr_py_lifecycle_001_failed_iterator_does_not_affect_peer(
    tmp_path,
):
    dataset = genoio.vcf(_write_block_vcf(tmp_path, missing=True))
    strict = dataset.iter_blocks(size=1, missing="raise")
    permissive = dataset.iter_blocks(size=1, missing="nan")

    strict_first = next(strict)
    permissive_first = next(permissive)
    with pytest.raises(genoio.MissingDataError, match="missing"):
        next(strict)
    permissive_second = next(permissive)

    np.testing.assert_array_equal(strict_first, permissive_first)
    np.testing.assert_array_equal(
        permissive_second,
        np.array([[1.0], [np.nan]], dtype=np.float32),
    )
    with pytest.raises(StopIteration):
        next(permissive)
