# pattern: Imperative Shell

import struct
from pathlib import Path

import polars as pl
import pytest

FIXTURE_ROOT = Path(__file__).parent / "fixtures"

_BGEN_FLAG_LAYOUT2 = 2 << 2
_BGEN_FLAG_SAMPLE_IDENTIFIERS = 1 << 31


def _write_bgen_header(
    buffer: bytearray,
    *,
    n_samples: int,
    n_variants: int,
    has_sample_ids: bool,
) -> None:
    flags = _BGEN_FLAG_LAYOUT2
    if has_sample_ids:
        flags |= _BGEN_FLAG_SAMPLE_IDENTIFIERS
    buffer.extend(struct.pack("<IIII4sI", 20, 20, n_variants, n_samples, b"bgen", flags))


def _write_sample_identifier_block(buffer: bytearray, sample_ids: list[str]) -> None:
    encoded_ids = [sample_id.encode() for sample_id in sample_ids]
    block_len = 8 + sum(2 + len(sample_id) for sample_id in encoded_ids)
    buffer.extend(struct.pack("<II", block_len, len(encoded_ids)))
    for sample_id in encoded_ids:
        buffer.extend(struct.pack("<H", len(sample_id)))
        buffer.extend(sample_id)


def _write_layout2_variant(
    buffer: bytearray,
    *,
    variant_id: str,
    rsid: str,
    chrom: str,
    pos: int,
    alleles: tuple[str, str],
    n_samples: int,
) -> None:
    for value in (variant_id, rsid, chrom):
        encoded = value.encode()
        buffer.extend(struct.pack("<H", len(encoded)))
        buffer.extend(encoded)
    buffer.extend(struct.pack("<IH", pos, len(alleles)))
    for allele in alleles:
        encoded = allele.encode()
        buffer.extend(struct.pack("<I", len(encoded)))
        buffer.extend(encoded)
    buffer.extend(struct.pack("<IIHBBBB", 10, n_samples, len(alleles), 2, 2, 0, 0))


def _write_tiny_bgen(path: Path, *, sample_ids: list[str] | None) -> None:
    buffer = bytearray()
    _write_bgen_header(
        buffer,
        n_samples=2,
        n_variants=2,
        has_sample_ids=sample_ids is not None,
    )
    if sample_ids is not None:
        _write_sample_identifier_block(buffer, sample_ids)
    variant_offset = len(buffer) - 4
    buffer[0:4] = struct.pack("<I", variant_offset)
    _write_layout2_variant(
        buffer,
        variant_id="var1",
        rsid="rs1",
        chrom="1",
        pos=10,
        alleles=("A", "G"),
        n_samples=2,
    )
    _write_layout2_variant(
        buffer,
        variant_id="var2",
        rsid="rs2",
        chrom="2",
        pos=20,
        alleles=("C", "T"),
        n_samples=2,
    )
    path.write_bytes(buffer)


def _write_sample_file(path: Path, sample_ids: list[str]) -> None:
    rows = ["ID_1 ID_2 missing", "0 0 0"]
    rows.extend(f"{sample_id} {sample_id} 0" for sample_id in sample_ids)
    path.write_text("\n".join(rows) + "\n")


def test_vcf_samples_and_variants_return_source_ordered_polars_frames():
    import genoio

    dataset = genoio.vcf(FIXTURE_ROOT / "vcf" / "tiny.vcf")

    samples = dataset.samples()
    variants = dataset.variants()

    assert isinstance(samples, pl.DataFrame)
    assert samples.columns == ["fid", "iid", "father", "mother", "sex", "phenotype"]
    assert samples["iid"].to_list() == ["S1", "S2", "S3"]
    assert samples["fid"].to_list() == [None, None, None]

    assert isinstance(variants, pl.DataFrame)
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    assert variants.rows() == [
        ("1", 10, "rs1", "A", "G"),
        ("1", 20, "rs2", "C", "T"),
        ("2", 30, "indel1", "AT", "A"),
    ]


def test_rust_metadata_payload_is_column_oriented():
    import genoio._api as api

    metadata = api._rust.read_metadata("vcf", {"vcf": str(FIXTURE_ROOT / "vcf" / "tiny.vcf")})

    assert isinstance(metadata["samples"], dict)
    assert isinstance(metadata["variants"], dict)
    assert metadata["samples"]["iid"] == ["S1", "S2", "S3"]
    assert metadata["variants"]["id"] == ["rs1", "rs2", "indel1"]


def test_plink1_samples_and_variants_normalize_metadata_without_decoding_bed():
    import genoio

    dataset = genoio.bfile(FIXTURE_ROOT / "plink1" / "tiny")

    samples = dataset.samples()
    variants = dataset.variants()

    assert samples.rows() == [
        ("F1", "S1", None, None, "1", "-9"),
        ("F1", "S2", "S1", None, "2", "1.5"),
        ("F2", "S3", None, None, "0", "2.0"),
    ]
    assert variants.rows() == [
        ("1", 10, "rs1", "A", "G"),
        ("1", 20, "rs2", "C", "T"),
        ("2", 30, "indel1", "AT", "A"),
    ]


def test_metadata_is_cached_after_first_load(monkeypatch):
    import genoio
    import genoio._api as api

    dataset = genoio.vcf(FIXTURE_ROOT / "vcf" / "tiny.vcf")
    calls = 0
    original = api._rust.read_metadata

    def counted_read_metadata(format, members):
        nonlocal calls
        calls += 1
        return original(format, members)

    monkeypatch.setattr(api._rust, "read_metadata", counted_read_metadata)

    dataset.samples()
    dataset.variants()
    dataset.samples()

    assert calls == 1


def test_metadata_frames_are_cached_after_first_assembly(monkeypatch):
    import genoio
    import genoio._api as api

    dataset = genoio.vcf(FIXTURE_ROOT / "vcf" / "tiny.vcf")
    sample_calls = 0
    variant_calls = 0
    original_samples_frame = api.samples_frame
    original_variants_frame = api.variants_frame

    def counted_samples_frame(columns):
        nonlocal sample_calls
        sample_calls += 1
        return original_samples_frame(columns)

    def counted_variants_frame(columns):
        nonlocal variant_calls
        variant_calls += 1
        return original_variants_frame(columns)

    monkeypatch.setattr(api, "samples_frame", counted_samples_frame)
    monkeypatch.setattr(api, "variants_frame", counted_variants_frame)

    samples = dataset.samples()
    variants = dataset.variants()

    assert dataset.samples() is samples
    assert dataset.variants() is variants
    assert sample_calls == 1
    assert variant_calls == 1


def test_variant_stats_are_rejected_until_genotype_statistics_exist():
    import genoio

    dataset = genoio.vcf(FIXTURE_ROOT / "vcf" / "tiny.vcf")

    with pytest.raises(genoio.InvalidOptionError, match="variant stats"):
        dataset.variants(stats=["maf"])


def test_plink1_haplotype_reads_are_rejected_by_capabilities():
    import genoio

    dataset = genoio.bfile(FIXTURE_ROOT / "plink1" / "tiny")

    assert dataset._metadata()["capabilities"]["supports_haplo"] is False
    with pytest.raises(genoio.UnsupportedRepresentation, match="haplo"):
        dataset.read(kind="haplo")


def test_malformed_plink1_fam_line_is_reported_as_invalid_source(tmp_path):
    import genoio

    prefix = tmp_path / "bad"
    (tmp_path / "bad.bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x00]))
    (tmp_path / "bad.bim").write_text("1 rs1 0 10 G A\n")
    (tmp_path / "bad.fam").write_text("F1 S1 0 0 1\n")

    dataset = genoio.bfile(prefix)
    with pytest.raises(genoio.InvalidSourceError, match="fam"):
        dataset.samples()


def test_malformed_plink1_bim_line_is_reported_as_invalid_source(tmp_path):
    import genoio

    prefix = tmp_path / "bad"
    (tmp_path / "bad.bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x00]))
    (tmp_path / "bad.bim").write_text("1 rs1 0 10 G\n")
    (tmp_path / "bad.fam").write_text("F1 S1 0 0 1 -9\n")

    dataset = genoio.bfile(prefix)
    with pytest.raises(genoio.InvalidSourceError, match="bim"):
        dataset.variants()


def test_bgen_samples_and_variants_return_embedded_metadata(tmp_path):
    import genoio

    bgen = tmp_path / "tiny.bgen"
    _write_tiny_bgen(bgen, sample_ids=["sample_1", "sample_2"])

    dataset = genoio.bgen(bgen)

    assert dataset.samples()["iid"].to_list() == ["sample_1", "sample_2"]
    assert dataset.variants().rows() == [
        ("1", 10, "rs1", "A", "G"),
        ("2", 20, "rs2", "C", "T"),
    ]


def test_bgen_samples_can_come_from_companion_sample_file(tmp_path):
    import genoio

    prefix = tmp_path / "tiny"
    _write_tiny_bgen(prefix.with_suffix(".bgen"), sample_ids=None)
    _write_sample_file(prefix.with_suffix(".sample"), ["sample_a", "sample_b"])

    dataset = genoio.bgen(prefix)

    assert dataset.samples()["iid"].to_list() == ["sample_a", "sample_b"]
    assert dataset.variants().rows() == [
        ("1", 10, "rs1", "A", "G"),
        ("2", 20, "rs2", "C", "T"),
    ]


def test_bgen_missing_sample_ids_raise_invalid_source_error(tmp_path):
    import genoio

    bgen = tmp_path / "tiny.bgen"
    _write_tiny_bgen(bgen, sample_ids=None)

    dataset = genoio.bgen(bgen)
    with pytest.raises(genoio.InvalidSourceError, match="sample"):
        dataset.samples()
