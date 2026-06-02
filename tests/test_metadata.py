from pathlib import Path

import polars as pl
import pytest

FIXTURE_ROOT = Path(__file__).parent / "fixtures"


def test_vcf_samples_and_variants_return_source_ordered_polars_frames():
    import genoio

    dataset = genoio.open(FIXTURE_ROOT / "vcf" / "tiny.vcf")

    samples = dataset.samples()
    variants = dataset.variants()

    assert isinstance(samples, pl.DataFrame)
    assert samples.columns == ["fid", "iid", "father", "mother", "sex", "phenotype"]
    assert samples["iid"].to_list() == ["S1", "S2", "S3"]
    assert samples["fid"].to_list() == [None, None, None]

    assert isinstance(variants, pl.DataFrame)
    assert variants.columns == [
        "chrom",
        "pos",
        "id",
        "a0",
        "a1",
        "ref",
        "alt",
        "source_a0",
        "source_a1",
        "flipped",
        "af",
        "maf",
        "mac",
        "missing_rate",
        "n_called",
    ]
    assert variants.select("chrom", "pos", "id", "a0", "a1", "ref", "alt").rows() == [
        ("1", 10, "rs1", "A", "G", "A", "G"),
        ("1", 20, "rs2", "C", "T", "C", "T,A"),
        ("2", 30, "indel1", "AT", "A", "AT", "A"),
    ]
    assert variants["flipped"].to_list() == [False, False, False]
    assert variants["af"].to_list() == [None, None, None]


def test_plink1_samples_and_variants_normalize_metadata_without_decoding_bed():
    import genoio

    dataset = genoio.open(FIXTURE_ROOT / "plink1" / "tiny")

    samples = dataset.samples()
    variants = dataset.variants()

    assert samples.rows() == [
        ("F1", "S1", None, None, "1", "-9"),
        ("F1", "S2", "S1", None, "2", "1.5"),
        ("F2", "S3", None, None, "0", "2.0"),
    ]
    assert variants.select("chrom", "pos", "id", "a0", "a1", "ref", "alt").rows() == [
        ("1", 10, "rs1", "A", "G", None, None),
        ("1", 20, "rs2", "C", "T", None, None),
        ("2", 30, "indel1", "AT", "A", None, None),
    ]


def test_metadata_is_cached_after_first_load(monkeypatch):
    import genoio
    import genoio._api as api

    dataset = genoio.open(FIXTURE_ROOT / "vcf" / "tiny.vcf")
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


def test_variant_stats_are_rejected_until_genotype_statistics_exist():
    import genoio

    dataset = genoio.open(FIXTURE_ROOT / "vcf" / "tiny.vcf")

    with pytest.raises(genoio.InvalidOptionError, match="variant stats"):
        dataset.variants(stats=["maf"])


def test_plink1_haplotype_reads_are_rejected_by_capabilities():
    import genoio

    dataset = genoio.open(FIXTURE_ROOT / "plink1" / "tiny")

    assert dataset._metadata()["capabilities"]["supports_haplo"] is False
    with pytest.raises(genoio.UnsupportedRepresentation, match="haplo"):
        dataset.read(kind="haplo")


def test_malformed_plink1_fam_line_is_reported_as_invalid_source(tmp_path):
    import genoio

    prefix = tmp_path / "bad"
    (tmp_path / "bad.bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x00]))
    (tmp_path / "bad.bim").write_text("1 rs1 0 10 G A\n")
    (tmp_path / "bad.fam").write_text("F1 S1 0 0 1\n")

    dataset = genoio.open(prefix)
    with pytest.raises(genoio.InvalidSourceError, match="fam"):
        dataset.samples()


def test_malformed_plink1_bim_line_is_reported_as_invalid_source(tmp_path):
    import genoio

    prefix = tmp_path / "bad"
    (tmp_path / "bad.bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x00]))
    (tmp_path / "bad.bim").write_text("1 rs1 0 10 G\n")
    (tmp_path / "bad.fam").write_text("F1 S1 0 0 1 -9\n")

    dataset = genoio.open(prefix)
    with pytest.raises(genoio.InvalidSourceError, match="bim"):
        dataset.variants()
