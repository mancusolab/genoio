import sys
from pathlib import Path

import pytest

FIXTURE_ROOT = Path(__file__).parent / "fixtures"


def test_import_exposes_public_names_without_reference_packages():
    import genoio

    expected_names = {
        "open",
        "read",
        "samples",
        "variants",
        "chrom",
        "region",
        "snp",
        "biallelic",
        "maf",
        "mac",
        "missing_rate",
        "polymorphic",
        "id_in",
        "GenoioError",
        "AmbiguousSourceError",
        "MissingCompanionFileError",
        "UnsupportedFormatError",
        "InvalidSourceError",
        "UnsupportedRepresentation",
        "InvalidOptionError",
    }

    assert expected_names <= set(dir(genoio))
    assert "jaxqtl" not in sys.modules
    assert "linear_dag" not in sys.modules


def test_open_returns_lightweight_dataset_for_vcf(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf"
    source_path.touch()

    dataset = genoio.open(source_path)

    assert dataset.source.path == source_path
    assert dataset.source.format.value == "vcf"


def test_dataset_read_accepts_approved_options_before_later_implementation(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.open(source_path)

    read_options = {
        "kind": "geno",
        "sparse": "csc",
        "variants": genoio.snp() & genoio.id_in(["rs1", "rs2"]),
        "samples": ["s2", "s1"],
        "missing": "nan",
        "dtype": "float32",
        "return_samples": True,
        "return_variants": True,
    }

    with pytest.raises(NotImplementedError, match="implemented in a later phase"):
        dataset.read(**read_options)


def test_dataset_blocks_accepts_read_options_and_validates_size(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.open(source_path)

    read_options = {
        "kind": "geno",
        "sparse": False,
        "variants": genoio.chrom("1"),
        "samples": ("s2", "s1"),
        "missing": "raise",
        "dtype": "uint8",
        "return_samples": False,
        "return_variants": True,
    }

    with pytest.raises(NotImplementedError, match="implemented in a later phase"):
        dataset.blocks(8192, **read_options)

    with pytest.raises(genoio.InvalidOptionError, match="block size"):
        dataset.blocks(0, **read_options)


def test_dataset_variants_accepts_documented_default_stats_keyword():
    import genoio

    dataset = genoio.open(FIXTURE_ROOT / "vcf" / "tiny.vcf")

    variants = dataset.variants(stats=None)

    assert variants["id"].to_list() == ["rs1", "rs2", "indel1"]


def test_top_level_variants_accepts_documented_default_stats_keyword():
    import genoio

    variants = genoio.variants(FIXTURE_ROOT / "vcf" / "tiny.vcf", stats=None)

    assert variants["id"].to_list() == ["rs1", "rs2", "indel1"]


def test_dataset_variants_rejects_stats_until_stat_metadata_is_implemented(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.open(source_path)

    with pytest.raises(genoio.InvalidOptionError, match="variant stats"):
        dataset.variants(stats=["maf"])


def test_dataset_read_rejects_unsupported_representation_options(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.open(source_path)

    with pytest.raises(genoio.UnsupportedRepresentation):
        dataset.read(kind="unsupported")

    with pytest.raises(genoio.UnsupportedRepresentation):
        dataset.read(sparse="unsupported")


def test_filter_helpers_build_serializable_expressions():
    import genoio

    expression = genoio.chrom("1") & genoio.region("1:10-20") & ~genoio.missing_rate(max=0.05)

    assert expression.to_ir() == {
        "op": "and",
        "args": [
            {
                "op": "and",
                "args": [
                    {"op": "chrom", "value": "1"},
                    {"op": "region", "value": "1:10-20"},
                ],
            },
            {"op": "not", "arg": {"op": "missing_rate", "max": 0.05}},
        ],
    }


def test_region_rejects_malformed_region_syntax():
    import genoio

    with pytest.raises(genoio.InvalidOptionError):
        genoio.region("not-a-region")


def test_snp_helper_is_zero_argument_snp_only_predicate():
    import genoio

    assert genoio.snp().to_ir() == {"op": "snp", "value": None}


def test_id_in_helper_is_variant_id_matching_predicate():
    import genoio

    assert genoio.id_in(["rs2", "rs1"]).to_ir() == {"op": "id_in", "value": ("rs2", "rs1")}
