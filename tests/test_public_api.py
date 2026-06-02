import sys

import pytest


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


def test_dataset_methods_validate_representation_before_later_implementation(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.open(source_path)

    with pytest.raises(genoio.UnsupportedRepresentation):
        dataset.read(representation="unsupported")

    with pytest.raises(NotImplementedError, match="implemented in a later phase"):
        dataset.read()


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
