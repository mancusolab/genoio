import sys
from inspect import signature
from pathlib import Path

import numpy as np
import pytest

FIXTURE_ROOT = Path(__file__).parent / "fixtures"


def write_blocks_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "blocks.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1
"""
    )
    return path


def test_import_exposes_public_names_without_reference_packages():
    import genoio

    expected_names = {
        "open",
        "read",
        "samples",
        "variants",
        "blocks",
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


def test_dataset_read_recognizes_sparse_options_and_validates_missing_policy(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.open(source_path)

    with pytest.raises(genoio.InvalidOptionError, match="sparse missing values"):
        dataset.read(sparse="csc", missing="nan")


def test_sparse_default_missing_signature_is_readable():
    import genoio

    assert "missing: 'Any' = DEFAULT_MISSING" in str(signature(genoio.read))


def test_dataset_blocks_accepts_read_options_and_validates_size(tmp_path):
    import genoio

    source_path = tmp_path / "cohort.vcf.gz"
    source_path.touch()
    dataset = genoio.open(source_path)

    read_options = {
        "kind": "geno",
        "sparse": False,
        "variants": None,
        "samples": ("s2", "s1"),
        "missing": "raise",
        "dtype": "uint8",
        "return_samples": False,
        "return_variants": True,
    }

    block_iterator = dataset.blocks(8192, **read_options)

    assert iter(block_iterator) is block_iterator

    with pytest.raises(genoio.InvalidOptionError, match="block size"):
        dataset.blocks(0, **read_options)


def test_top_level_blocks_matches_dataset_blocks(tmp_path):
    import genoio

    path = write_blocks_vcf(tmp_path)

    top_level_blocks = list(genoio.blocks(path, size=2, return_variants=True))
    dataset_blocks = list(genoio.open(path).blocks(size=2, return_variants=True))

    assert len(top_level_blocks) == len(dataset_blocks)
    for top_level_block, dataset_block in zip(top_level_blocks, dataset_blocks, strict=True):
        np.testing.assert_array_equal(top_level_block[0], dataset_block[0])
        assert top_level_block[1].equals(dataset_block[1])


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

    with pytest.raises(genoio.InvalidOptionError):
        dataset.read(sparse="unsupported")


def test_filter_helpers_build_serializable_expressions():
    import genoio

    expression = genoio.chrom("1") & genoio.region("1:10-20") & ~genoio.missing_rate(max=0.05)

    assert expression.to_ir() == {
        "op": "and",
        "left": {
            "op": "and",
            "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
            "right": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}},
        },
        "right": {
            "op": "not",
            "expr": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.05}},
        },
    }


def test_region_rejects_malformed_region_syntax():
    import genoio

    with pytest.raises(genoio.InvalidOptionError):
        genoio.region("not-a-region")


def test_snp_helper_is_zero_argument_snp_only_predicate():
    import genoio

    assert genoio.snp().to_ir() == {"op": "predicate", "name": "snp", "params": {}}


def test_id_in_helper_is_variant_id_matching_predicate():
    import genoio

    assert genoio.id_in(["rs2", "rs1"]).to_ir() == {
        "op": "predicate",
        "name": "id_in",
        "params": {"values": ["rs2", "rs1"]},
    }
