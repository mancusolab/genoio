from pathlib import Path

import numpy as np
import polars as pl
import pytest

FIXTURE_ROOT = Path(__file__).parent / "fixtures"


def write_biallelic_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "tiny.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t./.\t0/0
"""
    )
    return path


def write_multi_alt_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "multi_alt.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG,T\t.\tPASS\t.\tGT\t0/1
"""
    )
    return path


def write_all_missing_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "all_missing.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t./.\t./.
"""
    )
    return path


def test_dense_vcf_read_returns_sample_by_variant_numpy_array_and_metadata(tmp_path):
    import genoio

    dataset = genoio.open(write_biallelic_vcf(tmp_path))

    G, samples, variants = dataset.read(return_samples=True, return_variants=True)

    assert isinstance(G, np.ndarray)
    assert G.dtype == np.dtype("float32")
    np.testing.assert_array_equal(G, np.array([[0.0, 1.0], [1.0, np.nan], [2.0, 0.0]], dtype=np.float32))
    assert isinstance(samples, pl.DataFrame)
    assert samples["iid"].to_list() == ["S1", "S2", "S3"]
    assert isinstance(variants, pl.DataFrame)
    assert variants["id"].to_list() == ["rs1", "rs2"]


def test_dense_plink1_read_matches_fixture_matrix_and_return_tuple_shapes():
    import genoio

    dataset = genoio.open(FIXTURE_ROOT / "plink1" / "tiny")

    G = dataset.read()
    G_samples = dataset.read(return_samples=True)
    G_variants = dataset.read(return_variants=True)
    G_both = dataset.read(return_samples=True, return_variants=True)

    expected = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, np.nan], [2.0, np.nan, 0.0]], dtype=np.float32)
    np.testing.assert_array_equal(G, expected)
    assert isinstance(G_samples, tuple)
    assert isinstance(G_variants, tuple)
    assert isinstance(G_both, tuple)
    assert len(G_samples) == 2
    assert len(G_variants) == 2
    assert len(G_both) == 3
    np.testing.assert_array_equal(G_samples[0], expected)
    assert G_samples[1]["iid"].to_list() == ["S1", "S2", "S3"]
    np.testing.assert_array_equal(G_variants[0], expected)
    assert G_variants[1]["id"].to_list() == ["rs1", "rs2", "indel1"]
    np.testing.assert_array_equal(G_both[0], expected)
    assert G_both[1]["iid"].to_list() == ["S1", "S2", "S3"]
    assert G_both[2]["id"].to_list() == ["rs1", "rs2", "indel1"]


def test_missing_policies_nan_raise_and_impute(tmp_path):
    import genoio

    dataset = genoio.open(write_biallelic_vcf(tmp_path))

    with_nan = dataset.read(missing="nan", dtype="float64")
    with_impute = dataset.read(missing="impute")

    assert with_nan.dtype == np.dtype("float64")
    assert np.isnan(with_nan[1, 1])
    np.testing.assert_array_equal(with_impute, np.array([[0.0, 1.0], [1.0, 0.5], [2.0, 0.0]], dtype=np.float32))
    with pytest.raises(genoio.MissingDataError, match="missing genotype"):
        dataset.read(missing="raise")


def test_missing_policy_impute_rejects_all_missing_variant(tmp_path):
    import genoio

    dataset = genoio.open(write_all_missing_vcf(tmp_path))

    with pytest.raises(genoio.MissingDataError, match="all-missing variant"):
        dataset.read(missing="impute")


def test_missing_policy_rejects_integer_dtype_combinations(tmp_path):
    import genoio

    dataset = genoio.open(write_biallelic_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match='missing="nan"'):
        dataset.read(dtype=np.int16, missing="nan")
    with pytest.raises(genoio.InvalidOptionError, match='missing="impute"'):
        dataset.read(dtype=np.int16, missing="impute")
    assert dataset.read(dtype=np.int16, missing="raise", samples=["S1"]).dtype == np.dtype("int16")


def test_read_option_validation_prioritizes_dtype_and_missing_before_samples(tmp_path):
    import genoio

    dataset = genoio.open(write_biallelic_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match='missing="nan"'):
        dataset.read(dtype=np.int16, missing="nan", samples=["S1", "S1"])


def test_dense_vcf_read_rejects_multi_alt_records_even_when_gt_uses_first_alt(tmp_path):
    import genoio

    dataset = genoio.open(write_multi_alt_vcf(tmp_path))

    with pytest.raises(genoio.InvalidSourceError, match="multi-ALT"):
        dataset.read()


def test_unordered_sample_keep_list_returns_rows_in_source_order_and_metadata_matches(tmp_path):
    import genoio

    dataset = genoio.open(write_biallelic_vcf(tmp_path))

    G, samples = dataset.read(samples=["S3", "S1"], return_samples=True)

    np.testing.assert_array_equal(G, np.array([[0.0, 1.0], [2.0, 0.0]], dtype=np.float32))
    assert samples["iid"].to_list() == ["S1", "S3"]


def test_missing_requested_sample_raises_structured_error_with_counts(tmp_path):
    import genoio

    dataset = genoio.open(write_biallelic_vcf(tmp_path))

    with pytest.raises(genoio.SampleFilterError, match="requested=2 retained=1 missing=1"):
        dataset.read(samples=["S1", "S4"])
