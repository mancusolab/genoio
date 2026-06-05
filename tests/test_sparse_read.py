from pathlib import Path

import numpy as np
import pytest
from scipy import sparse as scipy_sparse


def write_sparse_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "sparse.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t0/1
"""
    )
    return path


def write_missing_vcf(tmp_path: Path) -> Path:
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


def test_sparse_true_returns_csc_and_preserves_tuple_metadata(tmp_path):
    import genoio

    dataset = genoio.vcf(write_sparse_vcf(tmp_path))

    G, samples, variants = dataset.read(sparse=True, return_samples=True, return_variants=True)

    assert scipy_sparse.isspmatrix_csc(G)
    assert G.shape == (3, 2)
    assert samples["iid"].to_list() == ["S1", "S2", "S3"]
    assert variants["id"].to_list() == ["rs1", "rs2"]


def test_sparse_csr_returns_csr_matrix(tmp_path):
    import genoio

    dataset = genoio.vcf(write_sparse_vcf(tmp_path))

    G = dataset.read(sparse="csr", missing="raise")

    assert scipy_sparse.isspmatrix_csr(G)
    assert G.shape == (3, 2)


def test_sparse_metadata_reports_counted_alleles_for_sparse_orientation(tmp_path):
    import genoio

    dataset = genoio.vcf(write_sparse_vcf(tmp_path))

    sparse_matrix, variants = dataset.read(sparse="csc", missing="raise", return_variants=True)

    assert variants.select("id", "a0", "a1").rows() == [("rs1", "A", "G"), ("rs2", "T", "C")]
    np.testing.assert_array_equal(
        sparse_matrix.toarray(),
        np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 1.0]], dtype=np.float32),
    )


def test_sparse_missing_data_raises_structured_error(tmp_path):
    import genoio

    dataset = genoio.vcf(write_missing_vcf(tmp_path))

    with pytest.raises(genoio.MissingDataError, match="sparse missing values"):
        dataset.read(sparse=True, missing="raise")


@pytest.mark.parametrize("missing", ["nan", "impute"])
def test_sparse_rejects_missing_policies_that_require_stored_missing_values(tmp_path, missing):
    import genoio

    dataset = genoio.vcf(write_sparse_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="does not store sparse missing values"):
        dataset.read(sparse=True, missing=missing)


def test_sparse_rejects_unknown_options_before_calling_rust(tmp_path):
    import genoio

    dataset = genoio.vcf(write_sparse_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="unsupported sparse option"):
        dataset.read(sparse="coo", missing="raise")

    with pytest.raises(genoio.InvalidOptionError, match="unsupported sparse option"):
        dataset.read(sparse=[], missing="raise")  # ty: ignore[invalid-argument-type]


def test_sparse_invalid_missing_policy_raises_structured_error(tmp_path):
    import genoio

    dataset = genoio.vcf(write_sparse_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="unsupported missing-data policy"):
        dataset.read(sparse=True, missing=[])
