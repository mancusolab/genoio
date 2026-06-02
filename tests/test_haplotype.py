from pathlib import Path

import numpy as np
import pytest
from scipy import sparse as scipy_sparse

FIXTURE_ROOT = Path(__file__).parent / "fixtures"


def write_unphased_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "unphased.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t1/1\t0/0
"""
    )
    return path


def write_phased_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "phased.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|0
"""
    )
    return path


def write_mixed_phase_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "mixed_phase.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t1/1\t0|0
"""
    )
    return path


def test_plink1_haplotype_reads_raise_unsupported_representation():
    import genoio

    dataset = genoio.open(FIXTURE_ROOT / "plink1" / "tiny")

    with pytest.raises(genoio.UnsupportedRepresentation, match="haplo"):
        dataset.read(kind="haplo")


def test_unphased_vcf_haplotype_reads_raise_unsupported_representation(tmp_path):
    import genoio

    dataset = genoio.open(write_unphased_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="haplo"):
        dataset.read(kind="haplo")


def test_default_read_is_dense_genotype_read(tmp_path):
    import genoio

    dataset = genoio.open(write_unphased_vcf(tmp_path))

    np.testing.assert_array_equal(dataset.read(), dataset.read(kind="geno"))


def test_phased_vcf_haplotype_dense_counts_a1_in_sample_haplotype_order(tmp_path):
    import genoio

    dataset = genoio.open(write_phased_vcf(tmp_path))

    H, samples, variants = dataset.read(kind="haplo", return_samples=True, return_variants=True)

    np.testing.assert_array_equal(
        H,
        np.array(
            [
                [0.0, 1.0],
                [1.0, 1.0],
                [1.0, 0.0],
                [0.0, 0.0],
            ],
            dtype=np.float32,
        ),
    )
    assert samples["iid"].to_list() == ["S1", "S1", "S2", "S2"]
    assert variants["id"].to_list() == ["rs1", "rs2"]


def test_phased_vcf_haplotype_sparse_uses_requested_sparse_format(tmp_path):
    import genoio

    dataset = genoio.open(write_phased_vcf(tmp_path))

    H_csc = dataset.read(kind="haplo", sparse=True)
    H_csr = dataset.read(kind="haplo", sparse="csr")

    assert scipy_sparse.isspmatrix_csc(H_csc)
    assert scipy_sparse.isspmatrix_csr(H_csr)
    np.testing.assert_array_equal(
        H_csc.toarray(),
        np.array(
            [
                [0.0, 1.0],
                [1.0, 1.0],
                [1.0, 0.0],
                [0.0, 0.0],
            ],
            dtype=np.float32,
        ),
    )
    np.testing.assert_array_equal(H_csr.toarray(), H_csc.toarray())


def test_haplotype_read_rejects_unphased_separator_in_retained_variant(tmp_path):
    import genoio

    dataset = genoio.open(write_mixed_phase_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="unphased"):
        dataset.read(kind="haplo", variants=["rs2"])
