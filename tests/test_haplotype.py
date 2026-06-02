from pathlib import Path

import numpy as np
import pytest

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
