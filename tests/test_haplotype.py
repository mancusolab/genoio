# pattern: Imperative Shell

from pathlib import Path

import numpy as np
import pytest
from scipy import sparse as scipy_sparse
from test_dense_read import write_bgen_dosage

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


def write_mixed_phase_stat_filter_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "mixed_phase_stat_filter.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs_phased\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0
1\t20\trs_unphased_monomorphic\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0
"""
    )
    return path


def write_common_a1_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "common_a1.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t1/1\t1/1\t0/1
"""
    )
    return path


def test_plink1_haplotype_reads_raise_unsupported_representation():
    import genoio

    dataset = genoio.bfile(FIXTURE_ROOT / "plink1" / "tiny")

    with pytest.raises(genoio.UnsupportedRepresentation, match="haplo"):
        dataset.read(kind="haplo")


def test_unphased_vcf_haplotype_reads_raise_unsupported_representation(tmp_path):
    import genoio

    dataset = genoio.vcf(write_unphased_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="haplo"):
        dataset.read(kind="haplo")


def test_default_read_is_dense_genotype_read(tmp_path):
    import genoio

    dataset = genoio.vcf(write_unphased_vcf(tmp_path))

    np.testing.assert_array_equal(dataset.read(), dataset.read(kind="geno"))


def test_phased_vcf_haplotype_dense_counts_a1_in_sample_haplotype_order(tmp_path):
    import genoio

    dataset = genoio.vcf(write_phased_vcf(tmp_path))

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
    assert samples["haplotype_index"].to_list() == [0, 1, 0, 1]
    assert samples["source_sample_index"].to_list() == [0, 0, 1, 1]
    assert variants["id"].to_list() == ["rs1", "rs2"]


def test_phased_vcf_haplotype_dosage_reads_raise_hardcall_gt_message(tmp_path):
    import genoio

    dataset = genoio.vcf(write_phased_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="VCF haplotype dosage.*hardcall GT"):
        dataset.read(kind="haplo", dosage="dosage")


def test_filtered_haplotype_read_preserves_source_sample_index(tmp_path):
    import genoio

    dataset = genoio.vcf(write_phased_vcf(tmp_path))

    H, samples = dataset.read(kind="haplo", samples=["S2"], return_samples=True)

    np.testing.assert_array_equal(
        H,
        np.array(
            [
                [1.0, 0.0],
                [0.0, 0.0],
            ],
            dtype=np.float32,
        ),
    )
    assert samples["iid"].to_list() == ["S2", "S2"]
    assert samples["source_sample_index"].to_list() == [1, 1]
    assert samples["haplotype_index"].to_list() == [0, 1]


def test_phased_vcf_haplotype_sparse_uses_requested_sparse_format(tmp_path):
    import genoio

    dataset = genoio.vcf(write_phased_vcf(tmp_path))

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

    dataset = genoio.vcf(write_mixed_phase_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="unphased"):
        dataset.read(kind="haplo", variants=["rs2"])


def test_haplotype_stat_filter_drops_unphased_variant_before_separator_check(tmp_path):
    import genoio

    dataset = genoio.vcf(write_mixed_phase_stat_filter_vcf(tmp_path))

    H, variants = dataset.read(kind="haplo", variants=genoio.maf(min=0.1), return_variants=True)

    np.testing.assert_array_equal(
        H,
        np.array([[0.0], [1.0], [1.0], [0.0]], dtype=np.float32),
    )
    assert variants["id"].to_list() == ["rs_phased"]


def test_sparse_haplotype_stat_filter_drops_unphased_variant_before_separator_check(tmp_path):
    import genoio

    dataset = genoio.vcf(write_mixed_phase_stat_filter_vcf(tmp_path))

    H, variants = dataset.read(kind="haplo", sparse=True, variants=genoio.maf(min=0.1), return_variants=True)

    assert scipy_sparse.isspmatrix_csc(H)
    np.testing.assert_array_equal(
        H.toarray(),
        np.array([[0.0], [1.0], [1.0], [0.0]], dtype=np.float32),
    )
    assert variants["id"].to_list() == ["rs_phased"]


def test_haplotype_blocks_stream_dense_haplotype_columns(tmp_path):
    import genoio

    dataset = genoio.vcf(write_phased_vcf(tmp_path))

    full, full_samples, full_variants = dataset.read(kind="haplo", return_samples=True, return_variants=True)
    blocks = list(dataset.iter_blocks(size=1, kind="haplo", return_samples=True, return_variants=True))

    assert len(blocks) == 2
    np.testing.assert_array_equal(np.concatenate([block[0] for block in blocks], axis=1), full)
    assert blocks[0][1].equals(full_samples)
    assert [block[2]["id"].to_list() for block in blocks] == [["rs1"], ["rs2"]]
    assert [variant_id for block in blocks for variant_id in block[2]["id"].to_list()] == full_variants["id"].to_list()


def test_filtered_haplotype_blocks_preserve_source_sample_index(tmp_path):
    import genoio

    dataset = genoio.vcf(write_phased_vcf(tmp_path))

    blocks = list(dataset.iter_blocks(size=1, kind="haplo", samples=["S2"], return_samples=True))

    assert len(blocks) == 2
    first_block_samples = blocks[0][1]
    np.testing.assert_array_equal(
        np.concatenate([block[0] for block in blocks], axis=1),
        np.array([[1.0, 0.0], [0.0, 0.0]], dtype=np.float32),
    )
    assert first_block_samples["iid"].to_list() == ["S2", "S2"]
    assert first_block_samples["source_sample_index"].to_list() == [1, 1]
    assert first_block_samples["haplotype_index"].to_list() == [0, 1]


def test_haplotype_blocks_stream_sparse_haplotype_columns(tmp_path):
    import genoio

    dataset = genoio.vcf(write_phased_vcf(tmp_path))

    full = dataset.read(kind="haplo", sparse=True)
    blocks = list(dataset.iter_blocks(size=1, kind="haplo", sparse=True))

    assert len(blocks) == 2
    assert all(scipy_sparse.isspmatrix_csc(block) for block in blocks)
    np.testing.assert_array_equal(scipy_sparse.hstack(blocks, format="csc").toarray(), full.toarray())


def test_bgen_haplotype_dosage_blocks_concatenate_to_full_matrix(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, phased=True))

    full = dataset.read(kind="haplo", dosage="dosage")
    blocks = list(dataset.iter_blocks(size=1, kind="haplo", dosage="dosage"))

    assert [block.shape for block in blocks] == [(4, 1), (4, 1)]
    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


def test_bgen_haplotype_dosage_iter_regions_yields_one_result_per_region(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, phased=True))
    regions = [genoio.region("1:1-30"), genoio.region("2:1-30")]

    region_reads = list(dataset.iter_regions(regions, kind="haplo", dosage="dosage", return_variants=True))

    assert [region for region, _ in region_reads] == regions
    assert [variants["id"].to_list() for _, (_, variants) in region_reads] == [["rs1"], ["rs2"]]
    assert [matrix.shape for _, (matrix, _) in region_reads] == [(4, 1), (4, 1)]


def test_bgen_haplotype_dosage_rejects_unphased_retained_records(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, phased=False))

    with pytest.raises(genoio.UnsupportedRepresentation, match="unphased"):
        dataset.read(kind="haplo", dosage="dosage")


def test_sparse_genotype_reads_still_minor_allele_flip_by_default(tmp_path):
    import genoio

    dataset = genoio.vcf(write_common_a1_vcf(tmp_path))

    G, variants = dataset.read(kind="geno", sparse=True, return_variants=True)

    assert scipy_sparse.isspmatrix_csc(G)
    np.testing.assert_array_equal(G.toarray(), np.array([[0.0], [0.0], [1.0]], dtype=np.float32))
    assert variants["a0"].to_list() == ["G"]
    assert variants["a1"].to_list() == ["A"]
