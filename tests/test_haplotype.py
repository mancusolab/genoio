# pattern: Imperative Shell

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from fixture_writers import (
    write_bgen_dosage,
    write_ld_phased_dosage_plink2,
    write_ld_phased_hardcall_plink2,
    write_phased_dosage_plink2,
    write_phased_hardcall_plink2,
    write_sample_filtered_unphased_hardcall_plink2,
)
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


def write_indexed_phased_vcf_with_outside_region(tmp_path: Path) -> Path:
    source = tmp_path / "indexed_phased.vcf"
    source.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0
1\t20\trs20\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|0
1\t40\toutside_region\tT\tC\t.\tPASS\t.\tGT\t0|0\t0|0
"""
    )
    compressed = tmp_path / "indexed_phased.vcf.gz"
    with compressed.open("wb") as output:
        subprocess.run(["bgzip", "-c", str(source)], stdout=output, check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(compressed)], check=True)
    return compressed


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


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="indexed VCF haplotype test requires bgzip and tabix",
)
@pytest.mark.parametrize("read_name", ["read_haplotypes_dense", "read_haplotypes_sparse"])
def test_indexed_vcf_region_pushdown_applies_to_haplotype_reads(tmp_path, read_name):
    import genoio
    from genoio import _rust
    from genoio._read_options import _variant_filter_ir

    path = write_indexed_phased_vcf_with_outside_region(tmp_path)
    options = {
        "samples": None,
        "variants": _variant_filter_ir(genoio.region("1:10-20")),
        "variant_window": None,
        "dosage": "hardcall",
        "return_samples": False,
        "return_variants": False,
        "matrix_only": False,
    }

    result = getattr(_rust, read_name)("vcf", {"vcf": str(path)}, options)

    assert result["diagnostics"]["candidate_variants"] == 2
    assert result["diagnostics"]["retained_variants"] == 2
    assert tuple(result["shape"]) == (4, 2)


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


def test_plink2_haplotype_hardcall_read_returns_haplotype_rows_and_metadata(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_hardcall_plink2(tmp_path))

    H, samples, variants = dataset.read(
        kind="haplo",
        dosage="hardcall",
        return_samples=True,
        return_variants=True,
    )

    np.testing.assert_array_equal(
        H,
        np.array(
            [
                [0.0, 1.0],
                [1.0, 0.0],
                [0.0, 0.0],
                [0.0, 1.0],
                [1.0, np.nan],
                [1.0, np.nan],
            ],
            dtype=np.float32,
        ),
    )
    assert samples["iid"].to_list() == ["S1", "S1", "S2", "S2", "S3", "S3"]
    assert samples["source_sample_index"].to_list() == [0, 0, 1, 1, 2, 2]
    assert samples["haplotype_index"].to_list() == [0, 1, 0, 1, 0, 1]
    assert variants["id"].to_list() == ["rs1", "rs2"]
    assert H.shape[1] == len(variants)


def test_plink2_haplotype_dosage_read_returns_haplotype_dosage_rows(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_dosage_plink2(tmp_path))

    H, samples, variants = dataset.read(
        kind="haplo",
        dosage="dosage",
        return_samples=True,
        return_variants=True,
    )

    np.testing.assert_allclose(
        H,
        np.array(
            [
                [0.25, 0.0],
                [0.75, 0.0],
                [0.0, 0.1],
                [0.5, 0.1],
                [1.0, 0.2],
                [1.0, 0.2],
            ],
            dtype=np.float32,
        ),
        rtol=0,
        atol=2.0 / 32768.0,
    )
    assert samples["iid"].to_list() == ["S1", "S1", "S2", "S2", "S3", "S3"]
    assert samples["source_sample_index"].to_list() == [0, 0, 1, 1, 2, 2]
    assert samples["haplotype_index"].to_list() == [0, 1, 0, 1, 0, 1]
    assert variants["id"].to_list() == ["rs1", "rs2"]
    assert H.shape[1] == len(variants)


def test_plink2_haplotype_sample_filter_preserves_source_order_and_haplotype_rows(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_hardcall_plink2(tmp_path))

    H, samples = dataset.read(kind="haplo", samples=["S3", "S1"], return_samples=True)

    np.testing.assert_array_equal(
        H,
        np.array(
            [
                [0.0, 1.0],
                [1.0, 0.0],
                [1.0, np.nan],
                [1.0, np.nan],
            ],
            dtype=np.float32,
        ),
    )
    assert samples["iid"].to_list() == ["S1", "S1", "S3", "S3"]
    assert samples["source_sample_index"].to_list() == [0, 0, 2, 2]
    assert samples["haplotype_index"].to_list() == [0, 1, 0, 1]


@pytest.mark.parametrize(
    ("writer", "read_options", "assert_matrix"),
    [
        (
            write_phased_hardcall_plink2,
            {"kind": "haplo", "dosage": "hardcall"},
            np.testing.assert_array_equal,
        ),
        (
            write_phased_dosage_plink2,
            {"kind": "haplo", "dosage": "dosage"},
            np.testing.assert_allclose,
        ),
    ],
)
def test_plink2_haplotype_blocks_concatenate_to_full_read(tmp_path, writer, read_options, assert_matrix):
    import genoio

    dataset = genoio.pfile(writer(tmp_path))

    full, full_variants = dataset.read(**read_options, return_variants=True)
    blocks = list(dataset.iter_blocks(size=1, **read_options, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1"], ["rs2"]]
    assert_matrix(np.concatenate([block for block, _ in blocks], axis=1), full)
    assert [variant_id for _, variants in blocks for variant_id in variants["id"].to_list()] == full_variants[
        "id"
    ].to_list()


@pytest.mark.parametrize(
    ("writer", "read_options", "assert_matrix"),
    [
        (
            write_ld_phased_hardcall_plink2,
            {"kind": "haplo", "dosage": "hardcall"},
            np.testing.assert_array_equal,
        ),
        (
            write_ld_phased_dosage_plink2,
            {"kind": "haplo", "dosage": "dosage"},
            np.testing.assert_allclose,
        ),
    ],
)
def test_plink2_haplotype_blocks_decode_ld_compressed_second_variant(tmp_path, writer, read_options, assert_matrix):
    import genoio

    dataset = genoio.pfile(writer(tmp_path))

    full = dataset.read(**read_options)
    blocks = list(dataset.iter_blocks(size=1, **read_options, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1"], ["rs2"]]
    assert_matrix(blocks[1][0], full[:, 1:2])


@pytest.mark.parametrize(
    ("writer", "read_options"),
    [
        (write_phased_hardcall_plink2, {"kind": "haplo", "dosage": "hardcall"}),
        (write_phased_dosage_plink2, {"kind": "haplo", "dosage": "dosage"}),
    ],
)
def test_plink2_haplotype_iter_regions_yields_one_result_per_region(tmp_path, writer, read_options):
    import genoio

    dataset = genoio.pfile(writer(tmp_path))
    regions = [genoio.region("1:1-15"), genoio.region("1:16-25")]

    region_reads = list(dataset.iter_regions(regions, **read_options, return_variants=True))

    assert [region for region, _ in region_reads] == regions
    assert [variants["id"].to_list() for _, (_, variants) in region_reads] == [["rs1"], ["rs2"]]
    assert [matrix.shape for _, (matrix, _) in region_reads] == [(6, 1), (6, 1)]


@pytest.mark.parametrize(
    ("writer", "read_options", "assert_matrix"),
    [
        (
            write_ld_phased_hardcall_plink2,
            {"kind": "haplo", "dosage": "hardcall"},
            np.testing.assert_array_equal,
        ),
        (
            write_ld_phased_dosage_plink2,
            {"kind": "haplo", "dosage": "dosage"},
            np.testing.assert_allclose,
        ),
    ],
)
def test_plink2_haplotype_region_decode_ld_compressed_second_variant(tmp_path, writer, read_options, assert_matrix):
    import genoio

    dataset = genoio.pfile(writer(tmp_path))
    full = dataset.read(**read_options)
    regions = [genoio.region("1:16-25")]

    region_reads = list(dataset.iter_regions(regions, **read_options, return_variants=True))

    assert [region for region, _ in region_reads] == regions
    matrix, variants = region_reads[0][1]
    assert variants["id"].to_list() == ["rs2"]
    assert_matrix(matrix, full[:, 1:2])


def test_plink2_haplotype_sample_filter_ignores_unselected_unphased_heterozygote(tmp_path):
    import genoio

    dataset = genoio.pfile(write_sample_filtered_unphased_hardcall_plink2(tmp_path))

    H, samples = dataset.read(kind="haplo", samples=["S2", "S3"], return_samples=True)

    np.testing.assert_array_equal(
        H,
        np.array([[0.0], [1.0], [1.0], [0.0]], dtype=np.float32),
    )
    assert samples["iid"].to_list() == ["S2", "S2", "S3", "S3"]


def test_plink2_haplotype_unphased_retained_record_fails(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_hardcall_plink2(tmp_path, unphased_second_variant=True))

    with pytest.raises(genoio.UnsupportedRepresentation, match="unphased"):
        dataset.read(kind="haplo", variants=["rs2"])


def test_plink2_haplotype_metadata_filter_skips_unsupported_record(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_dosage_plink2(tmp_path, unphased_second_variant=True))

    H, variants = dataset.read(kind="haplo", dosage="dosage", variants=["rs1"], return_variants=True)

    np.testing.assert_allclose(
        H,
        np.array([[0.25], [0.75], [0.0], [0.5], [1.0], [1.0]], dtype=np.float32),
        rtol=0,
        atol=2.0 / 32768.0,
    )
    assert variants["id"].to_list() == ["rs1"]


def test_plink2_sparse_hardcall_haplotypes_match_dense(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_hardcall_plink2(tmp_path))

    dense, dense_samples, dense_variants = dataset.read(
        kind="haplo",
        dosage="hardcall",
        samples=["S1", "S2"],
        return_samples=True,
        return_variants=True,
    )
    H, samples, variants = dataset.read(
        kind="haplo",
        dosage="hardcall",
        sparse=True,
        samples=["S1", "S2"],
        return_samples=True,
        return_variants=True,
    )

    assert scipy_sparse.isspmatrix_csc(H)
    np.testing.assert_array_equal(H.toarray(), dense)
    assert samples.equals(dense_samples)
    assert variants.equals(dense_variants)


def test_plink2_sparse_hardcall_haplotypes_reject_retained_missing(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_hardcall_plink2(tmp_path))

    with pytest.raises(genoio.MissingDataError, match="sparse missing values"):
        dataset.read(kind="haplo", dosage="hardcall", sparse=True)


def test_plink2_sparse_dosage_haplotypes_still_fail(tmp_path):
    import genoio

    dataset = genoio.pfile(write_phased_dosage_plink2(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="sparse haplotype reads.*dense"):
        dataset.read(kind="haplo", dosage="dosage", sparse=True)


def test_plink2_sparse_hardcall_haplotype_blocks_concatenate(tmp_path):
    import genoio

    dataset = genoio.pfile(write_ld_phased_hardcall_plink2(tmp_path))

    full = dataset.read(kind="haplo", dosage="hardcall", sparse=True, samples=["S1", "S2"])
    blocks = list(
        dataset.iter_blocks(
            size=1,
            kind="haplo",
            dosage="hardcall",
            sparse=True,
            samples=["S1", "S2"],
        )
    )

    assert len(blocks) == 2
    assert all(scipy_sparse.isspmatrix_csc(block) for block in blocks)
    np.testing.assert_array_equal(scipy_sparse.hstack(blocks, format="csc").toarray(), full.toarray())


def test_sparse_genotype_reads_still_minor_allele_flip_by_default(tmp_path):
    import genoio

    dataset = genoio.vcf(write_common_a1_vcf(tmp_path))

    G, variants = dataset.read(kind="geno", sparse=True, return_variants=True)

    assert scipy_sparse.isspmatrix_csc(G)
    np.testing.assert_array_equal(G.toarray(), np.array([[0.0], [0.0], [1.0]], dtype=np.float32))
    assert variants["a0"].to_list() == ["G"]
    assert variants["a1"].to_list() == ["A"]
