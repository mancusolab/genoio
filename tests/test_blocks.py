from pathlib import Path

import numpy as np
import pytest
from scipy import sparse as scipy_sparse


def write_blocks_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "blocks.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1
2\t30\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0
2\t40\trs4\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1
1\t50\trs5\tA\tC\t.\tPASS\t.\tGT\t1/1\t1/1\t0/1
"""
    )
    return path


def test_blocks_honor_size_and_concatenate_to_full_dense_read(tmp_path):
    import genoio

    dataset = genoio.open(write_blocks_vcf(tmp_path))

    full = dataset.read()
    blocks = list(dataset.blocks(size=2))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]
    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


def test_blocks_variant_metadata_aligns_with_each_block(tmp_path):
    import genoio

    dataset = genoio.open(write_blocks_vcf(tmp_path))

    blocks = list(dataset.blocks(size=2, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs3", "rs4"], ["rs5"]]
    for G_block, variants in blocks:
        assert G_block.shape[1] == len(variants)


def test_blocks_return_samples_keeps_sample_order_constant(tmp_path):
    import genoio

    dataset = genoio.open(write_blocks_vcf(tmp_path))

    blocks = list(dataset.blocks(size=3, samples=["S3", "S1"], return_samples=True, return_variants=True))

    assert len(blocks) == 2
    for G_block, samples, variants in blocks:
        assert G_block.shape[0] == 2
        assert G_block.shape[1] == len(variants)
        assert samples["iid"].to_list() == ["S1", "S3"]


def test_blocks_apply_filters_and_sample_keep_lists_like_full_reads(tmp_path):
    import genoio

    dataset = genoio.open(write_blocks_vcf(tmp_path))
    read_options = {"variants": genoio.chrom("1"), "samples": ["S3", "S1"]}

    full, full_variants = dataset.read(**read_options, return_variants=True)
    blocks = list(dataset.blocks(size=2, **read_options, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs5"]]
    np.testing.assert_array_equal(np.concatenate([block for block, _ in blocks], axis=1), full)
    assert full_variants["id"].to_list() == ["rs1", "rs2", "rs5"]


def test_blocks_validate_size(tmp_path):
    import genoio

    dataset = genoio.open(write_blocks_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="positive integer"):
        list(dataset.blocks(size=0))


def test_blocks_request_bounded_variant_windows_at_rust_boundary(tmp_path, monkeypatch):
    import genoio

    dataset = genoio.open(write_blocks_vcf(tmp_path))
    calls = []
    original = genoio.Dataset._read_block_from_rust

    def recording_read_block(self, *, variant_window, **kwargs):
        calls.append(dict(variant_window))
        return original(self, variant_window=variant_window, **kwargs)

    monkeypatch.setattr(genoio.Dataset, "_read_block_from_rust", recording_read_block)

    list(dataset.blocks(size=2))

    assert calls
    assert all(call["len"] <= 2 for call in calls)
    assert [call["start"] for call in calls] == [0, 2, 4]


def test_blocks_do_not_call_public_read_internally(tmp_path, monkeypatch):
    import genoio

    dataset = genoio.open(write_blocks_vcf(tmp_path))

    def fail_read(*args, **kwargs):
        raise AssertionError("blocks must use the bounded Rust call boundary, not public read()")

    monkeypatch.setattr(genoio.Dataset, "read", fail_read)

    blocks = list(dataset.blocks(size=2))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]


def test_sparse_blocks_work_with_default_missing_policy(tmp_path):
    import genoio

    dataset = genoio.open(write_blocks_vcf(tmp_path))

    blocks = list(dataset.blocks(size=2, sparse=True))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]
    assert all(scipy_sparse.isspmatrix_csc(block) for block in blocks)
