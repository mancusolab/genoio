# pattern: Imperative Shell

from pathlib import Path

import numpy as np
import pytest
from scipy import sparse as scipy_sparse
from test_dense_read import write_fixed_width_plink2


def write_bad_variable_width_block_offset_plink2(tmp_path: Path) -> Path:
    prefix = tmp_path / "bad_offset"
    record = bytes([0x00])
    header_len = 12 + 8 + 1 + 1
    bad_first_block_offset = header_len - 1
    prefix.with_suffix(".pgen").write_bytes(
        b"\x6c\x1b\x10"
        + (1).to_bytes(4, "little")
        + (4).to_bytes(4, "little")
        + bytes([0x04])
        + bad_first_block_offset.to_bytes(8, "little")
        + bytes([0x00])
        + bytes([len(record)])
        + record
    )
    prefix.with_suffix(".pvar").write_text(
        """\
#CHROM POS ID REF ALT
1 10 rs1 A G
"""
    )
    prefix.with_suffix(".psam").write_text(
        """\
#IID
S1
S2
S3
S4
"""
    )
    return prefix


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


def empty_dense_rust_result() -> dict:
    return {
        "values": [],
        "shape": (1, 0),
        "missing_mask": [],
        "samples": {
            "fid": [],
            "iid": [],
            "father": [],
            "mother": [],
            "sex": [],
            "phenotype": [],
            "source_sample_index": [],
            "haplotype_index": [],
        },
        "variants": {
            "chrom": [],
            "pos": [],
            "id": [],
            "a0": [],
            "a1": [],
            "ref_allele": [],
            "alt_allele": [],
            "source_a0": [],
            "source_a1": [],
            "flipped": [],
            "qual": [],
            "af": [],
            "maf": [],
            "mac": [],
            "missing_rate": [],
            "n_called": [],
        },
        "diagnostics": {},
    }


def test_blocks_honor_size_and_concatenate_to_full_dense_read(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    full = dataset.read()
    blocks = list(dataset.blocks(size=2))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]
    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


def test_blocks_variant_metadata_aligns_with_each_block(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    blocks = list(dataset.blocks(size=2, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs3", "rs4"], ["rs5"]]
    for G_block, variants in blocks:
        assert G_block.shape[1] == len(variants)


def test_blocks_return_samples_keeps_sample_order_constant(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    blocks = list(dataset.blocks(size=3, samples=["S3", "S1"], return_samples=True, return_variants=True))

    assert len(blocks) == 2
    for G_block, samples, variants in blocks:
        assert G_block.shape[0] == 2
        assert G_block.shape[1] == len(variants)
        assert samples["iid"].to_list() == ["S1", "S3"]


def test_blocks_apply_filters_and_sample_keep_lists_like_full_reads(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    read_options = {"variants": genoio.chrom("1"), "samples": ["S3", "S1"]}

    full, full_variants = dataset.read(**read_options, return_variants=True)
    blocks = list(dataset.blocks(size=2, **read_options, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs5"]]
    np.testing.assert_array_equal(np.concatenate([block for block, _ in blocks], axis=1), full)
    assert full_variants["id"].to_list() == ["rs1", "rs2", "rs5"]


def test_plink2_blocks_honor_size_and_concatenate_to_full_dense_read(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    full = dataset.read()
    blocks = list(dataset.blocks(size=2))

    assert [block.shape for block in blocks] == [(3, 2), (3, 1)]
    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


def test_plink2_blocks_variant_metadata_aligns_with_each_block(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    blocks = list(dataset.blocks(size=2, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs3"]]
    for G_block, variants in blocks:
        assert G_block.shape[1] == len(variants)


def test_plink2_blocks_return_samples_keeps_source_order_for_each_block(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    blocks = list(dataset.blocks(size=1, samples=["S3", "S1"], return_samples=True, return_variants=True))

    assert len(blocks) == 3
    for G_block, samples, variants in blocks:
        assert G_block.shape[0] == 2
        assert G_block.shape[1] == len(variants)
        assert samples["iid"].to_list() == ["S1", "S3"]


def test_plink2_matrix_only_blocks_pass_private_matrix_only_option(tmp_path, monkeypatch):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))
    calls = []

    def fake_read_dense_from_rust(self, members, options):
        calls.append(dict(options))
        return empty_dense_rust_result()

    monkeypatch.setattr(genoio.Dataset, "_read_dense_from_rust", fake_read_dense_from_rust)

    list(dataset.blocks(size=2))

    assert calls
    assert calls[0]["matrix_only"] is True


@pytest.mark.parametrize(
    "read_options",
    [
        {"return_samples": True},
        {"return_variants": True},
        {"samples": ["S1"]},
        {"variants": ["rs1"]},
    ],
)
def test_plink2_blocks_disable_matrix_only_when_metadata_or_filters_are_needed(tmp_path, monkeypatch, read_options):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))
    calls = []

    def fake_read_dense_from_rust(self, members, options):
        calls.append(dict(options))
        return empty_dense_rust_result()

    monkeypatch.setattr(genoio.Dataset, "_read_dense_from_rust", fake_read_dense_from_rust)

    list(dataset.blocks(size=2, **read_options))

    assert calls
    assert calls[0]["matrix_only"] is False


@pytest.mark.parametrize(
    ("read_options", "bad_member", "match"),
    [
        ({"return_samples": True}, ".psam", "too few fields"),
        ({"return_variants": True}, ".pvar", "invalid position"),
        ({"samples": ["S1"]}, ".psam", "too few fields"),
        ({"variants": "chrom"}, ".pvar", "invalid position"),
    ],
)
def test_plink2_blocks_metadata_required_paths_reject_malformed_companion_files(
    tmp_path, read_options, bad_member, match
):
    import genoio

    prefix = write_fixed_width_plink2(tmp_path)
    if bad_member == ".psam":
        prefix.with_suffix(".psam").write_text("#FID IID\nF1\n")
    else:
        prefix.with_suffix(".pvar").write_text("#CHROM POS ID REF ALT\n1 bad rs1 A G\n")
    if read_options.get("variants") == "chrom":
        read_options = {**read_options, "variants": genoio.chrom("1")}

    dataset = genoio.pfile(prefix)

    with pytest.raises(genoio.InvalidSourceError, match=match):
        list(dataset.blocks(size=1, **read_options))


@pytest.mark.parametrize(
    ("pvar_text", "match"),
    [
        (
            """\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 bad rs2 C T
2 30 rs3 G A
""",
            "invalid position",
        ),
        (
            """\
#CHROM POS ID REF ALT
1 10 rs1 A G
""",
            "pvar variant count 1",
        ),
    ],
)
def test_plink2_metadata_blocks_validate_full_pvar_before_first_block_return(tmp_path, pvar_text, match):
    import genoio

    prefix = write_fixed_width_plink2(tmp_path)
    prefix.with_suffix(".pvar").write_text(pvar_text)
    dataset = genoio.pfile(prefix)

    with pytest.raises(genoio.InvalidSourceError, match=match):
        next(dataset.blocks(size=1, return_variants=True))


def test_plink2_matrix_only_blocks_reject_bad_variable_width_block_offset(tmp_path):
    import genoio

    dataset = genoio.pfile(write_bad_variable_width_block_offset_plink2(tmp_path))

    with pytest.raises(genoio.InvalidSourceError, match="block offset|header length"):
        next(dataset.blocks(size=1))


def test_plink2_metadata_blocks_reject_bad_variable_width_block_offset(tmp_path):
    import genoio

    dataset = genoio.pfile(write_bad_variable_width_block_offset_plink2(tmp_path))

    with pytest.raises(genoio.InvalidSourceError, match="block offset|header length"):
        next(dataset.blocks(size=1, return_variants=True))


def test_blocks_validate_size(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="positive integer"):
        list(dataset.blocks(size=0))

    with pytest.raises(genoio.InvalidOptionError, match="unsupported sparse option"):
        list(dataset.blocks(size=2, sparse=[]))


def test_blocks_request_bounded_variant_windows_at_rust_boundary(tmp_path, monkeypatch):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    calls = []
    original = genoio.Dataset._read_validated

    def recording_read_validated(self, *, variant_window, **kwargs):
        calls.append(dict(variant_window))
        return original(self, variant_window=variant_window, **kwargs)

    monkeypatch.setattr(genoio.Dataset, "_read_validated", recording_read_validated)

    list(dataset.blocks(size=2))

    assert calls
    assert all(call["len"] <= 2 for call in calls)
    assert [call["start"] for call in calls] == [0, 2, 4]


def test_blocks_do_not_call_public_read_internally(tmp_path, monkeypatch):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    def fail_read(*args, **kwargs):
        raise AssertionError("blocks must use the bounded Rust call boundary, not public read()")

    monkeypatch.setattr(genoio.Dataset, "read", fail_read)

    blocks = list(dataset.blocks(size=2))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]


def test_sparse_blocks_work_with_default_missing_policy(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    blocks = list(dataset.blocks(size=2, sparse=True))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]
    assert all(scipy_sparse.isspmatrix_csc(block) for block in blocks)
