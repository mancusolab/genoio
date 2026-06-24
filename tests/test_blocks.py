# pattern: Imperative Shell

from pathlib import Path

import numpy as np
import pytest
from fixture_writers import write_bgen_dosage, write_fixed_width_plink2
from scipy import sparse as scipy_sparse


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


def test_iter_blocks_replaces_blocks_in_public_api(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    assert hasattr(dataset, "iter_blocks")
    assert not hasattr(dataset, "blocks")


def test_iter_blocks_honor_size_and_concatenate_to_full_dense_read(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    full = dataset.read()
    blocks = list(dataset.iter_blocks(size=2))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]
    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


def test_iter_regions_yields_region_and_read_result_for_each_region(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    regions = [genoio.region("1:1-25"), genoio.region("2:1-35")]

    region_reads = list(dataset.iter_regions(regions, return_variants=True))

    assert [region for region, _ in region_reads] == regions
    assert [variants["id"].to_list() for _, (_, variants) in region_reads] == [["rs1", "rs2"], ["rs3"]]


def test_iter_regions_rejects_variants_read_option(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="variants"):
        list(dataset.iter_regions([genoio.region("1:1-25")], variants=genoio.chrom("1")))


def test_blocks_variant_metadata_aligns_with_each_block(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    blocks = list(dataset.iter_blocks(size=2, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs3", "rs4"], ["rs5"]]
    for G_block, variants in blocks:
        assert G_block.shape[1] == len(variants)


def test_blocks_return_samples_keeps_sample_order_constant(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    blocks = list(dataset.iter_blocks(size=3, samples=["S3", "S1"], return_samples=True, return_variants=True))

    assert len(blocks) == 2
    for G_block, samples, variants in blocks:
        assert G_block.shape[0] == 2
        assert G_block.shape[1] == len(variants)
        assert samples["iid"].to_list() == ["S1", "S3"]


def test_blocks_apply_filters_and_sample_keep_lists_like_full_reads(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    read_options = {"variants": genoio.chrom("1"), "samples": ["S3", "S1"]}

    full, full_variants = dataset.read(variants=genoio.chrom("1"), samples=["S3", "S1"], return_variants=True)
    blocks = list(dataset.iter_blocks(size=2, **read_options, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs5"]]
    np.testing.assert_array_equal(np.concatenate([block for block, _ in blocks], axis=1), full)
    assert full_variants["id"].to_list() == ["rs1", "rs2", "rs5"]


def test_bgen_dosage_blocks_yield_no_blocks_for_empty_variant_filter(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    blocks = list(dataset.iter_blocks(size=1, dosage="dosage", variants=[], return_variants=True))

    assert blocks == []


def test_bgen_dosage_blocks_honor_size_and_concatenate_to_full_read(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    full = dataset.read(dosage="dosage")
    blocks = list(dataset.iter_blocks(1, dosage="dosage"))

    assert [block.shape for block in blocks] == [(2, 1), (2, 1)]
    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


def test_bgen_dosage_blocks_variant_metadata_aligns_with_each_block(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    blocks = list(dataset.iter_blocks(1, dosage="dosage", return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1"], ["rs2"]]
    for G_block, variants in blocks:
        assert G_block.shape == (2, len(variants))


def test_bgen_dosage_filtered_blocks_match_filtered_full_read(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))
    read_options = {"dosage": "dosage", "variants": genoio.chrom("2")}

    full, full_variants = dataset.read(dosage="dosage", variants=genoio.chrom("2"), return_variants=True)
    blocks = list(dataset.iter_blocks(1, **read_options, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs2"]]
    np.testing.assert_array_equal(np.concatenate([block for block, _ in blocks], axis=1), full)
    assert full_variants["id"].to_list() == ["rs2"]


def test_bgen_dosage_blocks_yield_no_blocks_for_nonmatching_metadata_filter(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    blocks = list(dataset.iter_blocks(size=1, dosage="dosage", variants=genoio.chrom("9")))

    assert blocks == []


def test_plink2_blocks_honor_size_and_concatenate_to_full_dense_read(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    full = dataset.read()
    blocks = list(dataset.iter_blocks(size=2))

    assert [block.shape for block in blocks] == [(3, 2), (3, 1)]
    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


def test_plink2_blocks_variant_metadata_aligns_with_each_block(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    blocks = list(dataset.iter_blocks(size=2, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs3"]]
    for G_block, variants in blocks:
        assert G_block.shape[1] == len(variants)


def test_plink2_blocks_return_samples_keeps_source_order_for_each_block(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    blocks = list(dataset.iter_blocks(size=1, samples=["S3", "S1"], return_samples=True, return_variants=True))

    assert len(blocks) == 3
    for G_block, samples, variants in blocks:
        assert G_block.shape[0] == 2
        assert G_block.shape[1] == len(variants)
        assert samples["iid"].to_list() == ["S1", "S3"]


@pytest.mark.parametrize(
    ("read_options", "expected_matrix_only"),
    [
        pytest.param({}, True, id="matrix-only-fast-path"),
        pytest.param({"return_samples": True}, False, id="sample-metadata"),
        pytest.param({"return_variants": True}, False, id="variant-metadata"),
        pytest.param({"samples": ["S1"]}, True, id="sample-filter"),
        pytest.param({"variants": ["rs1"]}, True, id="variant-filter"),
    ],
)
def test_plink2_blocks_set_matrix_only_by_metadata_needs(tmp_path, monkeypatch, read_options, expected_matrix_only):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))
    calls = []

    def fake_read_from_rust(self, kind, sparse, members, options):
        assert kind == "geno"
        assert sparse is False
        calls.append(dict(options))
        return empty_dense_rust_result()

    monkeypatch.setattr(genoio.Dataset, "_read_from_rust", fake_read_from_rust)

    list(dataset.iter_blocks(size=2, **read_options))

    assert calls
    assert calls[0]["matrix_only"] is expected_matrix_only


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
        list(dataset.iter_blocks(size=1, **read_options))


@pytest.mark.parametrize(
    "pvar_text",
    [
        pytest.param(
            """\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 bad rs2 C T
2 30 rs3 G A
""",
            id="malformed-later-row",
        ),
        pytest.param(
            """\
#CHROM POS ID REF ALT
1 10 rs1 A G
""",
            id="missing-later-row",
        ),
    ],
)
def test_plink2_metadata_blocks_skip_later_pvar_records_before_first_block_return(tmp_path, pvar_text):
    import genoio

    prefix = write_fixed_width_plink2(tmp_path)
    prefix.with_suffix(".pvar").write_text(pvar_text)
    dataset = genoio.pfile(prefix)

    _, variants = next(dataset.iter_blocks(size=1, return_variants=True))

    assert variants["id"].to_list() == ["rs1"]


def test_plink2_metadata_blocks_validate_requested_pvar_window(tmp_path):
    import genoio

    prefix = write_fixed_width_plink2(tmp_path)
    prefix.with_suffix(".pvar").write_text(
        """\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 bad rs2 C T
2 30 rs3 G A
"""
    )
    dataset = genoio.pfile(prefix)

    iterator = dataset.iter_blocks(size=1, return_variants=True)
    next(iterator)
    with pytest.raises(genoio.InvalidSourceError, match="invalid position"):
        next(iterator)


@pytest.mark.parametrize(
    "read_options",
    [
        pytest.param({}, id="matrix-only"),
        pytest.param({"return_variants": True}, id="metadata"),
    ],
)
def test_plink2_blocks_reject_bad_variable_width_block_offset(tmp_path, read_options):
    import genoio

    dataset = genoio.pfile(write_bad_variable_width_block_offset_plink2(tmp_path))

    with pytest.raises(genoio.InvalidSourceError, match="block offset|header length"):
        next(dataset.iter_blocks(size=1, **read_options))


def test_blocks_validate_size(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="positive integer"):
        list(dataset.iter_blocks(size=0))

    with pytest.raises(genoio.InvalidOptionError, match="unsupported sparse option"):
        list(dataset.iter_blocks(size=2, sparse=[]))


def test_blocks_request_bounded_variant_windows_at_rust_boundary(tmp_path, monkeypatch):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    calls = []
    original = genoio.Dataset._read_payload

    def recording_read_payload(self, *, variant_window, **kwargs):
        calls.append(dict(variant_window))
        return original(self, variant_window=variant_window, **kwargs)

    monkeypatch.setattr(genoio.Dataset, "_read_payload", recording_read_payload)

    list(dataset.iter_blocks(size=2))

    assert calls
    assert all(call["len"] <= 2 for call in calls)
    assert [call["start"] for call in calls] == [0, 2, 4]


def test_blocks_do_not_call_public_read_internally(tmp_path, monkeypatch):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    def fail_read(*args, **kwargs):
        raise AssertionError("iter_blocks must use the bounded Rust call boundary, not public read()")

    monkeypatch.setattr(genoio.Dataset, "read", fail_read)

    blocks = list(dataset.iter_blocks(size=2))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]


def test_sparse_blocks_work_with_default_missing_policy(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    blocks = list(dataset.iter_blocks(size=2, sparse=True))

    assert [block.shape for block in blocks] == [(3, 2), (3, 2), (3, 1)]
    assert all(scipy_sparse.isspmatrix_csc(block) for block in blocks)
