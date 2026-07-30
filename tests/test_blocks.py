# pattern: Imperative Shell

import sqlite3
from pathlib import Path
from typing import Any, cast

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


def write_bgen_index(path: Path) -> None:
    contents = path.read_bytes()
    starts = [
        contents.index(len(variant_id).to_bytes(2, "little") + variant_id.encode()) for variant_id in ("var1", "var2")
    ]
    sizes = [starts[1] - starts[0], len(contents) - starts[1]]
    rows = [
        ("1", 10, "rs1", "A", "G", starts[0], sizes[0]),
        ("2", 20, "rs2", "C", "T", starts[1], sizes[1]),
    ]
    with sqlite3.connect(f"{path}.bgi") as connection:
        connection.execute(
            """\
CREATE TABLE Variant (
    chromosome TEXT NOT NULL,
    position INT NOT NULL,
    rsid TEXT NOT NULL,
    number_of_alleles INT NOT NULL,
    allele1 TEXT NOT NULL,
    allele2 TEXT NULL,
    file_start_position INT NOT NULL,
    size_in_bytes INT NOT NULL,
    PRIMARY KEY (
        chromosome, position, rsid, allele1, allele2, file_start_position
    )
)"""
        )
        connection.executemany(
            """\
INSERT INTO Variant (
    chromosome, position, rsid, number_of_alleles, allele1, allele2,
    file_start_position, size_in_bytes
) VALUES (?, ?, ?, 2, ?, ?, ?, ?)""",
            rows,
        )


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


def test_pbr_py_excluded_001_read_remains_on_stateless_native_route(
    tmp_path,
    monkeypatch,
):
    import genoio
    from genoio import _rust

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    native_read_dense = _rust.read_dense
    stateless_calls = []

    def fail_if_block_reader_is_constructed(*args, **kwargs):
        raise AssertionError("Dataset.read() must not construct _BlockReader")

    def recording_read_dense(*args, **kwargs):
        stateless_calls.append((args, kwargs))
        return native_read_dense(*args, **kwargs)

    monkeypatch.setattr(_rust, "_BlockReader", fail_if_block_reader_is_constructed)
    monkeypatch.setattr(_rust, "read_dense", recording_read_dense)

    observed = dataset.read()

    np.testing.assert_array_equal(
        observed,
        np.array(
            [
                [0.0, 1.0, 2.0, 0.0, 2.0],
                [1.0, 0.0, 1.0, 0.0, 2.0],
                [2.0, 2.0, 0.0, 1.0, 1.0],
            ],
            dtype=np.float32,
        ),
    )
    assert len(stateless_calls) == 1


def test_pbr_py_excluded_001_iter_regions_retains_one_read_per_region(
    tmp_path,
    monkeypatch,
):
    import genoio
    from genoio import _rust

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    regions = [genoio.region("1:1-25"), genoio.region("2:1-35")]
    original_read = genoio.Dataset.read
    read_calls = []

    def fail_if_block_reader_is_constructed(*args, **kwargs):
        raise AssertionError("Dataset.iter_regions() must not construct _BlockReader")

    def recording_read(self, *args, **kwargs):
        read_calls.append((args, kwargs))
        return original_read(self, *args, **kwargs)

    monkeypatch.setattr(_rust, "_BlockReader", fail_if_block_reader_is_constructed)
    monkeypatch.setattr(genoio.Dataset, "read", recording_read)

    region_reads = list(dataset.iter_regions(regions, return_variants=True))

    assert [region for region, _ in region_reads] == regions
    assert [variants["id"].to_list() for _, (_, variants) in region_reads] == [["rs1", "rs2"], ["rs3"]]
    assert len(read_calls) == len(regions)
    assert all(args == () for args, _ in read_calls)
    assert all(options["variants"] is region for (_, options), region in zip(read_calls, regions, strict=True))


def test_pbr_py_excluded_001_iter_regions_empty_and_validation_remain_lazy(
    tmp_path,
    monkeypatch,
):
    import genoio
    from genoio import _rust

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    def fail_if_block_reader_is_constructed(*args, **kwargs):
        raise AssertionError("Dataset.iter_regions() must not construct _BlockReader")

    monkeypatch.setattr(_rust, "_BlockReader", fail_if_block_reader_is_constructed)

    empty = dataset.iter_regions([], sparse=[])
    invalid = dataset.iter_regions([genoio.region("1:1-25")], sparse=[])

    assert list(empty) == []
    with pytest.raises(genoio.InvalidOptionError, match="unsupported sparse option"):
        next(invalid)


def test_pbr_py_meta_001_variant_metadata_aligns_with_each_block(tmp_path):
    import genoio

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    blocks = list(dataset.iter_blocks(size=2, return_variants=True))

    assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "rs2"], ["rs3", "rs4"], ["rs5"]]
    for G_block, variants in blocks:
        assert G_block.shape[1] == len(variants)


def test_pbr_py_meta_001_sample_metadata_is_present_on_every_block(tmp_path):
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


@pytest.mark.parametrize(
    ("kind", "phased"),
    [
        pytest.param("geno", False, id="genotype-dosage"),
        pytest.param("haplo", True, id="haplotype-dosage"),
    ],
)
@pytest.mark.parametrize("indexed", [False, True], ids=["sequential", "indexed"])
def test_pbr_py_matrix_001_pbr_py_meta_001_bgen_dosage_blocks_match_sequential_oracle(
    tmp_path,
    kind,
    phased,
    indexed,
):
    import genoio

    path = write_bgen_dosage(tmp_path, phased=phased)
    read_options = {
        "kind": kind,
        "dosage": "dosage",
        "variants": genoio.region("1:1-15"),
        "dtype": "float64",
    }
    oracle, oracle_samples, oracle_variants = cast(Any, genoio.bgen(path).read)(
        **read_options,
        return_samples=True,
        return_variants=True,
    )
    if indexed:
        write_bgen_index(path)
    blocks = list(
        genoio.bgen(path).iter_blocks(
            size=1,
            **read_options,
            return_samples=True,
            return_variants=True,
        )
    )

    assert len(blocks) == 1
    block, samples, variants = blocks[0]
    np.testing.assert_array_equal(block, oracle)
    assert block.dtype == np.dtype("float64")
    assert samples.schema == oracle_samples.schema
    assert samples.equals(oracle_samples)
    assert variants.schema == oracle_variants.schema
    assert variants.equals(oracle_variants)


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
    from genoio import _rust

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))
    calls = []

    class RecordingReader:
        def __init__(self, format, members, kind, sparse, options, block_size):
            assert format == "plink2"
            assert members
            assert kind == "geno"
            assert sparse is False
            assert block_size == 2
            calls.append(dict(options))

        def next_block(self):
            return None

        def close(self):
            return None

    monkeypatch.setattr(_rust, "_BlockReader", RecordingReader)

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


def test_pbr_py_cutover_001_validates_eagerly_without_constructing_reader(
    tmp_path,
    monkeypatch,
):
    import genoio
    from genoio import _rust

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))

    def fail_if_constructed(*args, **kwargs):
        raise AssertionError("_BlockReader must not be constructed during eager validation")

    monkeypatch.setattr(_rust, "_BlockReader", fail_if_constructed)
    with pytest.raises(genoio.InvalidOptionError, match="positive integer"):
        dataset.iter_blocks(size=0)

    with pytest.raises(genoio.InvalidOptionError, match="unsupported sparse option"):
        dataset.iter_blocks(size=2, sparse=[])

    with pytest.raises(genoio.UnsupportedRepresentation):
        dataset.iter_blocks(size=2, dosage="dosage", sparse=True)


def test_pbr_py_cutover_001_constructs_one_lazy_session_and_avoids_stateless_reads(
    tmp_path,
    monkeypatch,
):
    import genoio
    from genoio import _rust

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    constructed = []

    class RecordingReader:
        def __init__(self, format, members, kind, sparse, options, block_size):
            constructed.append(self)
            self._reader = _rust_reader(format, members, kind, sparse, options, block_size)
            self.next_calls = 0
            self.close_calls = 0

        def next_block(self):
            self.next_calls += 1
            return self._reader.next_block()

        def close(self):
            self.close_calls += 1
            return self._reader.close()

    _rust_reader = _rust._BlockReader

    def fail_stateless(*args, **kwargs):
        raise AssertionError("iter_blocks must not call a stateless read path")

    monkeypatch.setattr(_rust, "_BlockReader", RecordingReader)
    monkeypatch.setattr(genoio.Dataset, "_read_payload", fail_stateless)
    monkeypatch.setattr(_rust, "read_dense", fail_stateless)
    monkeypatch.setattr(_rust, "read_sparse", fail_stateless)
    monkeypatch.setattr(_rust, "read_haplotypes_dense", fail_stateless)
    monkeypatch.setattr(_rust, "read_haplotypes_sparse", fail_stateless)

    iterator = dataset.iter_blocks(size=2)
    assert constructed == []

    first = next(iterator)
    assert first.shape == (3, 2)
    assert len(constructed) == 1

    remaining = list(iterator)

    assert [block.shape for block in remaining] == [(3, 2), (3, 1)]
    assert len(constructed) == 1
    assert constructed[0].next_calls == 4
    assert constructed[0].close_calls == 1


def test_pbr_py_iterator_001_interleaved_iterators_own_independent_sessions(
    tmp_path,
    monkeypatch,
):
    import genoio
    from genoio import _rust

    dataset = genoio.vcf(write_blocks_vcf(tmp_path))
    native_reader = _rust._BlockReader
    readers = []

    def recording_reader(*args, **kwargs):
        reader = native_reader(*args, **kwargs)
        readers.append(reader)
        return reader

    monkeypatch.setattr(_rust, "_BlockReader", recording_reader)
    left = dataset.iter_blocks(size=2)
    right = dataset.iter_blocks(size=2)

    left_first = next(left)
    right_first = next(right)
    left_second = next(left)
    cast(Any, left).close()
    right_second = next(right)
    right_tail = list(right)

    assert len(readers) == 2
    assert readers[0] is not readers[1]
    np.testing.assert_array_equal(left_first, right_first)
    np.testing.assert_array_equal(left_second, right_second)
    assert [block.shape for block in right_tail] == [(3, 1)]


def test_pbr_py_cutover_001_source_open_is_deferred_until_first_next(tmp_path):
    import genoio

    path = write_blocks_vcf(tmp_path)
    dataset = genoio.vcf(path)
    iterator = dataset.iter_blocks(size=2)
    path.unlink()

    with pytest.raises(genoio.InvalidSourceError):
        next(iterator)


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
