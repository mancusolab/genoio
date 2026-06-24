# pattern: Imperative Shell

from pathlib import Path
from typing import Any, cast

import numpy as np
import polars as pl
import pytest
from fixture_writers import (
    _bgen_sample_identifier_block,
    write_bgen_dosage,
    write_fixed_width_plink2,
    write_fixed_width_plink2_dosage,
)

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


def write_ds_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "dosage.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Expected alternate allele dosage">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DS\t0/0:0.2\t0/1:1.4\t1/1:1.8
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT:DS\t0/0:0\t0/0:.\t0/1:0.7
"""
    )
    return path


def write_wrong_magic_bgen(tmp_path: Path) -> Path:
    path = write_bgen_dosage(tmp_path)
    contents = bytearray(path.read_bytes())
    contents[16:20] = b"nope"
    path.write_bytes(contents)
    return path


def write_truncated_variant_bgen(tmp_path: Path) -> Path:
    path = tmp_path / "truncated_variant.bgen"
    contents = bytearray()
    flags = (2 << 2) | (1 << 31)
    contents.extend((20).to_bytes(4, "little"))
    contents.extend((20).to_bytes(4, "little"))
    contents.extend((1).to_bytes(4, "little"))
    contents.extend((2).to_bytes(4, "little"))
    contents.extend(b"bgen")
    contents.extend(flags.to_bytes(4, "little"))
    contents.extend(_bgen_sample_identifier_block(["sample_1", "sample_2"]))
    variant_offset = len(contents) - 4
    contents[0:4] = variant_offset.to_bytes(4, "little")
    contents.extend((4).to_bytes(2, "little"))
    contents.extend(b"rs")
    path.write_bytes(contents)
    return path


def write_gt_only_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "gt_only.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1
"""
    )
    return path


def test_dense_vcf_read_returns_sample_by_variant_numpy_array_and_metadata(tmp_path):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    G, samples, variants = dataset.read(return_samples=True, return_variants=True)

    assert isinstance(G, np.ndarray)
    assert G.dtype == np.dtype("float32")
    np.testing.assert_array_equal(G, np.array([[0.0, 1.0], [1.0, np.nan], [2.0, 0.0]], dtype=np.float32))
    assert isinstance(samples, pl.DataFrame)
    assert samples["iid"].to_list() == ["S1", "S2", "S3"]
    assert isinstance(variants, pl.DataFrame)
    assert variants["id"].to_list() == ["rs1", "rs2"]


def test_default_read_equals_explicit_genotype_read(tmp_path):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    np.testing.assert_array_equal(dataset.read(), dataset.read(kind="geno"))


def test_dense_genotype_reads_do_not_minor_allele_flip_by_default(tmp_path):
    import genoio

    dataset = genoio.vcf(write_common_a1_vcf(tmp_path))

    G, variants = dataset.read(return_variants=True)

    np.testing.assert_array_equal(G, np.array([[2.0], [2.0], [1.0]], dtype=np.float32))
    assert variants["a0"].to_list() == ["A"]
    assert variants["a1"].to_list() == ["G"]


def test_dense_vcf_dosage_reads_ds_values_without_gt_fallback(tmp_path):
    import genoio

    dataset = genoio.vcf(write_ds_vcf(tmp_path))

    G, variants = dataset.read(dosage="dosage", missing="nan", return_variants=True)

    np.testing.assert_array_equal(
        G,
        np.array([[0.2, 0.0], [1.4, np.nan], [1.8, 0.7]], dtype=np.float32),
    )
    assert variants["id"].to_list() == ["rs1", "rs2"]
    assert variants["a1"].to_list() == ["G", "T"]


def test_dense_bgen_dosage_read_returns_sample_by_variant_matrix(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    G = dataset.read(dosage="dosage")

    assert G.shape == (2, 2)
    assert G.dtype == np.dtype("float32")
    np.testing.assert_allclose(
        G,
        np.array(
            [
                [0.29803923, 1.0],
                [1.0980392, 0.8],
            ],
            dtype=np.float32,
        ),
        rtol=0,
        atol=2.0 / 255.0,
    )


def test_dense_bgen_phased_dosage_read_returns_collapsed_expected_a1_count(tmp_path):
    import genoio

    dataset = genoio.bgen(
        write_bgen_dosage(
            tmp_path,
            phased=True,
            variant_calls=[
                [(255, 0), (128, 64)],
                [(0, 255), None],
            ],
        )
    )

    G = dataset.read(dosage="dosage")

    assert G.shape == (2, 2)
    np.testing.assert_allclose(
        G,
        np.array(
            [
                [1.0, 1.0],
                [1.2470589, np.nan],
            ],
            dtype=np.float32,
        ),
        rtol=0,
        atol=2.0 / 255.0,
    )


def test_dense_bgen_haplotype_dosage_read_returns_haplotype_rows(tmp_path):
    import genoio

    dataset = genoio.bgen(
        write_bgen_dosage(
            tmp_path,
            phased=True,
            variant_calls=[
                [(255, 0), (128, 64)],
                [(0, 255), None],
            ],
        )
    )

    H, samples = dataset.read(kind="haplo", dosage="dosage", return_samples=True)

    assert H.shape == (4, 2)
    np.testing.assert_allclose(
        H,
        np.array(
            [
                [0.0, 1.0],
                [1.0, 0.0],
                [0.49803922, np.nan],
                [0.7490196, np.nan],
            ],
            dtype=np.float32,
        ),
        rtol=0,
        atol=1.0 / 255.0,
    )
    assert samples["iid"].to_list() == ["sample_1", "sample_1", "sample_2", "sample_2"]
    assert samples["source_sample_index"].to_list() == [0, 0, 1, 1]
    assert samples["haplotype_index"].to_list() == [0, 1, 0, 1]


def test_empty_bgen_haplotype_dosage_read_preserves_haplotype_sample_columns(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, phased=True))

    H, samples = dataset.read(kind="haplo", dosage="dosage", samples=[], return_samples=True)

    assert H.shape == (0, 2)
    assert samples.columns == [
        "fid",
        "iid",
        "father",
        "mother",
        "sex",
        "phenotype",
        "source_sample_index",
        "haplotype_index",
    ]
    assert samples.to_dict(as_series=False) == {
        "fid": [],
        "iid": [],
        "father": [],
        "mother": [],
        "sex": [],
        "phenotype": [],
        "source_sample_index": [],
        "haplotype_index": [],
    }


def test_empty_bgen_genotype_dosage_read_omits_haplotype_sample_columns(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    G, samples = dataset.read(dosage="dosage", samples=[], return_samples=True)

    assert G.shape == (0, 2)
    assert samples.columns == ["fid", "iid", "father", "mother", "sex", "phenotype"]


def test_dense_bgen_dosage_default_missing_policy_returns_nan(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, missing=True))

    G = dataset.read(dosage="dosage")

    assert G.shape == (2, 2)
    np.testing.assert_allclose(G[0], np.array([0.29803923, 1.0], dtype=np.float32), rtol=0, atol=2.0 / 255.0)
    assert np.isnan(G[1, 0])
    assert np.isclose(G[1, 1], 0.8, atol=2.0 / 255.0)


def test_dense_bgen_dosage_accepts_missing_samples_with_packed_zero_probabilities(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, missing=True, pack_missing_probabilities=True))

    G = dataset.read(dosage="dosage")

    assert G.shape == (2, 2)
    np.testing.assert_allclose(
        G[0],
        np.array([0.29803923, 1.0], dtype=np.float32),
        rtol=0,
        atol=2.0 / 255.0,
    )
    assert np.isnan(G[1, 0])
    assert np.isclose(G[1, 1], 0.8, atol=2.0 / 255.0)


def test_dense_bgen_phased_dosage_accepts_missing_samples_with_packed_zero_probabilities(tmp_path):
    import genoio

    dataset = genoio.bgen(
        write_bgen_dosage(
            tmp_path,
            phased=True,
            pack_missing_probabilities=True,
            variant_calls=[
                [(255, 0), None],
                [(0, 255), (128, 64)],
            ],
        )
    )

    G = dataset.read(dosage="dosage")

    assert G.shape == (2, 2)
    np.testing.assert_allclose(
        G,
        np.array(
            [
                [1.0, 1.0],
                [np.nan, 1.2470589],
            ],
            dtype=np.float32,
        ),
        rtol=0,
        atol=2.0 / 255.0,
    )


def test_dense_bgen_dosage_missing_raise_rejects_missing_calls(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, missing=True))

    with pytest.raises(genoio.MissingDataError, match="missing genotype"):
        dataset.read(dosage="dosage", missing="raise")


def test_dense_bgen_haplotype_dosage_missing_raise_rejects_missing_calls(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path, missing=True, phased=True))

    with pytest.raises(genoio.MissingDataError, match="missing genotype"):
        dataset.read(kind="haplo", dosage="dosage", missing="raise")


@pytest.mark.parametrize(
    ("writer", "match"),
    [
        (write_wrong_magic_bgen, "magic"),
        (write_truncated_variant_bgen, "failed to fill whole buffer"),
    ],
)
def test_dense_bgen_dosage_rejects_malformed_files_as_invalid_source(tmp_path, writer, match):
    import genoio

    dataset = genoio.bgen(writer(tmp_path))

    with pytest.raises(genoio.InvalidSourceError, match=match):
        dataset.read(dosage="dosage")


def test_dense_bgen_dosage_rejects_valid_unsupported_multiallelic_source(tmp_path):
    import genoio

    dataset = genoio.bgen(
        write_bgen_dosage(
            tmp_path,
            variants=[("var1", "rs1", "1", 10, ["A", "C", "G"])],
            variant_calls=[[(204, 26), (51, 128)]],
        )
    )

    with pytest.raises(genoio.UnsupportedRepresentation, match="biallelic"):
        dataset.read(dosage="dosage")


def test_dense_bgen_hardcall_error_points_to_supported_dosage_option(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match='dosage="dosage"'):
        dataset.read()


def test_dense_bgen_dosage_sample_filter_uses_source_order(tmp_path):
    import genoio

    dataset = genoio.bgen(
        write_bgen_dosage(
            tmp_path,
            sample_ids=["sample_1", "sample_2", "sample_3"],
            variant_calls=[
                [(204, 26), (51, 128), (0, 0)],
                [(0, 255), (102, 102), (255, 0)],
            ],
        )
    )

    G, samples = dataset.read(dosage="dosage", samples=["sample_3", "sample_1"], return_samples=True)

    np.testing.assert_allclose(
        G,
        np.array(
            [
                [0.29803923, 1.0],
                [2.0, 0.0],
            ],
            dtype=np.float32,
        ),
        rtol=0,
        atol=2.0 / 255.0,
    )
    assert samples["iid"].to_list() == ["sample_1", "sample_3"]


def test_dense_vcf_dosage_sample_filter_uses_selected_samples(tmp_path):
    import genoio

    dataset = genoio.vcf(write_ds_vcf(tmp_path))

    G, samples = dataset.read(dosage="dosage", missing="nan", samples=["S3", "S1"], return_samples=True)

    np.testing.assert_array_equal(G, np.array([[0.2, 0.0], [1.8, 0.7]], dtype=np.float32))
    assert samples["iid"].to_list() == ["S1", "S3"]


def test_dense_vcf_dosage_accepts_metadata_only_variant_filters(tmp_path):
    import genoio

    dataset = genoio.vcf(write_ds_vcf(tmp_path))

    G, variants = dataset.read(
        dosage="dosage",
        missing="nan",
        variants=genoio.id_in(["rs2"]),
        return_variants=True,
    )

    np.testing.assert_array_equal(G, np.array([[0.0], [np.nan], [0.7]], dtype=np.float32))
    assert variants["id"].to_list() == ["rs2"]


def test_dense_vcf_dosage_missing_policy_raise_rejects_missing_ds(tmp_path):
    import genoio

    dataset = genoio.vcf(write_ds_vcf(tmp_path))

    with pytest.raises(genoio.MissingDataError, match="missing genotype"):
        dataset.read(dosage="dosage", missing="raise")


def test_dense_vcf_dosage_requires_ds_field_without_gt_fallback(tmp_path):
    import genoio

    dataset = genoio.vcf(write_gt_only_vcf(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="FORMAT/DS"):
        dataset.read(dosage="dosage")


def test_dense_vcf_dosage_blocks_match_full_dosage_read(tmp_path):
    import genoio

    dataset = genoio.vcf(write_ds_vcf(tmp_path))

    full = dataset.read(dosage="dosage", missing="nan")
    blocks = list(dataset.iter_blocks(1, dosage="dosage", missing="nan"))

    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


def test_vcf_dosage_genotype_stat_filters_use_fractional_mac(tmp_path):
    import genoio

    dataset = genoio.vcf(write_ds_vcf(tmp_path))

    G, variants = dataset.read(
        dosage="dosage",
        missing="nan",
        variants=genoio.mac(max=2),
        return_variants=True,
    )

    np.testing.assert_array_equal(G, np.array([[0.0], [np.nan], [0.7]], dtype=np.float32))
    assert variants["id"].to_list() == ["rs2"]


def test_return_samples_and_variants_tuple_order_is_matrix_samples_variants(tmp_path):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    result = dataset.read(return_samples=True, return_variants=True)

    assert len(result) == 3
    assert isinstance(result[0], np.ndarray)
    assert isinstance(result[1], pl.DataFrame)
    assert isinstance(result[2], pl.DataFrame)
    assert result[1]["iid"].to_list() == ["S1", "S2", "S3"]
    assert result[2]["id"].to_list() == ["rs1", "rs2"]


def test_dense_plink1_read_matches_fixture_matrix_and_return_tuple_shapes():
    import genoio

    dataset = genoio.bfile(FIXTURE_ROOT / "plink1" / "tiny")

    G = dataset.read()
    G_samples = dataset.read(return_samples=True)
    G_variants = dataset.read(return_variants=True)
    G_both = dataset.read(return_samples=True, return_variants=True)

    expected = np.array([[0.0, np.nan, 2.0], [np.nan, 0.0, 1.0], [2.0, 1.0, 0.0]], dtype=np.float32)
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


def test_dense_plink2_read_matches_fixed_width_hardcall_fixture(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    G, samples, variants = dataset.read(return_samples=True, return_variants=True)

    assert G.dtype == np.dtype("float32")
    np.testing.assert_array_equal(
        G,
        np.array(
            [
                [0.0, 1.0, 2.0],
                [np.nan, 0.0, 1.0],
                [2.0, 1.0, 0.0],
            ],
            dtype=np.float32,
        ),
    )
    assert samples["iid"].to_list() == ["S1", "S2", "S3"]
    assert variants["id"].to_list() == ["rs1", "rs2", "rs3"]
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]


def test_dense_plink2_sample_filter_keeps_source_order_and_values(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    G, samples = dataset.read(samples=["S3", "S1"], return_samples=True)

    np.testing.assert_array_equal(
        G,
        np.array(
            [
                [0.0, 1.0, 2.0],
                [2.0, 1.0, 0.0],
            ],
            dtype=np.float32,
        ),
    )
    assert samples["iid"].to_list() == ["S1", "S3"]


def test_dense_plink2_dosage_requires_pgen_dosage_track(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    with pytest.raises(genoio.UnsupportedRepresentation, match="dosage"):
        dataset.read(dosage="dosage")


def test_dense_plink2_dosage_reads_stored_a1_dosage_values(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2_dosage(tmp_path))

    G, variants = dataset.read(dosage="dosage", missing="nan", return_variants=True)

    np.testing.assert_allclose(
        G,
        np.array([[0.2, 0.0], [1.4, np.nan], [1.8, 0.7]], dtype=np.float32),
        rtol=0,
        atol=2.0 / 32768.0,
    )
    assert variants["id"].to_list() == ["rs1", "rs2"]
    assert variants["a1"].to_list() == ["G", "T"]


def test_dense_plink2_dosage_sample_and_variant_filters_apply_before_decode(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2_dosage(tmp_path))

    G, samples, variants = dataset.read(
        dosage="dosage",
        missing="nan",
        samples=["S3", "S1"],
        variants=genoio.id_in(["rs2"]),
        return_samples=True,
        return_variants=True,
    )

    np.testing.assert_allclose(G, np.array([[0.0], [0.7]], dtype=np.float32), rtol=0, atol=2.0 / 32768.0)
    assert samples["iid"].to_list() == ["S1", "S3"]
    assert variants["id"].to_list() == ["rs2"]


def test_dense_plink2_dosage_file_default_read_uses_hardcalls(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2_dosage(tmp_path))

    G, samples, variants = dataset.read(missing="nan", return_samples=True, return_variants=True)

    np.testing.assert_array_equal(G, np.array([[0.0, 0.0], [1.0, np.nan], [2.0, 0.0]], dtype=np.float32))
    assert samples["iid"].to_list() == ["S1", "S2", "S3"]
    assert variants["id"].to_list() == ["rs1", "rs2"]


def test_dense_plink2_dosage_blocks_match_full_dosage_read(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2_dosage(tmp_path))

    full = dataset.read(dosage="dosage", missing="nan")
    blocks = list(dataset.iter_blocks(1, dosage="dosage", missing="nan"))

    np.testing.assert_array_equal(np.concatenate(blocks, axis=1), full)


@pytest.mark.parametrize(
    ("read_options", "bad_member", "match"),
    [
        ({"return_samples": True}, ".psam", "too few fields"),
        ({"return_variants": True}, ".pvar", "invalid position"),
        ({"samples": ["S1"]}, ".psam", "too few fields"),
        ({"variants": "chrom"}, ".pvar", "invalid position"),
    ],
)
def test_dense_plink2_metadata_required_paths_reject_malformed_companion_files(
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
        cast(Any, dataset.read)(**read_options)


def test_dense_plink2_rejects_psam_header_without_iid(tmp_path):
    import genoio

    prefix = write_fixed_width_plink2(tmp_path)
    prefix.with_suffix(".psam").write_text(
        """\
#FID SEX
F1 1
F2 2
F3 0
"""
    )
    dataset = genoio.pfile(prefix)

    with pytest.raises(genoio.InvalidSourceError, match="psam header missing IID"):
        dataset.read(return_samples=True)


@pytest.mark.parametrize(
    ("psam_text", "expected_iids"),
    [
        (
            """\
#IID
S1
S2
S3
""",
            ["S1", "S2", "S3"],
        ),
        (
            """\
#FID IID PAT MAT SEX PHENO
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
""",
            ["S1", "S2", "S3"],
        ),
        (
            """\
#SEX IID FID
1 S1 F1
2 S2 F1
0 S3 F2
""",
            ["S1", "S2", "S3"],
        ),
    ],
)
def test_dense_plink2_accepts_supported_psam_iid_header_forms(tmp_path, psam_text, expected_iids):
    import genoio

    prefix = write_fixed_width_plink2(tmp_path)
    prefix.with_suffix(".psam").write_text(psam_text)
    dataset = genoio.pfile(prefix)

    _, samples = dataset.read(return_samples=True)

    assert samples["iid"].to_list() == expected_iids


def test_missing_policies_nan_raise_and_impute(tmp_path):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    with_nan = dataset.read(missing="nan", dtype="float64")
    with_impute = dataset.read(missing="impute")

    assert with_nan.dtype == np.dtype("float64")
    assert np.isnan(with_nan[1, 1])
    np.testing.assert_array_equal(with_impute, np.array([[0.0, 1.0], [1.0, 0.5], [2.0, 0.0]], dtype=np.float32))
    with pytest.raises(genoio.MissingDataError, match="missing genotype"):
        dataset.read(missing="raise")


def test_missing_policy_impute_rejects_all_missing_variant(tmp_path):
    import genoio

    dataset = genoio.vcf(write_all_missing_vcf(tmp_path))

    with pytest.raises(genoio.MissingDataError, match="all-missing variant"):
        dataset.read(missing="impute")


@pytest.mark.parametrize(
    "missing",
    [
        pytest.param("nan", id="nan"),
        pytest.param("impute", id="impute"),
    ],
)
def test_missing_policy_rejects_integer_dtype_combinations(tmp_path, missing):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match=f'missing="{missing}"'):
        dataset.read(dtype=np.int16, missing=missing)


def test_missing_policy_raise_accepts_integer_dtype_when_no_missing_calls_are_retained(tmp_path):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    assert dataset.read(dtype=np.int16, missing="raise", samples=["S1"]).dtype == np.dtype("int16")


def test_read_option_validation_prioritizes_dtype_and_missing_before_samples(tmp_path):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match='missing="nan"'):
        dataset.read(dtype=np.int16, missing="nan", samples=["S1", "S1"])


def test_dense_vcf_read_rejects_multi_alt_records_even_when_gt_uses_first_alt(tmp_path):
    import genoio

    dataset = genoio.vcf(write_multi_alt_vcf(tmp_path))

    with pytest.raises(genoio.InvalidSourceError, match="multi-ALT"):
        dataset.read()


def test_unordered_sample_keep_list_returns_rows_in_source_order_and_metadata_matches(tmp_path):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    G, samples = dataset.read(samples=["S3", "S1"], return_samples=True)

    np.testing.assert_array_equal(G, np.array([[0.0, 1.0], [2.0, 0.0]], dtype=np.float32))
    assert samples["iid"].to_list() == ["S1", "S3"]


def test_missing_requested_sample_raises_structured_error_with_counts(tmp_path):
    import genoio

    dataset = genoio.vcf(write_biallelic_vcf(tmp_path))

    with pytest.raises(genoio.SampleFilterError, match="requested=2 retained=1 missing=1"):
        dataset.read(samples=["S1", "S4"])
