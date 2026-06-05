# pattern: Imperative Shell

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from test_dense_read import write_bgen_dosage, write_fixed_width_plink2


def write_filter_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "filters.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t30\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t10\tPASS\t.\tGT\t0/0\t./.\t0/0
2\t30\tindel1\tAT\tA\t40\tPASS\t.\tGT\t0/1\t0/0\t0/0
1\t40\trs4\tG\tA\t.\tPASS\t.\tGT\t./.\t./.\t./.
"""
    )
    return path


def write_indexed_pushdown_vcf(tmp_path: Path) -> Path:
    source = tmp_path / "indexed_pushdown.vcf"
    source.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0
1\t20\trs20\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t30\trs30\tG\tA\t.\tPASS\t.\tGT\t0/0\t./.\t0/0
1\t40\tbad_outside_region\tT\tC\t.\tPASS\t.\tGT\t0/3\t0/0\t0/0
"""
    )
    compressed = tmp_path / "indexed_pushdown.vcf.gz"
    with compressed.open("wb") as output:
        subprocess.run(["bgzip", "-c", str(source)], stdout=output, check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(compressed)], check=True)
    return compressed


def write_plink1_equivalent_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "plink1_equivalent.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tG\tA\t.\tPASS\t.\tGT\t0/0\t./.\t1/1
1\t20\trs2\tT\tC\t.\tPASS\t.\tGT\t./.\t0/0\t0/1
2\t30\tindel1\tA\tAT\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0
"""
    )
    return path


def assert_read_matches(dataset, expr, expected_ids, expected_matrix):
    G, variants = dataset.read(variants=expr, return_variants=True)

    assert variants["id"].to_list() == expected_ids
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    np.testing.assert_array_equal(G, np.array(expected_matrix, dtype=np.float32))


def test_filter_expressions_are_frozen_composable_and_json_serializable():
    import genoio

    expr = (genoio.chrom("1") & genoio.snp()) | ~genoio.id_in(["rs4"])

    assert expr.to_ir() == {
        "op": "or",
        "left": {
            "op": "and",
            "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
            "right": {"op": "predicate", "name": "snp", "params": {}},
        },
        "right": {
            "op": "not",
            "expr": {
                "op": "predicate",
                "name": "id_in",
                "params": {"values": ["rs4"]},
            },
        },
    }
    json.dumps(expr.to_ir())

    with pytest.raises(AttributeError):
        expr.left = genoio.chrom("2")  # pyright: ignore[reportAttributeAccessIssue]


@pytest.mark.parametrize(
    "factory",
    [
        lambda genoio: genoio.region("1:0-10"),
        lambda genoio: genoio.region("1:20-10"),
        lambda genoio: genoio.region("1"),
        lambda genoio: genoio.maf(min=-0.1),
        lambda genoio: genoio.maf(max=1.1),
        lambda genoio: genoio.maf(min=0.4, max=0.2),
        lambda genoio: genoio.maf(min=float("nan")),
        lambda genoio: genoio.maf(max=float("nan")),
        lambda genoio: genoio.qual(),
        lambda genoio: genoio.qual(min=-1),
        lambda genoio: genoio.qual(max=float("nan")),
        lambda genoio: genoio.qual(min=40, max=20),
        lambda genoio: genoio.mac(min=-1),
        lambda genoio: genoio.missing_rate(float("nan")),
        lambda genoio: genoio.missing_rate(1.1),
        lambda genoio: genoio.id_in(["rs1", "rs1"]),
        lambda genoio: genoio.id_in(["rs1", object()]),
    ],
)
def test_filter_constructors_reject_invalid_values(factory):
    import genoio

    with pytest.raises(genoio.InvalidOptionError):
        factory(genoio)


def test_variants_accepts_composed_filter_and_matches_polars_numpy_reference(tmp_path):
    import genoio

    path = write_filter_vcf(tmp_path)
    dataset = genoio.vcf(path)
    expr = genoio.snp() & genoio.maf(min=0.1) & genoio.missing_rate(0.5)

    G, variants = dataset.read(variants=expr, return_variants=True)

    source = pl.DataFrame(
        {
            "id": ["rs1", "rs2", "indel1", "rs4"],
            "ref": ["A", "C", "AT", "G"],
            "alt": ["G", "T", "A", "A"],
            "maf": [0.5, 0.0, 1.0 / 6.0, None],
            "missing_rate": [0.0, 1.0 / 3.0, 0.0, 1.0],
        }
    )
    expected_ids = (
        source.lazy()
        .filter(
            (pl.col("ref").str.len_chars() == 1)
            & (pl.col("alt").str.len_chars() == 1)
            & (pl.col("maf") >= 0.1)
            & (pl.col("missing_rate") <= 0.5)
        )
        .select("id")
        .collect()
        .get_column("id")
        .to_list()
    )

    assert variants["id"].to_list() == expected_ids == ["rs1"]
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    np.testing.assert_array_equal(G, np.array([[0.0], [1.0], [2.0]], dtype=np.float32))


@pytest.mark.parametrize(
    ("expr_factory", "expected_ids", "expected_matrix"),
    [
        (
            lambda genoio: genoio.chrom("2") | genoio.id_in(["rs2"]),
            ["rs2", "indel1"],
            [[0.0, 1.0], [np.nan, 0.0], [0.0, 0.0]],
        ),
        (
            lambda genoio: genoio.chrom("1") & ~genoio.id_in(["rs2", "rs4"]),
            ["rs1"],
            [[0.0], [1.0], [2.0]],
        ),
        (
            lambda genoio: (genoio.qual(min=20) | genoio.maf(min=0.4)) & genoio.missing_rate(0.5),
            ["rs1", "indel1"],
            [[0.0, 1.0], [1.0, 0.0], [2.0, 0.0]],
        ),
        (
            lambda genoio: genoio.region("1:1-25") & genoio.maf(min=0.1),
            ["rs1"],
            [[0.0], [1.0], [2.0]],
        ),
    ],
)
def test_vcf_filter_boolean_combinations_match_fixture_reference_results(
    tmp_path, expr_factory, expected_ids, expected_matrix
):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    assert_read_matches(dataset, expr_factory(genoio), expected_ids, expected_matrix)


def test_qual_filter_matches_fixture_reference_results(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    G, variants = dataset.read(variants=genoio.qual(min=20) & genoio.biallelic(), return_variants=True)

    assert variants["id"].to_list() == ["rs1", "indel1"]
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    np.testing.assert_array_equal(
        G,
        np.array(
            [
                [0.0, 1.0],
                [1.0, 0.0],
                [2.0, 0.0],
            ],
            dtype=np.float32,
        ),
    )


def test_biallelic_filter_matches_fixture_reference_results(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    G, variants = dataset.read(variants=genoio.biallelic(), return_variants=True)

    assert variants["id"].to_list() == ["rs1", "rs2", "indel1", "rs4"]
    np.testing.assert_array_equal(
        G,
        np.array(
            [
                [0.0, 0.0, 1.0, np.nan],
                [1.0, np.nan, 0.0, np.nan],
                [2.0, 0.0, 0.0, np.nan],
            ],
            dtype=np.float32,
        ),
    )


def test_mac_filter_matches_called_genotype_reference_results(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    G, variants = dataset.read(variants=genoio.mac(min=1, max=1), return_variants=True)

    assert variants["id"].to_list() == ["indel1"]
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    np.testing.assert_array_equal(G, np.array([[1.0], [0.0], [0.0]], dtype=np.float32))


def test_plink2_genotype_filters_return_core_variant_metadata(tmp_path):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    G, variants = dataset.read(
        variants=genoio.maf(min=0.2) & genoio.missing_rate(0.5),
        return_variants=True,
    )

    assert variants["id"].to_list() == ["rs1", "rs2", "rs3"]
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    np.testing.assert_array_equal(
        G,
        np.array([[0.0, 1.0, 2.0], [np.nan, 0.0, 1.0], [2.0, 1.0, 0.0]], dtype=np.float32),
    )


@pytest.mark.parametrize(
    ("expr_factory", "expected_ids", "expected_matrix"),
    [
        (
            lambda genoio: genoio.maf(min=0.2),
            ["rs2"],
            [[1.0], [np.nan]],
        ),
        (
            lambda genoio: genoio.mac(min=1, max=1),
            ["rs2"],
            [[1.0], [np.nan]],
        ),
        (
            lambda genoio: genoio.missing_rate(0.0),
            ["rs1"],
            [[0.0], [0.0]],
        ),
        (
            lambda genoio: genoio.polymorphic(),
            ["rs2"],
            [[1.0], [np.nan]],
        ),
    ],
)
def test_bgen_dosage_genotype_filters_match_dosage_reference_results(
    tmp_path, expr_factory, expected_ids, expected_matrix
):
    import genoio

    path = write_bgen_dosage(
        tmp_path,
        variant_calls=[
            [(255, 0), (255, 0)],
            [(0, 255), None],
        ],
    )
    dataset = genoio.bgen(path)

    G, variants = dataset.read(
        dosage="dosage",
        variants=expr_factory(genoio),
        return_variants=True,
    )

    assert variants["id"].to_list() == expected_ids
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    np.testing.assert_array_equal(G, np.array(expected_matrix, dtype=np.float32))


def test_bgen_dosage_empty_variant_filter_returns_empty_matrix_and_variant_schema(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    G, variants = dataset.read(dosage="dosage", variants=[], return_variants=True)

    assert G.shape == (2, 0)
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    assert variants.height == 0


def test_bgen_dosage_nonmatching_metadata_filter_returns_empty_matrix_and_variant_schema(tmp_path):
    import genoio

    dataset = genoio.bgen(write_bgen_dosage(tmp_path))

    G, variants = dataset.read(dosage="dosage", variants=genoio.chrom("9"), return_variants=True)

    assert G.shape == (2, 0)
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    assert variants.height == 0


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="indexed VCF public API test requires bgzip and tabix",
)
def test_indexed_vcf_region_pushdown_combines_with_genotype_stat_filter(tmp_path):
    import genoio

    dataset = genoio.vcf(write_indexed_pushdown_vcf(tmp_path))

    G, variants = dataset.read(
        variants=genoio.region("1:10-30") & genoio.maf(min=0.2),
        return_variants=True,
    )

    assert variants["id"].to_list() == ["rs20"]
    np.testing.assert_array_equal(G, np.array([[0.0], [1.0], [2.0]], dtype=np.float32))
    with pytest.raises(genoio.InvalidSourceError, match="multiallelic GT"):
        dataset.read(variants=genoio.maf(min=0.2))


@pytest.mark.parametrize(
    ("expr_factory", "expected_ids", "expected_matrix"),
    [
        (
            lambda genoio: genoio.maf(min=0.2) & genoio.missing_rate(0.5),
            ["rs1", "rs2", "indel1"],
            [[0.0, np.nan, 2.0], [np.nan, 0.0, 1.0], [2.0, 1.0, 0.0]],
        ),
        (
            lambda genoio: genoio.mac(min=1, max=1),
            ["rs2"],
            [[np.nan], [0.0], [1.0]],
        ),
        (
            lambda genoio: genoio.polymorphic() & genoio.missing_rate(0.0),
            ["indel1"],
            [[2.0], [1.0], [0.0]],
        ),
        (
            lambda genoio: genoio.chrom("1") & ~genoio.id_in(["rs2"]),
            ["rs1"],
            [[0.0], [np.nan], [2.0]],
        ),
    ],
)
def test_plink1_filter_combinations_match_fixture_reference_results(expr_factory, expected_ids, expected_matrix):
    import genoio

    dataset = genoio.bfile(Path(__file__).parent / "fixtures" / "plink1" / "tiny")

    assert_read_matches(dataset, expr_factory(genoio), expected_ids, expected_matrix)


@pytest.mark.parametrize(
    ("expr_factory", "expected_ids", "expected_matrix"),
    [
        (
            lambda genoio: genoio.chrom("1") & ~genoio.id_in(["rs2"]),
            ["rs1"],
            [[0.0], [np.nan], [2.0]],
        ),
        (
            lambda genoio: genoio.polymorphic() & genoio.missing_rate(0.0),
            ["rs2", "rs3"],
            [[1.0, 2.0], [0.0, 1.0], [1.0, 0.0]],
        ),
        (
            lambda genoio: genoio.chrom("2") | genoio.id_in(["rs1"]),
            ["rs1", "rs3"],
            [[0.0, 2.0], [np.nan, 1.0], [2.0, 0.0]],
        ),
    ],
)
def test_plink2_filter_boolean_combinations_match_fixture_reference_results(
    tmp_path, expr_factory, expected_ids, expected_matrix
):
    import genoio

    dataset = genoio.pfile(write_fixed_width_plink2(tmp_path))

    assert_read_matches(dataset, expr_factory(genoio), expected_ids, expected_matrix)


@pytest.mark.parametrize(
    "expr_factory",
    [
        lambda genoio: genoio.chrom("1") & genoio.missing_rate(0.5),
        lambda genoio: (genoio.chrom("1") | genoio.id_in(["indel1"])) & genoio.polymorphic(),
        lambda genoio: genoio.chrom("2") | genoio.mac(min=1, max=1),
    ],
)
def test_vcf_and_plink1_filters_retain_equivalent_variants_on_matching_sources(tmp_path, expr_factory):
    import genoio

    expr = expr_factory(genoio)
    vcf_dataset = genoio.vcf(write_plink1_equivalent_vcf(tmp_path))
    plink1_dataset = genoio.bfile(Path(__file__).parent / "fixtures" / "plink1" / "tiny")

    G_vcf, variants_vcf = vcf_dataset.read(variants=expr, return_variants=True)
    G_plink1, variants_plink1 = plink1_dataset.read(variants=expr, return_variants=True)

    assert variants_vcf["id"].to_list() == variants_plink1["id"].to_list()
    np.testing.assert_array_equal(G_vcf, G_plink1)


def test_polymorphic_filter_matches_called_genotype_reference_results(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    G, variants = dataset.read(variants=genoio.polymorphic(), return_variants=True)

    assert variants["id"].to_list() == ["rs1", "indel1"]
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    np.testing.assert_array_equal(
        G,
        np.array(
            [
                [0.0, 1.0],
                [1.0, 0.0],
                [2.0, 0.0],
            ],
            dtype=np.float32,
        ),
    )


def test_iterable_variant_id_selection_preserves_source_variant_order(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    G, variants = dataset.read(variants=["rs2", "rs1"], return_variants=True)

    assert variants["id"].to_list() == ["rs1", "rs2"]
    np.testing.assert_array_equal(
        G,
        np.array([[0.0, 0.0], [1.0, np.nan], [2.0, 0.0]], dtype=np.float32),
    )


def test_generator_variant_id_selection_is_consumed_once(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    G, variants = dataset.read(
        variants=(variant_id for variant_id in ["rs2", "rs1"]),
        return_variants=True,
    )

    assert variants["id"].to_list() == ["rs1", "rs2"]
    np.testing.assert_array_equal(
        G,
        np.array([[0.0, 0.0], [1.0, np.nan], [2.0, 0.0]], dtype=np.float32),
    )


def test_variants_rejects_python_callbacks(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    with pytest.raises(genoio.InvalidOptionError, match="serializable"):
        dataset.read(variants=lambda variant: True)
