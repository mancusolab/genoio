import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest


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
        expr.left = genoio.chrom("2")


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
    assert variants["maf"].to_list() == [0.5]
    assert variants["mac"].to_list() == [3]
    assert variants["missing_rate"].to_list() == [0.0]
    assert variants["n_called"].to_list() == [3]
    np.testing.assert_array_equal(G, np.array([[0.0], [1.0], [2.0]], dtype=np.float32))


def test_qual_filter_matches_fixture_reference_results(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    G, variants = dataset.read(variants=genoio.qual(min=20) & genoio.biallelic(), return_variants=True)

    assert variants["id"].to_list() == ["rs1", "indel1"]
    assert variants["qual"].to_list() == [30.0, 40.0]
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
    assert variants["mac"].to_list() == [1]
    assert variants["n_called"].to_list() == [3]
    np.testing.assert_array_equal(G, np.array([[1.0], [0.0], [0.0]], dtype=np.float32))


def test_polymorphic_filter_matches_called_genotype_reference_results(tmp_path):
    import genoio

    dataset = genoio.vcf(write_filter_vcf(tmp_path))

    G, variants = dataset.read(variants=genoio.polymorphic(), return_variants=True)

    assert variants["id"].to_list() == ["rs1", "indel1"]
    assert variants["mac"].to_list() == [3, 1]
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
