# pattern: Imperative Shell

from pathlib import Path

import numpy as np
import pytest
from scipy import sparse as scipy_sparse

EXPECTED_MATRIX = np.array(
    [
        [0.0, np.nan, 2.0, 0.0],
        [np.nan, 0.0, 1.0, 1.0],
        [2.0, 1.0, 0.0, 2.0],
        [1.0, 2.0, np.nan, 0.0],
    ],
    dtype=np.float32,
)
EXPECTED_SAMPLES = ["S1", "S2", "S3", "S4"]
EXPECTED_VARIANTS = ["rs1", "rs2", "indel1", "rs4"]
EXPECTED_VARIANT_ROWS = [
    ("1", 10, "rs1", "G", "A"),
    ("1", 20, "rs2", "T", "C"),
    ("2", 30, "indel1", "A", "AT"),
    ("2", 40, "rs4", "C", "T"),
]


def write_canonical_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "canonical.vcf"
    path.write_text(
        """\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
1\t10\trs1\tG\tA\t.\tPASS\t.\tGT\t0/0\t./.\t1/1\t0/1
1\t20\trs2\tT\tC\t.\tPASS\t.\tGT\t./.\t0/0\t0/1\t1/1
2\t30\tindel1\tA\tAT\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0\t./.
2\t40\trs4\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0
"""
    )
    return path


def write_canonical_plink1(tmp_path: Path) -> Path:
    prefix = tmp_path / "canonical_bed"
    prefix.with_suffix(".bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x87, 0x2D, 0x78, 0xCB]))
    prefix.with_suffix(".bim").write_text(
        """\
1 rs1 0 10 A G
1 rs2 0 20 C T
2 indel1 0 30 AT A
2 rs4 0 40 T C
"""
    )
    prefix.with_suffix(".fam").write_text(
        """\
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
F2 S4 0 0 2 -9
"""
    )
    return prefix


def write_canonical_plink2(tmp_path: Path) -> Path:
    prefix = tmp_path / "canonical_pgen"
    prefix.with_suffix(".pgen").write_bytes(
        bytes([0x6C, 0x1B, 0x02])
        + (4).to_bytes(4, "little")
        + (4).to_bytes(4, "little")
        + bytes([0x00, 0x6C, 0x93, 0xC6, 0x24])
    )
    prefix.with_suffix(".pvar").write_text(
        """\
#CHROM POS ID REF ALT
1 10 rs1 G A
1 20 rs2 T C
2 30 indel1 A AT
2 40 rs4 C T
"""
    )
    prefix.with_suffix(".psam").write_text(
        """\
#FID IID PAT MAT SEX PHENO
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
F2 S4 0 0 2 -9
"""
    )
    return prefix


@pytest.fixture
def canonical_datasets(tmp_path):
    import genoio

    return {
        "vcf": genoio.vcf(write_canonical_vcf(tmp_path)),
        "plink1": genoio.bfile(write_canonical_plink1(tmp_path)),
        "plink2": genoio.pfile(write_canonical_plink2(tmp_path)),
    }


def assert_core_metadata(samples, variants):
    assert samples["iid"].to_list() == EXPECTED_SAMPLES
    assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
    assert variants.rows() == EXPECTED_VARIANT_ROWS


def test_canonical_sources_read_same_matrix_and_metadata(canonical_datasets):
    for dataset in canonical_datasets.values():
        G, samples, variants = dataset.read(return_samples=True, return_variants=True)

        np.testing.assert_array_equal(G, EXPECTED_MATRIX)
        assert_core_metadata(samples, variants)


def test_canonical_sources_keep_sample_selection_in_source_order(canonical_datasets):
    expected = EXPECTED_MATRIX[[0, 2], :]

    for dataset in canonical_datasets.values():
        G, samples = dataset.read(samples=["S3", "S1"], return_samples=True)

        np.testing.assert_array_equal(G, expected)
        assert samples["iid"].to_list() == ["S1", "S3"]


@pytest.mark.parametrize(
    ("expr_factory", "expected_ids", "expected_matrix"),
    [
        (
            lambda genoio: genoio.chrom("1"),
            ["rs1", "rs2"],
            EXPECTED_MATRIX[:, [0, 1]],
        ),
        (
            lambda genoio: genoio.chrom("2") | genoio.id_in(["rs1"]),
            ["rs1", "indel1", "rs4"],
            EXPECTED_MATRIX[:, [0, 2, 3]],
        ),
        (
            lambda genoio: genoio.chrom("1") & ~genoio.id_in(["rs2"]),
            ["rs1"],
            EXPECTED_MATRIX[:, [0]],
        ),
        (
            lambda genoio: genoio.maf(min=0.2) & genoio.missing_rate(0.25),
            ["rs1", "rs2", "indel1", "rs4"],
            EXPECTED_MATRIX,
        ),
        (
            lambda genoio: genoio.mac(min=3, max=3),
            ["rs1", "rs2", "indel1", "rs4"],
            EXPECTED_MATRIX,
        ),
        (
            lambda genoio: genoio.polymorphic() & genoio.missing_rate(0.0),
            ["rs4"],
            EXPECTED_MATRIX[:, [3]],
        ),
    ],
)
def test_canonical_sources_apply_filters_equivalently(canonical_datasets, expr_factory, expected_ids, expected_matrix):
    import genoio

    expr = expr_factory(genoio)

    for dataset in canonical_datasets.values():
        G, variants = dataset.read(variants=expr, return_variants=True)

        assert variants["id"].to_list() == expected_ids
        np.testing.assert_array_equal(G, expected_matrix)


def test_canonical_sources_impute_missing_values_equivalently(canonical_datasets):
    expected = np.array(
        [
            [0.0, 1.0, 2.0, 0.0],
            [1.0, 0.0, 1.0, 1.0],
            [2.0, 1.0, 0.0, 2.0],
            [1.0, 2.0, 1.0, 0.0],
        ],
        dtype=np.float32,
    )

    for dataset in canonical_datasets.values():
        np.testing.assert_array_equal(dataset.read(missing="impute"), expected)


def test_canonical_sources_sparse_reads_match_dense_when_no_missing_calls_remain(canonical_datasets):
    import genoio

    expr = genoio.polymorphic() & genoio.missing_rate(0.0)

    for dataset in canonical_datasets.values():
        dense = dataset.read(variants=expr)
        sparse = dataset.read(variants=expr, sparse=True)

        assert scipy_sparse.isspmatrix_csc(sparse)
        np.testing.assert_array_equal(sparse.toarray(), dense)


def test_canonical_sources_blocks_concatenate_to_filtered_full_reads(canonical_datasets):
    import genoio

    read_options = {"variants": genoio.chrom("2") | genoio.id_in(["rs1"]), "samples": ["S3", "S1"]}

    for dataset in canonical_datasets.values():
        full, full_variants = dataset.read(**read_options, return_variants=True)
        blocks = list(dataset.iter_blocks(size=2, **read_options, return_variants=True))

        assert [variants["id"].to_list() for _, variants in blocks] == [["rs1", "indel1"], ["rs4"]]
        np.testing.assert_array_equal(np.concatenate([block for block, _ in blocks], axis=1), full)
        assert full_variants["id"].to_list() == ["rs1", "indel1", "rs4"]


def test_canonical_sources_allow_filters_that_retain_zero_variants(canonical_datasets):
    import genoio

    expr = genoio.chrom("21")

    for dataset in canonical_datasets.values():
        dense, samples, variants = dataset.read(variants=expr, return_samples=True, return_variants=True)
        sparse = dataset.read(variants=expr, sparse=True)
        imputed = dataset.read(variants=expr, missing="impute")
        blocks = list(dataset.iter_blocks(size=2, variants=expr, return_variants=True))

        assert dense.shape == (len(EXPECTED_SAMPLES), 0)
        assert samples["iid"].to_list() == EXPECTED_SAMPLES
        assert variants.shape == (0, 5)
        assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
        assert scipy_sparse.isspmatrix_csc(sparse)
        assert sparse.shape == dense.shape
        assert sparse.nnz == 0
        np.testing.assert_array_equal(imputed, dense)
        assert blocks == []
