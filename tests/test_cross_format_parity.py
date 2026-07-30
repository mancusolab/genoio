# pattern: Imperative Shell

import numpy as np
import pytest
from fixture_writers import (
    EXPECTED_MATRIX,
    EXPECTED_SAMPLES,
    EXPECTED_VARIANT_ROWS,
    write_canonical_plink1,
    write_canonical_plink2,
    write_canonical_vcf,
)
from scipy import sparse as scipy_sparse


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
        sparse_matrix = dataset.read(variants=expr, sparse=True)

        assert scipy_sparse.isspmatrix_csc(sparse_matrix)
        np.testing.assert_array_equal(sparse_matrix.toarray(), dense)


@pytest.mark.parametrize(
    ("sparse", "dtype", "size", "filter_mode"),
    [
        pytest.param(False, "float32", 1, "three-variants", id="dense-size-1"),
        pytest.param(False, "float64", 2, "three-variants", id="dense-partial-final"),
        pytest.param(False, "float32", 3, "three-variants", id="dense-exact"),
        pytest.param(False, "float64", 5, "three-variants", id="dense-oversized"),
        pytest.param("csc", "float32", 1, "complete-only", id="csc"),
        pytest.param("csr", "float64", 3, "complete-only", id="csr"),
    ],
)
def test_pbr_py_matrix_001_pbr_py_meta_001_canonical_blocks_match_filtered_reads(
    canonical_datasets,
    sparse,
    dtype,
    size,
    filter_mode,
):
    import genoio

    variants = (
        genoio.chrom("2") | genoio.id_in(["rs1"]) if filter_mode == "three-variants" else genoio.missing_rate(0.0)
    )
    read_options = {
        "variants": variants,
        "samples": ["S3", "S1"],
        "sparse": sparse,
        "dtype": dtype,
    }

    for dataset in canonical_datasets.values():
        full, full_samples, full_variants = dataset.read(
            **read_options,
            return_samples=True,
            return_variants=True,
        )
        blocks = list(
            dataset.iter_blocks(
                size=size,
                **read_options,
                return_samples=True,
                return_variants=True,
            )
        )

        assert blocks
        assert all(block.shape[1] <= size for block, _, _ in blocks)
        if sparse:
            assert all(getattr(scipy_sparse, f"isspmatrix_{sparse}")(block) for block, _, _ in blocks)
            combined = scipy_sparse.hstack([block for block, _, _ in blocks], format=sparse)
            np.testing.assert_array_equal(combined.toarray(), full.toarray())
        else:
            np.testing.assert_array_equal(
                np.concatenate([block for block, _, _ in blocks], axis=1),
                full,
            )
        assert full.dtype == np.dtype(dtype)
        assert all(block.dtype == np.dtype(dtype) for block, _, _ in blocks)
        assert all(samples.schema == full_samples.schema for _, samples, _ in blocks)
        assert all(samples.equals(full_samples) for _, samples, _ in blocks)
        assert all(block_variants.schema == full_variants.schema for _, _, block_variants in blocks)
        assert [row for _, _, block_variants in blocks for row in block_variants.rows()] == full_variants.rows()


@pytest.mark.parametrize("sparse", [False, "csr"])
def test_pbr_py_matrix_001_canonical_blocks_allow_filters_that_retain_zero_variants(
    canonical_datasets,
    sparse,
):
    import genoio

    expr = genoio.chrom("21")

    for dataset in canonical_datasets.values():
        dense, samples, variants = dataset.read(variants=expr, return_samples=True, return_variants=True)
        sparse_matrix = dataset.read(variants=expr, sparse=True)
        imputed = dataset.read(variants=expr, missing="impute")
        blocks = list(
            dataset.iter_blocks(
                size=2,
                variants=expr,
                sparse=sparse,
                dtype="float64",
                return_variants=True,
            )
        )

        assert dense.shape == (len(EXPECTED_SAMPLES), 0)
        assert samples["iid"].to_list() == EXPECTED_SAMPLES
        assert variants.shape == (0, 5)
        assert variants.columns == ["chrom", "pos", "id", "a0", "a1"]
        assert scipy_sparse.isspmatrix_csc(sparse_matrix)
        assert sparse_matrix.shape == dense.shape
        assert sparse_matrix.nnz == 0
        np.testing.assert_array_equal(imputed, dense)
        assert blocks == []
