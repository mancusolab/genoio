# pattern: Imperative Shell

import numpy as np
import polars as pl
import pytest
from test_cross_format_parity import EXPECTED_SAMPLES, write_canonical_plink2


def aligned_design(samples: pl.DataFrame, phenotypes: pl.DataFrame, covariates: pl.DataFrame) -> pl.DataFrame:
    design = (
        # Start from genoio sample order so y and C match genotype matrix rows.
        samples.select("iid")
        .join(phenotypes.select("iid", "expression"), on="iid", how="left")
        .join(covariates.select("iid", "age", "sex", "PC1", "PC2"), on="iid", how="left")
    )
    missing_input_rows = design.select(
        # Check every analysis column, but keep iid out of the missing-value test.
        pl.any_horizontal(pl.all().exclude("iid").is_null())
    ).to_series()
    if missing_input_rows.any():
        raise ValueError("phenotype or covariate table is missing retained samples")
    return design


def test_documented_association_workflow_aligns_phenotypes_to_genotype_sample_order(tmp_path):
    import genoio

    dataset = genoio.pfile(write_canonical_plink2(tmp_path))
    samples = dataset.samples()
    phenotypes = pl.DataFrame(
        {
            "iid": ["S3", "S1", "S_extra", "S4", "S2"],
            "expression": [3.0, 1.0, 99.0, 4.0, 2.0],
        }
    )
    covariates = pl.DataFrame(
        {
            "iid": ["S4", "S2", "S1", "S_extra", "S3"],
            "age": [44.0, 42.0, 41.0, 99.0, 43.0],
            "sex": [0.0, 1.0, 0.0, 1.0, 1.0],
            "PC1": [0.4, 0.2, 0.1, 9.9, 0.3],
            "PC2": [1.4, 1.2, 1.1, 9.9, 1.3],
        }
    )

    design = aligned_design(samples, phenotypes, covariates)
    y = design["expression"].to_numpy()
    C = design.select("age", "sex", "PC1", "PC2").to_numpy()
    scanned = []

    for X, variants in dataset.blocks(2, variants=genoio.chrom("1"), return_variants=True):
        scanned.append((X.shape, variants["id"].to_list(), y.copy(), C.copy()))

    assert design["iid"].to_list() == EXPECTED_SAMPLES
    np.testing.assert_array_equal(y, np.array([1.0, 2.0, 3.0, 4.0]))
    np.testing.assert_array_equal(
        C,
        np.array(
            [
                [41.0, 0.0, 0.1, 1.1],
                [42.0, 1.0, 0.2, 1.2],
                [43.0, 1.0, 0.3, 1.3],
                [44.0, 0.0, 0.4, 1.4],
            ]
        ),
    )
    assert [(shape, variant_ids) for shape, variant_ids, _, _ in scanned] == [
        ((4, 2), ["rs1", "rs2"])
    ]
    for _, _, block_y, block_covariates in scanned:
        np.testing.assert_array_equal(block_y, y)
        np.testing.assert_array_equal(block_covariates, C)


def test_documented_sample_filtered_association_workflow_uses_returned_sample_frame(tmp_path):
    import genoio

    dataset = genoio.pfile(write_canonical_plink2(tmp_path))
    phenotypes = pl.DataFrame(
        {
            "iid": ["S1", "S2", "S3", "S4"],
            "expression": [1.0, 2.0, 3.0, 4.0],
        }
    )
    covariates = pl.DataFrame(
        {
            "iid": ["S1", "S2", "S3", "S4"],
            "age": [41.0, 42.0, 43.0, 44.0],
            "sex": [0.0, 1.0, 1.0, 0.0],
            "PC1": [0.1, 0.2, 0.3, 0.4],
            "PC2": [1.1, 1.2, 1.3, 1.4],
        }
    )

    blocks = list(
        dataset.blocks(
            2,
            samples=["S3", "S1"],
            variants=genoio.chrom("1"),
            return_samples=True,
            return_variants=True,
        )
    )

    assert len(blocks) == 1
    X, block_samples, variants = blocks[0]
    design = aligned_design(block_samples, phenotypes, covariates)

    assert block_samples["iid"].to_list() == ["S1", "S3"]
    assert variants["id"].to_list() == ["rs1", "rs2"]
    assert X.shape == (2, 2)
    np.testing.assert_array_equal(design["expression"].to_numpy(), np.array([1.0, 3.0]))
    np.testing.assert_array_equal(
        design.select("age", "sex", "PC1", "PC2").to_numpy(),
        np.array(
            [
                [41.0, 0.0, 0.1, 1.1],
                [43.0, 1.0, 0.3, 1.3],
            ]
        ),
    )


def test_documented_association_workflow_rejects_missing_phenotype_rows(tmp_path):
    import genoio

    samples = genoio.pfile(write_canonical_plink2(tmp_path)).samples()
    phenotypes = pl.DataFrame(
        {
            "iid": ["S1", "S2", "S3"],
            "expression": [1.0, 2.0, 3.0],
        }
    )
    covariates = pl.DataFrame(
        {
            "iid": ["S1", "S2", "S3", "S4"],
            "age": [41.0, 42.0, 43.0, 44.0],
            "sex": [0.0, 1.0, 1.0, 0.0],
            "PC1": [0.1, 0.2, 0.3, 0.4],
            "PC2": [1.1, 1.2, 1.3, 1.4],
        }
    )

    with pytest.raises(ValueError, match="missing retained samples"):
        aligned_design(samples, phenotypes, covariates)
