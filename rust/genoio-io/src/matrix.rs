// pattern: Functional Core

use genoio_core::{
    transpose_variant_major_to_sample_major, DenseDiagnostics, DenseGenotypeMatrix, GenoioError,
    SampleRecord, SparseGenotypeMatrix, VariantRecord,
};

use crate::error::Result;

pub(crate) struct DenseMatrixParts {
    pub(crate) n_samples: usize,
    pub(crate) n_variants: usize,
    pub(crate) values: Vec<f32>,
    pub(crate) missing_mask: Vec<bool>,
    pub(crate) samples: Vec<SampleRecord>,
    pub(crate) variants: Vec<VariantRecord>,
    pub(crate) diagnostics: DenseDiagnostics,
}

pub(crate) struct VariantMajorDenseParts {
    pub(crate) n_samples: usize,
    pub(crate) n_variants: usize,
    pub(crate) variant_major_values: Vec<f32>,
    pub(crate) variant_major_missing: Vec<bool>,
    pub(crate) samples: Vec<SampleRecord>,
    pub(crate) variants: Vec<VariantRecord>,
    pub(crate) diagnostics: DenseDiagnostics,
}

pub(crate) fn finish_dense_matrix(
    parts: DenseMatrixParts,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let DenseMatrixParts {
        n_samples,
        n_variants,
        values,
        missing_mask,
        samples,
        variants,
        diagnostics,
    } = parts;
    if matrix_only {
        DenseGenotypeMatrix::new_matrix_only(
            n_samples,
            n_variants,
            values,
            missing_mask,
            diagnostics,
        )
    } else {
        DenseGenotypeMatrix::new(
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples,
            variants,
            diagnostics,
        )
    }
}

pub(crate) fn finish_variant_major_dense_matrix(
    parts: VariantMajorDenseParts,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let VariantMajorDenseParts {
        n_samples,
        n_variants,
        variant_major_values,
        variant_major_missing,
        samples,
        variants,
        diagnostics,
    } = parts;
    validate_variant_major_len("values", variant_major_values.len(), n_samples, n_variants)?;
    validate_variant_major_len(
        "missing mask",
        variant_major_missing.len(),
        n_samples,
        n_variants,
    )?;
    let values =
        transpose_variant_major_to_sample_major(&variant_major_values, n_samples, n_variants);
    let missing_mask =
        transpose_variant_major_to_sample_major(&variant_major_missing, n_samples, n_variants);

    finish_dense_matrix(
        DenseMatrixParts {
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

pub(crate) fn shrink_sample_major_width<T: Copy>(
    values: &mut Vec<T>,
    n_samples: usize,
    old_width: usize,
    new_width: usize,
) {
    debug_assert!(new_width <= old_width);
    if old_width == new_width {
        return;
    }
    for sample_index in 1..n_samples {
        let source_start = sample_index * old_width;
        let target_start = sample_index * new_width;
        values.copy_within(source_start..source_start + new_width, target_start);
    }
    values.truncate(n_samples * new_width);
}

fn validate_variant_major_len(
    field: &str,
    actual_len: usize,
    n_samples: usize,
    n_variants: usize,
) -> Result<()> {
    let expected_len = n_samples.checked_mul(n_variants).ok_or_else(|| {
        GenoioError::internal_contract("variant-major dense matrix shape is out of range")
    })?;
    if actual_len != expected_len {
        return Err(GenoioError::internal_contract(format!(
            "variant-major dense {field} length {actual_len} does not match shape {n_samples} x {n_variants}",
        )));
    }
    Ok(())
}

pub(crate) fn empty_sparse_matrix(
    samples: Vec<SampleRecord>,
    mut diagnostics: DenseDiagnostics,
) -> Result<SparseGenotypeMatrix> {
    diagnostics.retained_variants = 0;
    SparseGenotypeMatrix::new(
        samples.len(),
        0,
        vec![0],
        Vec::new(),
        Vec::new(),
        samples,
        Vec::new(),
        diagnostics,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_record(index: usize) -> SampleRecord {
        SampleRecord {
            fid: None,
            iid: format!("sample-{index}"),
            father: None,
            mother: None,
            sex: None,
            phenotype: None,
            source_sample_index: Some(index),
            haplotype_index: None,
        }
    }

    fn variant_record(index: usize) -> VariantRecord {
        VariantRecord {
            chrom: "1".to_string(),
            pos: u32::try_from(index + 1).expect("test variant position fits in u32"),
            id: format!("variant-{index}"),
            a0: "A".to_string(),
            a1: "C".to_string(),
            ref_allele: Some("A".to_string()),
            alt_allele: Some("C".to_string()),
            source_a0: "A".to_string(),
            source_a1: "C".to_string(),
            flipped: false,
            qual: None,
            af: None,
            maf: None,
            mac: None,
            missing_rate: None,
            n_called: None,
        }
    }

    fn variant_major_parts(
        variant_major_values: Vec<f32>,
        variant_major_missing: Vec<bool>,
    ) -> VariantMajorDenseParts {
        VariantMajorDenseParts {
            n_samples: 2,
            n_variants: 2,
            variant_major_values,
            variant_major_missing,
            samples: vec![sample_record(0), sample_record(1)],
            variants: vec![variant_record(0), variant_record(1)],
            diagnostics: DenseDiagnostics::default(),
        }
    }

    #[test]
    fn finish_variant_major_dense_matrix_transposes_valid_buffers() {
        let matrix = finish_variant_major_dense_matrix(
            variant_major_parts(vec![0.0, 1.0, 2.0, 3.0], vec![false, true, false, false]),
            false,
        )
        .expect("valid variant-major buffers should build a dense matrix");

        assert_eq!(matrix.values, vec![0.0, 2.0, 1.0, 3.0]);
        assert_eq!(matrix.missing_mask, vec![false, false, true, false]);
    }

    #[test]
    fn finish_variant_major_dense_matrix_rejects_short_values_before_transpose() {
        let error = finish_variant_major_dense_matrix(
            variant_major_parts(vec![0.0, 1.0, 2.0], vec![false, true, false, false]),
            false,
        )
        .expect_err("short values should fail validation");

        assert!(error
            .to_string()
            .contains("variant-major dense values length 3 does not match shape 2 x 2"));
    }

    #[test]
    fn finish_variant_major_dense_matrix_rejects_extra_missing_mask_values() {
        let error = finish_variant_major_dense_matrix(
            variant_major_parts(
                vec![0.0, 1.0, 2.0, 3.0],
                vec![false, true, false, false, true],
            ),
            false,
        )
        .expect_err("extra missing mask values should fail validation");

        assert!(error
            .to_string()
            .contains("variant-major dense missing mask length 5 does not match shape 2 x 2"));
    }
}
