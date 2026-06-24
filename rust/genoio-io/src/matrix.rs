// pattern: Functional Core
//! Matrix assembly helpers shared by format readers.
//!
//! Format modules decode variants in the order that is cheapest for their
//! source files. These helpers validate shapes, preserve explicit dense buffer
//! layouts, and construct the public dense or sparse core structs.

use genoio_core::{
    DenseDiagnostics, DenseGenotypeMatrix, DenseLayout, DenseMissingPolicy, GenoioError,
    SampleRecord, SparseGenotypeMatrix, VariantRecord,
};

use crate::error::Result;

/// Already sample-major dense matrix components.
///
/// Use this when a reader wrote directly into the public sample-by-variant
/// layout and only needs final shape validation and metadata elision.
pub(crate) struct DenseMatrixParts {
    pub(crate) n_samples: usize,
    pub(crate) n_variants: usize,
    pub(crate) values: Vec<f32>,
    pub(crate) samples: Vec<SampleRecord>,
    pub(crate) variants: Vec<VariantRecord>,
    pub(crate) diagnostics: DenseDiagnostics,
}

/// Variant-major dense matrix components.
///
/// Use this when a decoder naturally appends one complete variant at a time.
/// The layout tag lets Python expose the public sample-by-variant shape without
/// forcing Rust to materialize a second transposed buffer.
pub(crate) struct VariantMajorDenseParts {
    pub(crate) n_samples: usize,
    pub(crate) n_variants: usize,
    pub(crate) variant_major_values: Vec<f32>,
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
        samples,
        variants,
        diagnostics,
    } = parts;
    if matrix_only {
        DenseGenotypeMatrix::new_matrix_only(n_samples, n_variants, values, diagnostics)
    } else {
        DenseGenotypeMatrix::new(
            n_samples,
            n_variants,
            values,
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
        samples,
        variants,
        diagnostics,
    } = parts;
    validate_variant_major_len("values", variant_major_values.len(), n_samples, n_variants)?;

    // Preserve the decoder's natural order. The Python bridge uses this tag to
    // build a strided NumPy view in public sample-by-variant order.
    if matrix_only {
        DenseGenotypeMatrix::new_matrix_only_with_layout(
            n_samples,
            n_variants,
            variant_major_values,
            DenseLayout::VariantMajor,
            diagnostics,
        )
    } else {
        DenseGenotypeMatrix::new_with_layout(
            n_samples,
            n_variants,
            variant_major_values,
            DenseLayout::VariantMajor,
            samples,
            variants,
            diagnostics,
        )
    }
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

pub(crate) fn apply_dense_missing_policy_to_variant(
    values: &mut [f32],
    missing_indices: &[usize],
    policy: DenseMissingPolicy,
) -> Result<()> {
    if missing_indices.is_empty() {
        return Ok(());
    }
    validate_missing_indices(values.len(), missing_indices)?;
    match policy {
        DenseMissingPolicy::Raise => Err(GenoioError::missing_data(
            "missing genotype calls are present in retained data",
        )),
        DenseMissingPolicy::Nan => {
            for &index in missing_indices {
                values[index] = f32::NAN;
            }
            Ok(())
        }
        DenseMissingPolicy::Impute => {
            let called_count = values.len() - missing_indices.len();
            if called_count == 0 {
                return Err(GenoioError::missing_data(
                    "cannot impute all-missing variant",
                ));
            }
            let called_sum = sum_called_values(values, missing_indices);
            let mean = called_sum / called_count as f32;
            for &index in missing_indices {
                values[index] = mean;
            }
            Ok(())
        }
    }
}

fn sum_called_values(values: &[f32], missing_indices: &[usize]) -> f32 {
    let mut called_sum = 0.0_f32;
    let mut missing_cursor = 0_usize;
    for (index, &value) in values.iter().enumerate() {
        if missing_indices
            .get(missing_cursor)
            .is_some_and(|&missing_index| missing_index == index)
        {
            missing_cursor += 1;
        } else {
            called_sum += value;
        }
    }
    called_sum
}

fn validate_missing_indices(values_len: usize, missing_indices: &[usize]) -> Result<()> {
    let mut previous = None;
    for &index in missing_indices {
        if index >= values_len {
            return Err(GenoioError::internal_contract(
                "dense missing index is outside variant values",
            ));
        }
        if previous.is_some_and(|previous| index <= previous) {
            return Err(GenoioError::internal_contract(
                "dense missing indices must be sorted and unique",
            ));
        }
        previous = Some(index);
    }
    Ok(())
}

/// Write a variant-major vector into one column of preallocated sample-major
/// dense buffers.
pub(crate) fn write_sample_major_variant_slot(
    values: &mut [f32],
    n_samples: usize,
    row_width: usize,
    variant_index: usize,
    variant_values: &[f32],
) -> Result<()> {
    if variant_values.len() != n_samples {
        return Err(GenoioError::internal_contract(format!(
            "variant value count {} does not match sample count {n_samples}",
            variant_values.len(),
        )));
    }
    let expected_len = n_samples.checked_mul(row_width).ok_or_else(|| {
        GenoioError::internal_contract("sample-major dense matrix shape is out of range")
    })?;
    if values.len() != expected_len {
        return Err(GenoioError::internal_contract(
            "sample-major dense buffer does not match declared shape",
        ));
    }
    if variant_index >= row_width {
        return Err(GenoioError::internal_contract(
            "sample-major variant index is outside row width",
        ));
    }

    // Readers accept variants one at a time; this fills the corresponding
    // column in each preallocated sample row without a later full transpose.
    for (sample_index, &value) in variant_values.iter().enumerate().take(n_samples) {
        let target_index = sample_index * row_width + variant_index;
        values[target_index] = value;
    }
    Ok(())
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

    fn variant_major_parts(variant_major_values: Vec<f32>) -> VariantMajorDenseParts {
        VariantMajorDenseParts {
            n_samples: 2,
            n_variants: 2,
            variant_major_values,
            samples: vec![sample_record(0), sample_record(1)],
            variants: vec![variant_record(0), variant_record(1)],
            diagnostics: DenseDiagnostics::default(),
        }
    }

    #[test]
    fn finish_variant_major_dense_matrix_preserves_variant_major_layout() {
        let matrix =
            finish_variant_major_dense_matrix(variant_major_parts(vec![0.0, 1.0, 2.0, 3.0]), false)
                .expect("valid variant-major buffers should build a dense matrix");

        assert_eq!(matrix.values, vec![0.0, 1.0, 2.0, 3.0]);
        assert_eq!(matrix.layout, DenseLayout::VariantMajor);
    }

    #[test]
    fn finish_variant_major_dense_matrix_rejects_short_values_before_layout_tagging() {
        let error =
            finish_variant_major_dense_matrix(variant_major_parts(vec![0.0, 1.0, 2.0]), false)
                .expect_err("short values should fail validation");

        assert!(error
            .to_string()
            .contains("variant-major dense values length 3 does not match shape 2 x 2"));
    }

    #[test]
    fn apply_dense_missing_policy_to_variant_imputes_missing_indices() {
        let mut values = vec![2.0, f32::NAN, 6.0];

        apply_dense_missing_policy_to_variant(&mut values, &[1], DenseMissingPolicy::Impute)
            .expect("single missing value should be imputed");

        assert_eq!(values, vec![2.0, 4.0, 6.0]);
    }

    #[test]
    fn apply_dense_missing_policy_to_variant_rejects_unsorted_indices() {
        let mut values = vec![0.0, 1.0, 2.0];
        let error =
            apply_dense_missing_policy_to_variant(&mut values, &[2, 1], DenseMissingPolicy::Nan)
                .expect_err("unsorted indices should fail");

        assert!(error
            .to_string()
            .contains("dense missing indices must be sorted and unique"));
    }

    #[test]
    fn write_sample_major_variant_slot_writes_one_column_per_sample() {
        let mut values = vec![0.0; 6];

        write_sample_major_variant_slot(&mut values, 2, 3, 1, &[4.0, 5.0])
            .expect("column write should succeed");

        assert_eq!(values, vec![0.0, 4.0, 0.0, 0.0, 5.0, 0.0]);
    }
}
