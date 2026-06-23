//! Dense output staging for the text VCF backend.
//!
//! The preferred layout depends on the read shape. All-sample matrix-only reads
//! can write directly into the final sample-major layout, while filtered reads
//! keep variant-major staging for contiguous appends and a single final
//! transpose.

use genoio_core::{
    DenseGenotypeMatrix, DenseMissingPolicy, SampleRecord, VariantFilter, VariantRecord,
};

use crate::error::Result;
use crate::matrix::{
    apply_dense_missing_policy_to_variant, finish_dense_matrix, finish_variant_major_dense_matrix,
    missing_indices_from_mask, shrink_sample_major_width, write_sample_major_variant_slot,
    DenseMatrixParts, VariantMajorDenseParts,
};

pub(super) fn can_write_sample_major_directly(
    selection: &genoio_core::DenseSampleSelection,
    source_sample_count: usize,
    variant_filter: Option<&VariantFilter>,
) -> bool {
    // Strided sample-major writes avoid a final transpose, but profiling showed
    // they lose locality for sample subsets and genotype-stat filters. Keep the
    // direct path to the all-sample, metadata-only case where it pays.
    !variant_filter.is_some_and(VariantFilter::requires_genotype_stats)
        && selection.source_indices.len() == source_sample_count
        && selection
            .source_indices
            .iter()
            .copied()
            .eq(0..selection.source_indices.len())
}

pub(super) enum TextDenseOutput {
    SampleMajor {
        n_samples: usize,
        row_width: usize,
        values: Vec<f32>,
        variant_values: Vec<f32>,
        missing_indices: Vec<usize>,
    },
    VariantMajor {
        n_samples: usize,
        values: Vec<f32>,
        variant_values: Vec<f32>,
        missing_indices: Vec<usize>,
    },
}

impl TextDenseOutput {
    pub(super) fn new(n_samples: usize, variant_capacity: usize, sample_major: bool) -> Self {
        let len = n_samples * variant_capacity;
        if sample_major {
            Self::SampleMajor {
                n_samples,
                row_width: variant_capacity,
                values: vec![0.0; len],
                variant_values: Vec::with_capacity(n_samples),
                missing_indices: Vec::new(),
            }
        } else {
            Self::VariantMajor {
                n_samples,
                values: Vec::with_capacity(len),
                variant_values: Vec::with_capacity(n_samples),
                missing_indices: Vec::new(),
            }
        }
    }

    pub(super) fn write_variant(
        &mut self,
        variant_index: usize,
        decoded_values: &[f32],
        decoded_missing: &[bool],
        missing_policy: DenseMissingPolicy,
    ) -> Result<()> {
        match self {
            Self::SampleMajor {
                n_samples,
                row_width,
                values,
                variant_values,
                missing_indices,
            } => write_sample_major_variant_slot(
                values,
                *n_samples,
                *row_width,
                variant_index,
                finalize_variant_values(
                    variant_values,
                    missing_indices,
                    decoded_values,
                    decoded_missing,
                    missing_policy,
                )?,
            ),
            Self::VariantMajor {
                values: output_values,
                variant_values,
                missing_indices,
                ..
            } => {
                output_values.extend_from_slice(finalize_variant_values(
                    variant_values,
                    missing_indices,
                    decoded_values,
                    decoded_missing,
                    missing_policy,
                )?);
                Ok(())
            }
        }
    }

    pub(super) fn finish(
        self,
        n_variants: usize,
        samples: Vec<SampleRecord>,
        variants: Vec<VariantRecord>,
        diagnostics: genoio_core::DenseDiagnostics,
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::SampleMajor {
                n_samples,
                row_width,
                mut values,
                ..
            } => {
                shrink_sample_major_width(&mut values, n_samples, row_width, n_variants);
                finish_dense_matrix(
                    DenseMatrixParts {
                        n_samples,
                        n_variants,
                        values,
                        samples,
                        variants,
                        diagnostics,
                    },
                    matrix_only,
                )
            }
            Self::VariantMajor {
                n_samples, values, ..
            } => finish_variant_major_dense_matrix(
                VariantMajorDenseParts {
                    n_samples,
                    n_variants,
                    variant_major_values: values,
                    samples,
                    variants,
                    diagnostics,
                },
                matrix_only,
            ),
        }
    }
}

fn finalize_variant_values<'a>(
    scratch_values: &'a mut Vec<f32>,
    missing_indices: &mut Vec<usize>,
    decoded_values: &[f32],
    decoded_missing: &[bool],
    missing_policy: DenseMissingPolicy,
) -> Result<&'a [f32]> {
    scratch_values.clear();
    scratch_values.extend_from_slice(decoded_values);
    missing_indices_from_mask(decoded_missing, missing_indices);
    apply_dense_missing_policy_to_variant(scratch_values, missing_indices, missing_policy)?;
    Ok(scratch_values)
}
