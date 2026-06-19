//! Dense output staging for the text VCF backend.
//!
//! The preferred layout depends on the read shape. All-sample matrix-only reads
//! can write directly into the final sample-major layout, while filtered reads
//! keep variant-major staging for contiguous appends and a single final
//! transpose.

use genoio_core::{DenseGenotypeMatrix, SampleRecord, VariantFilter, VariantRecord};

use crate::error::Result;
use crate::matrix::{
    finish_dense_matrix, finish_variant_major_dense_matrix, shrink_sample_major_width,
    write_sample_major_variant_slot, DenseMatrixParts, VariantMajorDenseParts,
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
        missing_mask: Vec<bool>,
    },
    VariantMajor {
        n_samples: usize,
        values: Vec<f32>,
        missing: Vec<bool>,
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
                missing_mask: vec![false; len],
            }
        } else {
            Self::VariantMajor {
                n_samples,
                values: Vec::with_capacity(len),
                missing: Vec::with_capacity(len),
            }
        }
    }

    pub(super) fn write_variant(
        &mut self,
        variant_index: usize,
        decoded_values: &[f32],
        decoded_missing: &[bool],
    ) -> Result<()> {
        match self {
            Self::SampleMajor {
                n_samples,
                row_width,
                values,
                missing_mask,
            } => write_sample_major_variant_slot(
                values,
                missing_mask,
                *n_samples,
                *row_width,
                variant_index,
                decoded_values,
                decoded_missing,
            ),
            Self::VariantMajor {
                values: output_values,
                missing: output_missing,
                ..
            } => {
                output_values.extend_from_slice(decoded_values);
                output_missing.extend_from_slice(decoded_missing);
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
                mut missing_mask,
            } => {
                shrink_sample_major_width(&mut values, n_samples, row_width, n_variants);
                shrink_sample_major_width(&mut missing_mask, n_samples, row_width, n_variants);
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
            Self::VariantMajor {
                n_samples,
                values,
                missing,
            } => finish_variant_major_dense_matrix(
                VariantMajorDenseParts {
                    n_samples,
                    n_variants,
                    variant_major_values: values,
                    variant_major_missing: missing,
                    samples,
                    variants,
                    diagnostics,
                },
                matrix_only,
            ),
        }
    }
}
