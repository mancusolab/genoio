//! Dense output staging for the text VCF backend.
//!
//! The preferred layout depends on the read shape. All-sample matrix-only reads
//! can write directly into the final sample-major layout, while filtered reads
//! keep variant-major staging for contiguous appends and a single final
//! transpose.

use genoio_core::{
    DenseGenotypeMatrix, DenseLayout, DenseMissingPolicy, GenoioError, SampleMetadataBuffers,
    VariantFilter, VariantMetadataBuffers,
};

use crate::error::Result;
use crate::matrix::{
    apply_dense_missing_policy_to_variant, shrink_sample_major_width,
    write_sample_major_variant_slot,
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
        decoded_missing_indices: &[usize],
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
                    decoded_missing_indices,
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
                    decoded_missing_indices,
                    missing_policy,
                )?);
                Ok(())
            }
        }
    }

    /// Write a decoded variant that needs no missing-value policy work.
    ///
    /// This is the hot path for `missing="nan"` and `missing="raise"` records
    /// without missing calls. It avoids copying through the reusable scratch
    /// vector only to discover that the missing policy is a no-op.
    pub(super) fn write_variant_no_missing_direct(
        &mut self,
        variant_index: usize,
        decoded_values: &[f32],
    ) -> Result<()> {
        match self {
            Self::SampleMajor {
                n_samples,
                row_width,
                values,
                ..
            } => write_sample_major_variant_slot(
                values,
                *n_samples,
                *row_width,
                variant_index,
                decoded_values,
            ),
            Self::VariantMajor {
                n_samples,
                values: output_values,
                ..
            } => {
                if decoded_values.len() != *n_samples {
                    return Err(GenoioError::internal_contract(format!(
                        "variant value count {} does not match sample count {n_samples}",
                        decoded_values.len()
                    )));
                }
                output_values.extend_from_slice(decoded_values);
                Ok(())
            }
        }
    }

    pub(super) fn finish(
        self,
        n_variants: usize,
        samples: Option<SampleMetadataBuffers>,
        variants: Option<VariantMetadataBuffers>,
        diagnostics: genoio_core::DenseDiagnostics,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::SampleMajor {
                n_samples,
                row_width,
                mut values,
                ..
            } => {
                shrink_sample_major_width(&mut values, n_samples, row_width, n_variants);
                DenseGenotypeMatrix::new_with_layout(
                    n_samples,
                    n_variants,
                    values,
                    DenseLayout::SampleMajor,
                    samples,
                    variants,
                    diagnostics,
                )
            }
            Self::VariantMajor {
                n_samples, values, ..
            } => DenseGenotypeMatrix::new_with_layout(
                n_samples,
                n_variants,
                values,
                DenseLayout::VariantMajor,
                samples,
                variants,
                diagnostics,
            ),
        }
    }
}

fn finalize_variant_values<'a>(
    scratch_values: &'a mut Vec<f32>,
    missing_indices: &mut Vec<usize>,
    decoded_values: &[f32],
    decoded_missing_indices: &[usize],
    missing_policy: DenseMissingPolicy,
) -> Result<&'a [f32]> {
    scratch_values.clear();
    scratch_values.extend_from_slice(decoded_values);
    missing_indices.clear();
    missing_indices.extend_from_slice(decoded_missing_indices);
    apply_dense_missing_policy_to_variant(scratch_values, missing_indices, missing_policy)?;
    Ok(scratch_values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_missing_direct_write_preserves_variant_values_without_scratch_policy_work() {
        let mut output = TextDenseOutput::new(3, 2, false);

        output
            .write_variant_no_missing_direct(0, &[0.0, 1.0, 2.0])
            .expect("direct write should append variant-major values");
        let matrix = output
            .finish(1, None, None, genoio_core::DenseDiagnostics::default())
            .expect("matrix should finish");

        assert_eq!(matrix.values, vec![0.0, 1.0, 2.0]);
    }
}
