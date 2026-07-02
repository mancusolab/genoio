//! Dense output staging for the text VCF backend.
//!
//! Text VCF reads stage retained variants contiguously in variant-major order.
//! The Python boundary handles public sample-by-variant assembly from the
//! layout metadata, so the hot record loop avoids strided sample-major writes.

use genoio_core::{
    DenseGenotypeMatrix, DenseLayout, DenseMissingPolicy, GenoioError, SampleMetadataBuffers,
    VariantMetadataBuffers,
};

use crate::error::Result;
use crate::matrix::apply_dense_missing_policy_to_variant;

pub(super) struct TextDenseOutput {
    n_samples: usize,
    values: Vec<f32>,
    variant_values: Vec<f32>,
    missing_indices: Vec<usize>,
}

impl TextDenseOutput {
    pub(super) fn new(n_samples: usize, variant_capacity: usize) -> Self {
        let len = n_samples * variant_capacity;
        Self {
            n_samples,
            values: Vec::with_capacity(len),
            variant_values: Vec::with_capacity(n_samples),
            missing_indices: Vec::new(),
        }
    }

    pub(super) fn write_variant(
        &mut self,
        decoded_values: &[f32],
        decoded_missing_indices: &[usize],
        missing_policy: DenseMissingPolicy,
    ) -> Result<()> {
        self.values.extend_from_slice(finalize_variant_values(
            &mut self.variant_values,
            &mut self.missing_indices,
            decoded_values,
            decoded_missing_indices,
            missing_policy,
        )?);
        Ok(())
    }

    /// Write a decoded variant that needs no missing-value policy work.
    ///
    /// This is the hot path for `missing="nan"` and `missing="raise"` records
    /// without missing calls. It avoids copying through the reusable scratch
    /// vector only to discover that the missing policy is a no-op.
    pub(super) fn write_variant_no_missing_direct(&mut self, decoded_values: &[f32]) -> Result<()> {
        if decoded_values.len() != self.n_samples {
            return Err(GenoioError::internal_contract(format!(
                "variant value count {} does not match sample count {}",
                decoded_values.len(),
                self.n_samples
            )));
        }
        self.values.extend_from_slice(decoded_values);
        Ok(())
    }

    pub(super) fn finish(
        self,
        n_variants: usize,
        samples: Option<SampleMetadataBuffers>,
        variants: Option<VariantMetadataBuffers>,
        diagnostics: genoio_core::DenseDiagnostics,
    ) -> Result<DenseGenotypeMatrix> {
        DenseGenotypeMatrix::new_with_layout(
            self.n_samples,
            n_variants,
            self.values,
            DenseLayout::VariantMajor,
            samples,
            variants,
            diagnostics,
        )
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
        let mut output = TextDenseOutput::new(3, 2);

        output
            .write_variant_no_missing_direct(&[0.0, 1.0, 2.0])
            .expect("direct write should append variant-major values");
        let matrix = output
            .finish(1, None, None, genoio_core::DenseDiagnostics::default())
            .expect("matrix should finish");

        assert_eq!(matrix.values, vec![0.0, 1.0, 2.0]);
    }
}
