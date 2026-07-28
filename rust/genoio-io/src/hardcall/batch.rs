// pattern: Functional Core

use genoio_core::DenseMissingPolicy;

use crate::error::Result;
use crate::matrix::{apply_dense_missing_policy_to_variant, write_sample_major_variant_slot};

use super::packed::PackedHardcalls;
use super::HARDCALL_BATCH_SIZE;

/// Small batch of packed variants waiting for sample-major expansion.
///
/// Hard-call dense readers decode variants one at a time but write the final
/// matrix by sample rows. Batching keeps the transpose local without storing the
/// whole variant-major matrix.
#[derive(Debug, Clone)]
pub(crate) struct HardcallBatch {
    variants: Vec<PackedHardcalls>,
    sample_ct: usize,
}

impl HardcallBatch {
    pub(crate) fn new(sample_ct: usize) -> Self {
        Self {
            variants: Vec::with_capacity(HARDCALL_BATCH_SIZE),
            sample_ct,
        }
    }

    pub(crate) fn push(&mut self, packed: &PackedHardcalls) {
        debug_assert_eq!(packed.sample_ct(), self.sample_ct);
        let mut copy = PackedHardcalls::default();
        copy.copy_from(packed);
        self.variants.push(copy);
    }

    pub(crate) fn len(&self) -> usize {
        self.variants.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.variants.is_empty()
    }

    pub(crate) fn is_full(&self) -> bool {
        self.variants.len() == HARDCALL_BATCH_SIZE
    }

    pub(crate) fn clear(&mut self) {
        self.variants.clear();
    }
}

#[inline]
#[expect(
    clippy::too_many_arguments,
    reason = "batch flushing keeps output buffers and reusable scratch explicit on the hot path"
)]
pub(crate) fn flush_hardcall_batch_into_sample_major(
    batch: &mut HardcallBatch,
    source_indices: &[usize],
    batch_start: &mut usize,
    n_variants: usize,
    values: &mut [f32],
    missing_policy: DenseMissingPolicy,
    variant_values: &mut Vec<f32>,
    missing_indices: &mut Vec<usize>,
) -> Result<()> {
    if batch.is_empty() {
        return Ok(());
    }
    for (batch_variant_index, packed) in batch.variants.iter().enumerate() {
        let variant_index = *batch_start + batch_variant_index;
        packed.expand_selected(source_indices, variant_values, missing_indices);
        apply_dense_missing_policy_to_variant(variant_values, missing_indices, missing_policy)?;
        write_sample_major_variant_slot(
            values,
            source_indices.len(),
            n_variants,
            variant_index,
            variant_values,
        )?;
    }
    *batch_start += batch.len();
    batch.clear();
    Ok(())
}
