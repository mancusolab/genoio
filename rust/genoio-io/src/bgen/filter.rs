// pattern: Functional Core
//! BGEN dosage filter evaluation helpers.
//!
//! Matrix-only reads can often decide genotype-stat filters from dosage counts
//! without building full `VariantStats`. Metadata-returning reads still attach
//! complete stats to retained variants.

use std::path::Path;

use genoio_core::{
    attach_variant_stats, DenseDiagnostics, VariantFilter, VariantRecord, VariantStats,
};

use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::retention::{RetainedVariantState, RetentionAction};

use super::decode::{
    decode_buffered_dosage_values, try_decode_buffered_dosage_values_with_counts,
    DosageDecodeBuffers,
};

pub(super) fn apply_genotype_filter_result(
    retention: &mut RetainedVariantState,
    diagnostics: &mut DenseDiagnostics,
    variant: &mut VariantRecord,
    retain_variant: bool,
    stats: Option<VariantStats>,
) -> RetentionAction {
    // Genotype-filtered reads attach computed stats only to variants that
    // survive the retained-variant window and filter decision.
    let action = retention.genotype_decision(retain_variant, diagnostics);
    if matches!(action, RetentionAction::Include) {
        if let Some(stats) = stats {
            attach_variant_stats(variant, stats);
        }
    }
    action
}

pub(super) fn decode_and_evaluate_dosage_filter(
    bgen: &Path,
    sample_count: u32,
    source_indices: &[usize],
    buffers: &mut DosageDecodeBuffers,
    filter: &VariantFilter,
    variant: &VariantRecord,
    matrix_only: bool,
) -> Result<(bool, Option<VariantStats>)> {
    let fast_counts = if matrix_only {
        // Matrix-only genotype filters only need a retain/drop decision. For
        // common fully-called phased BGENs, decode once into scratch while
        // accumulating counts and reuse that scratch if the variant is kept.
        try_decode_buffered_dosage_values_with_counts(bgen, sample_count, source_indices, buffers)?
    } else {
        None
    };
    if fast_counts.is_none() {
        decode_buffered_dosage_values(bgen, sample_count, source_indices, buffers)?;
    }
    if let Some(counts) = fast_counts {
        if let Some(retain) = counts.evaluate_plan(filter.genotype_filter_plan())? {
            return Ok((retain, None));
        }
    }
    evaluate_dosage_filter(
        &buffers.selected_values,
        &buffers.selected_missing,
        filter,
        variant,
        !matrix_only,
    )
}
