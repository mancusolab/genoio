// pattern: Imperative Shell
//! Sparse hard-call PLINK2 read orchestration.
//!
//! The sparse path decodes retained hard calls, applies the same filters as the
//! dense path, flips common alternate columns when needed, and emits CSC columns
//! without materializing a dense matrix.

use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, flip_values_to_minor_allele, reject_sparse_missing,
    GenotypeFilterPlan, PartialFilterDecision, SampleMetadataArrowBuffers,
    SparseGenotypeMatrixArrowVariants, VariantFilter, VariantMetadataArrowBuffers, VariantRecord,
    VariantWindow,
};

use crate::error::Result;
use crate::hardcall::evaluate_packed_hardcall_filter;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::metadata::{parse_pvar_source_window, PvarRecordReader};
use super::pgen::{
    open_pgen_payload, read_plink2_variant_packed, read_plink2_variant_values, PgenDecoderState,
    PgenLayout,
};
use super::require_genotype_decision_filter;
use super::source::{
    empty_sparse_arrow_for_selection, require_pvar, variant_output_capacity, Plink2ReadContext,
};

/// Source-window row that may omit metadata for matrix-only reads.
struct SourceWindowVariant {
    source_index: usize,
    variant: Option<VariantRecord>,
}

#[inline]
fn append_decoded_sparse_column(
    decoder_state: &mut PgenDecoderState,
    variant: &mut VariantRecord,
    indptr: &mut Vec<i32>,
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
) -> Result<()> {
    reject_sparse_missing(!decoder_state.missing_indices.is_empty())?;
    flip_values_to_minor_allele(&mut decoder_state.values, variant);
    append_sparse_column(indptr, indices, data, &decoder_state.values)
}

#[inline]
fn append_decoded_sparse_column_without_variant_metadata(
    decoder_state: &mut PgenDecoderState,
    indptr: &mut Vec<i32>,
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
) -> Result<()> {
    reject_sparse_missing(!decoder_state.missing_indices.is_empty())?;
    // Matrix-only source windows must preserve sparse minor-allele orientation
    // without mutating metadata that the caller did not request.
    flip_values_to_minor_allele_without_metadata(&mut decoder_state.values);
    append_sparse_column(indptr, indices, data, &decoder_state.values)
}

fn flip_values_to_minor_allele_without_metadata(values: &mut [f32]) {
    let a1_count = values.iter().sum::<f32>();
    let a0_count = 2.0 * values.len() as f32 - a1_count;
    if a1_count <= a0_count {
        return;
    }
    for value in values {
        *value = 2.0 - *value;
    }
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors sparse read options plus metadata return choices"
)]
pub fn read_plink2_sparse_windowed_with_arrow_variants(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    // See the dense fast path: unfiltered windows can be interpreted directly
    // in source coordinates, but filtered windows cannot.
    if let (None, Some(window)) = (variant_filter, variant_window) {
        return read_plink2_sparse_source_window_arrow(
            pgen,
            pvar,
            psam,
            requested_samples,
            window,
            return_samples,
            return_variants,
        );
    }

    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected,
    } = Plink2ReadContext::new(pgen, psam, requested_samples)?;
    let mut diagnostics = selection.diagnostics.clone();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        require_pvar(pvar)?;
        return empty_sparse_arrow_for_selection(selection, return_samples, return_variants);
    }
    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let mut indptr = Vec::with_capacity(output_variant_capacity + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut output_variant_count = 0_usize;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let genotype_filter_plan = variant_filter.map_or(
        GenotypeFilterPlan::Generic,
        VariantFilter::genotype_filter_plan,
    );
    let requires_sequential_decode = matches!(header.layout, PgenLayout::VariableWidth);
    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let mut decoded_packed = false;
        if requires_sequential_decode {
            read_plink2_variant_packed(
                pgen,
                &mut file,
                &header,
                variant_index,
                &mut decoder_state,
            )?;
            decoded_packed = true;
        }
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => {
                stopped_after_window = true;
                break;
            }
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        if !decoded_packed {
            read_plink2_variant_packed(
                pgen,
                &mut file,
                &header,
                variant_index,
                &mut decoder_state,
            )?;
        }
        let mut stats = None;
        if needs_genotype_decision {
            let filter = require_genotype_decision_filter(variant_filter)?;
            let (retain_variant, computed_stats) = evaluate_packed_hardcall_filter(
                &decoder_state.packed,
                &selection.source_indices,
                all_samples_selected,
                filter,
                genotype_filter_plan,
                Some(&variant),
                return_variants,
            )?;
            stats = computed_stats;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        decoder_state.packed.expand_selected(
            &selection.source_indices,
            &mut decoder_state.values,
            &mut decoder_state.missing_indices,
        );
        append_decoded_sparse_column(
            &mut decoder_state,
            &mut variant,
            &mut indptr,
            &mut indices,
            &mut data,
        )?;
        output_variant_count += 1;
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    let samples = SampleMetadataArrowBuffers::optional_from_records(
        &selection.samples,
        return_samples,
        false,
    )?;
    SparseGenotypeMatrixArrowVariants::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        samples,
        variants,
        diagnostics,
    )
}

fn read_plink2_sparse_source_window_arrow(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    window: VariantWindow,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    let decode_variant_ct = window.start.saturating_add(window.len);
    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected: _,
    } = Plink2ReadContext::new_prefix(pgen, psam, requested_samples, decode_variant_ct)?;
    let mut diagnostics = selection.diagnostics;
    let n_variants = header
        .variant_ct
        .saturating_sub(window.start)
        .min(window.len);
    let window_variants: Vec<SourceWindowVariant> = if return_variants {
        parse_pvar_source_window(pvar, window, header.variant_ct)?
            .into_iter()
            .map(|(source_index, variant)| SourceWindowVariant {
                source_index,
                variant: Some(variant),
            })
            .collect()
    } else {
        require_pvar(pvar)?;
        (window.start..window.start + n_variants)
            .map(|source_index| SourceWindowVariant {
                source_index,
                variant: None,
            })
            .collect()
    };
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let output_variant_capacity = n_variants;
    let mut indptr = Vec::with_capacity(output_variant_capacity + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));

    match header.layout {
        PgenLayout::FixedWidth
        | PgenLayout::FixedWidthDosage
        | PgenLayout::FixedWidthPhasedDosage => {
            // Fixed-width records can be decoded by direct source index.
            for source_variant in window_variants {
                read_plink2_variant_values(
                    pgen,
                    &mut file,
                    &header,
                    source_variant.source_index,
                    &selection.source_indices,
                    &mut decoder_state,
                )?;
                if let Some(mut variant) = source_variant.variant {
                    append_decoded_sparse_column(
                        &mut decoder_state,
                        &mut variant,
                        &mut indptr,
                        &mut indices,
                        &mut data,
                    )?;
                    if let Some(variants) = variants.as_mut() {
                        variants.push_record(&variant)?;
                    }
                } else {
                    append_decoded_sparse_column_without_variant_metadata(
                        &mut decoder_state,
                        &mut indptr,
                        &mut indices,
                        &mut data,
                    )?;
                }
            }
        }
        PgenLayout::VariableWidth => {
            let mut window_iter = window_variants.into_iter().peekable();
            let prefix_end = header.record_types.len();
            // Preserve LD state exactly as dense reads do, then append only
            // requested variants to sparse columns.
            for variant_index in 0..prefix_end {
                read_plink2_variant_values(
                    pgen,
                    &mut file,
                    &header,
                    variant_index,
                    &selection.source_indices,
                    &mut decoder_state,
                )?;
                if window_iter
                    .peek()
                    .is_some_and(|source_variant| source_variant.source_index == variant_index)
                {
                    if let Some(source_variant) = window_iter.next() {
                        if let Some(mut variant) = source_variant.variant {
                            append_decoded_sparse_column(
                                &mut decoder_state,
                                &mut variant,
                                &mut indptr,
                                &mut indices,
                                &mut data,
                            )?;
                            if let Some(variants) = variants.as_mut() {
                                variants.push_record(&variant)?;
                            }
                        } else {
                            append_decoded_sparse_column_without_variant_metadata(
                                &mut decoder_state,
                                &mut indptr,
                                &mut indices,
                                &mut data,
                            )?;
                        }
                    }
                }
            }
        }
    }

    let n_variants = indptr.len().saturating_sub(1);
    diagnostics.candidate_variants = n_variants;
    diagnostics.retained_variants = n_variants;
    let samples = SampleMetadataArrowBuffers::optional_from_records(
        &selection.samples,
        return_samples,
        false,
    )?;
    SparseGenotypeMatrixArrowVariants::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        samples,
        variants,
        diagnostics,
    )
}
