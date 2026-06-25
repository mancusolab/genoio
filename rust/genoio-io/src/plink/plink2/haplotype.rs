// pattern: Imperative Shell
//! Haplotype PLINK2 read orchestration.
//!
//! Haplotype readers use PGEN phase or phased-dosage auxiliary tracks and emit
//! two rows per selected sample. Sparse hard-call output rejects retained
//! missing haplotypes because sparse CSC has no missing-value channel.

use std::path::Path;

use genoio_core::{
    attach_variant_stats, reject_sparse_missing, DenseGenotypeMatrixArrowVariants, DenseLayout,
    DenseMissingPolicy, GenotypeFilterPlan, PartialFilterDecision,
    SparseGenotypeMatrixArrowVariants, VariantFilter, VariantMetadataArrowBuffers, VariantWindow,
};

use crate::error::Result;
use crate::hardcall::evaluate_packed_hardcall_filter;
use crate::matrix::apply_dense_missing_policy_to_variant;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::evaluate_dosage_filter;
use super::metadata::PvarRecordReader;
use super::pgen::{
    decode_plink2_haplotype_dosage_aux, decode_plink2_haplotype_hardcall_aux, open_pgen_payload,
    read_plink2_variant_haplotype_dosage_track, read_plink2_variant_haplotype_main_track,
    PgenDecoderState, PgenHaplotypeDecodeState,
};
use super::require_genotype_decision_filter;
use super::source::{
    empty_dense_arrow_for_samples, empty_sparse_arrow_for_samples,
    expand_selected_samples_to_haplotypes, require_pvar, variant_output_capacity,
    Plink2ReadContext,
};

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors haplotype dense read options plus metadata return choices"
)]
pub fn read_plink2_haplotypes_dense_windowed_with_arrow_variants(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected,
    } = Plink2ReadContext::new(pgen, psam, requested_samples)?;
    let mut diagnostics = selection.diagnostics.clone();
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        require_pvar(pvar)?;
        return empty_dense_arrow_for_samples(
            haplotype_samples,
            diagnostics,
            return_samples,
            return_variants,
        );
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let mut haplotype_state = PgenHaplotypeDecodeState::default();
    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let mut output_variant_count = 0_usize;
    let genotype_filter_plan = variant_filter.map_or(
        GenotypeFilterPlan::Generic,
        VariantFilter::genotype_filter_plan,
    );

    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let main_track_cursor = read_plink2_variant_haplotype_main_track(
            pgen,
            &mut file,
            &header,
            variant_index,
            &mut decoder_state,
        )?;
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
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
            if let Some(stats) = computed_stats {
                attach_variant_stats(&mut variant, stats);
            }
        }
        decode_plink2_haplotype_hardcall_aux(
            pgen,
            &header,
            variant_index,
            main_track_cursor,
            &selection.source_indices,
            &decoder_state,
            &mut haplotype_state,
        )?;
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        apply_dense_missing_policy_to_variant(
            &mut haplotype_state.selected_haplotype_values,
            &haplotype_state.selected_haplotype_missing_indices,
            missing_policy,
        )?;
        variant_major_values.extend_from_slice(&haplotype_state.selected_haplotype_values);
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = n_haplotypes;
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    let samples = if return_samples {
        haplotype_samples
    } else {
        Vec::new()
    };
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        n_variants,
        variant_major_values,
        DenseLayout::VariantMajor,
        samples,
        variants,
        diagnostics,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors haplotype dosage read options plus metadata return choices"
)]
pub fn read_plink2_haplotypes_dosage_dense_windowed_with_arrow_variants(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected: _,
    } = Plink2ReadContext::new(pgen, psam, requested_samples)?;
    let mut diagnostics = selection.diagnostics.clone();
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        require_pvar(pvar)?;
        return empty_dense_arrow_for_samples(
            haplotype_samples,
            diagnostics,
            return_samples,
            return_variants,
        );
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let mut haplotype_state = PgenHaplotypeDecodeState::default();
    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let mut output_variant_count = 0_usize;

    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let main_track_cursor = read_plink2_variant_haplotype_dosage_track(
            pgen,
            &mut file,
            &header,
            variant_index,
            &mut decoder_state,
        )?;
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

        // Haplotype buffers are needed for every retained variant; genotype
        // filters additionally read the collapsed diploid dosage scratch.
        decode_plink2_haplotype_dosage_aux(
            pgen,
            &header,
            variant_index,
            main_track_cursor,
            &selection.source_indices,
            &decoder_state,
            &mut haplotype_state,
        )?;
        if needs_genotype_decision {
            let filter = require_genotype_decision_filter(variant_filter)?;
            let (retain_variant, stats) = evaluate_dosage_filter(
                &haplotype_state.selected_collapsed_values,
                &haplotype_state.selected_collapsed_missing_indices,
                filter,
                &variant,
                return_variants,
            )?;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
            if let Some(stats) = stats {
                attach_variant_stats(&mut variant, stats);
            }
        }
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        apply_dense_missing_policy_to_variant(
            &mut haplotype_state.selected_haplotype_values,
            &haplotype_state.selected_haplotype_missing_indices,
            missing_policy,
        )?;
        variant_major_values.extend_from_slice(&haplotype_state.selected_haplotype_values);
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = n_haplotypes;
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    let samples = if return_samples {
        haplotype_samples
    } else {
        Vec::new()
    };
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        n_variants,
        variant_major_values,
        DenseLayout::VariantMajor,
        samples,
        variants,
        diagnostics,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors haplotype sparse read options plus metadata return choices"
)]
pub fn read_plink2_haplotypes_sparse_windowed_with_arrow_variants(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected,
    } = Plink2ReadContext::new(pgen, psam, requested_samples)?;
    let mut diagnostics = selection.diagnostics.clone();
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    let n_rows = haplotype_samples.len();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        require_pvar(pvar)?;
        return empty_sparse_arrow_for_samples(
            haplotype_samples,
            diagnostics,
            return_samples,
            return_variants,
        );
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let mut haplotype_state = PgenHaplotypeDecodeState::default();
    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut output_variant_count = 0_usize;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let genotype_filter_plan = variant_filter.map_or(
        GenotypeFilterPlan::Generic,
        VariantFilter::genotype_filter_plan,
    );

    let mut indptr = Vec::with_capacity(output_variant_capacity + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();

    // This mirrors the dense hard-call loop, but emits retained haplotype
    // columns directly to CSC to avoid materializing a dense intermediate.
    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let main_track_cursor = read_plink2_variant_haplotype_main_track(
            pgen,
            &mut file,
            &header,
            variant_index,
            &mut decoder_state,
        )?;
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
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
            if let Some(stats) = computed_stats {
                attach_variant_stats(&mut variant, stats);
            }
        }

        decode_plink2_haplotype_hardcall_aux(
            pgen,
            &header,
            variant_index,
            main_track_cursor,
            &selection.source_indices,
            &decoder_state,
            &mut haplotype_state,
        )?;
        append_haplotype_sparse_column(
            &mut indptr,
            &mut indices,
            &mut data,
            &haplotype_state.selected_haplotype_values,
            &haplotype_state.selected_haplotype_missing_indices,
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

    let n_cols = output_variant_count;
    diagnostics.retained_variants = n_cols;

    let samples = if return_samples {
        haplotype_samples
    } else {
        Vec::new()
    };
    SparseGenotypeMatrixArrowVariants::new(
        n_rows,
        n_cols,
        indptr,
        indices,
        data,
        samples,
        variants,
        diagnostics,
    )
}

#[inline]
fn append_haplotype_sparse_column(
    indptr: &mut Vec<usize>,
    indices: &mut Vec<usize>,
    data: &mut Vec<f32>,
    values: &[f32],
    missing_indices: &[usize],
) -> Result<()> {
    reject_sparse_missing(!missing_indices.is_empty())?;
    for (row, &value) in values.iter().enumerate() {
        if value != 0.0 {
            indices.push(row);
            data.push(value);
        }
    }
    indptr.push(indices.len());
    Ok(())
}
