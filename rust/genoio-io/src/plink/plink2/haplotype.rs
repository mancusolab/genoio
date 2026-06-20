// pattern: Imperative Shell
//! Haplotype PLINK2 read orchestration.
//!
//! Haplotype readers use PGEN phase or phased-dosage auxiliary tracks and emit
//! two rows per selected sample. Sparse hard-call output rejects retained
//! missing haplotypes because sparse CSC has no missing-value channel.

use std::path::Path;

use genoio_core::{
    attach_variant_stats, reject_sparse_missing_values, DenseGenotypeMatrix, GenotypeFilterPlan,
    PartialFilterDecision, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};

use crate::error::Result;
use crate::hardcall::evaluate_packed_hardcall_filter;
use crate::matrix::{finish_variant_major_dense_matrix, VariantMajorDenseParts};
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
    empty_dense_for_samples, empty_sparse_for_samples, expand_selected_samples_to_haplotypes,
    require_pvar, variant_output_capacity, Plink2ReadContext,
};

/// Read retained explicit-phased PLINK2 hard calls as dense haplotype rows.
pub fn read_plink2_haplotypes_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected,
    } = Plink2ReadContext::new(pgen, psam, requested_samples)?;
    let mut diagnostics = selection.diagnostics.clone();
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        require_pvar(pvar)?;
        return empty_dense_for_samples(haplotype_samples, diagnostics, matrix_only);
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let mut haplotype_state = PgenHaplotypeDecodeState::default();
    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut variant_major_missing = Vec::with_capacity(n_haplotypes * output_variant_capacity);
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
                true,
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
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&haplotype_state.selected_haplotype_values);
        variant_major_missing.extend_from_slice(&haplotype_state.selected_haplotype_missing);
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
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: haplotype_samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

/// Read retained explicit-phased PLINK2 dosages as dense haplotype rows.
pub fn read_plink2_haplotypes_dosage_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected: _,
    } = Plink2ReadContext::new(pgen, psam, requested_samples)?;
    let mut diagnostics = selection.diagnostics.clone();
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        require_pvar(pvar)?;
        return empty_dense_for_samples(haplotype_samples, diagnostics, matrix_only);
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let mut haplotype_state = PgenHaplotypeDecodeState::default();
    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut variant_major_missing = Vec::with_capacity(n_haplotypes * output_variant_capacity);
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
                &haplotype_state.selected_collapsed_missing,
                filter,
                &variant,
                !matrix_only,
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
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&haplotype_state.selected_haplotype_values);
        variant_major_missing.extend_from_slice(&haplotype_state.selected_haplotype_missing);
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
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: haplotype_samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

/// Read all retained explicit-phased PLINK2 hard calls as sparse haplotype CSC.
pub fn read_plink2_haplotypes_sparse(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_plink2_haplotypes_sparse_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        None,
    )
}

/// Read retained explicit-phased PLINK2 hard calls as sparse haplotype CSC.
pub fn read_plink2_haplotypes_sparse_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
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
        return empty_sparse_for_samples(haplotype_samples, diagnostics);
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let mut haplotype_state = PgenHaplotypeDecodeState::default();
    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let mut variants = Vec::with_capacity(output_variant_capacity);
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
                true,
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
            &haplotype_state.selected_haplotype_missing,
        )?;
        variants.push(variant);
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_cols = variants.len();
    diagnostics.retained_variants = n_cols;

    SparseGenotypeMatrix::new(
        n_rows,
        n_cols,
        indptr,
        indices,
        data,
        haplotype_samples,
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
    missing: &[bool],
) -> Result<()> {
    reject_sparse_missing_values(missing)?;
    for (row, &value) in values.iter().enumerate() {
        if value != 0.0 {
            indices.push(row);
            data.push(value);
        }
    }
    indptr.push(indices.len());
    Ok(())
}
