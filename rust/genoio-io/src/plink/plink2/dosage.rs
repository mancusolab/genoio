// pattern: Imperative Shell
//! Dense dosage PLINK2 read orchestration.
//!
//! Dosage reads start from hard-call main-track values and apply PGEN dosage
//! overlays before filter evaluation. Output is staged variant-major because
//! dosage overlays naturally decode one retained variant at a time.

use std::path::Path;

use genoio_core::{
    attach_variant_stats, DenseGenotypeMatrixArrowVariants, DenseLayout, DenseMissingPolicy,
    PartialFilterDecision, SampleMetadataArrowBuffers, VariantFilter, VariantMetadataArrowBuffers,
    VariantWindow,
};

use crate::error::Result;
use crate::matrix::apply_dense_missing_policy_to_variant;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::evaluate_dosage_filter;
use super::metadata::PvarRecordReader;
use super::pgen::{open_pgen_payload, read_plink2_variant_dosage, PgenDecoderState};
use super::require_genotype_decision_filter;
use super::source::{
    empty_dense_arrow_for_samples, require_pvar, variant_output_capacity, Plink2ReadContext,
};

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors dense dosage read options plus metadata return choices"
)]
pub fn read_plink2_dosage_dense_windowed_with_arrow_variants(
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
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        require_pvar(pvar)?;
        return empty_dense_arrow_for_samples(
            selection.samples,
            diagnostics,
            return_samples,
            return_variants,
        );
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let mut output_variant_count = 0_usize;

    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
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

        read_plink2_variant_dosage(
            pgen,
            &mut file,
            &header,
            variant_index,
            &selection.source_indices,
            &mut decoder_state,
        )?;
        if matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
            let filter = require_genotype_decision_filter(variant_filter)?;
            let (retain_variant, stats) = evaluate_dosage_filter(
                &decoder_state.values,
                &decoder_state.missing_indices,
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
            &mut decoder_state.values,
            &decoder_state.missing_indices,
            missing_policy,
        )?;
        variant_major_values.extend_from_slice(&decoder_state.values);
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = selection.samples.len();
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    let samples = SampleMetadataArrowBuffers::optional_from_records(
        &selection.samples,
        return_samples,
        false,
    )?;
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
