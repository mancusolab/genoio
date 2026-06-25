// pattern: Imperative Shell
//! BGEN phased haplotype dosage dense read orchestration.
//!
//! The reader expands phased probabilities into two haplotype rows per selected
//! sample. Genotype-stat filters use collapsed diploid dosages, while retained
//! output preserves haplotype order.

use std::path::Path;

use genoio_core::{
    select_samples_source_order, DenseGenotypeMatrix, DenseGenotypeMatrixArrowVariants,
    DenseLayout, DenseMissingPolicy, DenseSampleSelection, GenoioError, PartialFilterDecision,
    SampleRecord, VariantFilter, VariantMetadataArrowBuffers, VariantWindow,
};

use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::matrix::apply_dense_missing_policy_to_variant;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::decode::{decode_buffered_haplotype_values, HaplotypeDecodeBuffers};
use super::filter::apply_genotype_filter_result;
use super::index::{indexed_region_records, BgenIndexRecord};
use super::session::{BgenIndexedReadContext, BgenReadSession, BgenVariantCursor};

fn dense_arrow_output_to_rows(
    output: DenseGenotypeMatrixArrowVariants,
    context: &'static str,
) -> Result<DenseGenotypeMatrix> {
    output.into_matrix().map_err(|error| {
        GenoioError::internal_contract(format!(
            "BGEN {context} Arrow-to-row compatibility conversion failed: {error}"
        ))
    })
}

fn empty_dense_arrow_for_samples(
    samples: Vec<SampleRecord>,
    mut diagnostics: genoio_core::DenseDiagnostics,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    diagnostics.retained_variants = 0;
    let n_samples = samples.len();
    let samples = if return_samples { samples } else { Vec::new() };
    let variants = return_variants.then(|| VariantMetadataArrowBuffers::with_capacity(0));
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        0,
        Vec::new(),
        DenseLayout::SampleMajor,
        samples,
        variants,
        diagnostics,
    )
}

pub fn read_bgen_haplotypes_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_bgen_haplotypes_dosage_dense_windowed_with_missing_policy(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        matrix_only,
    )
}

pub fn read_bgen_haplotypes_dosage_dense_windowed_with_missing_policy(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_bgen_haplotypes_dosage_dense_windowed_with_arrow_variants(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        !matrix_only,
        !matrix_only,
    )
    .and_then(|output| dense_arrow_output_to_rows(output, "haplotype dosage dense"))
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors haplotype dosage read options plus metadata return choices"
)]
pub fn read_bgen_haplotypes_dosage_dense_windowed_with_arrow_variants(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let mut session = BgenReadSession::open(bgen)?;
    let all_samples = session.read_samples(sample)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bgen)?;
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    let mut diagnostics = selection.diagnostics.clone();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        diagnostics.retained_variants = 0;
        return empty_dense_arrow_for_samples(
            haplotype_samples,
            diagnostics,
            return_samples,
            return_variants,
        );
    }
    if let Some(index_records) = indexed_region_records(bgen, variant_filter)? {
        let context = BgenIndexedReadContext {
            session: &mut session,
            selection,
            diagnostics,
            variant_filter,
            variant_window,
            missing_policy,
            return_samples,
            return_variants,
        };
        return read_bgen_haplotypes_dosage_dense_indexed(context, &index_records);
    }

    session.seek_to_variants()?;

    let header_variant_count = usize::try_from(session.header.variant_count)
        .map_err(|_| GenoioError::invalid_source(bgen, "bgen variant count is out of range"))?;
    let output_variant_capacity = variant_window.map_or(header_variant_count, |window| {
        window.len.min(header_variant_count)
    });
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut decode_buffers = HaplotypeDecodeBuffers::default();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;

    let variant_count = session.header.variant_count;
    let sample_count = session.header.sample_count;
    let mut cursor = BgenVariantCursor::sequential(variant_count);
    loop {
        if retention.window_is_satisfied() {
            break;
        }
        if cursor.next(&mut session)?.is_none() {
            break;
        }
        let mut variant = session.read_variant()?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Skip => {
                session.skip_payload()?;
                continue;
            }
            MetadataRetentionAction::Stop => break,
            MetadataRetentionAction::Include => {
                session.read_payload_into(&mut decode_buffers.probability)?;
            }
            MetadataRetentionAction::DecodeGenotypes => {
                session.read_payload_into(&mut decode_buffers.probability)?;
                decode_buffered_haplotype_values(
                    bgen,
                    sample_count,
                    &selection.source_indices,
                    &mut decode_buffers,
                )?;
                let (retain_variant, stats) = evaluate_dosage_filter(
                    &decode_buffers.selected_collapsed_values,
                    &decode_buffers.selected_collapsed_missing_indices,
                    variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?,
                    &variant,
                    return_variants,
                )?;
                match apply_genotype_filter_result(
                    &mut retention,
                    &mut diagnostics,
                    &mut variant,
                    retain_variant,
                    stats,
                ) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => {
                        break;
                    }
                }
            }
        }

        if !matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
            decode_buffered_haplotype_values(
                bgen,
                sample_count,
                &selection.source_indices,
                &mut decode_buffers,
            )?;
        }
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        apply_dense_missing_policy_to_variant(
            &mut decode_buffers.selected_haplotype_values,
            &decode_buffers.selected_haplotype_missing_indices,
            missing_policy,
        )?;
        variant_major_values.extend_from_slice(&decode_buffers.selected_haplotype_values);
        output_variant_count += 1;
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

pub(super) fn expand_selected_samples_to_haplotypes(
    selection: &DenseSampleSelection,
) -> Vec<SampleRecord> {
    let mut haplotype_samples = Vec::with_capacity(selection.samples.len() * 2);
    for (sample, &source_index) in selection.samples.iter().zip(&selection.source_indices) {
        for haplotype_index in 0..2 {
            let mut haplotype_sample = sample.clone();
            haplotype_sample.source_sample_index = Some(source_index);
            haplotype_sample.haplotype_index = Some(haplotype_index);
            haplotype_samples.push(haplotype_sample);
        }
    }
    haplotype_samples
}

fn read_bgen_haplotypes_dosage_dense_indexed(
    context: BgenIndexedReadContext<'_>,
    index_records: &[BgenIndexRecord],
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let BgenIndexedReadContext {
        session,
        selection,
        mut diagnostics,
        variant_filter,
        variant_window,
        missing_policy,
        return_samples,
        return_variants,
    } = context;
    let bgen = session.bgen;
    let sample_count = session.header.sample_count;
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    let output_variant_capacity = variant_window.map_or(index_records.len(), |window| {
        window.len.min(index_records.len())
    });
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut decode_buffers = HaplotypeDecodeBuffers::default();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;

    let mut cursor = BgenVariantCursor::indexed(index_records);
    loop {
        if retention.window_is_satisfied() {
            break;
        }
        let Some(position) = cursor.next(session)? else {
            break;
        };
        let mut variant = session.read_variant()?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Skip => {
                session.skip_payload()?;
                position.validate_if_indexed(session)?;
                continue;
            }
            MetadataRetentionAction::Stop => {
                session.skip_payload()?;
                position.validate_if_indexed(session)?;
                break;
            }
            MetadataRetentionAction::Include => {
                session.read_payload_into(&mut decode_buffers.probability)?;
            }
            MetadataRetentionAction::DecodeGenotypes => {
                session.read_payload_into(&mut decode_buffers.probability)?;
                decode_buffered_haplotype_values(
                    bgen,
                    sample_count,
                    &selection.source_indices,
                    &mut decode_buffers,
                )?;
                let (retain_variant, stats) = evaluate_dosage_filter(
                    &decode_buffers.selected_collapsed_values,
                    &decode_buffers.selected_collapsed_missing_indices,
                    variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?,
                    &variant,
                    return_variants,
                )?;
                match apply_genotype_filter_result(
                    &mut retention,
                    &mut diagnostics,
                    &mut variant,
                    retain_variant,
                    stats,
                ) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => {
                        position.validate_if_indexed(session)?;
                        continue;
                    }
                    RetentionAction::Stop => {
                        position.validate_if_indexed(session)?;
                        break;
                    }
                }
            }
        }

        if !matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
            decode_buffered_haplotype_values(
                bgen,
                sample_count,
                &selection.source_indices,
                &mut decode_buffers,
            )?;
        }
        position.validate_if_indexed(session)?;
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        apply_dense_missing_policy_to_variant(
            &mut decode_buffers.selected_haplotype_values,
            &decode_buffers.selected_haplotype_missing_indices,
            missing_policy,
        )?;
        variant_major_values.extend_from_slice(&decode_buffers.selected_haplotype_values);
        output_variant_count += 1;
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
