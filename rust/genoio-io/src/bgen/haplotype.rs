// pattern: Imperative Shell
//! BGEN phased haplotype dosage dense read orchestration.

use std::path::Path;

use genoio_core::{
    select_samples_source_order, DenseGenotypeMatrix, DenseSampleSelection, GenoioError,
    PartialFilterDecision, SampleRecord, VariantFilter, VariantWindow,
};

use crate::error::Result;
use crate::matrix::{
    finish_dense_matrix, finish_variant_major_dense_matrix, DenseMatrixParts,
    VariantMajorDenseParts,
};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::decode::{decode_buffered_haplotype_values, HaplotypeDecodeBuffers};
use super::filter::{apply_genotype_filter_result, evaluate_dosage_filter};
use super::index::{indexed_region_records, BgenIndexRecord};
use super::session::{BgenIndexedReadContext, BgenReadSession, BgenVariantCursor};

pub fn read_bgen_haplotypes_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let mut session = BgenReadSession::open(bgen)?;
    let all_samples = session.read_samples(sample)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bgen)?;
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    let mut diagnostics = selection.diagnostics.clone();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        diagnostics.retained_variants = 0;
        return finish_dense_matrix(
            DenseMatrixParts {
                n_samples: haplotype_samples.len(),
                n_variants: 0,
                values: Vec::new(),
                missing_mask: Vec::new(),
                samples: haplotype_samples,
                variants: Vec::new(),
                diagnostics,
            },
            matrix_only,
        );
    }
    if let Some(index_records) = indexed_region_records(bgen, variant_filter)? {
        let context = BgenIndexedReadContext {
            session: &mut session,
            selection,
            diagnostics,
            variant_filter,
            variant_window,
            matrix_only,
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
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut variant_major_missing = Vec::with_capacity(n_haplotypes * output_variant_capacity);
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
                    &decode_buffers.selected_collapsed_missing,
                    variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?,
                    &variant,
                    !matrix_only,
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
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&decode_buffers.selected_haplotype_values);
        variant_major_missing.extend_from_slice(&decode_buffers.selected_haplotype_missing);
        output_variant_count += 1;
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
) -> Result<DenseGenotypeMatrix> {
    let BgenIndexedReadContext {
        session,
        selection,
        mut diagnostics,
        variant_filter,
        variant_window,
        matrix_only,
    } = context;
    let bgen = session.bgen;
    let sample_count = session.header.sample_count;
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    let output_variant_capacity = variant_window.map_or(index_records.len(), |window| {
        window.len.min(index_records.len())
    });
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut variant_major_missing = Vec::with_capacity(n_haplotypes * output_variant_capacity);
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
                    &decode_buffers.selected_collapsed_missing,
                    variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?,
                    &variant,
                    !matrix_only,
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
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&decode_buffers.selected_haplotype_values);
        variant_major_missing.extend_from_slice(&decode_buffers.selected_haplotype_missing);
        output_variant_count += 1;
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
