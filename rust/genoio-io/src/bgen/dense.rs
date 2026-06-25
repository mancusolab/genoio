// pattern: Imperative Shell
//! BGEN diploid dosage dense read orchestration.
//!
//! The loop combines metadata filtering, optional `.bgi` region cursors, dosage
//! decode buffers, and retained-window accounting for sample-by-variant output.

use std::path::Path;

use genoio_core::{
    select_samples_source_order, DenseGenotypeMatrix, DenseGenotypeMatrixArrowVariants,
    DenseLayout, DenseMissingPolicy, DenseSampleSelection, GenoioError, PartialFilterDecision,
    VariantFilter, VariantMetadataArrowBuffers, VariantWindow,
};

use crate::error::Result;
use crate::matrix::{
    apply_dense_missing_policy_to_variant, shrink_sample_major_width,
    write_sample_major_variant_slot,
};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::decode::{
    decode_buffered_dosage_values, try_decode_buffered_dosage_values_into_sample_major_slot,
    DosageDecodeBuffers, SampleMajorSlotMut,
};
use super::filter::{apply_genotype_filter_result, decode_and_evaluate_dosage_filter};
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
    samples: Vec<genoio_core::SampleRecord>,
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

pub fn read_bgen_dosage_dense(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_bgen_dosage_dense_windowed_with_missing_policy(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        None,
        DenseMissingPolicy::Nan,
        false,
    )
}

/// Read retained BGEN biallelic diploid dosages as a dense matrix.
pub fn read_bgen_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_bgen_dosage_dense_windowed_with_missing_policy(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        matrix_only,
    )
}

pub fn read_bgen_dosage_dense_windowed_with_missing_policy(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_bgen_dosage_dense_windowed_with_arrow_variants(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        !matrix_only,
        !matrix_only,
    )
    .and_then(|output| dense_arrow_output_to_rows(output, "dosage dense"))
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors dense dosage read options plus metadata return choices"
)]
pub fn read_bgen_dosage_dense_windowed_with_arrow_variants(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let matrix_only = !return_samples && !return_variants;
    let mut session = BgenReadSession::open(bgen)?;
    let all_samples = session.read_samples(sample)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bgen)?;
    let mut diagnostics = selection.diagnostics.clone();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        diagnostics.retained_variants = 0;
        return empty_dense_arrow_for_samples(
            selection.samples,
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
        return read_bgen_dosage_dense_indexed(context, &index_records);
    }
    if matrix_only && variant_filter.is_none() {
        // Matrix-only reads do not expose variant strings or positions. Skipping
        // those bytes preserves the same matrix contract while avoiding string
        // allocation and UTF-8 validation on the hot path.
        return read_bgen_dosage_dense_matrix_only_unfiltered(
            &mut session,
            selection,
            diagnostics,
            variant_window,
            missing_policy,
        );
    }

    session.seek_to_variants()?;

    let header_variant_count = usize::try_from(session.header.variant_count)
        .map_err(|_| GenoioError::invalid_source(bgen, "bgen variant count is out of range"))?;
    let output_variant_capacity = variant_window.map_or(header_variant_count, |window| {
        window.len.min(header_variant_count)
    });
    let n_samples = selection.samples.len();
    let mut values = sample_major_buffer(n_samples, output_variant_capacity)?;
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut decode_buffers = DosageDecodeBuffers::default();
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
                let filter = variant_filter.ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = decode_and_evaluate_dosage_filter(
                    bgen,
                    sample_count,
                    &selection.source_indices,
                    &mut decode_buffers,
                    filter,
                    &variant,
                    !return_variants,
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

        write_dosage_slot(
            BgenDosageSlotWrite {
                bgen,
                sample_count,
                source_indices: &selection.source_indices,
                buffers: &mut decode_buffers,
                values: &mut values,
                row_width: output_variant_capacity,
                variant_index: output_variant_count,
                missing_policy,
            },
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes),
        )?;
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        output_variant_count += 1;
    }

    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    shrink_sample_major_width(&mut values, n_samples, output_variant_capacity, n_variants);
    let samples = if return_samples {
        selection.samples
    } else {
        Vec::new()
    };
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        n_variants,
        values,
        DenseLayout::SampleMajor,
        samples,
        variants,
        diagnostics,
    )
}

fn read_bgen_dosage_dense_matrix_only_unfiltered(
    session: &mut BgenReadSession<'_>,
    selection: DenseSampleSelection,
    mut diagnostics: genoio_core::DenseDiagnostics,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    session.seek_to_variants()?;

    let header_variant_count = usize::try_from(session.header.variant_count).map_err(|_| {
        GenoioError::invalid_source(session.bgen, "bgen variant count is out of range")
    })?;
    let output_variant_capacity = variant_window.map_or(header_variant_count, |window| {
        window.len.min(header_variant_count)
    });
    let n_samples = selection.samples.len();
    let mut values = sample_major_buffer(n_samples, output_variant_capacity)?;
    let mut decode_buffers = DosageDecodeBuffers::default();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;

    let variant_count = session.header.variant_count;
    let sample_count = session.header.sample_count;
    let mut cursor = BgenVariantCursor::sequential(variant_count);
    loop {
        if retention.window_is_satisfied() {
            break;
        }
        if cursor.next(session)?.is_none() {
            break;
        }
        match retention.metadata_decision(PartialFilterDecision::Accept, &mut diagnostics) {
            MetadataRetentionAction::Include => {
                session.skip_variant()?;
                session.read_payload_into(&mut decode_buffers.probability)?;
                write_dosage_slot(
                    BgenDosageSlotWrite {
                        bgen: session.bgen,
                        sample_count,
                        source_indices: &selection.source_indices,
                        buffers: &mut decode_buffers,
                        values: &mut values,
                        row_width: output_variant_capacity,
                        variant_index: output_variant_count,
                        missing_policy,
                    },
                    false,
                )?;
                output_variant_count += 1;
            }
            MetadataRetentionAction::Skip => {
                session.skip_variant()?;
                session.skip_payload()?;
            }
            MetadataRetentionAction::Stop => break,
            MetadataRetentionAction::DecodeGenotypes => {
                return Err(GenoioError::internal_contract(
                    "unfiltered bgen matrix-only path requested genotype filtering",
                ));
            }
        }
    }

    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    shrink_sample_major_width(&mut values, n_samples, output_variant_capacity, n_variants);
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        n_variants,
        values,
        DenseLayout::SampleMajor,
        Vec::new(),
        None,
        diagnostics,
    )
}

/// Read retained BGEN biallelic diploid phased dosages as dense haplotype rows.
fn sample_major_buffer(n_samples: usize, n_variants: usize) -> Result<Vec<f32>> {
    let len = n_samples.checked_mul(n_variants).ok_or_else(|| {
        GenoioError::internal_contract("sample-major dense matrix shape is out of range")
    })?;
    Ok(vec![0.0; len])
}

struct BgenDosageSlotWrite<'a> {
    bgen: &'a Path,
    sample_count: u32,
    source_indices: &'a [usize],
    buffers: &'a mut DosageDecodeBuffers,
    values: &'a mut [f32],
    row_width: usize,
    variant_index: usize,
    missing_policy: DenseMissingPolicy,
}

fn write_dosage_slot(
    request: BgenDosageSlotWrite<'_>,
    already_decoded_for_filter: bool,
) -> Result<()> {
    let BgenDosageSlotWrite {
        bgen,
        sample_count,
        source_indices,
        buffers,
        values,
        row_width,
        variant_index,
        missing_policy,
    } = request;

    if !already_decoded_for_filter {
        let mut slot = SampleMajorSlotMut {
            values,
            row_width,
            variant_index,
        };
        // Common byte-aligned BGEN records can fill the final matrix slot
        // directly. Other shapes fall back to the generic selected decode.
        if try_decode_buffered_dosage_values_into_sample_major_slot(
            bgen,
            sample_count,
            source_indices,
            buffers,
            &mut slot,
        )? {
            return Ok(());
        }
        decode_buffered_dosage_values(bgen, sample_count, source_indices, buffers)?;
    }

    apply_dense_missing_policy_to_variant(
        &mut buffers.selected_values,
        &buffers.selected_missing_indices,
        missing_policy,
    )?;
    write_sample_major_variant_slot(
        values,
        source_indices.len(),
        row_width,
        variant_index,
        &buffers.selected_values,
    )
}

fn read_bgen_dosage_dense_indexed(
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
    let output_variant_capacity = variant_window.map_or(index_records.len(), |window| {
        window.len.min(index_records.len())
    });
    let n_samples = selection.samples.len();
    let mut values = sample_major_buffer(n_samples, output_variant_capacity)?;
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut decode_buffers = DosageDecodeBuffers::default();
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
                let filter = variant_filter.ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = decode_and_evaluate_dosage_filter(
                    bgen,
                    sample_count,
                    &selection.source_indices,
                    &mut decode_buffers,
                    filter,
                    &variant,
                    !return_variants,
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

        write_dosage_slot(
            BgenDosageSlotWrite {
                bgen,
                sample_count,
                source_indices: &selection.source_indices,
                buffers: &mut decode_buffers,
                values: &mut values,
                row_width: output_variant_capacity,
                variant_index: output_variant_count,
                missing_policy,
            },
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes),
        )?;
        position.validate_if_indexed(session)?;
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        output_variant_count += 1;
    }

    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    shrink_sample_major_width(&mut values, n_samples, output_variant_capacity, n_variants);
    let samples = if return_samples {
        selection.samples
    } else {
        Vec::new()
    };
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        n_variants,
        values,
        DenseLayout::SampleMajor,
        samples,
        variants,
        diagnostics,
    )
}
