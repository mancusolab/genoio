// pattern: Imperative Shell
//! BGEN diploid dosage dense read orchestration.
//!
//! The loop combines metadata filtering, optional `.bgi` region cursors, dosage
//! decode buffers, and retained-window accounting for sample-by-variant output.

use std::path::Path;

use genoio_core::{
    select_samples_source_order, DenseGenotypeMatrix, DenseLayout, DenseMissingPolicy, GenoioError,
    PartialFilterDecision, SampleMetadataBuffers, VariantFilter, VariantMetadataBuffers,
    VariantWindow,
};

use crate::blocks::{
    block_diagnostics_snapshot, checked_dense_block_len, BlockReadOptions, DosageSource, MatrixKind,
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
use super::session::{
    BgenBlockSession, BgenIndexedReadContext, BgenReadSession, BgenVariantCursor,
};

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors dense dosage read options plus metadata return choices"
)]
pub fn read_bgen_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    if let Some(index_records) = indexed_region_records(bgen, variant_filter)? {
        let mut session = BgenReadSession::open(bgen)?;
        let all_samples = session.read_samples(sample)?;
        let selection = select_samples_source_order(&all_samples, requested_samples, bgen)?;
        let diagnostics = selection.diagnostics.clone();
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
    let retained_skip = variant_window.map_or(0, |window| window.start);
    let options = BlockReadOptions {
        matrix_kind: MatrixKind::Genotype,
        sparse: false,
        requested_samples: requested_samples.map(<[String]>::to_vec),
        variant_filter: variant_filter.cloned(),
        dosage_source: DosageSource::Dosage,
        missing_policy,
        return_samples,
        return_variants,
    };
    let mut session = BgenBlockSession::open_windowed(
        bgen.to_path_buf(),
        sample.map(Path::to_path_buf),
        options,
        retained_skip,
    )?;
    let block_size = match variant_window {
        Some(window) => window.len,
        None => usize::try_from(session.remaining_variants)
            .map_err(|_| GenoioError::invalid_source(bgen, "bgen variant count is out of range"))?,
    };
    match session.next_dosage_block(block_size)? {
        Some(matrix) => Ok(matrix),
        None => session.empty_genotype_output(),
    }
}

impl BgenBlockSession {
    pub(super) fn next_dosage_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        if self.eof || block_size == 0 {
            return Ok(None);
        }

        let n_samples = self.selection.samples.len();
        let output_len = checked_dense_block_len(n_samples, block_size)?;
        self.record_dense_allocation(output_len);
        let mut values = vec![0.0; output_len];
        let mut variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(block_size));
        let retained_skip = std::mem::take(&mut self.retained_skip);
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: retained_skip,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;
        let matrix_only = !self.return_samples && !self.return_variants;
        let sample_count = self.io.header.sample_count;

        while !retention.window_is_satisfied() {
            if self.remaining_variants == 0 {
                self.eof = true;
                break;
            }
            self.remaining_variants -= 1;
            self.record_candidate_visit();

            if matrix_only && self.variant_filter.is_none() {
                match retention
                    .metadata_decision(PartialFilterDecision::Accept, &mut self.diagnostics)
                {
                    MetadataRetentionAction::Include => {
                        self.io.skip_variant()?;
                        self.io
                            .read_payload_into(&mut self.decode_buffers.probability)?;
                        self.record_payload_decode();
                        write_dosage_slot(
                            BgenDosageSlotWrite {
                                bgen: &self.io.bgen,
                                sample_count,
                                source_indices: &self.selection.source_indices,
                                buffers: &mut self.decode_buffers,
                                values: &mut values,
                                row_width: block_size,
                                variant_index: output_variant_count,
                                missing_policy: self.missing_policy,
                            },
                            false,
                        )?;
                        output_variant_count += 1;
                    }
                    MetadataRetentionAction::Skip => {
                        self.io.skip_variant()?;
                        self.io.skip_payload()?;
                    }
                    MetadataRetentionAction::Stop => {
                        self.io.skip_variant()?;
                        self.io.skip_payload()?;
                        break;
                    }
                    MetadataRetentionAction::DecodeGenotypes => {
                        return Err(GenoioError::internal_contract(
                            "unfiltered bgen matrix-only path requested genotype filtering",
                        ));
                    }
                }
                continue;
            }

            let mut variant = self.io.read_variant()?;
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => {
                    self.io.skip_payload()?;
                    continue;
                }
                MetadataRetentionAction::Stop => {
                    self.io.skip_payload()?;
                    break;
                }
                MetadataRetentionAction::Include => {
                    self.io
                        .read_payload_into(&mut self.decode_buffers.probability)?;
                    self.record_payload_decode();
                }
                MetadataRetentionAction::DecodeGenotypes => {
                    self.io
                        .read_payload_into(&mut self.decode_buffers.probability)?;
                    self.record_payload_decode();
                    let filter = self.variant_filter.as_ref().ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?;
                    let (retain_variant, stats) = decode_and_evaluate_dosage_filter(
                        &self.io.bgen,
                        sample_count,
                        &self.selection.source_indices,
                        &mut self.decode_buffers,
                        filter,
                        &variant,
                        !self.return_variants,
                    )?;
                    match apply_genotype_filter_result(
                        &mut retention,
                        &mut self.diagnostics,
                        &mut variant,
                        retain_variant,
                        stats,
                    ) {
                        RetentionAction::Include => {}
                        RetentionAction::Skip => continue,
                        RetentionAction::Stop => break,
                    }
                }
            }

            write_dosage_slot(
                BgenDosageSlotWrite {
                    bgen: &self.io.bgen,
                    sample_count,
                    source_indices: &self.selection.source_indices,
                    buffers: &mut self.decode_buffers,
                    values: &mut values,
                    row_width: block_size,
                    variant_index: output_variant_count,
                    missing_policy: self.missing_policy,
                },
                matches!(partial_decision, PartialFilterDecision::NeedGenotypes),
            )?;
            if let Some(variants) = variants.as_mut() {
                variants.push_record(&variant)?;
            }
            output_variant_count += 1;
        }

        if output_variant_count == 0 {
            return Ok(None);
        }
        shrink_sample_major_width(&mut values, n_samples, block_size, output_variant_count);
        let samples = SampleMetadataBuffers::optional_from_records(
            &self.selection.samples,
            self.return_samples,
            false,
        )?;
        let diagnostics = block_diagnostics_snapshot(&self.diagnostics, output_variant_count);
        DenseGenotypeMatrix::new_with_layout(
            n_samples,
            output_variant_count,
            values,
            DenseLayout::SampleMajor,
            samples,
            variants,
            diagnostics,
        )
        .map(Some)
    }
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
) -> Result<DenseGenotypeMatrix> {
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
    let bgen = session.bgen.clone();
    let sample_count = session.header.sample_count;
    let output_variant_capacity = variant_window.map_or(index_records.len(), |window| {
        window.len.min(index_records.len())
    });
    let n_samples = selection.samples.len();
    let mut values = sample_major_buffer(n_samples, output_variant_capacity)?;
    let mut variants =
        return_variants.then(|| VariantMetadataBuffers::with_capacity(output_variant_capacity));
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
                    &bgen,
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
                bgen: &bgen,
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
    let samples =
        SampleMetadataBuffers::optional_from_records(&selection.samples, return_samples, false)?;
    DenseGenotypeMatrix::new_with_layout(
        n_samples,
        n_variants,
        values,
        DenseLayout::SampleMajor,
        samples,
        variants,
        diagnostics,
    )
}
