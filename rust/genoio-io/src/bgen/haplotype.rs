// pattern: Imperative Shell
//! BGEN phased haplotype dosage dense read orchestration.
//!
//! The reader expands phased probabilities into two haplotype rows per selected
//! sample. Genotype-stat filters use collapsed diploid dosages, while retained
//! output preserves haplotype order.

use std::path::Path;

use genoio_core::{
    DenseGenotypeMatrix, DenseLayout, DenseMissingPolicy, DenseSampleSelection, GenoioError,
    PartialFilterDecision, SampleMetadataBuffers, SampleRecord, VariantFilter,
    VariantMetadataBuffers, VariantWindow,
};

use crate::blocks::{
    block_diagnostics_snapshot, checked_dense_block_len, BlockReadOptions, DosageSource, MatrixKind,
};
use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::matrix::apply_dense_missing_policy_to_variant;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::decode::{decode_buffered_haplotype_values, HaplotypeDecodeBuffers};
use super::filter::apply_genotype_filter_result;
use super::session::BgenBlockSession;

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors haplotype dosage read options plus metadata return choices"
)]
pub fn read_bgen_haplotypes_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    let retained_skip = variant_window.map_or(0, |window| window.start);
    let options = BlockReadOptions {
        matrix_kind: MatrixKind::Haplotype,
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
        None => session.source_record_capacity(),
    };
    match session.next_haplotype_block(block_size)? {
        Some(matrix) => Ok(matrix),
        None => session.empty_haplotype_output(),
    }
}

impl BgenBlockSession {
    pub(super) fn next_haplotype_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        if self.eof || block_size == 0 {
            return Ok(None);
        }

        let n_haplotypes = self.selection.samples.len().checked_mul(2).ok_or_else(|| {
            GenoioError::internal_contract("bgen haplotype row count is out of range")
        })?;
        let output_len = checked_dense_block_len(n_haplotypes, block_size)?;
        self.record_dense_allocation(output_len);
        let mut variant_major_values = Vec::with_capacity(output_len);
        let mut variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(block_size));
        let retained_skip = std::mem::take(&mut self.retained_skip);
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: retained_skip,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;
        let sample_count = self.io.header.sample_count;

        while !retention.window_is_satisfied() {
            let Some(position) = self.next_position()? else {
                break;
            };
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
                    position.validate_if_indexed(&mut self.io)?;
                    continue;
                }
                MetadataRetentionAction::Stop => {
                    self.io.skip_payload()?;
                    position.validate_if_indexed(&mut self.io)?;
                    break;
                }
                MetadataRetentionAction::Include => {
                    self.io.read_payload_into(
                        &mut haplotype_buffers_mut(&mut self.haplotype_buffers)?.probability,
                    )?;
                    position.validate_if_indexed(&mut self.io)?;
                    self.record_payload_decode();
                }
                MetadataRetentionAction::DecodeGenotypes => {
                    self.io.read_payload_into(
                        &mut haplotype_buffers_mut(&mut self.haplotype_buffers)?.probability,
                    )?;
                    position.validate_if_indexed(&mut self.io)?;
                    self.record_payload_decode();
                    let buffers = haplotype_buffers_mut(&mut self.haplotype_buffers)?;
                    decode_buffered_haplotype_values(
                        &self.io.bgen,
                        sample_count,
                        &self.selection.source_indices,
                        buffers,
                    )?;
                    let filter = self.variant_filter.as_ref().ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?;
                    let (retain_variant, stats) = evaluate_dosage_filter(
                        &buffers.selected_collapsed_values,
                        &buffers.selected_collapsed_missing_indices,
                        filter,
                        &variant,
                        self.return_variants,
                    )?;
                    match apply_genotype_filter_result(
                        &mut retention,
                        &mut self.diagnostics,
                        &mut variant,
                        retain_variant,
                        stats,
                    ) {
                        RetentionAction::Include => {}
                        RetentionAction::Skip => {
                            continue;
                        }
                        RetentionAction::Stop => {
                            break;
                        }
                    }
                }
            }

            if !matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
                let buffers = haplotype_buffers_mut(&mut self.haplotype_buffers)?;
                decode_buffered_haplotype_values(
                    &self.io.bgen,
                    sample_count,
                    &self.selection.source_indices,
                    buffers,
                )?;
            }
            let buffers = haplotype_buffers_mut(&mut self.haplotype_buffers)?;
            apply_dense_missing_policy_to_variant(
                &mut buffers.selected_haplotype_values,
                &buffers.selected_haplotype_missing_indices,
                self.missing_policy,
            )?;
            variant_major_values.extend_from_slice(&buffers.selected_haplotype_values);
            if let Some(variants) = variants.as_mut() {
                variants.push_record(&variant)?;
            }
            output_variant_count += 1;
        }

        if output_variant_count == 0 {
            return Ok(None);
        }
        let haplotype_samples = expand_selected_samples_to_haplotypes(&self.selection);
        let samples = SampleMetadataBuffers::optional_from_records(
            &haplotype_samples,
            self.return_samples,
            true,
        )?;
        let diagnostics = block_diagnostics_snapshot(&self.diagnostics, output_variant_count);
        DenseGenotypeMatrix::new_with_layout(
            n_haplotypes,
            output_variant_count,
            variant_major_values,
            DenseLayout::VariantMajor,
            samples,
            variants,
            diagnostics,
        )
        .map(Some)
    }
}

fn haplotype_buffers_mut(
    buffers: &mut Option<HaplotypeDecodeBuffers>,
) -> Result<&mut HaplotypeDecodeBuffers> {
    buffers.as_mut().ok_or_else(|| {
        GenoioError::internal_contract("bgen haplotype session is missing haplotype decode buffers")
    })
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
