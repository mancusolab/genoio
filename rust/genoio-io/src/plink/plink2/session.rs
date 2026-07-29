// pattern: Imperative Shell
//! Persistent PLINK2 PGEN/PVAR/PSAM block-reader session.
//!
//! The session owns one open PGEN file, one streaming PVAR reader, selected
//! PSAM samples, cumulative diagnostics, and reusable PGEN decoder scratch.
//! Fixed-width records retain direct addressing, while variable-width records
//! advance one shared LD decoder state across block boundaries.

use std::fs::File;
use std::path::PathBuf;

use genoio_core::{
    append_sparse_column, attach_variant_stats, flip_values_to_minor_allele, reject_sparse_missing,
    select_samples_source_order, DenseDiagnostics, DenseGenotypeMatrix, DenseLayout,
    DenseMissingPolicy, DenseSampleSelection, GenoioError, GenotypeFilterPlan,
    PartialFilterDecision, SampleMetadataBuffers, SparseGenotypeMatrix, VariantFilter,
    VariantMetadataBuffers, VariantRecord, VariantWindow,
};

use crate::blocks::{
    block_diagnostics_snapshot, checked_dense_block_len, checked_sparse_indptr_len, BlockOutput,
    BlockReadOptions, DosageSource, MatrixKind,
};
use crate::error::Result;
use crate::hardcall::{
    evaluate_packed_hardcall_filter, flush_hardcall_batch_into_sample_major, HardcallBatch,
};
use crate::matrix::{apply_dense_missing_policy_to_variant, shrink_sample_major_width};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

#[cfg(not(test))]
use super::metadata::{parse_psam, PvarRecordReader};
#[cfg(test)]
use super::metadata::{parse_psam_with_probe, PvarRecordReader};
use super::pgen::{
    decode_plink2_haplotype_dosage_aux, decode_plink2_variant_dosage_aux,
    read_plink2_variant_dosage, read_plink2_variant_dosage_main_track,
    read_plink2_variant_haplotype_dosage_track, read_plink2_variant_packed,
    read_supported_pgen_header_from_file, validate_plink2_sample_count, PgenDecoderState,
    PgenHaplotypeDecodeState, PgenHeader, PgenLayout,
};
use super::source::expand_selected_samples_to_haplotypes;
use super::{evaluate_dosage_filter, require_genotype_decision_filter};

/// Persistent PLINK2 hard-call state over one PGEN/PVAR/PSAM source set.
pub(crate) struct Plink2BlockSession {
    pgen: PathBuf,
    pvar: PathBuf,
    pgen_reader: File,
    pvar_reader: PvarRecordReader,
    header: PgenHeader,
    selection: DenseSampleSelection,
    source_position: usize,
    diagnostics: DenseDiagnostics,
    variant_filter: Option<VariantFilter>,
    genotype_filter_plan: GenotypeFilterPlan,
    matrix_kind: MatrixKind,
    dosage_source: DosageSource,
    missing_policy: DenseMissingPolicy,
    sparse: bool,
    return_samples: bool,
    return_variants: bool,
    all_samples_selected: bool,
    decoder_state: PgenDecoderState,
    haplotype_state: Option<Box<PgenHaplotypeDecodeState>>,
    batch: HardcallBatch,
    selected_values: Vec<f32>,
    missing_indices: Vec<usize>,
    eof: bool,
    failed: bool,
    #[cfg(test)]
    probe: Option<Plink2WorkProbe>,
}

impl Plink2BlockSession {
    pub(crate) fn open(
        pgen: PathBuf,
        pvar: PathBuf,
        psam: PathBuf,
        options: BlockReadOptions,
    ) -> Result<Self> {
        #[cfg(test)]
        {
            Self::open_impl(pgen, pvar, psam, options, None)
        }
        #[cfg(not(test))]
        {
            Self::open_impl(pgen, pvar, psam, options)
        }
    }

    fn open_impl(
        pgen: PathBuf,
        pvar: PathBuf,
        psam: PathBuf,
        options: BlockReadOptions,
        #[cfg(test)] probe: Option<Plink2WorkProbe>,
    ) -> Result<Self> {
        validate_plink2_options(&options)?;

        let mut pgen_reader = File::open(&pgen).map_err(|source| GenoioError::Io {
            path: pgen.clone(),
            source,
        })?;
        #[cfg(test)]
        if let Some(probe) = probe.as_ref() {
            probe.record_pgen_open();
        }
        let header = read_supported_pgen_header_from_file(&pgen, &mut pgen_reader)?;

        #[cfg(test)]
        let all_samples = if let Some(probe) = probe.as_ref() {
            parse_psam_with_probe(&psam, probe)?
        } else {
            super::metadata::parse_psam(&psam)?
        };
        #[cfg(not(test))]
        let all_samples = parse_psam(&psam)?;
        validate_plink2_sample_count(&pgen, &header, all_samples.len())?;
        let selection =
            select_samples_source_order(&all_samples, options.requested_samples.as_deref(), &pgen)?;

        #[cfg(test)]
        let pvar_reader = if let Some(probe) = probe.as_ref() {
            PvarRecordReader::new_with_probe(&pvar, probe)?
        } else {
            PvarRecordReader::new(&pvar)?
        };
        #[cfg(not(test))]
        let pvar_reader = PvarRecordReader::new(&pvar)?;

        let decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
        let haplotype_state = (options.matrix_kind == MatrixKind::Haplotype)
            .then(|| Box::new(PgenHaplotypeDecodeState::default()));
        let batch = HardcallBatch::new(header.sample_ct);
        let selected_values = Vec::with_capacity(selection.samples.len());
        let genotype_filter_plan = options.variant_filter.as_ref().map_or(
            GenotypeFilterPlan::Generic,
            VariantFilter::genotype_filter_plan,
        );
        let eof = options
            .variant_filter
            .as_ref()
            .is_some_and(VariantFilter::is_always_false);
        let diagnostics = selection.diagnostics.clone();

        Ok(Self {
            pgen,
            pvar,
            pgen_reader,
            pvar_reader,
            header,
            selection,
            source_position: 0,
            diagnostics,
            variant_filter: options.variant_filter,
            genotype_filter_plan,
            matrix_kind: options.matrix_kind,
            dosage_source: options.dosage_source,
            missing_policy: options.missing_policy,
            sparse: options.sparse,
            return_samples: options.return_samples,
            return_variants: options.return_variants,
            all_samples_selected: options.requested_samples.is_none(),
            decoder_state,
            haplotype_state,
            batch,
            selected_values,
            missing_indices: Vec::new(),
            eof,
            failed: false,
            #[cfg(test)]
            probe,
        })
    }

    #[cfg(test)]
    fn open_with_probe(
        pgen: PathBuf,
        pvar: PathBuf,
        psam: PathBuf,
        options: BlockReadOptions,
        probe: Plink2WorkProbe,
    ) -> Result<Self> {
        Self::open_impl(pgen, pvar, psam, options, Some(probe))
    }

    pub(crate) fn next_block(&mut self, block_size: usize) -> Result<Option<BlockOutput>> {
        if self.eof || self.failed {
            return Ok(None);
        }
        let result = match (self.matrix_kind, self.dosage_source, self.sparse) {
            (MatrixKind::Genotype, DosageSource::Hardcall, false) => self
                .next_dense_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            (MatrixKind::Genotype, DosageSource::Hardcall, true) => self
                .next_sparse_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Sparse)),
            (MatrixKind::Genotype, DosageSource::Dosage, false) => self
                .next_dense_dosage_block(block_size, false)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            (MatrixKind::Haplotype, DosageSource::Dosage, false) => self
                .next_dense_dosage_block(block_size, true)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            (MatrixKind::Haplotype, DosageSource::Hardcall, _) => Err(GenoioError::unsupported(
                "plink2 block session does not support hardcall haplotypes yet",
            )),
            (_, DosageSource::Dosage, true) => Err(GenoioError::unsupported(
                "plink2 dosage block sessions do not support sparse matrices",
            )),
        };
        if result.is_err() {
            self.failed = true;
        }
        result
    }

    fn next_dense_block(&mut self, block_size: usize) -> Result<Option<DenseGenotypeMatrix>> {
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
        self.batch.clear();
        let mut batch_start = 0_usize;
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            let Some((variant_index, mut variant)) = self.next_source_variant()? else {
                break;
            };
            let decoded_for_state = self.decode_for_future_ld_state(variant_index)?;
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }

            if !decoded_for_state {
                self.decode_variant(variant_index)?;
            }
            if matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
                let filter = self.variant_filter.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = evaluate_packed_hardcall_filter(
                    &self.decoder_state.packed,
                    &self.selection.source_indices,
                    self.all_samples_selected,
                    filter,
                    self.genotype_filter_plan,
                    Some(&variant),
                    self.return_variants,
                )?;
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                if let Some(stats) = stats {
                    attach_variant_stats(&mut variant, stats);
                }
            }

            if let Some(variants) = variants.as_mut() {
                variants.push_record(&variant)?;
            }
            self.batch.push(&self.decoder_state.packed);
            output_variant_count += 1;
            if self.batch.is_full() {
                flush_hardcall_batch_into_sample_major(
                    &mut self.batch,
                    &self.selection.source_indices,
                    &mut batch_start,
                    block_size,
                    &mut values,
                    self.missing_policy,
                    &mut self.selected_values,
                    &mut self.missing_indices,
                )?;
            }
        }

        self.finish_source_if_at_pgen_end()?;
        flush_hardcall_batch_into_sample_major(
            &mut self.batch,
            &self.selection.source_indices,
            &mut batch_start,
            block_size,
            &mut values,
            self.missing_policy,
            &mut self.selected_values,
            &mut self.missing_indices,
        )?;
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

    fn next_sparse_block(&mut self, block_size: usize) -> Result<Option<SparseGenotypeMatrix>> {
        if self.eof || block_size == 0 {
            return Ok(None);
        }

        let indptr_len = checked_sparse_indptr_len(block_size)?;
        self.record_sparse_allocation(indptr_len);
        let mut indptr = Vec::with_capacity(indptr_len);
        indptr.push(0);
        let mut indices = Vec::new();
        let mut data = Vec::new();
        let mut variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(block_size));
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            let Some((variant_index, mut variant)) = self.next_source_variant()? else {
                break;
            };
            let decoded_for_state = self.decode_for_future_ld_state(variant_index)?;
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }

            if !decoded_for_state {
                self.decode_variant(variant_index)?;
            }
            if matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
                let filter = self.variant_filter.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = evaluate_packed_hardcall_filter(
                    &self.decoder_state.packed,
                    &self.selection.source_indices,
                    self.all_samples_selected,
                    filter,
                    self.genotype_filter_plan,
                    Some(&variant),
                    self.return_variants,
                )?;
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                if let Some(stats) = stats {
                    attach_variant_stats(&mut variant, stats);
                }
            }

            self.decoder_state.packed.expand_selected(
                &self.selection.source_indices,
                &mut self.decoder_state.values,
                &mut self.decoder_state.missing_indices,
            );
            reject_sparse_missing(!self.decoder_state.missing_indices.is_empty())?;
            flip_values_to_minor_allele(&mut self.decoder_state.values, &mut variant);
            append_sparse_column(
                &mut indptr,
                &mut indices,
                &mut data,
                &self.decoder_state.values,
            )?;
            output_variant_count += 1;
            if let Some(variants) = variants.as_mut() {
                variants.push_record(&variant)?;
            }
        }

        self.finish_source_if_at_pgen_end()?;
        if output_variant_count == 0 {
            return Ok(None);
        }
        let samples = SampleMetadataBuffers::optional_from_records(
            &self.selection.samples,
            self.return_samples,
            false,
        )?;
        let diagnostics = block_diagnostics_snapshot(&self.diagnostics, output_variant_count);
        SparseGenotypeMatrix::new(
            self.selection.samples.len(),
            output_variant_count,
            indptr,
            indices,
            data,
            samples,
            variants,
            diagnostics,
        )
        .map(Some)
    }

    fn next_dense_dosage_block(
        &mut self,
        block_size: usize,
        haplotype: bool,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        if self.eof || block_size == 0 {
            return Ok(None);
        }

        let n_rows = if haplotype {
            self.selection.samples.len().checked_mul(2).ok_or_else(|| {
                GenoioError::internal_contract("haplotype row count is out of range")
            })?
        } else {
            self.selection.samples.len()
        };
        let output_len = checked_dense_block_len(n_rows, block_size)?;
        self.record_dense_allocation(output_len);
        let mut variant_major_values = Vec::with_capacity(output_len);
        let mut variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(block_size));
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            let Some((variant_index, mut variant)) = self.next_source_variant()? else {
                break;
            };
            let main_track_cursor = if matches!(self.header.layout, PgenLayout::VariableWidth) {
                let cursor = read_plink2_variant_dosage_main_track(
                    &self.pgen,
                    &mut self.pgen_reader,
                    &self.header,
                    variant_index,
                    &mut self.decoder_state,
                )?;
                self.record_main_decode();
                Some(cursor)
            } else {
                None
            };
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }

            if haplotype {
                let cursor = if let Some(cursor) = main_track_cursor {
                    cursor
                } else {
                    let cursor = read_plink2_variant_haplotype_dosage_track(
                        &self.pgen,
                        &mut self.pgen_reader,
                        &self.header,
                        variant_index,
                        &mut self.decoder_state,
                    )?;
                    self.record_main_decode();
                    cursor
                };
                decode_plink2_haplotype_dosage_aux(
                    &self.pgen,
                    &self.header,
                    variant_index,
                    cursor,
                    &self.selection.source_indices,
                    &self.decoder_state,
                    self.haplotype_state.as_deref_mut().ok_or_else(|| {
                        GenoioError::internal_contract(
                            "haplotype dosage session is missing decode state",
                        )
                    })?,
                )?;
            } else if let Some(cursor) = main_track_cursor {
                decode_plink2_variant_dosage_aux(
                    &self.pgen,
                    &self.header,
                    variant_index,
                    cursor,
                    &self.selection.source_indices,
                    &mut self.decoder_state,
                )?;
            } else {
                read_plink2_variant_dosage(
                    &self.pgen,
                    &mut self.pgen_reader,
                    &self.header,
                    variant_index,
                    &self.selection.source_indices,
                    &mut self.decoder_state,
                )?;
                self.record_main_decode();
            }
            self.record_auxiliary_decode();

            if matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
                let filter = require_genotype_decision_filter(self.variant_filter.as_ref())?;
                let (values, missing_indices) = if haplotype {
                    let haplotype_state = self.haplotype_state.as_deref().ok_or_else(|| {
                        GenoioError::internal_contract(
                            "haplotype dosage session is missing decode state",
                        )
                    })?;
                    (
                        haplotype_state.selected_collapsed_values.as_slice(),
                        haplotype_state
                            .selected_collapsed_missing_indices
                            .as_slice(),
                    )
                } else {
                    (
                        self.decoder_state.values.as_slice(),
                        self.decoder_state.missing_indices.as_slice(),
                    )
                };
                let (retain_variant, stats) = evaluate_dosage_filter(
                    values,
                    missing_indices,
                    filter,
                    &variant,
                    self.return_variants,
                )?;
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                if let Some(stats) = stats {
                    attach_variant_stats(&mut variant, stats);
                }
            }

            if let Some(variants) = variants.as_mut() {
                variants.push_record(&variant)?;
            }
            let (values, missing_indices) = if haplotype {
                let haplotype_state = self.haplotype_state.as_deref_mut().ok_or_else(|| {
                    GenoioError::internal_contract(
                        "haplotype dosage session is missing decode state",
                    )
                })?;
                (
                    &mut haplotype_state.selected_haplotype_values,
                    haplotype_state
                        .selected_haplotype_missing_indices
                        .as_slice(),
                )
            } else {
                (
                    &mut self.decoder_state.values,
                    self.decoder_state.missing_indices.as_slice(),
                )
            };
            apply_dense_missing_policy_to_variant(values, missing_indices, self.missing_policy)?;
            variant_major_values.extend_from_slice(values);
            output_variant_count += 1;
        }

        self.finish_source_if_at_pgen_end()?;
        if output_variant_count == 0 {
            return Ok(None);
        }
        let samples = if haplotype {
            let samples = expand_selected_samples_to_haplotypes(&self.selection);
            SampleMetadataBuffers::optional_from_records(&samples, self.return_samples, true)?
        } else {
            SampleMetadataBuffers::optional_from_records(
                &self.selection.samples,
                self.return_samples,
                false,
            )?
        };
        let diagnostics = block_diagnostics_snapshot(&self.diagnostics, output_variant_count);
        DenseGenotypeMatrix::new_with_layout(
            n_rows,
            output_variant_count,
            variant_major_values,
            DenseLayout::VariantMajor,
            samples,
            variants,
            diagnostics,
        )
        .map(Some)
    }

    fn next_source_variant(&mut self) -> Result<Option<(usize, VariantRecord)>> {
        if self.source_position >= self.header.variant_ct {
            self.finish_source_if_at_pgen_end()?;
            return Ok(None);
        }
        let Some((variant_index, variant)) = self.pvar_reader.next_record()? else {
            return Err(GenoioError::invalid_source(
                &self.pvar,
                "pvar contains fewer variants than pgen",
            ));
        };
        if variant_index != self.source_position {
            return Err(GenoioError::internal_contract(
                "PLINK2 PVAR source position is not monotonic",
            ));
        }
        self.source_position += 1;
        self.record_candidate_visit();
        Ok(Some((variant_index, variant)))
    }

    fn finish_source_if_at_pgen_end(&mut self) -> Result<()> {
        if self.source_position != self.header.variant_ct || self.eof {
            return Ok(());
        }
        if self.pvar_reader.next_record()?.is_some() {
            return Err(GenoioError::invalid_source(
                &self.pvar,
                "pvar variant count exceeds pgen variant count",
            ));
        }
        self.pvar_reader.validate_count(self.header.variant_ct)?;
        self.eof = true;
        Ok(())
    }

    fn decode_for_future_ld_state(&mut self, variant_index: usize) -> Result<bool> {
        if !matches!(self.header.layout, PgenLayout::VariableWidth) {
            return Ok(false);
        }
        self.decode_variant(variant_index)?;
        Ok(true)
    }

    fn decode_variant(&mut self, variant_index: usize) -> Result<()> {
        read_plink2_variant_packed(
            &self.pgen,
            &mut self.pgen_reader,
            &self.header,
            variant_index,
            &mut self.decoder_state,
        )?;
        self.record_main_decode();
        Ok(())
    }

    #[cfg(test)]
    fn record_candidate_visit(&self) {
        if let Some(probe) = &self.probe {
            probe.record_candidate_visit();
        }
    }

    #[cfg(not(test))]
    fn record_candidate_visit(&self) {}

    #[cfg(test)]
    fn record_main_decode(&self) {
        if let Some(probe) = &self.probe {
            probe.record_main_decode();
        }
    }

    #[cfg(not(test))]
    fn record_main_decode(&self) {}

    #[cfg(test)]
    fn record_auxiliary_decode(&self) {
        if let Some(probe) = &self.probe {
            probe.record_auxiliary_decode();
        }
    }

    #[cfg(not(test))]
    fn record_auxiliary_decode(&self) {}

    #[cfg(test)]
    fn record_dense_allocation(&self, len: usize) {
        if let Some(probe) = &self.probe {
            probe.record_dense_allocation(len);
        }
    }

    #[cfg(not(test))]
    fn record_dense_allocation(&self, _len: usize) {}

    #[cfg(test)]
    fn record_sparse_allocation(&self, indptr_len: usize) {
        if let Some(probe) = &self.probe {
            probe.record_sparse_allocation(indptr_len);
        }
    }

    #[cfg(not(test))]
    fn record_sparse_allocation(&self, _indptr_len: usize) {}
}

impl Drop for Plink2BlockSession {
    fn drop(&mut self) {
        #[cfg(test)]
        if let Some(probe) = &self.probe {
            probe.record_drop();
        }
    }
}

fn validate_plink2_options(options: &BlockReadOptions) -> Result<()> {
    if options.sparse && options.dosage_source == DosageSource::Dosage {
        return Err(GenoioError::unsupported(
            "plink2 dosage block sessions do not support sparse matrices",
        ));
    }
    Ok(())
}

#[cfg(test)]
#[derive(Debug, Clone, Default)]
pub(super) struct Plink2WorkProbe {
    counts: std::sync::Arc<std::sync::Mutex<Plink2WorkCounts>>,
}

#[cfg(test)]
impl Plink2WorkProbe {
    fn snapshot(&self) -> Plink2WorkCounts {
        self.counts
            .lock()
            .expect("plink2 work probe lock should not be poisoned")
            .clone()
    }

    fn update(&self, update: impl FnOnce(&mut Plink2WorkCounts)) {
        update(
            &mut self
                .counts
                .lock()
                .expect("plink2 work probe lock should not be poisoned"),
        );
    }

    fn record_pgen_open(&self) {
        self.update(|counts| counts.pgen_opens += 1);
    }

    pub(super) fn record_pvar_open(&self) {
        self.update(|counts| counts.pvar_opens += 1);
    }

    pub(super) fn record_psam_open(&self) {
        self.update(|counts| counts.psam_opens += 1);
    }

    fn record_candidate_visit(&self) {
        self.update(|counts| counts.candidate_visits += 1);
    }

    fn record_main_decode(&self) {
        self.update(|counts| counts.main_decodes += 1);
    }

    fn record_auxiliary_decode(&self) {
        self.update(|counts| counts.auxiliary_decodes += 1);
    }

    fn record_dense_allocation(&self, len: usize) {
        self.update(|counts| counts.max_dense_output_len = counts.max_dense_output_len.max(len));
    }

    fn record_sparse_allocation(&self, indptr_len: usize) {
        self.update(|counts| {
            counts.max_sparse_indptr_len = counts.max_sparse_indptr_len.max(indptr_len);
        });
    }

    fn record_drop(&self) {
        self.update(|counts| counts.drops += 1);
    }
}

#[cfg(test)]
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct Plink2WorkCounts {
    pgen_opens: usize,
    pvar_opens: usize,
    psam_opens: usize,
    candidate_visits: usize,
    main_decodes: usize,
    auxiliary_decodes: usize,
    max_dense_output_len: usize,
    max_sparse_indptr_len: usize,
    drops: usize,
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::{Path, PathBuf};

    use genoio_core::{DenseMissingPolicy, VariantFilter};
    use serde_json::json;

    use crate::blocks::{BlockReadOptions, DosageSource, MatrixKind};

    use super::{Plink2BlockSession, Plink2WorkProbe};

    fn write_text(path: &Path, contents: &str) {
        fs::write(path, contents).expect("test fixture should be written");
    }

    fn write_fixed_fixture(dir: &Path, records: &[u8]) -> (PathBuf, PathBuf, PathBuf) {
        let pgen = dir.join("tiny.pgen");
        let pvar = dir.join("tiny.pvar");
        let psam = dir.join("tiny.psam");
        let mut pgen_bytes = vec![0x6c, 0x1b, 0x02];
        pgen_bytes.extend(
            u32::try_from(records.len())
                .expect("test variant count fits u32")
                .to_le_bytes(),
        );
        pgen_bytes.extend(2_u32.to_le_bytes());
        pgen_bytes.push(0);
        pgen_bytes.extend(records);
        fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
        let pvar_body = (0..records.len())
            .map(|index| format!("1 {} rs{} A G\n", index + 1, index + 1))
            .collect::<String>();
        write_text(&pvar, &format!("#CHROM POS ID REF ALT\n{pvar_body}"));
        write_text(&psam, "#IID\nS1\nS2\n");
        (pgen, pvar, psam)
    }

    fn write_variable_fixture(
        dir: &Path,
        record_types: &[u8],
        records: &[&[u8]],
        variants: &str,
    ) -> (PathBuf, PathBuf, PathBuf) {
        let pgen = dir.join("variable.pgen");
        let pvar = dir.join("variable.pvar");
        let psam = dir.join("variable.psam");
        let header_len = 12 + 8 + record_types.len() + records.len();
        let mut pgen_bytes = vec![0x6c, 0x1b, 0x10];
        pgen_bytes.extend(
            u32::try_from(records.len())
                .expect("test variant count fits u32")
                .to_le_bytes(),
        );
        pgen_bytes.extend(2_u32.to_le_bytes());
        pgen_bytes.push(0x04);
        pgen_bytes.extend(
            u64::try_from(header_len)
                .expect("test header length fits u64")
                .to_le_bytes(),
        );
        pgen_bytes.extend(record_types);
        pgen_bytes.extend(
            records.iter().map(|record| {
                u8::try_from(record.len()).expect("test record length fits one byte")
            }),
        );
        for record in records {
            pgen_bytes.extend(*record);
        }
        fs::write(&pgen, pgen_bytes).expect("variable pgen fixture should be written");
        write_text(&pvar, variants);
        write_text(&psam, "#IID\nS1\nS2\n");
        (pgen, pvar, psam)
    }

    fn options(filter: Option<VariantFilter>) -> BlockReadOptions {
        BlockReadOptions {
            matrix_kind: MatrixKind::Genotype,
            sparse: false,
            requested_samples: None,
            variant_filter: filter,
            dosage_source: DosageSource::Hardcall,
            missing_policy: DenseMissingPolicy::Nan,
            return_samples: true,
            return_variants: true,
        }
    }

    fn dosage_options(matrix_kind: MatrixKind, filter: Option<VariantFilter>) -> BlockReadOptions {
        BlockReadOptions {
            matrix_kind,
            dosage_source: DosageSource::Dosage,
            ..options(filter)
        }
    }

    fn chrom_filter(chrom: &str) -> VariantFilter {
        VariantFilter::from_json_value(json!({
            "op": "predicate",
            "name": "chrom",
            "params": {"value": chrom}
        }))
        .expect("chromosome filter should parse")
    }

    #[test]
    fn pbr_rust_plink2_001_work_probe_counts_authoritative_opens_and_fixed_width_work() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (pgen, pvar, psam) = write_fixed_fixture(dir.path(), &[0x04, 0x08, 0x00]);
        let probe = Plink2WorkProbe::default();

        {
            let mut session = Plink2BlockSession::open_with_probe(
                pgen,
                pvar,
                psam,
                options(Some(chrom_filter("1"))),
                probe.clone(),
            )
            .expect("persistent plink2 session should open");
            while session
                .next_block(1)
                .expect("persistent plink2 block should decode")
                .is_some()
            {}
            let at_eof = probe.snapshot();
            assert!(session
                .next_block(1)
                .expect("persistent plink2 EOF should be sticky")
                .is_none());
            assert_eq!(probe.snapshot(), at_eof);
        }

        let counts = probe.snapshot();
        assert_eq!(counts.pgen_opens, 1);
        assert_eq!(counts.pvar_opens, 1);
        assert_eq!(counts.psam_opens, 1);
        assert_eq!(counts.candidate_visits, 3);
        assert_eq!(counts.main_decodes, 3);
        assert_eq!(counts.auxiliary_decodes, 0);
        assert_eq!(counts.max_dense_output_len, 2);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_plink2_001_fixed_width_metadata_rejects_decode_no_payload() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (pgen, pvar, psam) = write_fixed_fixture(dir.path(), &[0x04, 0x08, 0x00]);
        write_text(
            &pvar,
            "#CHROM POS ID REF ALT\n1 1 rs1 A G\n2 2 rs2 A G\n2 3 rs3 A G\n",
        );
        let probe = Plink2WorkProbe::default();
        let mut session = Plink2BlockSession::open_with_probe(
            pgen,
            pvar,
            psam,
            options(Some(chrom_filter("1"))),
            probe.clone(),
        )
        .expect("persistent plink2 session should open");

        while session
            .next_block(1)
            .expect("persistent plink2 block should decode")
            .is_some()
        {}

        let counts = probe.snapshot();
        assert_eq!(counts.candidate_visits, 3);
        assert_eq!(counts.main_decodes, 1);
    }

    #[test]
    fn pbr_rust_plink2_001_probe_counts_only_successful_authoritative_opens() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (pgen, pvar, psam) = write_fixed_fixture(dir.path(), &[0x04]);
        fs::remove_file(&pvar).expect("pvar fixture should be removed before opening");
        let probe = Plink2WorkProbe::default();

        assert!(Plink2BlockSession::open_with_probe(
            pgen,
            pvar,
            psam,
            options(None),
            probe.clone(),
        )
        .is_err());

        let counts = probe.snapshot();
        assert_eq!(counts.pgen_opens, 1);
        assert_eq!(counts.psam_opens, 1);
        assert_eq!(counts.pvar_opens, 0);
    }

    #[test]
    fn pbr_rust_plink2_001_sparse_probe_records_linear_work_and_checked_allocation() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (pgen, pvar, psam) = write_fixed_fixture(dir.path(), &[0x04, 0x08, 0x00]);
        let probe = Plink2WorkProbe::default();
        let mut sparse_options = options(None);
        sparse_options.sparse = true;
        sparse_options.missing_policy = DenseMissingPolicy::Raise;
        let mut session =
            Plink2BlockSession::open_with_probe(pgen, pvar, psam, sparse_options, probe.clone())
                .expect("persistent sparse plink2 session should open");

        while session
            .next_block(1)
            .expect("persistent sparse plink2 block should decode")
            .is_some()
        {}

        let counts = probe.snapshot();
        assert_eq!(counts.candidate_visits, 3);
        assert_eq!(counts.main_decodes, 3);
        assert_eq!(counts.max_sparse_indptr_len, 2);
    }

    #[test]
    fn pbr_rust_plink2_001_pbr_rust_alloc_001_early_drop_is_block_bounded() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (pgen, pvar, psam) = write_fixed_fixture(dir.path(), &[0x04, 0x08, 0x00]);
        let probe = Plink2WorkProbe::default();
        {
            let mut session =
                Plink2BlockSession::open_with_probe(pgen, pvar, psam, options(None), probe.clone())
                    .expect("persistent plink2 session should open");
            assert!(session
                .next_block(1)
                .expect("first persistent plink2 block should decode")
                .is_some());
        }
        let counts = probe.snapshot();
        assert_eq!(counts.candidate_visits, 1);
        assert_eq!(counts.main_decodes, 1);
        assert_eq!(counts.max_dense_output_len, 2);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_plink2_001_session_is_send() {
        fn assert_send<T: Send>() {}
        assert_send::<Plink2BlockSession>();
        assert_send::<crate::blocks::BlockReader>();
    }

    #[test]
    fn pbr_rust_plink2_002_probe_tracks_rejected_non_ld_base_without_prefetch() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (pgen, pvar, psam) = write_variable_fixture(
            dir.path(),
            &[0, 0, 2],
            &[&[0x04], &[0x02], &[0]],
            "#CHROM POS ID REF ALT\n1 1 rs1 A G\n2 2 base A G\n1 3 rs3 A G\n",
        );
        let probe = Plink2WorkProbe::default();
        let mut session = Plink2BlockSession::open_with_probe(
            pgen,
            pvar,
            psam,
            options(Some(chrom_filter("1"))),
            probe.clone(),
        )
        .expect("persistent variable-width plink2 session should open");

        assert!(session
            .next_block(1)
            .expect("first variable-width block should decode")
            .is_some());
        let after_first = probe.snapshot();
        assert_eq!(after_first.candidate_visits, 1);
        assert_eq!(after_first.main_decodes, 1);

        let second = session
            .next_block(1)
            .expect("LD-dependent variable-width block should decode")
            .expect("second retained block should exist");
        let crate::blocks::BlockOutput::Dense(second) = second else {
            panic!("dense session should return a dense block");
        };
        assert_eq!(second.values, vec![2.0, 0.0]);
        let after_second = probe.snapshot();
        assert_eq!(after_second.candidate_visits, 3);
        assert_eq!(after_second.main_decodes, 3);
        assert_eq!(after_second.auxiliary_decodes, 0);
    }

    #[test]
    fn pbr_rust_plink2_002_failed_decode_preserves_ld_base_and_terminates_without_prefetch() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (pgen, pvar, psam) = write_variable_fixture(
            dir.path(),
            &[0, 4, 2],
            &[&[0x04], &[1], &[0]],
            "#CHROM POS ID REF ALT\n1 1 rs1 A G\n1 2 bad A G\n1 3 rs3 A G\n",
        );
        let probe = Plink2WorkProbe::default();
        let mut session =
            Plink2BlockSession::open_with_probe(pgen, pvar, psam, options(None), probe.clone())
                .expect("persistent variable-width plink2 session should open");

        assert!(session
            .next_block(1)
            .expect("first variable-width block should decode")
            .is_some());
        let error = session
            .next_block(1)
            .expect_err("malformed later main track should fail");
        assert!(matches!(
            error,
            genoio_core::GenoioError::InvalidSource { .. }
        ));
        let at_error = probe.snapshot();
        assert_eq!(at_error.candidate_visits, 2);
        assert_eq!(at_error.main_decodes, 1);
        assert!(session
            .next_block(1)
            .expect("failed session should be terminal")
            .is_none());
        assert_eq!(probe.snapshot(), at_error);

        // Resume only inside this private test to prove the failed non-LD
        // decode did not replace the prior valid LD base.
        session.failed = false;
        let third = session
            .next_block(1)
            .expect("test-only resumed session should decode prior-base LD record")
            .expect("third retained block should exist");
        let crate::blocks::BlockOutput::Dense(third) = third else {
            panic!("dense session should return a dense block");
        };
        assert_eq!(third.values, vec![0.0, 1.0]);
    }

    #[test]
    fn pbr_rust_plink2_003_probe_separates_main_and_auxiliary_dosage_work() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let mut rejected_base = vec![0x04];
        rejected_base.extend(16_384_u16.to_le_bytes());
        rejected_base.extend(8_192_u16.to_le_bytes());
        let mut retained_ld = vec![0];
        retained_ld.extend(0_u16.to_le_bytes());
        retained_ld.extend(3_277_u16.to_le_bytes());
        let (pgen, pvar, psam) = write_variable_fixture(
            dir.path(),
            &[0x40, 0x42],
            &[&rejected_base, &retained_ld],
            "#CHROM POS ID REF ALT\n2 1 rejected_base A G\n1 2 retained_ld A G\n",
        );
        let probe = Plink2WorkProbe::default();
        let mut session = Plink2BlockSession::open_with_probe(
            pgen,
            pvar,
            psam,
            dosage_options(MatrixKind::Genotype, Some(chrom_filter("1"))),
            probe.clone(),
        )
        .expect("persistent dosage session should open");

        let block = session
            .next_block(1)
            .expect("retained LD dosage block should decode")
            .expect("retained dosage block should exist");
        let crate::blocks::BlockOutput::Dense(block) = block else {
            panic!("dosage session should return a dense block");
        };
        assert_eq!(block.n_variants, 1);
        let counts = probe.snapshot();
        assert_eq!(counts.candidate_visits, 2);
        assert_eq!(counts.main_decodes, 2);
        assert_eq!(counts.auxiliary_decodes, 1);
        assert_eq!(counts.max_dense_output_len, 2);
    }
}
