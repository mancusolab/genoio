// pattern: Imperative Shell
//! Persistent sequential BCF block-reader sessions.

use std::path::PathBuf;

use genoio_core::{
    append_sparse_column, reject_sparse_missing, DenseDiagnostics, DenseGenotypeMatrix,
    DenseLayout, DenseSampleSelection, GenoioError, SampleMetadataBuffers, SparseGenotypeMatrix,
    VariantFilter, VariantMetadataBuffers, VariantWindow,
};
use noodles_bcf as bcf;
use noodles_vcf as noodles;

use crate::blocks::{
    block_diagnostics_snapshot, checked_dense_block_len, checked_sparse_indptr_len, BlockOutput,
    BlockReadOptions, DosageSource, MatrixKind,
};
use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::matrix::apply_dense_missing_policy_to_variant;
use crate::retention::{RetainedVariantState, RetentionAction};

use super::super::haplotype_sample_records;
use super::decode::{decode_ds_record, decode_gt_record, BcfDenseDecodeBuffers, BcfStatsMode};
use super::haplotype::{decode_phased_haplotype_record, BcfHaplotypeDecodeBuffers};
use super::record::{prepare_bcf_candidate, push_bcf_variant_row, BcfCandidateAction};
#[cfg(not(test))]
use super::source::open_bcf_input;
#[cfg(test)]
use super::source::open_bcf_input_with_hooks;
use super::source::{
    evaluate_bcf_gt_filter, flip_haplotype_values_to_minor_allele, flip_values_to_minor_allele,
    BcfInput,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BcfMode {
    DenseGenotype,
    DenseDosage,
    SparseGenotype,
    DenseHaplotype,
    SparseHaplotype,
}

/// One forward-only BCF reader with reusable record and decode buffers.
pub(crate) struct BcfBlockSession {
    path: PathBuf,
    reader: bcf::io::Reader<noodles_bgzf::io::Reader<std::fs::File>>,
    header: noodles::Header,
    selection: DenseSampleSelection,
    diagnostics: DenseDiagnostics,
    variant_filter: Option<VariantFilter>,
    missing_policy: genoio_core::DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    mode: BcfMode,
    record: bcf::Record,
    decoded: BcfDenseDecodeBuffers,
    haplotype_decoded: BcfHaplotypeDecodeBuffers,
    eof: bool,
    #[cfg(test)]
    probe: Option<BcfWorkProbe>,
}

impl BcfBlockSession {
    pub(crate) fn open(path: PathBuf, options: BlockReadOptions) -> Result<Self> {
        #[cfg(test)]
        {
            Self::open_impl(path, options, None)
        }
        #[cfg(not(test))]
        {
            Self::open_impl(path, options)
        }
    }

    fn open_impl(
        path: PathBuf,
        options: BlockReadOptions,
        #[cfg(test)] probe: Option<BcfWorkProbe>,
    ) -> Result<Self> {
        let mode = validate_bcf_options(&options)?;
        #[cfg(test)]
        let input = {
            let source_probe = probe.clone();
            let header_probe = probe.clone();
            open_bcf_input_with_hooks(
                &path,
                options.requested_samples.as_deref(),
                move || {
                    if let Some(probe) = source_probe {
                        probe.record_source_open();
                    }
                },
                move || {
                    if let Some(probe) = header_probe {
                        probe.record_header_parse();
                    }
                },
            )?
        };
        #[cfg(not(test))]
        let input = open_bcf_input(&path, options.requested_samples.as_deref())?;
        let BcfInput {
            reader,
            header,
            selection,
        } = input;
        let n_samples = selection.source_indices.len();
        let eof = options
            .variant_filter
            .as_ref()
            .is_some_and(VariantFilter::is_always_false);
        let diagnostics = selection.diagnostics.clone();

        Ok(Self {
            path,
            reader,
            header,
            selection,
            diagnostics,
            variant_filter: options.variant_filter,
            missing_policy: options.missing_policy,
            return_samples: options.return_samples,
            return_variants: options.return_variants,
            mode,
            record: bcf::Record::default(),
            decoded: BcfDenseDecodeBuffers::with_capacity(n_samples),
            haplotype_decoded: BcfHaplotypeDecodeBuffers::with_capacity(n_samples),
            eof,
            #[cfg(test)]
            probe,
        })
    }

    #[cfg(test)]
    fn open_with_probe(
        path: PathBuf,
        options: BlockReadOptions,
        probe: BcfWorkProbe,
    ) -> Result<Self> {
        Self::open_impl(path, options, Some(probe))
    }

    pub(crate) fn next_block(&mut self, block_size: usize) -> Result<Option<BlockOutput>> {
        if self.eof || block_size == 0 {
            return Ok(None);
        }
        let result = match self.mode {
            BcfMode::DenseGenotype => self
                .next_dense_block(block_size, DenseField::Gt)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            BcfMode::DenseDosage => self
                .next_dense_block(block_size, DenseField::Ds)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            BcfMode::SparseGenotype => self
                .next_sparse_genotype_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Sparse)),
            BcfMode::DenseHaplotype => self
                .next_dense_haplotype_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            BcfMode::SparseHaplotype => self
                .next_sparse_haplotype_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Sparse)),
        };
        if result.is_err() {
            self.eof = true;
        }
        result
    }

    fn next_dense_block(
        &mut self,
        block_size: usize,
        field: DenseField,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        let n_samples = self.selection.samples.len();
        let allocation_len = checked_dense_block_len(n_samples, block_size)?;
        self.record_dense_allocation(allocation_len);
        let mut variant_major_values = Vec::with_capacity(allocation_len);
        let mut variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(block_size));
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let prepared = match prepare_bcf_candidate(
                &self.path,
                &self.header,
                &self.record,
                self.variant_filter.as_ref(),
                &mut retention,
                &mut self.diagnostics,
            )? {
                BcfCandidateAction::Skip => continue,
                BcfCandidateAction::Stop => break,
                BcfCandidateAction::Decode(prepared) => prepared,
            };
            let variant = prepared.variant;
            let needs_genotype_decision = prepared.needs_genotype_decision;
            match field {
                DenseField::Gt => {
                    decode_gt_record(
                        &self.path,
                        &self.header,
                        &self.record,
                        &self.selection.source_indices,
                        match (needs_genotype_decision, self.return_variants) {
                            (true, false) => BcfStatsMode::Counts,
                            (true, true) => BcfStatsMode::Compute,
                            (false, _) => BcfStatsMode::Skip,
                        },
                        &mut self.decoded,
                    )?;
                    self.record_gt_decode();
                }
                DenseField::Ds => {
                    decode_ds_record(
                        &self.path,
                        &self.header,
                        &self.record,
                        &self.selection.source_indices,
                        false,
                        &mut self.decoded,
                    )?;
                    self.record_ds_decode();
                }
            }

            let mut stats_to_attach = None;
            if needs_genotype_decision {
                let filter = self.variant_filter.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = match field {
                    DenseField::Gt => evaluate_bcf_gt_filter(
                        &self.decoded,
                        filter,
                        &variant,
                        self.return_variants,
                        "GT",
                    )?,
                    DenseField::Ds => evaluate_dosage_filter(
                        &self.decoded.values,
                        &self.decoded.missing_indices,
                        filter,
                        &variant,
                        self.return_variants,
                    )?,
                };
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                stats_to_attach = stats;
            }

            push_bcf_variant_row(&mut variants, &variant, stats_to_attach, false)?;
            apply_dense_missing_policy_to_variant(
                &mut self.decoded.values,
                &self.decoded.missing_indices,
                self.missing_policy,
            )?;
            variant_major_values.extend_from_slice(&self.decoded.values);
            output_variant_count += 1;
        }

        self.finish_dense_output(variant_major_values, variants, output_variant_count)
    }

    fn next_sparse_genotype_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<SparseGenotypeMatrix>> {
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

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let prepared = match prepare_bcf_candidate(
                &self.path,
                &self.header,
                &self.record,
                self.variant_filter.as_ref(),
                &mut retention,
                &mut self.diagnostics,
            )? {
                BcfCandidateAction::Skip => continue,
                BcfCandidateAction::Stop => break,
                BcfCandidateAction::Decode(prepared) => prepared,
            };
            let variant = prepared.variant;
            let needs_genotype_decision = prepared.needs_genotype_decision;
            decode_gt_record(
                &self.path,
                &self.header,
                &self.record,
                &self.selection.source_indices,
                BcfStatsMode::from_needed(needs_genotype_decision),
                &mut self.decoded,
            )?;
            self.record_gt_decode();
            let mut stats_to_attach = None;
            if needs_genotype_decision {
                let stats = self.decoded.stats;
                let retain_variant = self
                    .variant_filter
                    .as_ref()
                    .is_none_or(|filter| filter.evaluate_view(&variant, stats.as_ref()));
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                stats_to_attach = stats;
            }

            reject_sparse_missing(!self.decoded.missing_indices.is_empty())?;
            let flipped = flip_values_to_minor_allele(self.decoded.values.as_mut_slice());
            append_sparse_column(&mut indptr, &mut indices, &mut data, &self.decoded.values)?;
            push_bcf_variant_row(&mut variants, &variant, stats_to_attach, flipped)?;
        }

        self.finish_sparse_output(indptr, indices, data, variants)
    }

    fn next_dense_haplotype_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        let n_rows = self.selection.samples.len().checked_mul(2).ok_or_else(|| {
            GenoioError::internal_contract("BCF haplotype row count is out of range")
        })?;
        let allocation_len = checked_dense_block_len(n_rows, block_size)?;
        self.record_dense_allocation(allocation_len);
        let mut variant_major_values = Vec::with_capacity(allocation_len);
        let mut variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(block_size));
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let prepared = match prepare_bcf_candidate(
                &self.path,
                &self.header,
                &self.record,
                self.variant_filter.as_ref(),
                &mut retention,
                &mut self.diagnostics,
            )? {
                BcfCandidateAction::Skip => continue,
                BcfCandidateAction::Stop => break,
                BcfCandidateAction::Decode(prepared) => prepared,
            };
            let variant = prepared.variant;
            let needs_genotype_decision = prepared.needs_genotype_decision;
            let mut stats_to_attach = None;
            if needs_genotype_decision {
                decode_gt_record(
                    &self.path,
                    &self.header,
                    &self.record,
                    &self.selection.source_indices,
                    if self.return_variants {
                        BcfStatsMode::Compute
                    } else {
                        BcfStatsMode::Counts
                    },
                    &mut self.decoded,
                )?;
                self.record_gt_decode();
                let filter = self.variant_filter.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = evaluate_bcf_gt_filter(
                    &self.decoded,
                    filter,
                    &variant,
                    self.return_variants,
                    "haplotype",
                )?;
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                stats_to_attach = stats;
            }

            decode_phased_haplotype_record(
                &self.path,
                &self.header,
                &self.record,
                &self.selection.source_indices,
                &mut self.haplotype_decoded,
            )?;
            self.record_phase_decode();
            apply_dense_missing_policy_to_variant(
                &mut self.haplotype_decoded.values,
                &self.haplotype_decoded.missing_indices,
                self.missing_policy,
            )?;
            push_bcf_variant_row(&mut variants, &variant, stats_to_attach, false)?;
            variant_major_values.extend_from_slice(&self.haplotype_decoded.values);
            output_variant_count += 1;
        }

        self.finish_dense_haplotype_output(variant_major_values, variants, output_variant_count)
    }

    fn next_sparse_haplotype_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<SparseGenotypeMatrix>> {
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

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let prepared = match prepare_bcf_candidate(
                &self.path,
                &self.header,
                &self.record,
                self.variant_filter.as_ref(),
                &mut retention,
                &mut self.diagnostics,
            )? {
                BcfCandidateAction::Skip => continue,
                BcfCandidateAction::Stop => break,
                BcfCandidateAction::Decode(prepared) => prepared,
            };
            let variant = prepared.variant;
            let needs_genotype_decision = prepared.needs_genotype_decision;
            let mut stats_to_attach = None;
            if needs_genotype_decision {
                decode_gt_record(
                    &self.path,
                    &self.header,
                    &self.record,
                    &self.selection.source_indices,
                    BcfStatsMode::Compute,
                    &mut self.decoded,
                )?;
                self.record_gt_decode();
                let stats = self.decoded.stats;
                let retain_variant = self
                    .variant_filter
                    .as_ref()
                    .is_none_or(|filter| filter.evaluate_view(&variant, stats.as_ref()));
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                stats_to_attach = stats;
            }

            decode_phased_haplotype_record(
                &self.path,
                &self.header,
                &self.record,
                &self.selection.source_indices,
                &mut self.haplotype_decoded,
            )?;
            self.record_phase_decode();
            reject_sparse_missing(!self.haplotype_decoded.missing_indices.is_empty())?;
            let flipped =
                flip_haplotype_values_to_minor_allele(self.haplotype_decoded.values.as_mut_slice());
            append_sparse_column(
                &mut indptr,
                &mut indices,
                &mut data,
                &self.haplotype_decoded.values,
            )?;
            push_bcf_variant_row(&mut variants, &variant, stats_to_attach, flipped)?;
        }

        self.finish_sparse_haplotype_output(indptr, indices, data, variants)
    }

    fn read_next_record(&mut self) -> Result<bool> {
        let read = self.reader.read_record(&mut self.record).map_err(|error| {
            GenoioError::invalid_source(&self.path, format!("bcf record error: {error}"))
        })?;
        self.record_read_call();
        if read == 0 {
            self.eof = true;
            return Ok(false);
        }
        self.record_candidate_visit();
        Ok(true)
    }

    fn finish_dense_output(
        &self,
        variant_major_values: Vec<f32>,
        variants: Option<VariantMetadataBuffers>,
        output_variant_count: usize,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        if output_variant_count == 0 {
            return Ok(None);
        }
        let samples = SampleMetadataBuffers::optional_from_records(
            &self.selection.samples,
            self.return_samples,
            false,
        )?;
        DenseGenotypeMatrix::new_with_layout(
            self.selection.samples.len(),
            output_variant_count,
            variant_major_values,
            DenseLayout::VariantMajor,
            samples,
            variants,
            block_diagnostics_snapshot(&self.diagnostics, output_variant_count),
        )
        .map(Some)
    }

    fn finish_sparse_output(
        &self,
        indptr: Vec<i32>,
        indices: Vec<i32>,
        data: Vec<f32>,
        variants: Option<VariantMetadataBuffers>,
    ) -> Result<Option<SparseGenotypeMatrix>> {
        let output_variant_count = indptr.len().saturating_sub(1);
        if output_variant_count == 0 {
            return Ok(None);
        }
        let samples = SampleMetadataBuffers::optional_from_records(
            &self.selection.samples,
            self.return_samples,
            false,
        )?;
        SparseGenotypeMatrix::new(
            self.selection.samples.len(),
            output_variant_count,
            indptr,
            indices,
            data,
            samples,
            variants,
            block_diagnostics_snapshot(&self.diagnostics, output_variant_count),
        )
        .map(Some)
    }

    fn finish_dense_haplotype_output(
        &self,
        variant_major_values: Vec<f32>,
        variants: Option<VariantMetadataBuffers>,
        output_variant_count: usize,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        if output_variant_count == 0 {
            return Ok(None);
        }
        let haplotype_samples =
            haplotype_sample_records(&self.selection.samples, &self.selection.source_indices);
        let samples = SampleMetadataBuffers::optional_from_records(
            &haplotype_samples,
            self.return_samples,
            true,
        )?;
        DenseGenotypeMatrix::new_with_layout(
            haplotype_samples.len(),
            output_variant_count,
            variant_major_values,
            DenseLayout::VariantMajor,
            samples,
            variants,
            block_diagnostics_snapshot(&self.diagnostics, output_variant_count),
        )
        .map(Some)
    }

    fn finish_sparse_haplotype_output(
        &self,
        indptr: Vec<i32>,
        indices: Vec<i32>,
        data: Vec<f32>,
        variants: Option<VariantMetadataBuffers>,
    ) -> Result<Option<SparseGenotypeMatrix>> {
        let output_variant_count = indptr.len().saturating_sub(1);
        if output_variant_count == 0 {
            return Ok(None);
        }
        let haplotype_samples =
            haplotype_sample_records(&self.selection.samples, &self.selection.source_indices);
        let samples = SampleMetadataBuffers::optional_from_records(
            &haplotype_samples,
            self.return_samples,
            true,
        )?;
        SparseGenotypeMatrix::new(
            haplotype_samples.len(),
            output_variant_count,
            indptr,
            indices,
            data,
            samples,
            variants,
            block_diagnostics_snapshot(&self.diagnostics, output_variant_count),
        )
        .map(Some)
    }

    #[cfg(test)]
    fn record_read_call(&self) {
        if let Some(probe) = &self.probe {
            probe.record_read_call();
        }
    }

    #[cfg(not(test))]
    fn record_read_call(&self) {}

    #[cfg(test)]
    fn record_candidate_visit(&self) {
        if let Some(probe) = &self.probe {
            probe.record_candidate_visit();
        }
    }

    #[cfg(not(test))]
    fn record_candidate_visit(&self) {}

    #[cfg(test)]
    fn record_gt_decode(&self) {
        if let Some(probe) = &self.probe {
            probe.record_gt_decode();
        }
    }

    #[cfg(not(test))]
    fn record_gt_decode(&self) {}

    #[cfg(test)]
    fn record_ds_decode(&self) {
        if let Some(probe) = &self.probe {
            probe.record_ds_decode();
        }
    }

    #[cfg(not(test))]
    fn record_ds_decode(&self) {}

    #[cfg(test)]
    fn record_phase_decode(&self) {
        if let Some(probe) = &self.probe {
            probe.record_phase_decode();
        }
    }

    #[cfg(not(test))]
    fn record_phase_decode(&self) {}

    #[cfg(test)]
    fn record_dense_allocation(&self, len: usize) {
        if let Some(probe) = &self.probe {
            probe.record_dense_allocation(len);
        }
    }

    #[cfg(not(test))]
    fn record_dense_allocation(&self, _len: usize) {}

    #[cfg(test)]
    fn record_sparse_allocation(&self, len: usize) {
        if let Some(probe) = &self.probe {
            probe.record_sparse_allocation(len);
        }
    }

    #[cfg(not(test))]
    fn record_sparse_allocation(&self, _len: usize) {}
}

impl Drop for BcfBlockSession {
    fn drop(&mut self) {
        #[cfg(test)]
        if let Some(probe) = &self.probe {
            probe.record_drop();
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DenseField {
    Gt,
    Ds,
}

fn validate_bcf_options(options: &BlockReadOptions) -> Result<BcfMode> {
    match (options.matrix_kind, options.sparse, options.dosage_source) {
        (MatrixKind::Genotype, false, DosageSource::Hardcall) => Ok(BcfMode::DenseGenotype),
        (MatrixKind::Genotype, false, DosageSource::Dosage) => Ok(BcfMode::DenseDosage),
        (MatrixKind::Genotype, true, DosageSource::Hardcall) => Ok(BcfMode::SparseGenotype),
        (MatrixKind::Genotype, true, DosageSource::Dosage) => Err(GenoioError::unsupported(
            "BCF dosage blocks support dense genotype matrices only",
        )),
        (MatrixKind::Haplotype, false, DosageSource::Hardcall) => Ok(BcfMode::DenseHaplotype),
        (MatrixKind::Haplotype, true, DosageSource::Hardcall) => Ok(BcfMode::SparseHaplotype),
        (MatrixKind::Haplotype, _, DosageSource::Dosage) => Err(GenoioError::unsupported(
            "BCF dosage blocks support dense genotype matrices only",
        )),
    }
}

#[cfg(test)]
#[derive(Debug, Clone, Default)]
struct BcfWorkProbe {
    counts: std::sync::Arc<std::sync::Mutex<BcfWorkCounts>>,
}

#[cfg(test)]
impl BcfWorkProbe {
    fn snapshot(&self) -> BcfWorkCounts {
        self.counts
            .lock()
            .expect("BCF probe lock should not be poisoned")
            .clone()
    }

    fn update(&self, update: impl FnOnce(&mut BcfWorkCounts)) {
        update(
            &mut self
                .counts
                .lock()
                .expect("BCF probe lock should not be poisoned"),
        );
    }

    fn record_source_open(&self) {
        self.update(|counts| counts.source_opens += 1);
    }

    fn record_header_parse(&self) {
        self.update(|counts| counts.header_parses += 1);
    }

    fn record_read_call(&self) {
        self.update(|counts| counts.read_record_calls += 1);
    }

    fn record_candidate_visit(&self) {
        self.update(|counts| counts.candidate_visits += 1);
    }

    fn record_gt_decode(&self) {
        self.update(|counts| counts.gt_decodes += 1);
    }

    fn record_ds_decode(&self) {
        self.update(|counts| counts.ds_decodes += 1);
    }

    fn record_phase_decode(&self) {
        self.update(|counts| counts.phase_decodes += 1);
    }

    fn record_dense_allocation(&self, len: usize) {
        self.update(|counts| counts.max_dense_output_len = counts.max_dense_output_len.max(len));
    }

    fn record_sparse_allocation(&self, len: usize) {
        self.update(|counts| counts.max_sparse_indptr_len = counts.max_sparse_indptr_len.max(len));
    }

    fn record_drop(&self) {
        self.update(|counts| counts.drops += 1);
    }
}

#[cfg(test)]
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct BcfWorkCounts {
    source_opens: usize,
    header_parses: usize,
    read_record_calls: usize,
    candidate_visits: usize,
    gt_decodes: usize,
    ds_decodes: usize,
    phase_decodes: usize,
    max_dense_output_len: usize,
    max_sparse_indptr_len: usize,
    drops: usize,
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::Path;

    use genoio_core::DenseMissingPolicy;
    use noodles_core::Position;
    use noodles_vcf::{
        header::record::value::{
            map::{
                format::{Number, Type},
                Contig, Format,
            },
            Map,
        },
        variant::{
            io::Write as _,
            record::samples::keys::key,
            record_buf::{samples::sample::Value, samples::Keys, AlternateBases, Ids, Samples},
        },
    };

    use super::*;

    fn record(
        chrom: &str,
        id: &str,
        pos: usize,
        genotypes: [&str; 2],
    ) -> noodles_vcf::variant::RecordBuf {
        let ids: Ids = [id.to_owned()].into_iter().collect();
        let keys: Keys = [String::from(key::GENOTYPE)].into_iter().collect();
        let samples = Samples::new(
            keys,
            genotypes
                .into_iter()
                .map(|gt| vec![Some(Value::from(gt))])
                .collect(),
        );
        noodles_vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name(chrom)
            .set_variant_start(Position::try_from(pos).expect("position should be valid"))
            .set_ids(ids)
            .set_reference_bases("A")
            .set_alternate_bases(AlternateBases::from(vec!["G".to_owned()]))
            .set_samples(samples)
            .build()
    }

    fn record_with_ds(
        chrom: &str,
        id: &str,
        pos: usize,
        genotypes: [&str; 2],
        dosages: [Option<f32>; 2],
    ) -> noodles_vcf::variant::RecordBuf {
        let ids: Ids = [id.to_owned()].into_iter().collect();
        let keys: Keys = [String::from(key::GENOTYPE), "DS".to_owned()]
            .into_iter()
            .collect();
        let samples = Samples::new(
            keys,
            genotypes
                .into_iter()
                .zip(dosages)
                .map(|(gt, ds)| vec![Some(Value::from(gt)), ds.map(Value::from)])
                .collect(),
        );
        noodles_vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name(chrom)
            .set_variant_start(Position::try_from(pos).expect("position should be valid"))
            .set_ids(ids)
            .set_reference_bases("A")
            .set_alternate_bases(AlternateBases::from(vec!["G".to_owned()]))
            .set_samples(samples)
            .build()
    }

    fn write_fixture(path: &Path, records: &[noodles_vcf::variant::RecordBuf]) {
        let file = fs::File::create(path).expect("test BCF should be created");
        let mut writer = noodles_bcf::io::Writer::new(file);
        let ds_format = Map::<Format>::builder()
            .set_number(Number::Count(1))
            .set_type(Type::Float)
            .set_description("Expected alternate allele dosage")
            .build()
            .expect("DS format should build");
        let header = noodles_vcf::Header::builder()
            .add_contig("1", Map::<Contig>::new())
            .add_contig("2", Map::<Contig>::new())
            .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
            .add_format("DS", ds_format)
            .add_sample_name("s1")
            .add_sample_name("s2")
            .build();
        writer
            .write_header(&header)
            .expect("test BCF header should be written");
        for record in records {
            writer
                .write_variant_record(&header, record)
                .expect("test BCF record should be written");
        }
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

    fn chrom_filter(chrom: &str) -> VariantFilter {
        VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "chrom",
            "params": {"value": chrom}
        }))
        .expect("chromosome filter should parse")
    }

    fn stateless_block(
        path: &Path,
        options: &BlockReadOptions,
        start: usize,
        len: usize,
    ) -> BlockOutput {
        let window = Some(VariantWindow { start, len });
        match (options.matrix_kind, options.sparse, options.dosage_source) {
            (MatrixKind::Genotype, false, DosageSource::Hardcall) => BlockOutput::Dense(
                crate::vcf::bcf::source::read_dense_windowed(
                    path,
                    options.requested_samples.as_deref(),
                    options.variant_filter.as_ref(),
                    window,
                    options.missing_policy,
                    options.return_samples,
                    options.return_variants,
                )
                .expect("stateless dense GT window should decode"),
            ),
            (MatrixKind::Genotype, false, DosageSource::Dosage) => BlockOutput::Dense(
                crate::vcf::bcf::source::read_dosage_dense_windowed(
                    path,
                    options.requested_samples.as_deref(),
                    options.variant_filter.as_ref(),
                    window,
                    options.missing_policy,
                    options.return_samples,
                    options.return_variants,
                )
                .expect("stateless dense DS window should decode"),
            ),
            (MatrixKind::Genotype, true, DosageSource::Hardcall) => BlockOutput::Sparse(
                crate::vcf::bcf::source::read_sparse_windowed(
                    path,
                    options.requested_samples.as_deref(),
                    options.variant_filter.as_ref(),
                    window,
                    options.return_samples,
                    options.return_variants,
                )
                .expect("stateless sparse GT window should decode"),
            ),
            (MatrixKind::Haplotype, false, DosageSource::Hardcall) => BlockOutput::Dense(
                crate::vcf::bcf::source::read_haplotypes_dense_windowed(
                    path,
                    options.requested_samples.as_deref(),
                    options.variant_filter.as_ref(),
                    window,
                    options.missing_policy,
                    options.return_samples,
                    options.return_variants,
                )
                .expect("stateless dense haplotype window should decode"),
            ),
            (MatrixKind::Haplotype, true, DosageSource::Hardcall) => BlockOutput::Sparse(
                crate::vcf::bcf::source::read_haplotypes_sparse_windowed(
                    path,
                    options.requested_samples.as_deref(),
                    options.variant_filter.as_ref(),
                    window,
                    options.return_samples,
                    options.return_variants,
                )
                .expect("stateless sparse haplotype window should decode"),
            ),
            _ => panic!("unsupported stateless BCF test mode"),
        }
    }

    fn output_width(output: &BlockOutput) -> usize {
        match output {
            BlockOutput::Dense(matrix) => matrix.n_variants,
            BlockOutput::Sparse(matrix) => matrix.n_cols,
        }
    }

    fn assert_f32_values(actual: &[f32], expected: &[f32], context: &str) {
        assert_eq!(actual.len(), expected.len(), "{context}: value length");
        for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (actual.is_nan() && expected.is_nan()) || actual == expected,
                "{context}: value {index} differs: actual={actual:?}, expected={expected:?}"
            );
        }
    }

    #[test]
    fn pbr_rust_bcf_001_concrete_session_is_send() {
        fn assert_send<T: Send>() {}
        assert_send::<BcfBlockSession>();
    }

    #[test]
    fn pbr_rust_bcf_001_002_five_mode_blocks_match_stateless_windows_and_probes() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("matrix.bcf");
        write_fixture(
            &path,
            &[
                record_with_ds("1", "rs1", 10, ["0|0", "0|1"], [Some(0.1), Some(0.9)]),
                record_with_ds(
                    "2",
                    "metadata_drop",
                    20,
                    ["0|1", "1|1"],
                    [Some(1.2), Some(1.8)],
                ),
                record_with_ds("1", "rs3", 30, ["1|1", "0|0"], [Some(1.9), Some(0.2)]),
                record_with_ds("1", "rs4", 40, ["0|1", "0|0"], [Some(0.8), Some(0.1)]),
            ],
        );
        let modes = [
            (MatrixKind::Genotype, false, DosageSource::Hardcall),
            (MatrixKind::Genotype, false, DosageSource::Dosage),
            (MatrixKind::Genotype, true, DosageSource::Hardcall),
            (MatrixKind::Haplotype, false, DosageSource::Hardcall),
            (MatrixKind::Haplotype, true, DosageSource::Hardcall),
        ];

        for (matrix_kind, sparse, dosage_source) in modes {
            for block_size in [1, 2, 3] {
                let probe = BcfWorkProbe::default();
                let mut mode_options = options(Some(chrom_filter("1")));
                mode_options.matrix_kind = matrix_kind;
                mode_options.sparse = sparse;
                mode_options.dosage_source = dosage_source;
                mode_options.missing_policy = if sparse {
                    DenseMissingPolicy::Raise
                } else {
                    DenseMissingPolicy::Impute
                };
                mode_options.requested_samples = Some(vec!["s2".to_owned(), "s1".to_owned()]);
                let mut session = BcfBlockSession::open_with_probe(
                    path.clone(),
                    mode_options.clone(),
                    probe.clone(),
                )
                .expect("BCF mode session should open");
                let mut start = 0;
                while let Some(actual) = session
                    .next_block(block_size)
                    .expect("BCF mode block should decode")
                {
                    let expected = stateless_block(&path, &mode_options, start, block_size);
                    assert_eq!(
                        actual, expected,
                        "persistent BCF output must equal the stateless window oracle"
                    );
                    start += output_width(&actual);
                }
                assert_eq!(start, 3);
                let at_eof = probe.snapshot();
                assert!(session
                    .next_block(block_size)
                    .expect("BCF mode EOF should be sticky")
                    .is_none());
                assert_eq!(probe.snapshot(), at_eof);
                drop(session);

                let counts = probe.snapshot();
                assert_eq!(counts.source_opens, 1);
                assert_eq!(counts.header_parses, 1);
                assert_eq!(counts.read_record_calls, 5);
                assert_eq!(counts.candidate_visits, 4);
                match (matrix_kind, dosage_source) {
                    (MatrixKind::Genotype, DosageSource::Hardcall) => {
                        assert_eq!(counts.gt_decodes, 3);
                        assert_eq!(counts.phase_decodes, 0);
                    }
                    (MatrixKind::Genotype, DosageSource::Dosage) => {
                        assert_eq!(counts.ds_decodes, 3);
                    }
                    (MatrixKind::Haplotype, DosageSource::Hardcall) => {
                        assert_eq!(counts.gt_decodes, 0);
                        assert_eq!(counts.phase_decodes, 3);
                    }
                    (_, DosageSource::Dosage) => unreachable!(),
                }
                assert_eq!(
                    counts.max_dense_output_len,
                    if sparse {
                        0
                    } else {
                        let rows = if matrix_kind == MatrixKind::Haplotype {
                            4
                        } else {
                            2
                        };
                        rows * block_size
                    }
                );
                assert_eq!(
                    counts.max_sparse_indptr_len,
                    if sparse { block_size + 1 } else { 0 }
                );
                assert_eq!(counts.drops, 1);
            }
        }
    }

    #[test]
    fn pbr_rust_bcf_001_002_missing_policy_and_sparse_rejection_are_delayed() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("missing.bcf");
        write_fixture(
            &path,
            &[
                record_with_ds("1", "valid", 10, ["0|0", "0|1"], [Some(0.1), Some(0.9)]),
                record_with_ds("1", "missing", 20, ["0|.", "0|1"], [None, Some(1.0)]),
                record_with_ds("1", "tail", 30, ["1|1", "0|0"], [Some(1.9), Some(0.2)]),
            ],
        );

        for (matrix_kind, dosage_source) in [
            (MatrixKind::Genotype, DosageSource::Hardcall),
            (MatrixKind::Genotype, DosageSource::Dosage),
        ] {
            for missing_policy in [DenseMissingPolicy::Nan, DenseMissingPolicy::Impute] {
                let mut mode_options = options(None);
                mode_options.matrix_kind = matrix_kind;
                mode_options.dosage_source = dosage_source;
                mode_options.missing_policy = missing_policy;
                let mut session = BcfBlockSession::open(path.clone(), mode_options)
                    .expect("dense missing-policy session should open");
                let BlockOutput::Dense(matrix) = session
                    .next_block(3)
                    .expect("dense missing-policy block should decode")
                    .expect("dense missing-policy block should exist")
                else {
                    panic!("dense missing-policy route should return dense output");
                };
                let missing_value = if missing_policy == DenseMissingPolicy::Nan {
                    f32::NAN
                } else {
                    1.0
                };
                let expected = if dosage_source == DosageSource::Hardcall {
                    vec![0.0, 1.0, missing_value, 1.0, 2.0, 0.0]
                } else {
                    vec![0.1, 0.9, missing_value, 1.0, 1.9, 0.2]
                };
                let context = format!("{dosage_source:?}/{missing_policy:?}");
                assert_f32_values(&matrix.values, &expected, &context);
                assert_eq!(
                    matrix
                        .values
                        .iter()
                        .map(|value| value.is_nan())
                        .collect::<Vec<_>>(),
                    if missing_policy == DenseMissingPolicy::Nan {
                        vec![false, false, true, false, false, false]
                    } else {
                        vec![false; 6]
                    },
                    "{context}: missing mask"
                );
                assert_eq!(
                    matrix
                        .samples
                        .as_ref()
                        .expect("missing-policy BCF output should include samples")
                        .iter()
                        .map(|sample| sample.iid.as_str())
                        .collect::<Vec<_>>(),
                    vec!["s1", "s2"],
                    "{context}: sample metadata"
                );
                assert_eq!(
                    matrix
                        .variants
                        .as_ref()
                        .expect("missing-policy BCF output should include variants")
                        .positions,
                    vec![10, 20, 30],
                    "{context}: variant metadata"
                );
                assert_eq!(
                    matrix.diagnostics,
                    DenseDiagnostics {
                        requested_samples: 2,
                        retained_samples: 2,
                        missing_samples: 0,
                        candidate_variants: 3,
                        retained_variants: 3,
                        dropped_metadata_variants: 0,
                        dropped_genotype_variants: 0,
                    },
                    "{context}: diagnostics"
                );
            }
        }

        for (matrix_kind, sparse, dosage_source) in [
            (MatrixKind::Genotype, false, DosageSource::Hardcall),
            (MatrixKind::Genotype, false, DosageSource::Dosage),
            (MatrixKind::Genotype, true, DosageSource::Hardcall),
            (MatrixKind::Haplotype, false, DosageSource::Hardcall),
            (MatrixKind::Haplotype, true, DosageSource::Hardcall),
        ] {
            let probe = BcfWorkProbe::default();
            let mut mode_options = options(None);
            mode_options.matrix_kind = matrix_kind;
            mode_options.sparse = sparse;
            mode_options.dosage_source = dosage_source;
            mode_options.missing_policy = DenseMissingPolicy::Raise;
            let mut session =
                BcfBlockSession::open_with_probe(path.clone(), mode_options, probe.clone())
                    .expect("delayed missing-error session should open");
            assert!(session
                .next_block(1)
                .expect("first valid BCF block should decode")
                .is_some());
            let error = session
                .next_block(1)
                .expect_err("second BCF block should expose retained missing data");
            if matrix_kind == MatrixKind::Haplotype {
                assert!(error.to_string().contains("unphased"));
            } else {
                assert!(error.to_string().contains("missing"));
            }
            let after_error = probe.snapshot();
            assert!(session
                .next_block(1)
                .expect("failed BCF session should be sticky")
                .is_none());
            assert_eq!(probe.snapshot(), after_error);
            assert_eq!(after_error.read_record_calls, 2);
            assert_eq!(after_error.candidate_visits, 2);
        }
    }

    #[test]
    fn pbr_rust_bcf_001_002_matrix_only_all_filtered_early_drop_and_setup_errors() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("lifecycle.bcf");
        write_fixture(
            &path,
            &[
                record_with_ds("1", "rs1", 10, ["0|0", "0|1"], [Some(0.1), Some(0.9)]),
                record_with_ds("1", "rs2", 20, ["0|1", "1|1"], [Some(1.2), Some(1.8)]),
                record_with_ds("1", "rs3", 30, ["1|1", "0|0"], [Some(1.9), Some(0.2)]),
            ],
        );
        let modes = [
            (MatrixKind::Genotype, false, DosageSource::Hardcall),
            (MatrixKind::Genotype, false, DosageSource::Dosage),
            (MatrixKind::Genotype, true, DosageSource::Hardcall),
            (MatrixKind::Haplotype, false, DosageSource::Hardcall),
            (MatrixKind::Haplotype, true, DosageSource::Hardcall),
        ];

        for (matrix_kind, sparse, dosage_source) in modes {
            let early_probe = BcfWorkProbe::default();
            let mut matrix_only = options(None);
            matrix_only.matrix_kind = matrix_kind;
            matrix_only.sparse = sparse;
            matrix_only.dosage_source = dosage_source;
            matrix_only.return_samples = false;
            matrix_only.return_variants = false;
            let mut early =
                BcfBlockSession::open_with_probe(path.clone(), matrix_only, early_probe.clone())
                    .expect("matrix-only BCF session should open");
            let first = early
                .next_block(1)
                .expect("matrix-only BCF block should decode")
                .expect("matrix-only BCF block should exist");
            match first {
                BlockOutput::Dense(matrix) => {
                    assert!(matrix.samples.is_none());
                    assert!(matrix.variants.is_none());
                }
                BlockOutput::Sparse(matrix) => {
                    assert!(matrix.samples.is_none());
                    assert!(matrix.variants.is_none());
                }
            }
            drop(early);
            let early_counts = early_probe.snapshot();
            assert_eq!(early_counts.read_record_calls, 1);
            assert_eq!(early_counts.candidate_visits, 1);
            assert_eq!(early_counts.drops, 1);

            let filtered_probe = BcfWorkProbe::default();
            let mut all_filtered = options(Some(chrom_filter("absent")));
            all_filtered.matrix_kind = matrix_kind;
            all_filtered.sparse = sparse;
            all_filtered.dosage_source = dosage_source;
            let mut filtered = BcfBlockSession::open_with_probe(
                path.clone(),
                all_filtered,
                filtered_probe.clone(),
            )
            .expect("all-filtered BCF session should open");
            assert!(filtered
                .next_block(1)
                .expect("all-filtered BCF session should reach EOF")
                .is_none());
            let counts = filtered_probe.snapshot();
            assert_eq!(counts.read_record_calls, 4);
            assert_eq!(counts.candidate_visits, 3);
            assert_eq!(counts.gt_decodes, 0);
            assert_eq!(counts.ds_decodes, 0);
            assert_eq!(counts.phase_decodes, 0);
        }

        let missing_probe = BcfWorkProbe::default();
        let missing_path = dir.path().join("missing-source.bcf");
        assert!(BcfBlockSession::open_with_probe(
            missing_path,
            options(None),
            missing_probe.clone()
        )
        .is_err());
        assert_eq!(missing_probe.snapshot(), BcfWorkCounts::default());

        let invalid_probe = BcfWorkProbe::default();
        let invalid_path = dir.path().join("invalid-header.bcf");
        fs::write(&invalid_path, b"not a BCF header\n")
            .expect("invalid BCF fixture should be written");
        assert!(BcfBlockSession::open_with_probe(
            invalid_path,
            options(None),
            invalid_probe.clone()
        )
        .is_err());
        assert_eq!(
            invalid_probe.snapshot(),
            BcfWorkCounts {
                source_opens: 1,
                ..BcfWorkCounts::default()
            }
        );
    }

    #[test]
    fn pbr_rust_bcf_001_probe_counts_one_setup_linear_work_and_bounded_allocation() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("probe.bcf");
        write_fixture(
            &path,
            &[
                record("1", "rs1", 10, ["0/0", "0/1"]),
                record("2", "metadata_drop", 20, ["0/2", "0/2"]),
                record("1", "rs3", 30, ["1/1", "0/0"]),
            ],
        );
        let probe = BcfWorkProbe::default();

        {
            let mut session = BcfBlockSession::open_with_probe(
                path,
                options(Some(chrom_filter("1"))),
                probe.clone(),
            )
            .expect("BCF session should open");
            assert!(session
                .next_block(1)
                .expect("first BCF block should decode")
                .is_some());
            assert!(session
                .next_block(1)
                .expect("second BCF block should decode")
                .is_some());
            assert!(session
                .next_block(1)
                .expect("BCF session should reach EOF")
                .is_none());
            let at_eof = probe.snapshot();
            assert!(session
                .next_block(1)
                .expect("BCF EOF should be sticky")
                .is_none());
            assert_eq!(probe.snapshot(), at_eof);
        }

        assert_eq!(
            probe.snapshot(),
            BcfWorkCounts {
                source_opens: 1,
                header_parses: 1,
                read_record_calls: 4,
                candidate_visits: 3,
                gt_decodes: 2,
                ds_decodes: 0,
                phase_decodes: 0,
                max_dense_output_len: 2,
                max_sparse_indptr_len: 0,
                drops: 1,
            }
        );
    }

    #[test]
    fn pbr_rust_bcf_001_later_gt_error_is_delayed_and_stops_further_work() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("delayed.bcf");
        write_fixture(
            &path,
            &[
                record("1", "rs1", 10, ["0/0", "0/1"]),
                record("1", "bad", 20, ["0/2", "1/1"]),
                record("1", "unreached", 30, ["0/0", "0/0"]),
            ],
        );
        let probe = BcfWorkProbe::default();
        let mut session = BcfBlockSession::open_with_probe(path, options(None), probe.clone())
            .expect("BCF session should open");

        assert!(session
            .next_block(1)
            .expect("first BCF block should decode")
            .is_some());
        let error = session
            .next_block(1)
            .expect_err("second BCF block should expose malformed GT");
        assert!(error.to_string().contains("multiallelic GT allele index"));
        let after_error = probe.snapshot();
        assert!(session
            .next_block(1)
            .expect("failed BCF session should not do more work")
            .is_none());
        assert_eq!(probe.snapshot(), after_error);
        assert_eq!(after_error.read_record_calls, 2);
        assert_eq!(after_error.candidate_visits, 2);
        assert_eq!(after_error.gt_decodes, 1);
    }

    #[test]
    fn pbr_rust_bcf_002_haplotype_probe_filters_before_phase_decode() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("filter-order.bcf");
        write_fixture(
            &path,
            &[
                record("1", "retained", 10, ["0|1", "1|0"]),
                record("1", "unphased_drop", 20, ["0/0", "0/0"]),
            ],
        );
        let probe = BcfWorkProbe::default();
        let mut haplotype_options = options(Some(
            VariantFilter::from_json_value(serde_json::json!({
                "op": "predicate",
                "name": "maf",
                "params": {"min": 0.1}
            }))
            .expect("MAF filter should parse"),
        ));
        haplotype_options.matrix_kind = MatrixKind::Haplotype;
        let mut session = BcfBlockSession::open_with_probe(path, haplotype_options, probe.clone())
            .expect("BCF haplotype session should open");

        assert!(session
            .next_block(1)
            .expect("retained phased block should decode")
            .is_some());
        assert!(session
            .next_block(1)
            .expect("unphased genotype-stat reject should be skipped")
            .is_none());

        let counts = probe.snapshot();
        assert_eq!(counts.source_opens, 1);
        assert_eq!(counts.header_parses, 1);
        assert_eq!(counts.read_record_calls, 3);
        assert_eq!(counts.candidate_visits, 2);
        assert_eq!(counts.gt_decodes, 2);
        assert_eq!(counts.phase_decodes, 1);
        assert_eq!(counts.max_dense_output_len, 4);
    }

    #[test]
    fn pbr_rust_bcf_002_later_unphased_error_is_delayed_and_sticky() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("delayed-phase.bcf");
        write_fixture(
            &path,
            &[
                record("1", "rs1", 10, ["0|1", "1|0"]),
                record("1", "bad", 20, ["0/1", "1|1"]),
                record("1", "unreached", 30, ["0|0", "0|0"]),
            ],
        );
        let probe = BcfWorkProbe::default();
        let mut haplotype_options = options(None);
        haplotype_options.matrix_kind = MatrixKind::Haplotype;
        let mut session = BcfBlockSession::open_with_probe(path, haplotype_options, probe.clone())
            .expect("BCF haplotype session should open");

        assert!(session
            .next_block(1)
            .expect("first phased block should decode")
            .is_some());
        let error = session
            .next_block(1)
            .expect_err("second block should expose unphased GT");
        assert!(error.to_string().contains("unphased"));
        let after_error = probe.snapshot();
        assert!(session
            .next_block(1)
            .expect("failed BCF haplotype session should be sticky")
            .is_none());
        assert_eq!(probe.snapshot(), after_error);
        assert_eq!(after_error.read_record_calls, 2);
        assert_eq!(after_error.candidate_visits, 2);
        assert_eq!(after_error.phase_decodes, 1);
    }

    #[test]
    fn pbr_rust_bcf_002_unsupported_dosage_modes_retain_structured_errors() {
        let mut sparse_dosage = options(None);
        sparse_dosage.sparse = true;
        sparse_dosage.dosage_source = DosageSource::Dosage;
        let sparse_error = validate_bcf_options(&sparse_dosage)
            .expect_err("sparse BCF dosage should remain unsupported");
        assert!(matches!(
            sparse_error,
            GenoioError::UnsupportedRepresentation { .. }
        ));

        let mut haplotype_dosage = options(None);
        haplotype_dosage.matrix_kind = MatrixKind::Haplotype;
        haplotype_dosage.dosage_source = DosageSource::Dosage;
        let haplotype_error = validate_bcf_options(&haplotype_dosage)
            .expect_err("BCF haplotype dosage should remain unsupported");
        assert!(matches!(
            haplotype_error,
            GenoioError::UnsupportedRepresentation { .. }
        ));
    }
}
