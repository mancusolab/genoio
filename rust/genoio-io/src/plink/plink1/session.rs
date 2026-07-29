// pattern: Imperative Shell
//! Persistent PLINK1 BED/BIM/FAM block-reader session.

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
use crate::matrix::shrink_sample_major_width;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::bed::{
    infer_bed_variant_count, open_bed_file, read_plink1_variant_packed, Plink1DecoderState,
};
use super::metadata::{parse_fam, BimRecordReader};

/// Persistent PLINK1 hard-call state over one BED/BIM/FAM source set.
pub(crate) struct Plink1BlockSession {
    bed: PathBuf,
    bim: PathBuf,
    bed_reader: File,
    bim_reader: BimRecordReader,
    selection: DenseSampleSelection,
    n_source_samples: usize,
    n_source_variants: usize,
    bytes_per_variant: usize,
    source_position: usize,
    retained_skip: usize,
    diagnostics: DenseDiagnostics,
    variant_filter: Option<VariantFilter>,
    genotype_filter_plan: GenotypeFilterPlan,
    missing_policy: DenseMissingPolicy,
    sparse: bool,
    return_samples: bool,
    return_variants: bool,
    all_samples_selected: bool,
    read_bim_records: bool,
    decoder_state: Plink1DecoderState,
    batch: HardcallBatch,
    selected_values: Vec<f32>,
    missing_indices: Vec<usize>,
    eof: bool,
    #[cfg(test)]
    probe: Option<Plink1WorkProbe>,
}

impl Plink1BlockSession {
    pub(crate) fn open(
        bed: PathBuf,
        bim: PathBuf,
        fam: PathBuf,
        options: BlockReadOptions,
    ) -> Result<Self> {
        Self::open_windowed(bed, bim, fam, options, 0, true)
    }

    pub(super) fn open_windowed(
        bed: PathBuf,
        bim: PathBuf,
        fam: PathBuf,
        options: BlockReadOptions,
        retained_skip: usize,
        read_bim_records: bool,
    ) -> Result<Self> {
        validate_plink1_options(&options)?;
        let bed_reader = open_bed_file(&bed)?;
        let all_samples = parse_fam(&fam)?;
        let selection =
            select_samples_source_order(&all_samples, options.requested_samples.as_deref(), &bed)?;
        let n_source_samples = all_samples.len();
        let bytes_per_variant = n_source_samples.div_ceil(4);
        let n_source_variants =
            infer_bed_variant_count(&bed, &bed_reader, n_source_samples, bytes_per_variant)?;
        let bim_reader = BimRecordReader::new(&bim)?;
        let decoder_state =
            Plink1DecoderState::new(n_source_samples, bytes_per_variant, selection.samples.len());
        let batch = HardcallBatch::new(n_source_samples);
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
            bed,
            bim,
            bed_reader,
            bim_reader,
            selection,
            n_source_samples,
            n_source_variants,
            bytes_per_variant,
            source_position: 0,
            retained_skip,
            diagnostics,
            variant_filter: options.variant_filter,
            genotype_filter_plan,
            missing_policy: options.missing_policy,
            sparse: options.sparse,
            return_samples: options.return_samples,
            return_variants: options.return_variants,
            all_samples_selected: options.requested_samples.is_none(),
            read_bim_records,
            decoder_state,
            batch,
            selected_values,
            missing_indices: Vec::new(),
            eof,
            #[cfg(test)]
            probe: None,
        })
    }

    #[cfg(test)]
    fn open_with_probe(
        bed: PathBuf,
        bim: PathBuf,
        fam: PathBuf,
        options: BlockReadOptions,
        probe: Plink1WorkProbe,
    ) -> Result<Self> {
        let mut session = Self::open(bed, bim, fam, options)?;
        probe.record_bed_open();
        probe.record_bim_open();
        probe.record_fam_open();
        session.probe = Some(probe);
        Ok(session)
    }

    pub(crate) fn next_block(&mut self, block_size: usize) -> Result<Option<BlockOutput>> {
        if self.sparse {
            self.next_sparse_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Sparse))
        } else {
            self.next_dense_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Dense))
        }
    }

    pub(super) fn source_record_capacity(&self) -> usize {
        self.n_source_variants
    }

    pub(super) fn empty_dense_output(&self) -> Result<DenseGenotypeMatrix> {
        let samples = SampleMetadataBuffers::optional_from_records(
            &self.selection.samples,
            self.return_samples,
            false,
        )?;
        let variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(0));
        DenseGenotypeMatrix::new_with_layout(
            self.selection.samples.len(),
            0,
            Vec::new(),
            DenseLayout::SampleMajor,
            samples,
            variants,
            block_diagnostics_snapshot(&self.diagnostics, 0),
        )
    }

    pub(super) fn empty_sparse_output(&self) -> Result<SparseGenotypeMatrix> {
        let samples = SampleMetadataBuffers::optional_from_records(
            &self.selection.samples,
            self.return_samples,
            false,
        )?;
        let variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(0));
        SparseGenotypeMatrix::new(
            self.selection.samples.len(),
            0,
            vec![0],
            Vec::new(),
            Vec::new(),
            samples,
            variants,
            block_diagnostics_snapshot(&self.diagnostics, 0),
        )
    }

    pub(super) fn next_dense_block(
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
        self.batch.clear();
        let mut batch_start = 0_usize;
        let retained_skip = std::mem::take(&mut self.retained_skip);
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: retained_skip,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            let Some((variant_index, mut variant)) = self.next_source_variant()? else {
                break;
            };
            let partial_decision = match variant.as_ref() {
                Some(variant) => self
                    .variant_filter
                    .as_ref()
                    .map_or(PartialFilterDecision::Accept, |filter| {
                        filter.partial_decision(variant)
                    }),
                None if self.variant_filter.is_none() => PartialFilterDecision::Accept,
                None => PartialFilterDecision::NeedGenotypes,
            };
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }

            self.decode_variant(variant_index)?;
            if matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
                let filter = self.variant_filter.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = match variant.as_ref() {
                    Some(variant) => evaluate_packed_hardcall_filter(
                        &self.decoder_state.packed,
                        &self.selection.source_indices,
                        self.all_samples_selected,
                        filter,
                        self.genotype_filter_plan,
                        Some(variant),
                        self.return_variants,
                    )?,
                    None => evaluate_packed_hardcall_filter::<VariantRecord>(
                        &self.decoder_state.packed,
                        &self.selection.source_indices,
                        self.all_samples_selected,
                        filter,
                        self.genotype_filter_plan,
                        None,
                        false,
                    )?,
                };
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                if let (Some(variant), Some(stats)) = (variant.as_mut(), stats) {
                    attach_variant_stats(variant, stats);
                }
            }

            if let Some(variants) = variants.as_mut() {
                variants.push_record(variant.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract(
                        "requested PLINK1 variant metadata was not parsed",
                    )
                })?)?;
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

        self.finish_source_if_at_bed_end()?;
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

    pub(super) fn next_sparse_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<SparseGenotypeMatrix>> {
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
        let retained_skip = std::mem::take(&mut self.retained_skip);
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: retained_skip,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            let Some((variant_index, variant)) = self.next_source_variant()? else {
                break;
            };
            let mut variant = variant.ok_or_else(|| {
                GenoioError::internal_contract("sparse PLINK1 reads require BIM metadata")
            })?;
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

            self.decode_variant(variant_index)?;
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

        self.finish_source_if_at_bed_end()?;
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

    fn next_source_variant(&mut self) -> Result<Option<(usize, Option<VariantRecord>)>> {
        if self.source_position >= self.n_source_variants {
            self.finish_source_if_at_bed_end()?;
            return Ok(None);
        }

        let variant = if self.read_bim_records {
            let Some((variant_index, variant)) = self.bim_reader.next_record()? else {
                return Err(GenoioError::invalid_source(
                    &self.bim,
                    "bim contains fewer variants than bed",
                ));
            };
            if variant_index != self.source_position {
                return Err(GenoioError::internal_contract(
                    "PLINK1 BIM source position is not monotonic",
                ));
            }
            Some(variant)
        } else {
            None
        };
        let variant_index = self.source_position;
        self.source_position += 1;
        self.record_candidate_visit();
        Ok(Some((variant_index, variant)))
    }

    fn finish_source_if_at_bed_end(&mut self) -> Result<()> {
        if self.source_position != self.n_source_variants || self.eof {
            return Ok(());
        }
        if self.read_bim_records && self.bim_reader.next_record()?.is_some() {
            return Err(GenoioError::invalid_source(
                &self.bim,
                "bim variant count exceeds bed variant count",
            ));
        }
        self.eof = true;
        Ok(())
    }

    fn decode_variant(&mut self, variant_index: usize) -> Result<()> {
        read_plink1_variant_packed(
            &self.bed,
            &mut self.bed_reader,
            variant_index,
            self.bytes_per_variant,
            self.n_source_samples,
            &mut self.decoder_state,
        )?;
        self.record_payload_decode();
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
    fn record_payload_decode(&self) {
        if let Some(probe) = &self.probe {
            probe.record_payload_decode();
        }
    }

    #[cfg(not(test))]
    fn record_payload_decode(&self) {}

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

impl Drop for Plink1BlockSession {
    fn drop(&mut self) {
        #[cfg(test)]
        if let Some(probe) = &self.probe {
            probe.record_drop();
        }
    }
}

fn validate_plink1_options(options: &BlockReadOptions) -> Result<()> {
    if options.matrix_kind != MatrixKind::Genotype {
        return Err(GenoioError::unsupported(
            "plink1 block reads support genotype matrices only",
        ));
    }
    if options.dosage_source != DosageSource::Hardcall {
        return Err(GenoioError::unsupported(
            "plink1 block reads support hardcall values only",
        ));
    }
    Ok(())
}

#[cfg(test)]
#[derive(Debug, Clone, Default)]
struct Plink1WorkProbe {
    counts: std::sync::Arc<std::sync::Mutex<Plink1WorkCounts>>,
}

#[cfg(test)]
impl Plink1WorkProbe {
    fn snapshot(&self) -> Plink1WorkCounts {
        self.counts
            .lock()
            .expect("plink1 work probe lock should not be poisoned")
            .clone()
    }

    fn update(&self, update: impl FnOnce(&mut Plink1WorkCounts)) {
        update(
            &mut self
                .counts
                .lock()
                .expect("plink1 work probe lock should not be poisoned"),
        );
    }

    fn record_bed_open(&self) {
        self.update(|counts| counts.bed_opens += 1);
    }

    fn record_bim_open(&self) {
        self.update(|counts| counts.bim_opens += 1);
    }

    fn record_fam_open(&self) {
        self.update(|counts| counts.fam_opens += 1);
    }

    fn record_candidate_visit(&self) {
        self.update(|counts| counts.candidate_visits += 1);
    }

    fn record_payload_decode(&self) {
        self.update(|counts| counts.payload_decodes += 1);
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
struct Plink1WorkCounts {
    bed_opens: usize,
    bim_opens: usize,
    fam_opens: usize,
    candidate_visits: usize,
    payload_decodes: usize,
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

    use super::{Plink1BlockSession, Plink1WorkProbe};

    fn write_text(path: &Path, contents: &str) {
        fs::write(path, contents).expect("test fixture should be written");
    }

    fn write_fixture(
        dir: &Path,
        bed_payload: &[u8],
        bim_contents: &str,
    ) -> (PathBuf, PathBuf, PathBuf) {
        let bed = dir.join("tiny.bed");
        let bim = dir.join("tiny.bim");
        let fam = dir.join("tiny.fam");
        let mut bed_bytes = vec![0x6c, 0x1b, 0x01];
        bed_bytes.extend_from_slice(bed_payload);
        fs::write(&bed, bed_bytes).expect("bed fixture should be written");
        write_text(&bim, bim_contents);
        write_text(&fam, "F1 S1 0 0 1 -9\nF1 S2 0 0 2 -9\n");
        (bed, bim, fam)
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

    fn sparse_options(filter: Option<VariantFilter>) -> BlockReadOptions {
        BlockReadOptions {
            sparse: true,
            missing_policy: DenseMissingPolicy::Raise,
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
    fn pbr_rust_plink1_001_work_probe_counts_opens_visits_decodes_eof_and_drop() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (bed, bim, fam) = write_fixture(
            dir.path(),
            &[0x0b, 0x2c, 0x08],
            "1 rs1 0 10 G A\n1 rs2 0 20 T C\n2 rs3 0 30 A G\n",
        );
        let probe = Plink1WorkProbe::default();

        {
            let mut session = Plink1BlockSession::open_with_probe(
                bed,
                bim,
                fam,
                options(Some(chrom_filter("2"))),
                probe.clone(),
            )
            .expect("persistent plink1 session should open");
            assert!(session
                .next_block(1)
                .expect("persistent plink1 block should decode")
                .is_some());
            assert!(session
                .next_block(1)
                .expect("persistent plink1 session should reach EOF")
                .is_none());
            let at_eof = probe.snapshot();
            assert!(session
                .next_block(1)
                .expect("persistent plink1 EOF should be sticky")
                .is_none());
            assert_eq!(probe.snapshot(), at_eof);
        }

        let counts = probe.snapshot();
        assert_eq!(counts.bed_opens, 1);
        assert_eq!(counts.bim_opens, 1);
        assert_eq!(counts.fam_opens, 1);
        assert_eq!(counts.candidate_visits, 3);
        assert_eq!(counts.payload_decodes, 1);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_plink1_001_early_drop_stops_dense_work() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (bed, bim, fam) = write_fixture(
            dir.path(),
            &[0x04, 0x0d, 0x03],
            "1 rs1 0 10 G A\n1 rs2 0 20 T C\n2 rs3 0 30 A G\n",
        );
        let probe = Plink1WorkProbe::default();

        {
            let mut session =
                Plink1BlockSession::open_with_probe(bed, bim, fam, options(None), probe.clone())
                    .expect("persistent plink1 session should open");
            assert!(session
                .next_block(1)
                .expect("first persistent plink1 block should decode")
                .is_some());
        }

        let counts = probe.snapshot();
        assert_eq!(counts.candidate_visits, 1);
        assert_eq!(counts.payload_decodes, 1);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_plink1_001_rejects_shorter_and_longer_bim_companions_while_advancing() {
        let short_dir = tempfile::tempdir().expect("test directory should be created");
        let (short_bed, short_bim, short_fam) =
            write_fixture(short_dir.path(), &[0x04, 0x0d], "1 rs1 0 10 G A\n");
        let mut short_session =
            Plink1BlockSession::open(short_bed, short_bim, short_fam, options(None))
                .expect("short-bim session should open");
        let short_error = short_session
            .next_block(2)
            .expect_err("shorter bim companion should fail while advancing");

        let long_dir = tempfile::tempdir().expect("test directory should be created");
        let (long_bed, long_bim, long_fam) = write_fixture(
            long_dir.path(),
            &[0x04, 0x0d],
            "1 rs1 0 10 G A\n1 rs2 0 20 T C\n2 rs3 0 30 A G\n",
        );
        let mut long_session =
            Plink1BlockSession::open(long_bed, long_bim, long_fam, options(None))
                .expect("long-bim session should open");
        let long_error = long_session
            .next_block(2)
            .expect_err("longer bim companion should fail while advancing");

        assert!(short_error.to_string().contains("fewer"));
        assert!(long_error.to_string().contains("exceeds"));
    }

    #[test]
    fn pbr_rust_plink1_001_pbr_rust_alloc_001_dense_allocation_is_block_bounded() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (bed, bim, fam) = write_fixture(
            dir.path(),
            &[0x04, 0x0d, 0x03],
            "1 rs1 0 10 G A\n1 rs2 0 20 T C\n2 rs3 0 30 A G\n",
        );
        let probe = Plink1WorkProbe::default();
        let mut session =
            Plink1BlockSession::open_with_probe(bed, bim, fam, options(None), probe.clone())
                .expect("persistent plink1 session should open");

        while session
            .next_block(1)
            .expect("persistent plink1 block should decode")
            .is_some()
        {}

        assert_eq!(probe.snapshot().max_dense_output_len, 2);
    }

    #[test]
    fn pbr_rust_plink1_002_sparse_probe_counts_linear_work_and_checked_block_allocation() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let (bed, bim, fam) = write_fixture(
            dir.path(),
            &[0x0b, 0x2c, 0x08],
            "1 rs1 0 10 G A\n1 rs2 0 20 T C\n2 rs3 0 30 A G\n",
        );
        let probe = Plink1WorkProbe::default();

        {
            let mut session = Plink1BlockSession::open_with_probe(
                bed,
                bim,
                fam,
                sparse_options(None),
                probe.clone(),
            )
            .expect("persistent sparse plink1 session should open");
            while session
                .next_block(1)
                .expect("persistent sparse plink1 block should decode")
                .is_some()
            {}
        }

        let counts = probe.snapshot();
        assert_eq!(counts.bed_opens, 1);
        assert_eq!(counts.bim_opens, 1);
        assert_eq!(counts.fam_opens, 1);
        assert_eq!(counts.candidate_visits, 3);
        assert_eq!(counts.payload_decodes, 3);
        assert_eq!(counts.max_sparse_indptr_len, 2);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_plink1_002_separate_sessions_keep_independent_cursors_and_probes() {
        let first_dir = tempfile::tempdir().expect("test directory should be created");
        let second_dir = tempfile::tempdir().expect("test directory should be created");
        let first_paths = write_fixture(
            first_dir.path(),
            &[0x0b, 0x2c, 0x08],
            "1 rs1 0 10 G A\n1 rs2 0 20 T C\n2 rs3 0 30 A G\n",
        );
        let second_paths = write_fixture(
            second_dir.path(),
            &[0x0b, 0x2c, 0x08],
            "1 rs1 0 10 G A\n1 rs2 0 20 T C\n2 rs3 0 30 A G\n",
        );
        let first_probe = Plink1WorkProbe::default();
        let second_probe = Plink1WorkProbe::default();
        let mut first = Plink1BlockSession::open_with_probe(
            first_paths.0,
            first_paths.1,
            first_paths.2,
            options(None),
            first_probe.clone(),
        )
        .expect("first persistent plink1 session should open");
        let mut second = Plink1BlockSession::open_with_probe(
            second_paths.0,
            second_paths.1,
            second_paths.2,
            sparse_options(None),
            second_probe.clone(),
        )
        .expect("second persistent sparse plink1 session should open");

        assert!(first
            .next_block(1)
            .expect("first session should advance")
            .is_some());
        assert!(second
            .next_block(1)
            .expect("second session should advance")
            .is_some());
        drop(first);
        assert!(second
            .next_block(1)
            .expect("second session should remain usable")
            .is_some());

        assert_eq!(first_probe.snapshot().candidate_visits, 1);
        assert_eq!(first_probe.snapshot().drops, 1);
        assert_eq!(second_probe.snapshot().candidate_visits, 2);
        assert_eq!(second_probe.snapshot().drops, 0);
    }
}
