// pattern: Imperative Shell
//! Persistent sequential and indexed text-VCF block-reader sessions.

use std::io::{self, BufRead, Read};
use std::path::PathBuf;

use genoio_core::{
    append_sparse_column, reject_sparse_missing, DenseDiagnostics, DenseGenotypeMatrix,
    DenseSampleSelection, GenoioError, PartialFilterDecision, RegionPredicate,
    SampleMetadataBuffers, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};
use noodles_bgzf as bgzf;
use noodles_vcf as noodles;

use crate::blocks::{
    block_diagnostics_snapshot, checked_dense_block_len, checked_sparse_indptr_len, BlockOutput,
    BlockReadOptions, DosageSource, MatrixKind,
};
use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::ds::{decode_ds_record, DsDecodeBuffers};
use super::gt::{
    decode_gt_record, decode_phased_gt_dense_record, decode_phased_gt_sparse_record,
    GtDecodeBuffers, GtStatsMode, HaplotypeDenseDecodeBuffers, HaplotypeSparseDecodeBuffers,
};
use super::record::{
    skip_variant_for_region, text_variant_view_from_text_record, validate_biallelic_variant,
};
use super::source::{
    index_chunks_for_region, open_bgzf_reader, open_text_vcf_input,
    open_text_vcf_input_from_reader, IndexChunk, TextVcfInput, TextVcfSource,
};
use super::sparse::{append_haplotype_minor_sparse_column, flip_values_to_minor_allele};
use super::{
    evaluate_text_gt_filter, haplotype_sample_records, write_dense_text_variant, TextDenseOutput,
    VariantMetadataSink, VariantMetadataSinkKind, VcfMetadataReturn,
};
use crate::vcf::is_compressed_vcf;
#[cfg(test)]
use crate::vcf::policy::companion_index_path;
use crate::vcf::policy::{has_vcf_index, reject_unindexed_compressed_region};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TextVcfMode {
    DenseGenotype,
    DenseDosage,
    SparseGenotype,
    DenseHaplotype,
    SparseHaplotype,
}

/// Owned equivalent of the borrowing `noodles_csi::io::Query`.
///
/// Chunks have already been normalized by `BinningIndex::query`. This adapter
/// retains the BGZF reader and executes each inclusive-start/exclusive-end
/// interval without storing a self-borrowing query.
pub(crate) struct OwnedIndexedBgzfReader {
    reader: bgzf::io::Reader<std::fs::File>,
    chunks: Vec<IndexChunk>,
    next_chunk: usize,
    active_end: Option<bgzf::VirtualPosition>,
}

impl OwnedIndexedBgzfReader {
    fn new(reader: bgzf::io::Reader<std::fs::File>, chunks: Vec<IndexChunk>) -> Self {
        Self {
            reader,
            chunks,
            next_chunk: 0,
            active_end: None,
        }
    }
}

impl Read for OwnedIndexedBgzfReader {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        let source = self.fill_buf()?;
        let amount = source.len().min(output.len());
        output[..amount].copy_from_slice(&source[..amount]);
        self.consume(amount);
        Ok(amount)
    }
}

impl BufRead for OwnedIndexedBgzfReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            let end = match self.active_end {
                Some(end) => end,
                None => {
                    let Some(chunk) = self.chunks.get(self.next_chunk).copied() else {
                        return Ok(&[]);
                    };
                    self.next_chunk += 1;
                    self.reader.seek(chunk.start())?;
                    self.active_end = Some(chunk.end());
                    chunk.end()
                }
            };

            let current = self.reader.virtual_position();
            if current >= end {
                self.active_end = None;
                continue;
            }
            let max_len = if current.compressed() == end.compressed() {
                usize::from(end.uncompressed() - current.uncompressed())
            } else {
                usize::MAX
            };
            let buffer = self.reader.fill_buf()?;
            return Ok(&buffer[..buffer.len().min(max_len)]);
        }
    }

    fn consume(&mut self, amount: usize) {
        self.reader.consume(amount);
    }
}

/// Persistent text-VCF state over one plain or gzip/multimember source.
pub(crate) enum TextVcfBlockSession {
    Plain(SequentialTextVcfSession<std::io::BufReader<std::fs::File>>),
    Compressed(
        SequentialTextVcfSession<std::io::BufReader<flate2::read::MultiGzDecoder<std::fs::File>>>,
    ),
    Indexed(SequentialTextVcfSession<OwnedIndexedBgzfReader>),
}

impl TextVcfBlockSession {
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
        #[cfg(test)] probe: Option<TextVcfWorkProbe>,
    ) -> Result<Self> {
        let mode = validate_text_options(&options)?;
        let indexed_region = options
            .variant_filter
            .as_ref()
            .and_then(VariantFilter::concrete_region_pushdown)
            .filter(|_| is_compressed_vcf(&path) && has_vcf_index(&path));
        reject_unindexed_compressed_region(&path, options.variant_filter.as_ref())?;
        if let Some(region) = indexed_region {
            #[cfg(test)]
            return Self::open_indexed(path, options, mode, region, probe);
            #[cfg(not(test))]
            return Self::open_indexed(path, options, mode, region);
        }
        let source = open_text_vcf_input(&path, options.requested_samples.as_deref(), None)?;
        #[cfg(test)]
        if let Some(probe) = &probe {
            probe.record_source_open();
            probe.record_header_parse();
        }

        match source {
            TextVcfSource::Plain(input) => {
                #[cfg(test)]
                let session =
                    SequentialTextVcfSession::new(path, input, options, mode, None, probe);
                #[cfg(not(test))]
                let session = SequentialTextVcfSession::new(path, input, options, mode, None);
                Ok(Self::Plain(session))
            }
            TextVcfSource::Compressed(input) => {
                #[cfg(test)]
                let session =
                    SequentialTextVcfSession::new(path, input, options, mode, None, probe);
                #[cfg(not(test))]
                let session = SequentialTextVcfSession::new(path, input, options, mode, None);
                Ok(Self::Compressed(session))
            }
            TextVcfSource::ThreadedCompressed(_) => Err(GenoioError::internal_contract(
                "persistent text VCF blocks do not configure threaded input",
            )),
        }
    }

    fn open_indexed(
        path: PathBuf,
        options: BlockReadOptions,
        mode: TextVcfMode,
        region: RegionPredicate,
        #[cfg(test)] probe: Option<TextVcfWorkProbe>,
    ) -> Result<Self> {
        let chunks = index_chunks_for_region(&path, &region)?.unwrap_or_default();
        #[cfg(test)]
        if let Some(probe) = &probe {
            if companion_index_path(&path, "tbi").exists() {
                probe.record_tabix_index_open();
            } else {
                probe.record_csi_index_open();
            }
        }
        let bgzf_reader = open_bgzf_reader(&path)?;
        let input = open_text_vcf_input_from_reader(
            &path,
            options.requested_samples.as_deref(),
            noodles::io::Reader::new(bgzf_reader),
        )?;
        let TextVcfInput { reader, selection } = input;
        let input = TextVcfInput {
            reader: noodles::io::Reader::new(OwnedIndexedBgzfReader::new(
                reader.into_inner(),
                chunks,
            )),
            selection,
        };
        #[cfg(test)]
        if let Some(probe) = &probe {
            probe.record_source_open();
            probe.record_header_parse();
        }
        #[cfg(test)]
        let session =
            SequentialTextVcfSession::new(path, input, options, mode, Some(region), probe);
        #[cfg(not(test))]
        let session = SequentialTextVcfSession::new(path, input, options, mode, Some(region));
        Ok(Self::Indexed(session))
    }

    #[cfg(test)]
    fn open_with_probe(
        path: PathBuf,
        options: BlockReadOptions,
        probe: TextVcfWorkProbe,
    ) -> Result<Self> {
        Self::open_impl(path, options, Some(probe))
    }

    pub(crate) fn next_block(&mut self, block_size: usize) -> Result<Option<BlockOutput>> {
        match self {
            Self::Plain(session) => session.next_block(block_size),
            Self::Compressed(session) => session.next_block(block_size),
            Self::Indexed(session) => session.next_block(block_size),
        }
    }
}

pub(crate) struct SequentialTextVcfSession<R> {
    path: PathBuf,
    reader: noodles::io::Reader<R>,
    selection: DenseSampleSelection,
    diagnostics: DenseDiagnostics,
    variant_filter: Option<VariantFilter>,
    missing_policy: genoio_core::DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    mode: TextVcfMode,
    region: Option<RegionPredicate>,
    record: noodles::Record,
    gt_decoded: GtDecodeBuffers,
    ds_decoded: DsDecodeBuffers,
    haplotype_dense_decoded: HaplotypeDenseDecodeBuffers,
    haplotype_sparse_decoded: HaplotypeSparseDecodeBuffers,
    eof: bool,
    #[cfg(test)]
    probe: Option<TextVcfWorkProbe>,
}

impl<R: BufRead> SequentialTextVcfSession<R> {
    fn new(
        path: PathBuf,
        input: TextVcfInput<R>,
        options: BlockReadOptions,
        mode: TextVcfMode,
        region: Option<RegionPredicate>,
        #[cfg(test)] probe: Option<TextVcfWorkProbe>,
    ) -> Self {
        let TextVcfInput { reader, selection } = input;
        let n_samples = selection.source_indices.len();
        let eof = options
            .variant_filter
            .as_ref()
            .is_some_and(VariantFilter::is_always_false);
        let diagnostics = selection.diagnostics.clone();
        Self {
            path,
            reader,
            selection,
            diagnostics,
            variant_filter: options.variant_filter,
            missing_policy: options.missing_policy,
            return_samples: options.return_samples,
            return_variants: options.return_variants,
            mode,
            region,
            record: noodles::Record::default(),
            gt_decoded: GtDecodeBuffers::with_capacity(n_samples),
            ds_decoded: DsDecodeBuffers::with_capacity(n_samples),
            haplotype_dense_decoded: HaplotypeDenseDecodeBuffers::with_capacity(n_samples),
            haplotype_sparse_decoded: HaplotypeSparseDecodeBuffers::with_capacity(n_samples),
            eof,
            #[cfg(test)]
            probe,
        }
    }

    fn next_block(&mut self, block_size: usize) -> Result<Option<BlockOutput>> {
        if self.eof || block_size == 0 {
            return Ok(None);
        }
        let result = match self.mode {
            TextVcfMode::DenseGenotype => self
                .next_dense_genotype_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            TextVcfMode::DenseDosage => self
                .next_dense_dosage_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            TextVcfMode::SparseGenotype => self
                .next_sparse_genotype_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Sparse)),
            TextVcfMode::DenseHaplotype => self
                .next_dense_haplotype_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            TextVcfMode::SparseHaplotype => self
                .next_sparse_haplotype_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Sparse)),
        };
        if result.is_err() {
            self.eof = true;
        }
        result
    }

    fn next_dense_genotype_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        let n_samples = self.selection.samples.len();
        self.record_dense_allocation(checked_dense_block_len(n_samples, block_size)?);
        let mut output = TextDenseOutput::new(n_samples, block_size);
        let metadata_return = self.metadata_return();
        let mut variants = VariantMetadataSink::new(
            VariantMetadataSinkKind::for_output(metadata_return),
            block_size,
        );
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let variant = text_variant_view_from_text_record(&self.path, &self.record)?;
            if skip_variant_for_region(&variant, self.region.as_ref()) {
                continue;
            }
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision_view(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }
            validate_biallelic_variant(&self.path, &variant)?;

            let needs_genotype_decision =
                matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
            let stats_mode = match (needs_genotype_decision, metadata_return.matrix_only()) {
                (true, true) => GtStatsMode::Counts,
                (true, false) => GtStatsMode::Compute,
                (false, _) => GtStatsMode::Skip,
            };
            decode_gt_record(
                &self.path,
                &self.record,
                &self.selection.source_indices,
                stats_mode,
                &mut self.gt_decoded,
            )?;
            self.record_gt_decode();

            if needs_genotype_decision {
                let filter = self.variant_filter.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = evaluate_text_gt_filter(
                    &self.gt_decoded,
                    filter,
                    &variant,
                    metadata_return.matrix_only(),
                    "GT",
                )?;
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                if let Some(stats) = stats {
                    variants.push_view_with_stats(&variant, stats)?;
                } else {
                    variants.push_view(&variant)?;
                }
            } else {
                variants.push_view(&variant)?;
            }

            write_dense_text_variant(
                &mut output,
                self.gt_decoded.values(),
                self.gt_decoded.missing_indices(),
                self.missing_policy,
            )?;
            output_variant_count += 1;
        }

        self.finish_dense_output(output, variants, output_variant_count)
    }

    fn next_dense_dosage_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        let n_samples = self.selection.samples.len();
        self.record_dense_allocation(checked_dense_block_len(n_samples, block_size)?);
        let mut output = TextDenseOutput::new(n_samples, block_size);
        let metadata_return = self.metadata_return();
        let mut variants = VariantMetadataSink::new(
            VariantMetadataSinkKind::for_output(metadata_return),
            block_size,
        );
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let variant = text_variant_view_from_text_record(&self.path, &self.record)?;
            if skip_variant_for_region(&variant, self.region.as_ref()) {
                continue;
            }
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision_view(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }
            validate_biallelic_variant(&self.path, &variant)?;

            decode_ds_record(
                &self.path,
                &self.record,
                &self.selection.source_indices,
                &mut self.ds_decoded,
            )?;
            self.record_ds_decode();
            let needs_genotype_decision =
                matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
            if needs_genotype_decision {
                let filter = self.variant_filter.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = evaluate_dosage_filter(
                    self.ds_decoded.values(),
                    self.ds_decoded.missing_indices(),
                    filter,
                    &variant,
                    !metadata_return.matrix_only(),
                )?;
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                if let Some(stats) = stats {
                    variants.push_view_with_stats(&variant, stats)?;
                } else {
                    variants.push_view(&variant)?;
                }
            } else {
                variants.push_view(&variant)?;
            }

            write_dense_text_variant(
                &mut output,
                self.ds_decoded.values(),
                self.ds_decoded.missing_indices(),
                self.missing_policy,
            )?;
            output_variant_count += 1;
        }

        self.finish_dense_output(output, variants, output_variant_count)
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
        let metadata_return = self.metadata_return();
        let mut variants = VariantMetadataSink::new(
            VariantMetadataSinkKind::for_output(metadata_return),
            block_size,
        );
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let variant = text_variant_view_from_text_record(&self.path, &self.record)?;
            if skip_variant_for_region(&variant, self.region.as_ref()) {
                continue;
            }
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision_view(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }
            validate_biallelic_variant(&self.path, &variant)?;

            let needs_genotype_decision =
                matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
            decode_gt_record(
                &self.path,
                &self.record,
                &self.selection.source_indices,
                GtStatsMode::from_needed(needs_genotype_decision),
                &mut self.gt_decoded,
            )?;
            self.record_gt_decode();
            let mut stats_to_attach = None;
            if needs_genotype_decision {
                let stats = self.gt_decoded.stats();
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

            reject_sparse_missing(!self.gt_decoded.missing_indices().is_empty())?;
            let flipped = flip_values_to_minor_allele(self.gt_decoded.values_mut());
            append_sparse_column(
                &mut indptr,
                &mut indices,
                &mut data,
                self.gt_decoded.values(),
            )?;
            variants.push_view_with_optional_stats_and_orientation(
                &variant,
                stats_to_attach,
                flipped,
            )?;
        }

        self.finish_sparse_output(indptr, indices, data, variants, false)
    }

    fn next_dense_haplotype_block(
        &mut self,
        block_size: usize,
    ) -> Result<Option<DenseGenotypeMatrix>> {
        let n_rows = self.selection.samples.len().checked_mul(2).ok_or_else(|| {
            GenoioError::internal_contract("text VCF haplotype row count is out of range")
        })?;
        self.record_dense_allocation(checked_dense_block_len(n_rows, block_size)?);
        let mut output = TextDenseOutput::new(n_rows, block_size);
        let metadata_return = self.metadata_return();
        let mut variants = VariantMetadataSink::new(
            VariantMetadataSinkKind::for_output(metadata_return),
            block_size,
        );
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));
        let mut output_variant_count = 0_usize;

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let variant = text_variant_view_from_text_record(&self.path, &self.record)?;
            if skip_variant_for_region(&variant, self.region.as_ref()) {
                continue;
            }
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision_view(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }
            validate_biallelic_variant(&self.path, &variant)?;

            let needs_genotype_decision =
                matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
            if needs_genotype_decision {
                let stats_mode = if metadata_return.matrix_only() {
                    GtStatsMode::Counts
                } else {
                    GtStatsMode::Compute
                };
                decode_gt_record(
                    &self.path,
                    &self.record,
                    &self.selection.source_indices,
                    stats_mode,
                    &mut self.gt_decoded,
                )?;
                self.record_gt_decode();
                let filter = self.variant_filter.as_ref().ok_or_else(|| {
                    GenoioError::internal_contract("genotype decision requires a variant filter")
                })?;
                let (retain_variant, stats) = evaluate_text_gt_filter(
                    &self.gt_decoded,
                    filter,
                    &variant,
                    metadata_return.matrix_only(),
                    "haplotype",
                )?;
                match retention.genotype_decision(retain_variant, &mut self.diagnostics) {
                    RetentionAction::Include => {}
                    RetentionAction::Skip => continue,
                    RetentionAction::Stop => break,
                }
                if let Some(stats) = stats {
                    variants.push_view_with_stats(&variant, stats)?;
                } else {
                    variants.push_view(&variant)?;
                }
            } else {
                variants.push_view(&variant)?;
            }

            decode_phased_gt_dense_record(
                &self.path,
                &self.record,
                &self.selection.source_indices,
                GtStatsMode::Skip,
                &mut self.haplotype_dense_decoded,
            )?;
            self.record_phase_decode();
            write_dense_text_variant(
                &mut output,
                self.haplotype_dense_decoded.values(),
                self.haplotype_dense_decoded.missing_indices(),
                self.missing_policy,
            )?;
            output_variant_count += 1;
        }

        self.finish_haplotype_dense_output(output, variants, output_variant_count)
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
        let metadata_return = self.metadata_return();
        let mut variants = VariantMetadataSink::new(
            VariantMetadataSinkKind::for_output(metadata_return),
            block_size,
        );
        let mut retention = RetainedVariantState::new(Some(VariantWindow {
            start: 0,
            len: block_size,
        }));

        while !retention.window_is_satisfied() {
            if !self.read_next_record()? {
                break;
            }
            let variant = text_variant_view_from_text_record(&self.path, &self.record)?;
            if skip_variant_for_region(&variant, self.region.as_ref()) {
                continue;
            }
            let partial_decision = self
                .variant_filter
                .as_ref()
                .map_or(PartialFilterDecision::Accept, |filter| {
                    filter.partial_decision_view(&variant)
                });
            match retention.metadata_decision(partial_decision, &mut self.diagnostics) {
                MetadataRetentionAction::Skip => continue,
                MetadataRetentionAction::Stop => break,
                MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            }
            validate_biallelic_variant(&self.path, &variant)?;

            let needs_genotype_decision =
                matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
            let mut stats_to_attach = None;
            if needs_genotype_decision {
                decode_gt_record(
                    &self.path,
                    &self.record,
                    &self.selection.source_indices,
                    GtStatsMode::Compute,
                    &mut self.gt_decoded,
                )?;
                self.record_gt_decode();
                let stats = self.gt_decoded.stats();
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

            decode_phased_gt_sparse_record(
                &self.path,
                &self.record,
                &self.selection.source_indices,
                GtStatsMode::Skip,
                &mut self.haplotype_sparse_decoded,
            )?;
            self.record_phase_decode();
            reject_sparse_missing(self.haplotype_sparse_decoded.has_missing())?;
            let flipped = append_haplotype_minor_sparse_column(
                &mut indptr,
                &mut indices,
                &mut data,
                &self.haplotype_sparse_decoded,
            )?;
            variants.push_view_with_optional_stats_and_orientation(
                &variant,
                stats_to_attach,
                flipped,
            )?;
        }

        self.finish_sparse_output(indptr, indices, data, variants, true)
    }

    fn read_next_record(&mut self) -> Result<bool> {
        let read = self.reader.read_record(&mut self.record).map_err(|error| {
            GenoioError::invalid_source(&self.path, format!("text VCF record error: {error}"))
        })?;
        if read == 0 {
            self.eof = true;
            return Ok(false);
        }
        self.record_candidate_visit();
        Ok(true)
    }

    fn finish_dense_output(
        &self,
        output: TextDenseOutput,
        variants: VariantMetadataSink,
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
        output
            .finish(
                output_variant_count,
                samples,
                variants.into_output()?,
                block_diagnostics_snapshot(&self.diagnostics, output_variant_count),
            )
            .map(Some)
    }

    fn finish_haplotype_dense_output(
        &self,
        output: TextDenseOutput,
        variants: VariantMetadataSink,
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
        output
            .finish(
                output_variant_count,
                samples,
                variants.into_output()?,
                block_diagnostics_snapshot(&self.diagnostics, output_variant_count),
            )
            .map(Some)
    }

    fn finish_sparse_output(
        &self,
        indptr: Vec<i32>,
        indices: Vec<i32>,
        data: Vec<f32>,
        variants: VariantMetadataSink,
        haplotype: bool,
    ) -> Result<Option<SparseGenotypeMatrix>> {
        let output_variant_count = indptr.len().saturating_sub(1);
        if output_variant_count == 0 {
            return Ok(None);
        }
        let (n_rows, samples) = if haplotype {
            let haplotype_samples =
                haplotype_sample_records(&self.selection.samples, &self.selection.source_indices);
            (
                haplotype_samples.len(),
                SampleMetadataBuffers::optional_from_records(
                    &haplotype_samples,
                    self.return_samples,
                    true,
                )?,
            )
        } else {
            (
                self.selection.samples.len(),
                SampleMetadataBuffers::optional_from_records(
                    &self.selection.samples,
                    self.return_samples,
                    false,
                )?,
            )
        };
        SparseGenotypeMatrix::new(
            n_rows,
            output_variant_count,
            indptr,
            indices,
            data,
            samples,
            variants.into_output()?,
            block_diagnostics_snapshot(&self.diagnostics, output_variant_count),
        )
        .map(Some)
    }

    fn metadata_return(&self) -> VcfMetadataReturn {
        VcfMetadataReturn {
            samples: self.return_samples,
            variants: self.return_variants,
        }
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

impl<R> Drop for SequentialTextVcfSession<R> {
    fn drop(&mut self) {
        #[cfg(test)]
        if let Some(probe) = &self.probe {
            probe.record_drop();
        }
    }
}

fn validate_text_options(options: &BlockReadOptions) -> Result<TextVcfMode> {
    if options.dosage_source == DosageSource::Dosage
        && (options.matrix_kind == MatrixKind::Haplotype || options.sparse)
    {
        return Err(GenoioError::unsupported(
            "text VCF dosage blocks support dense genotype matrices only",
        ));
    }
    Ok(
        match (options.matrix_kind, options.sparse, options.dosage_source) {
            (MatrixKind::Genotype, false, DosageSource::Hardcall) => TextVcfMode::DenseGenotype,
            (MatrixKind::Genotype, false, DosageSource::Dosage) => TextVcfMode::DenseDosage,
            (MatrixKind::Genotype, true, DosageSource::Hardcall) => TextVcfMode::SparseGenotype,
            (MatrixKind::Haplotype, false, DosageSource::Hardcall) => TextVcfMode::DenseHaplotype,
            (MatrixKind::Haplotype, true, DosageSource::Hardcall) => TextVcfMode::SparseHaplotype,
            (_, _, DosageSource::Dosage) => {
                return Err(GenoioError::unsupported(
                    "text VCF dosage blocks support dense genotype matrices only",
                ));
            }
        },
    )
}

#[cfg(test)]
#[derive(Debug, Clone, Default)]
struct TextVcfWorkProbe {
    counts: std::sync::Arc<std::sync::Mutex<TextVcfWorkCounts>>,
}

#[cfg(test)]
impl TextVcfWorkProbe {
    fn snapshot(&self) -> TextVcfWorkCounts {
        self.counts
            .lock()
            .expect("text VCF probe lock should not be poisoned")
            .clone()
    }

    fn update(&self, update: impl FnOnce(&mut TextVcfWorkCounts)) {
        update(
            &mut self
                .counts
                .lock()
                .expect("text VCF probe lock should not be poisoned"),
        );
    }

    fn record_source_open(&self) {
        self.update(|counts| counts.source_opens += 1);
    }

    fn record_header_parse(&self) {
        self.update(|counts| counts.header_parses += 1);
    }

    fn record_tabix_index_open(&self) {
        self.update(|counts| counts.tabix_index_opens += 1);
    }

    fn record_csi_index_open(&self) {
        self.update(|counts| counts.csi_index_opens += 1);
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
struct TextVcfWorkCounts {
    source_opens: usize,
    header_parses: usize,
    tabix_index_opens: usize,
    csi_index_opens: usize,
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
    use std::io::{BufRead, Read, Write};
    use std::path::{Path, PathBuf};

    use genoio_core::{DenseMissingPolicy, VariantFilter};
    use noodles_core::Position;
    use noodles_csi::binning_index::index::header::ReferenceSequenceNames;
    use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles_csi::binning_index::index::reference_sequence::index::BinnedIndex;
    use noodles_csi::binning_index::index::{header::Format, Header};
    use noodles_csi::{self as csi, binning_index::Indexer};
    use noodles_tabix as tabix;

    use crate::blocks::{BlockReadOptions, DosageSource, MatrixKind};

    use super::{OwnedIndexedBgzfReader, TextVcfBlockSession, TextVcfWorkCounts, TextVcfWorkProbe};

    fn write_fixture(path: &Path, records: &str) {
        fs::write(
            path,
            format!(
                "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
{records}"
            ),
        )
        .expect("text VCF fixture should be written");
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

    fn region_filter(region: &str) -> VariantFilter {
        VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "region",
            "params": {"value": region}
        }))
        .expect("region filter should parse")
    }

    fn build_tabix_index(path: &Path) {
        let file = fs::File::open(path).expect("BGZF fixture should open");
        let mut reader = noodles_bgzf::io::Reader::new(file);
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(Header::builder().set_format(Format::Vcf).build());
        let mut line = Vec::new();
        loop {
            line.clear();
            let chunk_start = reader.virtual_position();
            let len = reader
                .read_until(b'\n', &mut line)
                .expect("BGZF fixture should read");
            if len == 0 {
                break;
            }
            if line.starts_with(b"#") {
                continue;
            }
            let chunk_end = reader.virtual_position();
            let mut fields = line.split(|byte| *byte == b'\t');
            let chrom = std::str::from_utf8(fields.next().expect("record chrom"))
                .expect("record chrom should be UTF-8");
            let pos = std::str::from_utf8(fields.next().expect("record position"))
                .expect("record position should be UTF-8")
                .parse::<usize>()
                .expect("record position should parse");
            let pos = Position::try_from(pos).expect("record position should be valid");
            indexer
                .add_record(chrom, pos, pos, Chunk::new(chunk_start, chunk_end))
                .expect("tabix record should index");
        }
        let index_path = PathBuf::from(format!("{}.tbi", path.to_string_lossy()));
        tabix::fs::write(index_path, &indexer.build()).expect("tabix index should be written");
    }

    fn write_indexed_fixture(path: &Path) {
        let file = fs::File::create(path).expect("BGZF fixture should be created");
        let mut writer = noodles_bgzf::io::Writer::new(file);
        writer
            .write_all(
                b"##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n\
1\t20\trs20\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\n\
1\t30\trs30\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/0\n",
            )
            .expect("BGZF fixture should be written");
        writer
            .try_finish()
            .expect("BGZF fixture should be finished");
        build_tabix_index(path);
    }

    fn build_csi_index(path: &Path) {
        let reference_sequence_names: ReferenceSequenceNames =
            std::iter::once("1".into()).collect();
        let header = Header::builder()
            .set_format(Format::Vcf)
            .set_reference_sequence_names(reference_sequence_names)
            .build();
        let mut indexer = Indexer::<BinnedIndex>::default().set_header(header);
        let file = fs::File::open(path).expect("BGZF fixture should open");
        let mut reader = noodles_bgzf::io::Reader::new(file);
        let mut line = Vec::new();
        loop {
            line.clear();
            let chunk_start = reader.virtual_position();
            let len = reader
                .read_until(b'\n', &mut line)
                .expect("BGZF fixture should read");
            if len == 0 {
                break;
            }
            if line.starts_with(b"#") {
                continue;
            }
            let chunk_end = reader.virtual_position();
            let mut fields = line.split(|byte| *byte == b'\t');
            let chrom = fields.next().expect("record chrom");
            assert_eq!(chrom, b"1");
            let pos = std::str::from_utf8(fields.next().expect("record position"))
                .expect("record position should be UTF-8")
                .parse::<usize>()
                .expect("record position should parse");
            let pos = Position::try_from(pos).expect("record position should be valid");
            indexer
                .add_record(
                    Some((0, pos, pos, true)),
                    Chunk::new(chunk_start, chunk_end),
                )
                .expect("CSI record should index");
        }
        let index_path = PathBuf::from(format!("{}.csi", path.to_string_lossy()));
        csi::fs::write(index_path, &indexer.build(1)).expect("CSI index should be written");
    }

    #[test]
    fn pbr_rust_textvcf_001_concrete_sequential_sessions_are_send() {
        fn assert_send<T: Send>() {}
        assert_send::<TextVcfBlockSession>();
    }

    #[test]
    fn pbr_rust_textvcf_002_owned_chunk_reader_enforces_exact_same_block_ranges() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("chunks.bgz");
        let file = fs::File::create(&path).expect("BGZF fixture should be created");
        let mut writer = noodles_bgzf::io::Writer::new(file);
        writer
            .write_all(b"prefix\n")
            .expect("prefix should be written");
        let first_start = writer.virtual_position();
        writer
            .write_all(b"first\n")
            .expect("first row should be written");
        let first_end = writer.virtual_position();
        writer
            .write_all(b"middle\n")
            .expect("middle row should be written");
        let last_start = writer.virtual_position();
        writer
            .write_all(b"last\n")
            .expect("last row should be written");
        let last_end = writer.virtual_position();
        writer.try_finish().expect("BGZF fixture should finish");
        let reader = noodles_bgzf::io::Reader::new(
            fs::File::open(&path).expect("BGZF fixture should reopen"),
        );
        let mut query = OwnedIndexedBgzfReader::new(
            reader,
            vec![
                Chunk::new(first_start, first_end),
                Chunk::new(last_start, last_end),
            ],
        );
        let mut actual = String::new();
        query
            .read_to_string(&mut actual)
            .expect("owned chunk reader should decode");

        assert_eq!(actual, "first\nlast\n");
    }

    #[test]
    fn pbr_rust_textvcf_002_tabix_probe_counts_one_source_and_only_resolved_index() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("indexed.vcf.gz");
        write_indexed_fixture(&path);
        let probe = TextVcfWorkProbe::default();

        {
            let mut session = TextVcfBlockSession::open_with_probe(
                path,
                options(Some(region_filter("1:20-30"))),
                probe.clone(),
            )
            .expect("indexed text VCF session should open");
            assert!(session
                .next_block(1)
                .expect("first indexed block should decode")
                .is_some());
            assert!(session
                .next_block(1)
                .expect("second indexed block should decode")
                .is_some());
            assert!(session
                .next_block(1)
                .expect("indexed session should reach EOF")
                .is_none());
            let at_eof = probe.snapshot();
            assert!(session
                .next_block(1)
                .expect("indexed EOF should be sticky")
                .is_none());
            assert_eq!(probe.snapshot(), at_eof);
        }
        let counts = probe.snapshot();
        assert_eq!(counts.source_opens, 1);
        assert_eq!(counts.header_parses, 1);
        assert_eq!(counts.tabix_index_opens, 1);
        assert_eq!(counts.csi_index_opens, 0);
        assert_eq!(counts.gt_decodes, 2);
        assert_eq!(counts.max_dense_output_len, 2);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_textvcf_002_csi_probe_counts_one_source_and_only_resolved_index() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("indexed.vcf.gz");
        write_indexed_fixture(&path);
        build_csi_index(&path);
        fs::remove_file(PathBuf::from(format!("{}.tbi", path.to_string_lossy())))
            .expect("tabix companion should be removed");
        let probe = TextVcfWorkProbe::default();

        {
            let mut session = TextVcfBlockSession::open_with_probe(
                path,
                options(Some(region_filter("1:20-30"))),
                probe.clone(),
            )
            .expect("CSI text VCF session should open");
            assert!(session
                .next_block(2)
                .expect("CSI block should decode")
                .is_some());
            assert!(session
                .next_block(2)
                .expect("CSI session should reach EOF")
                .is_none());
        }
        let counts = probe.snapshot();
        assert_eq!(counts.source_opens, 1);
        assert_eq!(counts.header_parses, 1);
        assert_eq!(counts.tabix_index_opens, 0);
        assert_eq!(counts.csi_index_opens, 1);
        assert_eq!(counts.gt_decodes, 2);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_textvcf_001_probe_counts_one_setup_linear_work_bounded_allocation_and_drop() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("probe.vcf");
        write_fixture(
            &path,
            "\
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DS\t0/0:0.1\t0/1:0.9
2\t20\tdrop\tC\tT\t.\tPASS\t.\tGT:DS\tbad:bad\tbad:bad
1\t30\trs3\tG\tA\t.\tPASS\t.\tGT:DS\t1/1:1.9\t0/0:0.2
",
        );
        let probe = TextVcfWorkProbe::default();

        {
            let mut session = TextVcfBlockSession::open_with_probe(
                path,
                options(Some(chrom_filter("1"))),
                probe.clone(),
            )
            .expect("text VCF session should open");
            assert!(session
                .next_block(1)
                .expect("first block should decode")
                .is_some());
            assert!(session
                .next_block(1)
                .expect("second block should decode")
                .is_some());
            assert!(session
                .next_block(1)
                .expect("session should reach EOF")
                .is_none());
            let at_eof = probe.snapshot();
            assert!(session
                .next_block(1)
                .expect("EOF should be sticky")
                .is_none());
            assert_eq!(probe.snapshot(), at_eof);
        }

        assert_eq!(
            probe.snapshot(),
            TextVcfWorkCounts {
                source_opens: 1,
                header_parses: 1,
                tabix_index_opens: 0,
                csi_index_opens: 0,
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
    fn pbr_rust_textvcf_001_later_gt_error_is_delayed_and_stops_further_work() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let path = dir.path().join("delayed.vcf");
        write_fixture(
            &path,
            "\
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t20\tbad\tC\tT\t.\tPASS\t.\tGT\t0/2\t1/1
1\t30\tunreached\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0
",
        );
        let probe = TextVcfWorkProbe::default();
        let mut session = TextVcfBlockSession::open_with_probe(path, options(None), probe.clone())
            .expect("text VCF session should open");

        assert!(session
            .next_block(1)
            .expect("first block should decode")
            .is_some());
        let error = session
            .next_block(1)
            .expect_err("second block should expose malformed GT");
        assert!(error.to_string().contains("multiallelic GT"));
        let after_error = probe.snapshot();
        assert!(session
            .next_block(1)
            .expect("failed session should not do more work")
            .is_none());
        assert_eq!(probe.snapshot(), after_error);
        assert_eq!(after_error.candidate_visits, 2);
        assert_eq!(after_error.gt_decodes, 1);
    }
}
