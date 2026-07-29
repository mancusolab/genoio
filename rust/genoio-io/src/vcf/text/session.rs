// pattern: Imperative Shell
//! Persistent sequential and indexed text-VCF block-reader sessions.

use std::io::BufRead;
use std::path::PathBuf;

use genoio_core::{
    DenseDiagnostics, DenseGenotypeMatrix, DenseSampleSelection, GenoioError,
    PartialFilterDecision, SampleMetadataBuffers, VariantFilter, VariantWindow,
};
use noodles_vcf as noodles;

use crate::blocks::{
    block_diagnostics_snapshot, checked_dense_block_len, BlockOutput, BlockReadOptions,
    DosageSource, MatrixKind,
};
use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::ds::{decode_ds_record, DsDecodeBuffers};
use super::gt::{decode_gt_record, GtDecodeBuffers, GtStatsMode};
use super::record::{text_variant_view_from_text_record, validate_biallelic_variant};
use super::source::{open_text_vcf_input, TextVcfInput, TextVcfSource};
use super::{
    evaluate_text_gt_filter, write_dense_text_variant, TextDenseOutput, VariantMetadataSink,
    VariantMetadataSinkKind, VcfMetadataReturn,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TextVcfMode {
    DenseGenotype,
    DenseDosage,
}

/// Persistent text-VCF state over one plain or gzip/multimember source.
pub(crate) enum TextVcfBlockSession {
    Plain(SequentialTextVcfSession<std::io::BufReader<std::fs::File>>),
    Compressed(
        SequentialTextVcfSession<std::io::BufReader<flate2::read::MultiGzDecoder<std::fs::File>>>,
    ),
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
        let mode = validate_sequential_text_options(&options)?;
        let source = open_text_vcf_input(&path, options.requested_samples.as_deref(), None)?;
        #[cfg(test)]
        if let Some(probe) = &probe {
            probe.record_source_open();
            probe.record_header_parse();
        }

        match source {
            TextVcfSource::Plain(input) => {
                #[cfg(test)]
                let session = SequentialTextVcfSession::new(path, input, options, mode, probe);
                #[cfg(not(test))]
                let session = SequentialTextVcfSession::new(path, input, options, mode);
                Ok(Self::Plain(session))
            }
            TextVcfSource::Compressed(input) => {
                #[cfg(test)]
                let session = SequentialTextVcfSession::new(path, input, options, mode, probe);
                #[cfg(not(test))]
                let session = SequentialTextVcfSession::new(path, input, options, mode);
                Ok(Self::Compressed(session))
            }
            TextVcfSource::ThreadedCompressed(_) => Err(GenoioError::internal_contract(
                "persistent text VCF blocks do not configure threaded input",
            )),
        }
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
    record: noodles::Record,
    gt_decoded: GtDecodeBuffers,
    ds_decoded: DsDecodeBuffers,
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
            record: noodles::Record::default(),
            gt_decoded: GtDecodeBuffers::with_capacity(n_samples),
            ds_decoded: DsDecodeBuffers::with_capacity(n_samples),
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
            TextVcfMode::DenseGenotype => self.next_dense_genotype_block(block_size),
            TextVcfMode::DenseDosage => self.next_dense_dosage_block(block_size),
        };
        if result.is_err() {
            self.eof = true;
        }
        result.map(|matrix| matrix.map(BlockOutput::Dense))
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
    fn record_dense_allocation(&self, len: usize) {
        if let Some(probe) = &self.probe {
            probe.record_dense_allocation(len);
        }
    }

    #[cfg(not(test))]
    fn record_dense_allocation(&self, _len: usize) {}
}

impl<R> Drop for SequentialTextVcfSession<R> {
    fn drop(&mut self) {
        #[cfg(test)]
        if let Some(probe) = &self.probe {
            probe.record_drop();
        }
    }
}

fn validate_sequential_text_options(options: &BlockReadOptions) -> Result<TextVcfMode> {
    if options.matrix_kind != MatrixKind::Genotype {
        return Err(GenoioError::unsupported(
            "text VCF haplotype block reads are not implemented yet",
        ));
    }
    if options.sparse {
        return Err(GenoioError::unsupported(
            "text VCF sparse block reads are not implemented yet",
        ));
    }
    Ok(match options.dosage_source {
        DosageSource::Hardcall => TextVcfMode::DenseGenotype,
        DosageSource::Dosage => TextVcfMode::DenseDosage,
    })
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

    fn record_candidate_visit(&self) {
        self.update(|counts| counts.candidate_visits += 1);
    }

    fn record_gt_decode(&self) {
        self.update(|counts| counts.gt_decodes += 1);
    }

    fn record_ds_decode(&self) {
        self.update(|counts| counts.ds_decodes += 1);
    }

    fn record_dense_allocation(&self, len: usize) {
        self.update(|counts| counts.max_dense_output_len = counts.max_dense_output_len.max(len));
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
    candidate_visits: usize,
    gt_decodes: usize,
    ds_decodes: usize,
    max_dense_output_len: usize,
    drops: usize,
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::Path;

    use genoio_core::{DenseMissingPolicy, VariantFilter};

    use crate::blocks::{BlockReadOptions, DosageSource, MatrixKind};

    use super::{TextVcfBlockSession, TextVcfWorkCounts, TextVcfWorkProbe};

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

    #[test]
    fn pbr_rust_textvcf_001_concrete_sequential_sessions_are_send() {
        fn assert_send<T: Send>() {}
        assert_send::<TextVcfBlockSession>();
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
                candidate_visits: 3,
                gt_decodes: 2,
                ds_decodes: 0,
                max_dense_output_len: 2,
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
