// pattern: Imperative Shell
//! BGEN reader session and variant cursor helpers.
//!
//! A session owns the buffered file, parsed header, compression mode, and
//! reusable probability buffers for one read call. Cursor types hide the
//! difference between sequential scans and indexed byte-range reads.

use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use genoio_core::{
    select_samples_source_order, DenseDiagnostics, DenseMissingPolicy, DenseSampleSelection,
    GenoioError, SampleMetadataBuffers, SampleRecord, VariantFilter, VariantMetadataBuffers,
    VariantRecord, VariantWindow,
};

use crate::blocks::{BlockOutput, BlockReadOptions, DosageSource, MatrixKind};
use crate::error::Result;

use super::decode::{
    read_layout2_probability_payload_into, skip_layout2_probability_payload_raw,
    DosageDecodeBuffers, ProbabilityPayloadBuffers,
};
use super::header::{
    read_bgen_samples, read_layout2_variant_identifying_data, read_layout2_variant_metadata,
    skip_layout2_variant_identifying_data, BgenHeader,
};
use super::index::{validate_index_record_consumed, BgenIndexRecord};

const BGEN_READER_BUFFER_SIZE: usize = 1 << 20;

fn open_bgen_reader(bgen: &Path) -> Result<BufReader<File>> {
    let file = File::open(bgen).map_err(|source| GenoioError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
    Ok(BufReader::with_capacity(BGEN_READER_BUFFER_SIZE, file))
}

/// Shared BGEN reader state for one read call.
///
/// The session keeps path/header/compression context with the buffered reader so
/// sequential and indexed paths can share validation and payload I/O helpers.
pub(super) struct BgenReadSession {
    pub(super) reader: BufReader<File>,
    pub(super) bgen: PathBuf,
    pub(super) header: BgenHeader,
}

impl BgenReadSession {
    pub(super) fn open(bgen: &Path) -> Result<Self> {
        Self::open_owned(bgen.to_path_buf())
    }

    fn open_owned(bgen: PathBuf) -> Result<Self> {
        let mut reader = open_bgen_reader(&bgen)?;
        let header = BgenHeader::read_from(&mut reader, &bgen)?;
        header.validate(&bgen)?;
        Ok(Self {
            reader,
            bgen,
            header,
        })
    }

    pub(super) fn read_samples(&mut self, sample: Option<&Path>) -> Result<Vec<SampleRecord>> {
        read_bgen_samples(&mut self.reader, &self.bgen, sample, &self.header)
    }

    pub(super) fn read_all_variant_metadata(&mut self) -> Result<VariantMetadataBuffers> {
        read_layout2_variant_metadata(
            &mut self.reader,
            &self.bgen,
            self.header.variant_count,
            self.header.flags.compression,
        )
    }

    pub(super) fn seek_to_variants(&mut self) -> Result<()> {
        self.reader
            .seek(SeekFrom::Start(u64::from(self.header.offset) + 4))
            .map_err(|source| GenoioError::Io {
                path: self.bgen.clone(),
                source,
            })?;
        Ok(())
    }

    fn seek_to_index_record(&mut self, index_record: &BgenIndexRecord) -> Result<()> {
        self.reader
            .seek(SeekFrom::Start(index_record.file_start_position))
            .map_err(|source| GenoioError::Io {
                path: self.bgen.clone(),
                source,
            })?;
        Ok(())
    }

    pub(super) fn read_variant(&mut self) -> Result<VariantRecord> {
        read_layout2_variant_identifying_data(&mut self.reader, &self.bgen)
    }

    pub(super) fn skip_variant(&mut self) -> Result<()> {
        skip_layout2_variant_identifying_data(&mut self.reader, &self.bgen)
    }

    pub(super) fn read_payload_into(
        &mut self,
        buffers: &mut ProbabilityPayloadBuffers,
    ) -> Result<()> {
        read_layout2_probability_payload_into(
            &mut self.reader,
            &self.bgen,
            self.header.flags.compression,
            buffers,
        )
    }

    pub(super) fn skip_payload(&mut self) -> Result<()> {
        skip_layout2_probability_payload_raw(
            &mut self.reader,
            &self.bgen,
            self.header.flags.compression,
        )
    }

    fn validate_index_record_consumed(&mut self, index_record: &BgenIndexRecord) -> Result<()> {
        validate_index_record_consumed(&mut self.reader, &self.bgen, index_record)
    }
}

/// Persistent sequential BGEN dosage state.
pub(crate) struct BgenBlockSession {
    pub(super) io: BgenReadSession,
    pub(super) selection: DenseSampleSelection,
    pub(super) diagnostics: DenseDiagnostics,
    pub(super) variant_filter: Option<VariantFilter>,
    pub(super) missing_policy: DenseMissingPolicy,
    pub(super) return_samples: bool,
    pub(super) return_variants: bool,
    pub(super) remaining_variants: u32,
    pub(super) retained_skip: usize,
    pub(super) decode_buffers: DosageDecodeBuffers,
    pub(super) eof: bool,
    #[cfg(test)]
    probe: Option<BgenWorkProbe>,
}

impl BgenBlockSession {
    pub(crate) fn open(
        bgen: PathBuf,
        sample: Option<PathBuf>,
        options: BlockReadOptions,
    ) -> Result<Self> {
        Self::open_windowed(bgen, sample, options, 0)
    }

    pub(super) fn open_windowed(
        bgen: PathBuf,
        sample: Option<PathBuf>,
        options: BlockReadOptions,
        retained_skip: usize,
    ) -> Result<Self> {
        validate_bgen_options(&options)?;
        let mut io = BgenReadSession::open_owned(bgen)?;
        let all_samples = io.read_samples(sample.as_deref())?;
        let selection = select_samples_source_order(
            &all_samples,
            options.requested_samples.as_deref(),
            &io.bgen,
        )?;
        io.seek_to_variants()?;
        let remaining_variants = io.header.variant_count;
        let eof = options
            .variant_filter
            .as_ref()
            .is_some_and(VariantFilter::is_always_false);
        let diagnostics = selection.diagnostics.clone();

        Ok(Self {
            io,
            selection,
            diagnostics,
            variant_filter: options.variant_filter,
            missing_policy: options.missing_policy,
            return_samples: options.return_samples,
            return_variants: options.return_variants,
            remaining_variants,
            retained_skip,
            decode_buffers: DosageDecodeBuffers::default(),
            eof,
            #[cfg(test)]
            probe: None,
        })
    }

    #[cfg(test)]
    fn open_with_probe(
        bgen: PathBuf,
        sample: Option<PathBuf>,
        options: BlockReadOptions,
        probe: BgenWorkProbe,
    ) -> Result<Self> {
        let sample_was_supplied = sample.is_some();
        let mut session = Self::open(bgen, sample, options)?;
        probe.record_bgen_open();
        if sample_was_supplied {
            probe.record_sample_open();
        }
        session.probe = Some(probe);
        Ok(session)
    }

    pub(crate) fn next_block(&mut self, block_size: usize) -> Result<Option<BlockOutput>> {
        self.next_dosage_block(block_size)
            .map(|matrix| matrix.map(BlockOutput::Dense))
    }

    pub(super) fn empty_genotype_output(&self) -> Result<genoio_core::DenseGenotypeMatrix> {
        let samples = SampleMetadataBuffers::optional_from_records(
            &self.selection.samples,
            self.return_samples,
            false,
        )?;
        let variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(0));
        genoio_core::DenseGenotypeMatrix::new_with_layout(
            self.selection.samples.len(),
            0,
            Vec::new(),
            genoio_core::DenseLayout::SampleMajor,
            samples,
            variants,
            crate::blocks::block_diagnostics_snapshot(&self.diagnostics, 0),
        )
    }

    #[cfg(test)]
    pub(super) fn record_candidate_visit(&self) {
        if let Some(probe) = &self.probe {
            probe.record_candidate_visit();
        }
    }

    #[cfg(not(test))]
    pub(super) fn record_candidate_visit(&self) {}

    #[cfg(test)]
    pub(super) fn record_payload_decode(&self) {
        if let Some(probe) = &self.probe {
            probe.record_payload_decode();
        }
    }

    #[cfg(not(test))]
    pub(super) fn record_payload_decode(&self) {}

    #[cfg(test)]
    pub(super) fn record_dense_allocation(&self, len: usize) {
        if let Some(probe) = &self.probe {
            probe.record_dense_allocation(len);
        }
    }

    #[cfg(not(test))]
    pub(super) fn record_dense_allocation(&self, _len: usize) {}
}

impl Drop for BgenBlockSession {
    fn drop(&mut self) {
        #[cfg(test)]
        if let Some(probe) = &self.probe {
            probe.record_drop();
        }
    }
}

fn validate_bgen_options(options: &BlockReadOptions) -> Result<()> {
    if options.sparse {
        return Err(GenoioError::unsupported(
            "bgen block reads do not support sparse output",
        ));
    }
    if options.dosage_source != DosageSource::Dosage {
        return Err(GenoioError::unsupported(
            "bgen block reads support dosage values only",
        ));
    }
    if options.matrix_kind != MatrixKind::Genotype {
        return Err(GenoioError::unsupported(
            "bgen haplotype block reads are not available yet",
        ));
    }
    Ok(())
}

#[cfg(test)]
#[derive(Debug, Clone, Default)]
struct BgenWorkProbe {
    counts: std::sync::Arc<std::sync::Mutex<BgenWorkCounts>>,
}

#[cfg(test)]
impl BgenWorkProbe {
    fn snapshot(&self) -> BgenWorkCounts {
        self.counts
            .lock()
            .expect("bgen work probe lock should not be poisoned")
            .clone()
    }

    fn update(&self, update: impl FnOnce(&mut BgenWorkCounts)) {
        update(
            &mut self
                .counts
                .lock()
                .expect("bgen work probe lock should not be poisoned"),
        );
    }

    fn record_bgen_open(&self) {
        self.update(|counts| counts.bgen_opens += 1);
    }

    fn record_sample_open(&self) {
        self.update(|counts| counts.sample_opens += 1);
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

    fn record_drop(&self) {
        self.update(|counts| counts.drops += 1);
    }
}

#[cfg(test)]
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct BgenWorkCounts {
    bgen_opens: usize,
    sample_opens: usize,
    index_opens: usize,
    candidate_visits: usize,
    payload_decodes: usize,
    max_dense_output_len: usize,
    drops: usize,
}

/// Position of the current BGEN record in sequential or indexed reads.
///
/// Indexed positions carry the expected byte range so callers can verify that
/// variant and payload reads consumed the index record exactly.
#[derive(Debug, Clone, Copy)]
pub(super) enum BgenRecordPosition<'a> {
    Sequential,
    Indexed(&'a BgenIndexRecord),
}

impl BgenRecordPosition<'_> {
    pub(super) fn validate_if_indexed(self, session: &mut BgenReadSession) -> Result<()> {
        if let Self::Indexed(index_record) = self {
            session.validate_index_record_consumed(index_record)?;
        }
        Ok(())
    }
}

/// Static-dispatch cursor over either the full BGEN stream or `.bgi` byte ranges.
///
/// Indexed reads seek before yielding each record; sequential reads only count
/// down records already positioned after the header/sample block.
pub(super) enum BgenVariantCursor<'a> {
    Sequential {
        remaining: u32,
    },
    Indexed {
        records: &'a [BgenIndexRecord],
        next_index: usize,
    },
}

impl<'a> BgenVariantCursor<'a> {
    pub(super) fn sequential(variant_count: u32) -> Self {
        Self::Sequential {
            remaining: variant_count,
        }
    }

    pub(super) fn indexed(records: &'a [BgenIndexRecord]) -> Self {
        Self::Indexed {
            records,
            next_index: 0,
        }
    }

    pub(super) fn next(
        &mut self,
        session: &mut BgenReadSession,
    ) -> Result<Option<BgenRecordPosition<'a>>> {
        match self {
            Self::Sequential { remaining } => {
                if *remaining == 0 {
                    return Ok(None);
                }
                *remaining -= 1;
                Ok(Some(BgenRecordPosition::Sequential))
            }
            Self::Indexed {
                records,
                next_index,
            } => {
                let Some(index_record) = records.get(*next_index) else {
                    return Ok(None);
                };
                *next_index += 1;
                session.seek_to_index_record(index_record)?;
                Ok(Some(BgenRecordPosition::Indexed(index_record)))
            }
        }
    }
}

/// Indexed read state shared by BGEN dense output loops.
pub(super) struct BgenIndexedReadContext<'a> {
    pub(super) session: &'a mut BgenReadSession,
    pub(super) selection: DenseSampleSelection,
    pub(super) diagnostics: DenseDiagnostics,
    pub(super) variant_filter: Option<&'a VariantFilter>,
    pub(super) variant_window: Option<VariantWindow>,
    pub(super) missing_policy: DenseMissingPolicy,
    pub(super) return_samples: bool,
    pub(super) return_variants: bool,
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Write;
    use std::path::Path;

    use genoio_core::{DenseMissingPolicy, VariantFilter};
    use serde_json::json;

    use crate::blocks::{BlockOutput, BlockReadOptions, DosageSource, MatrixKind};

    use super::{BgenBlockSession, BgenWorkProbe};

    const FLAG_LAYOUT2: u32 = 2 << 2;

    fn push_u16_string(bytes: &mut Vec<u8>, value: &str) {
        bytes.extend_from_slice(
            &u16::try_from(value.len())
                .expect("test string length should fit u16")
                .to_le_bytes(),
        );
        bytes.extend_from_slice(value.as_bytes());
    }

    fn push_u32_string(bytes: &mut Vec<u8>, value: &str) {
        bytes.extend_from_slice(
            &u32::try_from(value.len())
                .expect("test string length should fit u32")
                .to_le_bytes(),
        );
        bytes.extend_from_slice(value.as_bytes());
    }

    fn push_variant(bytes: &mut Vec<u8>, id: &str, chrom: &str, pos: u32, dosage: (u8, u8)) {
        push_u16_string(bytes, id);
        push_u16_string(bytes, id);
        push_u16_string(bytes, chrom);
        bytes.extend_from_slice(&pos.to_le_bytes());
        bytes.extend_from_slice(&2_u16.to_le_bytes());
        push_u32_string(bytes, "A");
        push_u32_string(bytes, "G");

        let mut probability = Vec::new();
        probability.extend_from_slice(&1_u32.to_le_bytes());
        probability.extend_from_slice(&2_u16.to_le_bytes());
        probability.extend_from_slice(&[2, 2, 2, 0, 8, dosage.0, dosage.1]);
        bytes.extend_from_slice(
            &u32::try_from(probability.len())
                .expect("test probability length should fit u32")
                .to_le_bytes(),
        );
        bytes.extend_from_slice(&probability);
    }

    fn write_two_variant_bgen_without_embedded_samples(bgen: &Path, sample: &Path) {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&2_u32.to_le_bytes());
        bytes.extend_from_slice(&1_u32.to_le_bytes());
        bytes.extend_from_slice(b"bgen");
        bytes.extend_from_slice(&FLAG_LAYOUT2.to_le_bytes());
        push_variant(&mut bytes, "rs1", "1", 10, (255, 0));
        push_variant(&mut bytes, "rs2", "2", 20, (0, 255));
        fs::write(bgen, bytes).expect("test bgen should be written");

        let mut sample_file = fs::File::create(sample).expect("test sample file should open");
        sample_file
            .write_all(b"ID_1 ID_2 missing\n0 0 0\nsample_1 sample_1 0\n")
            .expect("test sample file should be written");
    }

    fn options(filter: Option<VariantFilter>) -> BlockReadOptions {
        BlockReadOptions {
            matrix_kind: MatrixKind::Genotype,
            sparse: false,
            requested_samples: None,
            variant_filter: filter,
            dosage_source: DosageSource::Dosage,
            missing_policy: DenseMissingPolicy::Nan,
            return_samples: true,
            return_variants: true,
        }
    }

    #[test]
    fn pbr_rust_bgen_001_work_probe_counts_opens_visits_decodes_eof_and_drop() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
        let filter = VariantFilter::from_json_value(json!({
            "op": "predicate",
            "name": "chrom",
            "params": {"value": "2"}
        }))
        .expect("test filter should parse");
        let probe = BgenWorkProbe::default();

        {
            let mut session = BgenBlockSession::open_with_probe(
                bgen,
                Some(sample),
                options(Some(filter)),
                probe.clone(),
            )
            .expect("persistent bgen session should open");
            let output = session
                .next_block(1)
                .expect("persistent bgen block should decode")
                .expect("one variant should be retained");
            let BlockOutput::Dense(matrix) = output else {
                panic!("bgen dosage session should return dense output");
            };
            assert_eq!(matrix.n_variants, 1);
            assert!(session
                .next_block(1)
                .expect("persistent bgen session should reach EOF")
                .is_none());
            let at_eof = probe.snapshot();
            assert!(session
                .next_block(1)
                .expect("persistent bgen EOF should be sticky")
                .is_none());
            assert_eq!(probe.snapshot(), at_eof);
        }

        let counts = probe.snapshot();
        assert_eq!(counts.bgen_opens, 1);
        assert_eq!(counts.sample_opens, 1);
        assert_eq!(counts.index_opens, 0);
        assert_eq!(counts.candidate_visits, 2);
        assert_eq!(counts.payload_decodes, 1);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_bgen_001_early_drop_stops_sequential_work() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
        let probe = BgenWorkProbe::default();

        {
            let mut session =
                BgenBlockSession::open_with_probe(bgen, Some(sample), options(None), probe.clone())
                    .expect("persistent bgen session should open");
            assert!(session
                .next_block(1)
                .expect("first persistent bgen block should decode")
                .is_some());
        }

        let counts = probe.snapshot();
        assert_eq!(counts.candidate_visits, 1);
        assert_eq!(counts.payload_decodes, 1);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_bgen_001_pbr_rust_alloc_001_dense_allocation_is_block_bounded() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
        let probe = BgenWorkProbe::default();
        let mut session =
            BgenBlockSession::open_with_probe(bgen, Some(sample), options(None), probe.clone())
                .expect("persistent bgen session should open");

        while session
            .next_block(1)
            .expect("persistent bgen block should decode")
            .is_some()
        {}

        assert_eq!(probe.snapshot().max_dense_output_len, 1);
    }
}
