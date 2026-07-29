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
    VariantRecord,
};

use crate::blocks::{BlockOutput, BlockReadOptions, DosageSource, MatrixKind};
use crate::error::Result;

use super::decode::{
    read_layout2_probability_payload_into, skip_layout2_probability_payload_raw,
    DosageDecodeBuffers, HaplotypeDecodeBuffers, ProbabilityPayloadBuffers,
};
use super::header::{
    read_bgen_samples, read_layout2_variant_identifying_data, read_layout2_variant_metadata,
    skip_layout2_variant_identifying_data, BgenHeader,
};
use super::index::{indexed_region_records, validate_index_record_consumed, BgenIndexRecord};

const BGEN_READER_BUFFER_SIZE: usize = 1 << 20;

fn open_bgen_reader(
    bgen: &Path,
    #[cfg(test)] probe: Option<&BgenWorkProbe>,
) -> Result<BufReader<File>> {
    let file = File::open(bgen).map_err(|source| GenoioError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
    #[cfg(test)]
    if let Some(probe) = probe {
        probe.record_bgen_open();
    }
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
        #[cfg(test)]
        {
            Self::open_owned_inner(bgen, None)
        }
        #[cfg(not(test))]
        {
            Self::open_owned_inner(bgen)
        }
    }

    #[cfg(test)]
    fn open_owned_with_probe(bgen: PathBuf, probe: &BgenWorkProbe) -> Result<Self> {
        Self::open_owned_inner(bgen, Some(probe))
    }

    fn open_owned_inner(bgen: PathBuf, #[cfg(test)] probe: Option<&BgenWorkProbe>) -> Result<Self> {
        #[cfg(test)]
        let mut reader = open_bgen_reader(&bgen, probe)?;
        #[cfg(not(test))]
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
        #[cfg(test)]
        {
            read_bgen_samples(&mut self.reader, &self.bgen, sample, &self.header)
        }
        #[cfg(not(test))]
        {
            read_bgen_samples(&mut self.reader, &self.bgen, sample, &self.header)
        }
    }

    #[cfg(test)]
    fn read_samples_with_probe(
        &mut self,
        sample: Option<&Path>,
        probe: &BgenWorkProbe,
    ) -> Result<Vec<SampleRecord>> {
        super::header::read_bgen_samples_with_probe(
            &mut self.reader,
            &self.bgen,
            sample,
            &self.header,
            probe,
        )
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

/// Persistent sequential or indexed BGEN genotype or haplotype dosage state.
pub(crate) struct BgenBlockSession {
    pub(super) io: BgenReadSession,
    pub(super) selection: DenseSampleSelection,
    pub(super) diagnostics: DenseDiagnostics,
    pub(super) variant_filter: Option<VariantFilter>,
    pub(super) missing_policy: DenseMissingPolicy,
    pub(super) return_samples: bool,
    pub(super) return_variants: bool,
    pub(super) matrix_kind: MatrixKind,
    pub(super) cursor: BgenVariantCursor,
    pub(super) retained_skip: usize,
    pub(super) dosage_buffers: Option<DosageDecodeBuffers>,
    pub(super) haplotype_buffers: Option<HaplotypeDecodeBuffers>,
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
        #[cfg(test)]
        {
            Self::open_windowed_impl(bgen, sample, options, retained_skip, None)
        }
        #[cfg(not(test))]
        {
            Self::open_windowed_impl(bgen, sample, options, retained_skip)
        }
    }

    fn open_windowed_impl(
        bgen: PathBuf,
        sample: Option<PathBuf>,
        options: BlockReadOptions,
        retained_skip: usize,
        #[cfg(test)] probe: Option<BgenWorkProbe>,
    ) -> Result<Self> {
        validate_bgen_options(&options)?;
        #[cfg(test)]
        let mut io = if let Some(probe) = probe.as_ref() {
            BgenReadSession::open_owned_with_probe(bgen, probe)?
        } else {
            BgenReadSession::open_owned(bgen)?
        };
        #[cfg(not(test))]
        let mut io = BgenReadSession::open_owned(bgen)?;
        #[cfg(test)]
        let all_samples = if let Some(probe) = probe.as_ref() {
            io.read_samples_with_probe(sample.as_deref(), probe)?
        } else {
            io.read_samples(sample.as_deref())?
        };
        #[cfg(not(test))]
        let all_samples = io.read_samples(sample.as_deref())?;
        let selection = select_samples_source_order(
            &all_samples,
            options.requested_samples.as_deref(),
            &io.bgen,
        )?;
        #[cfg(test)]
        let indexed_records = if let Some(probe) = probe.as_ref() {
            super::index::indexed_region_records_with_probe(
                &io.bgen,
                options.variant_filter.as_ref(),
                probe,
            )?
        } else {
            indexed_region_records(&io.bgen, options.variant_filter.as_ref())?
        };
        #[cfg(not(test))]
        let indexed_records = indexed_region_records(&io.bgen, options.variant_filter.as_ref())?;
        let cursor = match indexed_records {
            Some(records) => BgenVariantCursor::indexed(records),
            None => {
                io.seek_to_variants()?;
                BgenVariantCursor::sequential(io.header.variant_count)
            }
        };
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
            matrix_kind: options.matrix_kind,
            cursor,
            retained_skip,
            dosage_buffers: (options.matrix_kind == MatrixKind::Genotype)
                .then(DosageDecodeBuffers::default),
            haplotype_buffers: (options.matrix_kind == MatrixKind::Haplotype)
                .then(HaplotypeDecodeBuffers::default),
            eof,
            #[cfg(test)]
            probe,
        })
    }

    #[cfg(test)]
    fn open_with_probe(
        bgen: PathBuf,
        sample: Option<PathBuf>,
        options: BlockReadOptions,
        probe: BgenWorkProbe,
    ) -> Result<Self> {
        Self::open_windowed_impl(bgen, sample, options, 0, Some(probe))
    }

    pub(crate) fn next_block(&mut self, block_size: usize) -> Result<Option<BlockOutput>> {
        match self.matrix_kind {
            MatrixKind::Genotype => self
                .next_dosage_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
            MatrixKind::Haplotype => self
                .next_haplotype_block(block_size)
                .map(|matrix| matrix.map(BlockOutput::Dense)),
        }
    }

    pub(super) fn source_record_capacity(&self) -> usize {
        self.cursor.source_record_capacity()
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

    pub(super) fn empty_haplotype_output(&self) -> Result<genoio_core::DenseGenotypeMatrix> {
        let samples = super::haplotype::expand_selected_samples_to_haplotypes(&self.selection);
        let samples =
            SampleMetadataBuffers::optional_from_records(&samples, self.return_samples, true)?;
        let variants = self
            .return_variants
            .then(|| VariantMetadataBuffers::with_capacity(0));
        genoio_core::DenseGenotypeMatrix::new_with_layout(
            self.selection.samples.len() * 2,
            0,
            Vec::new(),
            genoio_core::DenseLayout::VariantMajor,
            samples,
            variants,
            crate::blocks::block_diagnostics_snapshot(&self.diagnostics, 0),
        )
    }

    pub(super) fn next_position(&mut self) -> Result<Option<BgenRecordPosition>> {
        let position = self.cursor.next(&mut self.io)?;
        if position.is_some() {
            self.record_candidate_visit();
        } else {
            self.eof = true;
        }
        Ok(position)
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
    Ok(())
}

#[cfg(test)]
#[derive(Debug, Clone, Default)]
pub(super) struct BgenWorkProbe {
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

    pub(super) fn record_bgen_open(&self) {
        self.update(|counts| counts.bgen_opens += 1);
    }

    pub(super) fn record_sample_open(&self) {
        self.update(|counts| counts.sample_opens += 1);
    }

    pub(super) fn record_index_open(&self) {
        self.update(|counts| counts.index_opens += 1);
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
pub(super) enum BgenRecordPosition {
    Sequential,
    Indexed(BgenIndexRecord),
}

impl BgenRecordPosition {
    pub(super) fn validate_if_indexed(self, session: &mut BgenReadSession) -> Result<()> {
        if let Self::Indexed(index_record) = self {
            session.validate_index_record_consumed(&index_record)?;
        }
        Ok(())
    }
}

/// Static-dispatch cursor over either the full BGEN stream or `.bgi` byte ranges.
///
/// Indexed reads seek before yielding each record; sequential reads only count
/// down records already positioned after the header/sample block.
pub(super) enum BgenVariantCursor {
    Sequential {
        remaining: u32,
    },
    Indexed {
        records: Vec<BgenIndexRecord>,
        next_index: usize,
    },
}

impl BgenVariantCursor {
    pub(super) fn sequential(variant_count: u32) -> Self {
        Self::Sequential {
            remaining: variant_count,
        }
    }

    pub(super) fn indexed(records: Vec<BgenIndexRecord>) -> Self {
        Self::Indexed {
            records,
            next_index: 0,
        }
    }

    pub(super) fn next(
        &mut self,
        session: &mut BgenReadSession,
    ) -> Result<Option<BgenRecordPosition>> {
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
                let Some(&index_record) = records.get(*next_index) else {
                    return Ok(None);
                };
                *next_index += 1;
                session.seek_to_index_record(&index_record)?;
                Ok(Some(BgenRecordPosition::Indexed(index_record)))
            }
        }
    }

    fn source_record_capacity(&self) -> usize {
        match self {
            Self::Sequential { remaining } => usize::try_from(*remaining).unwrap_or(usize::MAX),
            Self::Indexed { records, .. } => records.len(),
        }
    }
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
    const FLAG_SAMPLE_IDENTIFIERS: u32 = 1 << 31;

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
        push_variant_with_phase(bytes, id, chrom, pos, dosage, 0);
    }

    fn push_variant_with_phase(
        bytes: &mut Vec<u8>,
        id: &str,
        chrom: &str,
        pos: u32,
        dosage: (u8, u8),
        phased: u8,
    ) {
        push_variant_with_phase_and_ploidy(bytes, id, chrom, pos, dosage, phased, 2);
    }

    fn push_variant_with_phase_and_ploidy(
        bytes: &mut Vec<u8>,
        id: &str,
        chrom: &str,
        pos: u32,
        dosage: (u8, u8),
        phased: u8,
        ploidy: u8,
    ) {
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
        probability.extend_from_slice(&[2, 2, ploidy, phased, 8, dosage.0, dosage.1]);
        bytes.extend_from_slice(
            &u32::try_from(probability.len())
                .expect("test probability length should fit u32")
                .to_le_bytes(),
        );
        bytes.extend_from_slice(&probability);
    }

    fn write_two_variant_bgen_without_embedded_samples(
        bgen: &Path,
        sample: &Path,
    ) -> Vec<(u64, u64)> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&2_u32.to_le_bytes());
        bytes.extend_from_slice(&1_u32.to_le_bytes());
        bytes.extend_from_slice(b"bgen");
        bytes.extend_from_slice(&FLAG_LAYOUT2.to_le_bytes());
        let first_start = u64::try_from(bytes.len()).expect("first offset should fit u64");
        push_variant(&mut bytes, "rs1", "1", 10, (255, 0));
        let first_end = u64::try_from(bytes.len()).expect("first end should fit u64");
        let second_start = u64::try_from(bytes.len()).expect("second offset should fit u64");
        push_variant(&mut bytes, "rs2", "2", 20, (0, 255));
        let second_end = u64::try_from(bytes.len()).expect("second end should fit u64");
        fs::write(bgen, bytes).expect("test bgen should be written");

        let mut sample_file = fs::File::create(sample).expect("test sample file should open");
        sample_file
            .write_all(b"ID_1 ID_2 missing\n0 0 0\nsample_1 sample_1 0\n")
            .expect("test sample file should be written");

        vec![
            (first_start, first_end - first_start),
            (second_start, second_end - second_start),
        ]
    }

    fn write_two_variant_phased_bgen_without_embedded_samples(
        bgen: &Path,
        sample: &Path,
    ) -> Vec<(u64, u64)> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&2_u32.to_le_bytes());
        bytes.extend_from_slice(&1_u32.to_le_bytes());
        bytes.extend_from_slice(b"bgen");
        bytes.extend_from_slice(&FLAG_LAYOUT2.to_le_bytes());
        let first_start = u64::try_from(bytes.len()).expect("first offset should fit u64");
        push_variant_with_phase(&mut bytes, "rs1", "1", 10, (255, 0), 1);
        let first_end = u64::try_from(bytes.len()).expect("first end should fit u64");
        let second_start = u64::try_from(bytes.len()).expect("second offset should fit u64");
        push_variant_with_phase(&mut bytes, "rs2", "2", 20, (0, 255), 1);
        let second_end = u64::try_from(bytes.len()).expect("second end should fit u64");
        fs::write(bgen, bytes).expect("test bgen should be written");

        let mut sample_file = fs::File::create(sample).expect("test sample file should open");
        sample_file
            .write_all(b"ID_1 ID_2 missing\n0 0 0\nsample_1 sample_1 0\n")
            .expect("test sample file should be written");

        vec![
            (first_start, first_end - first_start),
            (second_start, second_end - second_start),
        ]
    }

    fn write_two_variant_missing_bgen_without_embedded_samples(
        bgen: &Path,
        sample: &Path,
        phased: u8,
    ) -> Vec<(u64, u64)> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&2_u32.to_le_bytes());
        bytes.extend_from_slice(&1_u32.to_le_bytes());
        bytes.extend_from_slice(b"bgen");
        bytes.extend_from_slice(&FLAG_LAYOUT2.to_le_bytes());
        let first_start = u64::try_from(bytes.len()).expect("first offset should fit u64");
        push_variant_with_phase(&mut bytes, "rs1", "1", 10, (255, 0), phased);
        let first_end = u64::try_from(bytes.len()).expect("first end should fit u64");
        let second_start = u64::try_from(bytes.len()).expect("second offset should fit u64");
        push_variant_with_phase_and_ploidy(&mut bytes, "rs2", "2", 20, (0, 255), phased, 0x82);
        let second_end = u64::try_from(bytes.len()).expect("second end should fit u64");
        fs::write(bgen, bytes).expect("test bgen should be written");

        let mut sample_file = fs::File::create(sample).expect("test sample file should open");
        sample_file
            .write_all(b"ID_1 ID_2 missing\n0 0 0\nsample_1 sample_1 0\n")
            .expect("test sample file should be written");

        vec![
            (first_start, first_end - first_start),
            (second_start, second_end - second_start),
        ]
    }

    fn write_one_variant_bgen_with_embedded_samples(bgen: &Path) {
        let mut sample_block = Vec::new();
        let sample_id = b"embedded_1";
        let sample_block_len = u32::try_from(8 + 2 + sample_id.len())
            .expect("test sample block length should fit u32");
        sample_block.extend_from_slice(&sample_block_len.to_le_bytes());
        sample_block.extend_from_slice(&1_u32.to_le_bytes());
        sample_block.extend_from_slice(
            &u16::try_from(sample_id.len())
                .expect("test sample identifier length should fit u16")
                .to_le_bytes(),
        );
        sample_block.extend_from_slice(sample_id);

        let mut bytes = Vec::new();
        let offset =
            u32::try_from(20 + sample_block.len()).expect("test bgen offset should fit u32");
        bytes.extend_from_slice(&offset.to_le_bytes());
        bytes.extend_from_slice(&20_u32.to_le_bytes());
        bytes.extend_from_slice(&1_u32.to_le_bytes());
        bytes.extend_from_slice(&1_u32.to_le_bytes());
        bytes.extend_from_slice(b"bgen");
        bytes.extend_from_slice(&(FLAG_LAYOUT2 | FLAG_SAMPLE_IDENTIFIERS).to_le_bytes());
        bytes.extend_from_slice(&sample_block);
        push_variant(&mut bytes, "rs1", "1", 10, (255, 0));
        fs::write(bgen, bytes).expect("test bgen should be written");
    }

    fn bgen_index_path(bgen: &Path) -> std::path::PathBuf {
        let mut path = bgen.as_os_str().to_os_string();
        path.push(".bgi");
        path.into()
    }

    fn write_bgen_index(bgen: &Path, ranges: &[(u64, u64)], second_size_delta: u64) {
        let connection =
            rusqlite::Connection::open(bgen_index_path(bgen)).expect("test bgen index should open");
        connection
            .execute_batch(
                "CREATE TABLE Variant (
                    chromosome TEXT NOT NULL,
                    position INT NOT NULL,
                    file_start_position INT NOT NULL,
                    size_in_bytes INT NOT NULL
                );",
            )
            .expect("test bgen index schema should be created");
        for (index, &(start, size)) in ranges.iter().enumerate() {
            let chromosome = if index == 0 { "1" } else { "2" };
            let position = i64::try_from((index + 1) * 10).expect("test position should fit i64");
            let size = if index == 1 {
                size + second_size_delta
            } else {
                size
            };
            connection
                .execute(
                    "INSERT INTO Variant (
                        chromosome, position, file_start_position, size_in_bytes
                     ) VALUES (?1, ?2, ?3, ?4)",
                    rusqlite::params![
                        chromosome,
                        position,
                        i64::try_from(start).expect("test start should fit i64"),
                        i64::try_from(size).expect("test size should fit i64")
                    ],
                )
                .expect("test bgen index row should insert");
        }
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
        let _ranges = write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
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
    fn pbr_rust_bgen_001_probe_counts_source_open_before_later_header_error() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("truncated.bgen");
        fs::write(&bgen, b"bgen").expect("truncated bgen fixture should be written");
        let probe = BgenWorkProbe::default();

        assert!(
            BgenBlockSession::open_with_probe(bgen, None, options(None), probe.clone()).is_err(),
            "truncated bgen header should fail after the source opens"
        );

        let counts = probe.snapshot();
        assert_eq!(counts.bgen_opens, 1);
        assert_eq!(counts.sample_opens, 0);
        assert_eq!(counts.index_opens, 0);
    }

    #[test]
    fn pbr_rust_bgen_001_embedded_samples_do_not_open_supplied_companion() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("embedded.bgen");
        let unused_sample = dir.path().join("does-not-exist.sample");
        write_one_variant_bgen_with_embedded_samples(&bgen);
        let probe = BgenWorkProbe::default();

        let _session = BgenBlockSession::open_with_probe(
            bgen,
            Some(unused_sample),
            options(None),
            probe.clone(),
        )
        .expect("embedded samples should make the supplied companion unnecessary");

        let counts = probe.snapshot();
        assert_eq!(counts.bgen_opens, 1);
        assert_eq!(counts.sample_opens, 0);
    }

    #[test]
    fn pbr_rust_bgen_001_early_drop_stops_sequential_work() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let _ranges = write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
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
        let _ranges = write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
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

    fn region_filter(value: &str) -> VariantFilter {
        VariantFilter::from_json_value(json!({
            "op": "predicate",
            "name": "region",
            "params": {"value": value}
        }))
        .expect("test region filter should parse")
    }

    #[test]
    fn pbr_rust_bgen_002_indexed_probe_owns_records_and_stops_at_eof() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let ranges = write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
        write_bgen_index(&bgen, &ranges, 0);
        let probe = BgenWorkProbe::default();

        {
            let mut session = BgenBlockSession::open_with_probe(
                bgen.clone(),
                Some(sample),
                options(Some(region_filter("2:1-30"))),
                probe.clone(),
            )
            .expect("indexed bgen session should open");
            fs::remove_file(bgen_index_path(&bgen))
                .expect("session should not retain the sqlite index handle");
            assert!(session
                .next_block(1)
                .expect("indexed bgen block should decode from owned records")
                .is_some());
            assert!(session
                .next_block(1)
                .expect("indexed bgen session should reach EOF")
                .is_none());
            let at_eof = probe.snapshot();
            assert!(session
                .next_block(1)
                .expect("indexed bgen EOF should be sticky")
                .is_none());
            assert_eq!(probe.snapshot(), at_eof);
        }

        let counts = probe.snapshot();
        assert_eq!(counts.bgen_opens, 1);
        assert_eq!(counts.sample_opens, 1);
        assert_eq!(counts.index_opens, 1);
        assert_eq!(counts.candidate_visits, 1);
        assert_eq!(counts.payload_decodes, 1);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_bgen_002_probe_counts_sqlite_open_before_later_query_error() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let _ranges = write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
        let connection = rusqlite::Connection::open(bgen_index_path(&bgen))
            .expect("empty test index should open");
        drop(connection);
        let probe = BgenWorkProbe::default();

        assert!(
            BgenBlockSession::open_with_probe(
                bgen,
                Some(sample),
                options(Some(region_filter("2:1-30"))),
                probe.clone(),
            )
            .is_err(),
            "an index without the Variant table should fail after sqlite opens"
        );

        let counts = probe.snapshot();
        assert_eq!(counts.bgen_opens, 1);
        assert_eq!(counts.sample_opens, 1);
        assert_eq!(counts.index_opens, 1);
    }

    #[test]
    fn pbr_rust_bgen_002_indexed_cursor_validates_every_consumed_byte_range() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let ranges = write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
        write_bgen_index(&bgen, &ranges, 1);
        let mut session = BgenBlockSession::open_with_probe(
            bgen,
            Some(sample),
            options(Some(region_filter("2:1-30"))),
            BgenWorkProbe::default(),
        )
        .expect("indexed bgen session should open");

        let error = session
            .next_block(1)
            .expect_err("mismatched indexed byte range should fail");

        assert!(
            error
                .to_string()
                .contains("byte range does not match decoded variant record"),
            "unexpected indexed byte-range error: {error}"
        );
    }

    #[test]
    fn pbr_rust_bgen_002_indexed_genotype_range_error_precedes_missing_policy_error() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let ranges = write_two_variant_missing_bgen_without_embedded_samples(&bgen, &sample, 0);
        write_bgen_index(&bgen, &ranges, 1);
        let mut read_options = options(Some(region_filter("2:1-30")));
        read_options.missing_policy = DenseMissingPolicy::Raise;
        let mut session = BgenBlockSession::open(bgen, Some(sample), read_options)
            .expect("indexed genotype session should open");

        let error = session
            .next_block(1)
            .expect_err("the consumed malformed index range must be validated first");

        assert!(
            error
                .to_string()
                .contains("byte range does not match decoded variant record"),
            "unexpected indexed byte-range error: {error}"
        );
    }

    #[test]
    fn pbr_rust_bgen_002_indexed_haplotype_range_error_precedes_missing_policy_error() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let ranges = write_two_variant_missing_bgen_without_embedded_samples(&bgen, &sample, 1);
        write_bgen_index(&bgen, &ranges, 1);
        let mut read_options = options(Some(region_filter("2:1-30")));
        read_options.matrix_kind = MatrixKind::Haplotype;
        read_options.missing_policy = DenseMissingPolicy::Raise;
        let mut session = BgenBlockSession::open(bgen, Some(sample), read_options)
            .expect("indexed haplotype session should open");

        let error = session
            .next_block(1)
            .expect_err("the consumed malformed index range must be validated first");

        assert!(
            error
                .to_string()
                .contains("byte range does not match decoded variant record"),
            "unexpected indexed byte-range error: {error}"
        );
    }

    #[test]
    fn pbr_rust_bgen_002_indexed_early_drop_stops_work() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let ranges = write_two_variant_bgen_without_embedded_samples(&bgen, &sample);
        write_bgen_index(&bgen, &ranges, 0);
        let probe = BgenWorkProbe::default();

        {
            let mut session = BgenBlockSession::open_with_probe(
                bgen,
                Some(sample),
                options(Some(region_filter("1:1-30"))),
                probe.clone(),
            )
            .expect("indexed bgen session should open");
            assert!(session
                .next_block(1)
                .expect("first indexed bgen block should decode")
                .is_some());
        }

        let counts = probe.snapshot();
        assert_eq!(counts.index_opens, 1);
        assert_eq!(counts.candidate_visits, 1);
        assert_eq!(counts.payload_decodes, 1);
        assert_eq!(counts.drops, 1);
    }

    #[test]
    fn pbr_rust_bgen_002_pbr_rust_alloc_001_haplotype_session_reuses_bounded_scratch() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let _ranges = write_two_variant_phased_bgen_without_embedded_samples(&bgen, &sample);
        let probe = BgenWorkProbe::default();
        let mut haplotype_options = options(None);
        haplotype_options.matrix_kind = MatrixKind::Haplotype;
        let mut session =
            BgenBlockSession::open_with_probe(bgen, Some(sample), haplotype_options, probe.clone())
                .expect("persistent bgen haplotype session should open");

        while session
            .next_block(1)
            .expect("persistent bgen haplotype block should decode")
            .is_some()
        {}

        let counts = probe.snapshot();
        assert_eq!(counts.candidate_visits, 2);
        assert_eq!(counts.payload_decodes, 2);
        assert_eq!(counts.max_dense_output_len, 2);
    }

    #[test]
    fn pbr_rust_bgen_002_indexed_haplotype_probe_counts_one_index_and_selected_range() {
        let dir = tempfile::tempdir().expect("test directory should be created");
        let bgen = dir.path().join("tiny.bgen");
        let sample = dir.path().join("tiny.sample");
        let ranges = write_two_variant_phased_bgen_without_embedded_samples(&bgen, &sample);
        write_bgen_index(&bgen, &ranges, 0);
        let probe = BgenWorkProbe::default();
        let mut haplotype_options = options(Some(region_filter("2:1-30")));
        haplotype_options.matrix_kind = MatrixKind::Haplotype;
        let mut session =
            BgenBlockSession::open_with_probe(bgen, Some(sample), haplotype_options, probe.clone())
                .expect("persistent indexed bgen haplotype session should open");

        assert!(session
            .next_block(1)
            .expect("indexed bgen haplotype block should decode")
            .is_some());
        assert!(session
            .next_block(1)
            .expect("indexed bgen haplotype session should reach EOF")
            .is_none());

        let counts = probe.snapshot();
        assert_eq!(counts.index_opens, 1);
        assert_eq!(counts.candidate_visits, 1);
        assert_eq!(counts.payload_decodes, 1);
    }
}
