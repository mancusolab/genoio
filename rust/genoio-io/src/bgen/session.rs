// pattern: Imperative Shell
//! BGEN reader session and variant cursor helpers.
//!
//! A session owns the buffered file, parsed header, compression mode, and
//! reusable probability buffers for one read call. Cursor types hide the
//! difference between sequential scans and indexed byte-range reads.

use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::Path;

use genoio_core::{
    DenseDiagnostics, DenseMissingPolicy, DenseSampleSelection, GenoioError, SampleRecord,
    VariantFilter, VariantMetadataArrowBuffers, VariantRecord, VariantWindow,
};

use crate::error::Result;

use super::decode::{
    read_layout2_probability_payload_into, skip_layout2_probability_payload_raw,
    ProbabilityPayloadBuffers,
};
use super::header::{
    read_bgen_samples, read_layout2_variant_identifying_data, read_layout2_variant_metadata_arrow,
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
pub(super) struct BgenReadSession<'a> {
    pub(super) reader: BufReader<File>,
    pub(super) bgen: &'a Path,
    pub(super) header: BgenHeader,
}

impl<'a> BgenReadSession<'a> {
    pub(super) fn open(bgen: &'a Path) -> Result<Self> {
        let mut reader = open_bgen_reader(bgen)?;
        let header = BgenHeader::read_from(&mut reader, bgen)?;
        header.validate(bgen)?;
        Ok(Self {
            reader,
            bgen,
            header,
        })
    }

    pub(super) fn read_samples(&mut self, sample: Option<&Path>) -> Result<Vec<SampleRecord>> {
        read_bgen_samples(&mut self.reader, self.bgen, sample, &self.header)
    }

    pub(super) fn read_all_variant_metadata_arrow(
        &mut self,
    ) -> Result<VariantMetadataArrowBuffers> {
        read_layout2_variant_metadata_arrow(
            &mut self.reader,
            self.bgen,
            self.header.variant_count,
            self.header.flags.compression,
        )
    }

    pub(super) fn seek_to_variants(&mut self) -> Result<()> {
        self.reader
            .seek(SeekFrom::Start(u64::from(self.header.offset) + 4))
            .map_err(|source| GenoioError::Io {
                path: self.bgen.to_path_buf(),
                source,
            })?;
        Ok(())
    }

    fn seek_to_index_record(&mut self, index_record: &BgenIndexRecord) -> Result<()> {
        self.reader
            .seek(SeekFrom::Start(index_record.file_start_position))
            .map_err(|source| GenoioError::Io {
                path: self.bgen.to_path_buf(),
                source,
            })?;
        Ok(())
    }

    pub(super) fn read_variant(&mut self) -> Result<VariantRecord> {
        read_layout2_variant_identifying_data(&mut self.reader, self.bgen)
    }

    pub(super) fn skip_variant(&mut self) -> Result<()> {
        skip_layout2_variant_identifying_data(&mut self.reader, self.bgen)
    }

    pub(super) fn read_payload_into(
        &mut self,
        buffers: &mut ProbabilityPayloadBuffers,
    ) -> Result<()> {
        read_layout2_probability_payload_into(
            &mut self.reader,
            self.bgen,
            self.header.flags.compression,
            buffers,
        )
    }

    pub(super) fn skip_payload(&mut self) -> Result<()> {
        skip_layout2_probability_payload_raw(
            &mut self.reader,
            self.bgen,
            self.header.flags.compression,
        )
    }

    fn validate_index_record_consumed(&mut self, index_record: &BgenIndexRecord) -> Result<()> {
        validate_index_record_consumed(&mut self.reader, self.bgen, index_record)
    }
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
    pub(super) fn validate_if_indexed(self, session: &mut BgenReadSession<'_>) -> Result<()> {
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
        session: &mut BgenReadSession<'_>,
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
    pub(super) session: &'a mut BgenReadSession<'a>,
    pub(super) selection: DenseSampleSelection,
    pub(super) diagnostics: DenseDiagnostics,
    pub(super) variant_filter: Option<&'a VariantFilter>,
    pub(super) variant_window: Option<VariantWindow>,
    pub(super) missing_policy: DenseMissingPolicy,
    pub(super) return_samples: bool,
    pub(super) return_variants: bool,
}
