// pattern: Imperative Shell
//! PGEN header dispatch and variant decode entry points.
//!
//! This module owns format-mode routing and decoder state. Submodules handle
//! header parsing, record I/O, main-track hard calls, dosage overlays, haplotype
//! auxiliary tracks, and small bit-packing primitives.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;
use crate::hardcall::PackedHardcalls as PackedGenotypes;

mod bitpack;
mod dosage_track;
mod haplotype_track;
mod header;
mod io;
mod main_track;

use self::dosage_track::{
    overlay_fixed_width_dosages, overlay_variable_width_dosages, DosageOverlayTarget,
};
pub(super) use self::haplotype_track::{
    decode_plink2_haplotype_dosage_aux, decode_plink2_haplotype_hardcall_aux,
    read_plink2_variant_haplotype_dosage_track, read_plink2_variant_haplotype_main_track,
};
use self::header::fixed_width_dosage_record_len;
pub(super) use self::header::{
    open_pgen_payload, read_supported_pgen_header, read_supported_pgen_header_from_file,
    read_supported_pgen_header_prefix, validate_plink2_dimensions, validate_plink2_sample_count,
};
use self::io::read_fixed_width_phased_dosage_variant_record;
pub(super) use self::io::{
    read_fixed_width_variant_packed_sequential, read_plink2_variant_packed,
    seek_fixed_width_variant_record,
};
use self::main_track::decode_variable_width_main_track;

const PGEN_MAGIC: [u8; 2] = [0x6c, 0x1b];
const PGEN_MODE_FIXED_WIDTH_HARDCALLS: u8 = 0x02;
const PGEN_MODE_FIXED_WIDTH_DOSAGE: u8 = 0x03;
const PGEN_MODE_FIXED_WIDTH_PHASED_DOSAGE: u8 = 0x04;
const PGEN_MODE_VARIABLE_WIDTH: u8 = 0x10;
const PGEN_HEADER_LEN: u64 = 12;
const PGEN_VARIANT_BLOCK_SIZE: usize = 65_536;
const PGEN_DOSAGE_SCALE: f32 = 2.0 / 32768.0;
const PGEN_PHASE_SCALE: f32 = 1.0 / 16384.0;
const PGEN_MAX_DOSAGE_RAW: u16 = 32_768;
const PGEN_MIN_PHASE_RAW: i16 = -16_384;
const PGEN_MAX_PHASE_RAW: i16 = 16_384;

/// Parsed PGEN header fields needed by all read paths.
///
/// Fixed-width layouts compute offsets from shape. Variable-width layouts carry
/// record types and block offsets so callers can seek directly to records.
#[derive(Debug, Clone)]
pub(super) struct PgenHeader {
    pub(super) layout: PgenLayout,
    pub(super) variant_ct: usize,
    pub(super) sample_ct: usize,
    bytes_per_variant: usize,
    pub(super) record_types: Vec<u8>,
    record_offsets: Vec<u64>,
}

/// Reused PGEN decode scratch for one read loop.
///
/// The state keeps variable-width record bytes, packed hard-call output, dense
/// selected values, sparse missing positions, and LD-compressed main-track
/// history. Missing indices are relative to the selected output vector.
#[derive(Debug, Clone)]
pub(super) struct PgenDecoderState {
    previous_non_ld_packed: PackedGenotypes,
    has_previous_non_ld: bool,
    record: Vec<u8>,
    pub(super) packed: PackedGenotypes,
    pub(super) values: Vec<f32>,
    pub(super) missing_indices: Vec<usize>,
}

/// Reused selected-output buffers for PLINK2 haplotype reads.
///
/// Collapsed diploid buffers are populated alongside haplotype rows when
/// genotype-stat filters need dosage semantics. Missing indices are sparse
/// positions into their corresponding selected-value buffers.
#[derive(Default)]
pub(super) struct PgenHaplotypeDecodeState {
    pub(super) selected_haplotype_values: Vec<f32>,
    pub(super) selected_haplotype_missing_indices: Vec<usize>,
    pub(super) selected_collapsed_values: Vec<f32>,
    pub(super) selected_collapsed_missing_indices: Vec<usize>,
}

struct SelectedSampleCursor<'a> {
    source_indices: &'a [usize],
    selected_index: usize,
}

impl<'a> SelectedSampleCursor<'a> {
    fn new(source_indices: &'a [usize]) -> Self {
        Self {
            source_indices,
            selected_index: 0,
        }
    }

    fn selected_index_for(&mut self, source_index: usize) -> Option<usize> {
        // source_indices are stored in PGEN source order, so a forward-only
        // cursor avoids a search for every stored dosage.
        while self
            .source_indices
            .get(self.selected_index)
            .is_some_and(|selected_source_index| *selected_source_index < source_index)
        {
            self.selected_index += 1;
        }
        self.source_indices
            .get(self.selected_index)
            .copied()
            .filter(|selected_source_index| *selected_source_index == source_index)
            .map(|_| self.selected_index)
    }
}

fn insert_sorted_unique_index(indices: &mut Vec<usize>, index: usize) {
    // PGEN dosage and phase tracks can override hard-call missingness after the
    // base genotype expansion, so keep the sparse index list sorted in place.
    match indices.binary_search(&index) {
        Ok(_) => {}
        Err(position) => indices.insert(position, index),
    }
}

fn remove_sorted_index(indices: &mut Vec<usize>, index: usize) {
    if let Ok(position) = indices.binary_search(&index) {
        indices.remove(position);
    }
}

impl PgenDecoderState {
    pub(super) fn new(sample_ct: usize, selected_sample_ct: usize) -> Self {
        Self {
            previous_non_ld_packed: PackedGenotypes::default(),
            has_previous_non_ld: false,
            record: Vec::with_capacity(sample_ct.div_ceil(4)),
            packed: PackedGenotypes::default(),
            values: Vec::with_capacity(selected_sample_ct),
            missing_indices: Vec::new(),
        }
    }
}

/// Supported PGEN storage layouts after header validation.
#[derive(Debug, Clone)]
pub(super) enum PgenLayout {
    FixedWidth,
    FixedWidthDosage,
    FixedWidthPhasedDosage,
    VariableWidth,
}

pub(super) fn read_plink2_variant_values(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    read_plink2_variant_packed(path, file, header, variant_index, decoder_state)?;
    decoder_state.packed.expand_selected(
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing_indices,
    );
    Ok(())
}

pub(super) fn read_plink2_variant_dosage(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    match header.layout {
        PgenLayout::FixedWidthDosage => read_fixed_width_dosage_variant_values(
            path,
            file,
            header,
            variant_index,
            source_indices,
            decoder_state,
        ),
        PgenLayout::FixedWidthPhasedDosage => read_fixed_width_phased_dosage_variant_values(
            path,
            file,
            header,
            variant_index,
            source_indices,
            decoder_state,
        ),
        PgenLayout::VariableWidth => read_variable_width_dosage_variant_values(
            path,
            file,
            header,
            variant_index,
            source_indices,
            decoder_state,
        ),
        PgenLayout::FixedWidth => Err(GenoioError::unsupported(
            "pgen does not contain dosage values",
        )),
    }
}

fn read_fixed_width_dosage_variant_values(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    let record_len = fixed_width_dosage_record_len(header.sample_ct);
    let payload_offset = variant_index
        .checked_mul(record_len)
        .ok_or_else(|| GenoioError::invalid_source(path, "pgen variant offset is out of range"))?;
    let offset = PGEN_HEADER_LEN
        .checked_add(u64::try_from(payload_offset).map_err(|_| {
            GenoioError::invalid_source(path, "pgen variant offset is out of range")
        })?)
        .ok_or_else(|| GenoioError::invalid_source(path, "pgen variant offset is out of range"))?;
    file.seek(SeekFrom::Start(offset))
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    decoder_state.record.resize(record_len, 0);
    file.read_exact(&mut decoder_state.record)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    decoder_state.packed.load_pgen_payload(
        &decoder_state.record[..header.bytes_per_variant],
        header.sample_ct,
    );
    decoder_state.packed.expand_selected(
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing_indices,
    );
    overlay_fixed_width_dosages(
        path,
        &decoder_state.record[header.bytes_per_variant..],
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing_indices,
    )
}

fn read_fixed_width_phased_dosage_variant_values(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    let cursor = read_fixed_width_phased_dosage_variant_record(
        path,
        file,
        header,
        variant_index,
        decoder_state,
    )?;
    decoder_state.packed.expand_selected(
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing_indices,
    );
    let dosage_end = cursor
        .checked_add(header.sample_ct.checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source(path, "pgen dosage byte count is out of range")
        })?)
        .ok_or_else(|| {
            GenoioError::invalid_source(path, "pgen dosage byte count is out of range")
        })?;
    overlay_fixed_width_dosages(
        path,
        &decoder_state.record[cursor..dosage_end],
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing_indices,
    )
}

fn read_variable_width_dosage_variant_values(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    let start = header.record_offsets[variant_index];
    let end = header.record_offsets[variant_index + 1];
    let record_len =
        usize::try_from(end.checked_sub(start).ok_or_else(|| {
            GenoioError::invalid_source(path, "pgen record length is out of range")
        })?)
        .map_err(|_| GenoioError::invalid_source(path, "pgen record length is out of range"))?;
    file.seek(SeekFrom::Start(start))
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    decoder_state.record.resize(record_len, 0);
    file.read_exact(&mut decoder_state.record)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    let record = decoder_state.record.as_slice();
    let record_type = header.record_types[variant_index];
    let dosage_bits = (record_type >> 5) & 0x03;
    if dosage_bits == 0 {
        return Err(GenoioError::unsupported(
            "pgen record does not contain dosage values",
        ));
    }
    let cursor = decode_variable_width_main_track(
        path,
        record,
        record_type,
        header.sample_ct,
        &decoder_state.previous_non_ld_packed,
        decoder_state.has_previous_non_ld,
        &mut decoder_state.packed,
    )?;
    if !matches!(record_type & 0x07, 2 | 3) {
        decoder_state
            .previous_non_ld_packed
            .copy_from(&decoder_state.packed);
        decoder_state.has_previous_non_ld = true;
    }

    decoder_state.packed.expand_selected(
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing_indices,
    );
    overlay_variable_width_dosages(
        path,
        record,
        cursor,
        dosage_bits,
        header.sample_ct,
        DosageOverlayTarget {
            source_indices,
            values: &mut decoder_state.values,
            missing_indices: &mut decoder_state.missing_indices,
        },
    )
}
