// pattern: Imperative Shell

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;
use crate::hardcall::PackedHardcalls as PackedGenotypes;

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

#[derive(Debug, Clone)]
pub(super) struct PgenHeader {
    pub(super) layout: PgenLayout,
    pub(super) variant_ct: usize,
    pub(super) sample_ct: usize,
    bytes_per_variant: usize,
    pub(super) record_types: Vec<u8>,
    record_offsets: Vec<u64>,
}

#[derive(Debug, Clone)]
pub(super) struct PgenDecoderState {
    previous_non_ld_packed: PackedGenotypes,
    has_previous_non_ld: bool,
    record: Vec<u8>,
    pub(super) packed: PackedGenotypes,
    pub(super) values: Vec<f32>,
    pub(super) missing: Vec<bool>,
}

struct DosageOverlayTarget<'a> {
    source_indices: &'a [usize],
    values: &'a mut [f32],
    missing: &'a mut [bool],
}

#[derive(Default)]
pub(super) struct PgenHaplotypeDecodeState {
    pub(super) selected_haplotype_values: Vec<f32>,
    pub(super) selected_haplotype_missing: Vec<bool>,
    pub(super) selected_collapsed_values: Vec<f32>,
    pub(super) selected_collapsed_missing: Vec<bool>,
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

impl PgenDecoderState {
    pub(super) fn new(sample_ct: usize, selected_sample_ct: usize) -> Self {
        Self {
            previous_non_ld_packed: PackedGenotypes::default(),
            has_previous_non_ld: false,
            record: Vec::with_capacity(sample_ct.div_ceil(4)),
            packed: PackedGenotypes::default(),
            values: Vec::with_capacity(selected_sample_ct),
            missing: Vec::with_capacity(selected_sample_ct),
        }
    }
}

#[derive(Debug, Clone)]
pub(super) enum PgenLayout {
    FixedWidth,
    FixedWidthDosage,
    FixedWidthPhasedDosage,
    VariableWidth,
}

pub(super) fn read_supported_pgen_header(path: &Path) -> Result<PgenHeader> {
    let mut file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = [0_u8; PGEN_HEADER_LEN as usize];
    file.read_exact(&mut header)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    if header[0..2] != PGEN_MAGIC {
        return Err(GenoioError::invalid_source(
            path,
            "invalid pgen magic bytes",
        ));
    }
    let (variant_ct, sample_ct) = parse_pgen_header_counts(path, &header)?;
    let bytes_per_variant = sample_ct.div_ceil(4);
    match header[2] {
        PGEN_MODE_FIXED_WIDTH_HARDCALLS => {
            if header[11] != 0 {
                return Err(GenoioError::invalid_source(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic hardcalls without header extensions are supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(path, &file, variant_ct, bytes_per_variant)?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidth,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_FIXED_WIDTH_DOSAGE => {
            if header[11] != 0 {
                return Err(GenoioError::invalid_source(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic dosage without header extensions is supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(
                path,
                &file,
                variant_ct,
                fixed_width_dosage_record_len(sample_ct),
            )?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidthDosage,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_FIXED_WIDTH_PHASED_DOSAGE => {
            if header[11] != 0 {
                return Err(GenoioError::invalid_source(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic phased dosage without header extensions is supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(
                path,
                &file,
                variant_ct,
                fixed_width_phased_dosage_record_len(sample_ct),
            )?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidthPhasedDosage,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_VARIABLE_WIDTH => {
            let (record_types, record_offsets) =
                read_variable_width_header_body(path, &mut file, variant_ct, header[11])?;
            Ok(PgenHeader {
                layout: PgenLayout::VariableWidth,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types,
                record_offsets,
            })
        }
        mode => Err(GenoioError::invalid_source(
            path,
            format!(
                "unsupported pgen mode 0x{mode:02x}; only fixed-width and variable-width biallelic hardcalls are supported"
            ),
        )),
    }
}

pub(super) fn read_supported_pgen_header_prefix(
    path: &Path,
    requested_variant_ct: usize,
) -> Result<PgenHeader> {
    let mut file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = [0_u8; PGEN_HEADER_LEN as usize];
    file.read_exact(&mut header)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    if header[0..2] != PGEN_MAGIC {
        return Err(GenoioError::invalid_source(
            path,
            "invalid pgen magic bytes",
        ));
    }
    let (variant_ct, sample_ct) = parse_pgen_header_counts(path, &header)?;
    let bytes_per_variant = sample_ct.div_ceil(4);
    let prefix_variant_ct = requested_variant_ct.min(variant_ct);
    match header[2] {
        PGEN_MODE_FIXED_WIDTH_HARDCALLS => {
            if header[11] != 0 {
                return Err(GenoioError::invalid_source(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic hardcalls without header extensions are supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(path, &file, variant_ct, bytes_per_variant)?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidth,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_FIXED_WIDTH_DOSAGE => {
            if header[11] != 0 {
                return Err(GenoioError::invalid_source(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic dosage without header extensions is supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(
                path,
                &file,
                variant_ct,
                fixed_width_dosage_record_len(sample_ct),
            )?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidthDosage,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_FIXED_WIDTH_PHASED_DOSAGE => {
            if header[11] != 0 {
                return Err(GenoioError::invalid_source(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic phased dosage without header extensions is supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(
                path,
                &file,
                variant_ct,
                fixed_width_phased_dosage_record_len(sample_ct),
            )?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidthPhasedDosage,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_VARIABLE_WIDTH => {
            let (record_types, record_offsets) = read_variable_width_header_body_prefix(
                path,
                &mut file,
                variant_ct,
                header[11],
                prefix_variant_ct,
            )?;
            Ok(PgenHeader {
                layout: PgenLayout::VariableWidth,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types,
                record_offsets,
            })
        }
        mode => Err(GenoioError::invalid_source(
            path,
            format!(
                "unsupported pgen mode 0x{mode:02x}; only fixed-width and variable-width biallelic hardcalls or unphased dosages are supported"
            ),
        )),
    }
}

fn parse_pgen_header_counts(
    path: &Path,
    header: &[u8; PGEN_HEADER_LEN as usize],
) -> Result<(usize, usize)> {
    let variant_ct = usize::try_from(u32::from_le_bytes([
        header[3], header[4], header[5], header[6],
    ]))
    .map_err(|_| GenoioError::invalid_source(path, "pgen variant count is out of range"))?;
    let sample_ct = usize::try_from(u32::from_le_bytes([
        header[7], header[8], header[9], header[10],
    ]))
    .map_err(|_| GenoioError::invalid_source(path, "pgen sample count is out of range"))?;
    Ok((variant_ct, sample_ct))
}

fn read_variable_width_header_body(
    path: &Path,
    file: &mut File,
    variant_ct: usize,
    header_format: u8,
) -> Result<(Vec<u8>, Vec<u64>)> {
    let type_length_format = header_format & 0x0f;
    let type_width_bits = match type_length_format {
        0..=3 => 4,
        4..=7 => 8,
        other => {
            return Err(GenoioError::invalid_source(
                path,
                format!("unsupported pgen variant-record type/length format {other}"),
            ));
        }
    };
    let length_width = usize::from((type_length_format & 0x03) + 1);
    let allele_count_format = (header_format >> 4) & 0x03;
    if allele_count_format != 0 {
        return Err(GenoioError::invalid_source(
            path,
            "unsupported pgen allele-count table; multiallelic PGEN decode is not implemented",
        ));
    }

    let block_ct = variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let mut block_offsets = Vec::with_capacity(block_ct);
    for _ in 0..block_ct {
        let mut bytes = [0_u8; 8];
        file.read_exact(&mut bytes)
            .map_err(|source| GenoioError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        block_offsets.push(u64::from_le_bytes(bytes));
    }

    let mut record_types = Vec::with_capacity(variant_ct);
    let mut record_lengths = Vec::with_capacity(variant_ct);
    for block_index in 0..block_ct {
        let block_variant_ct = block_variant_count(variant_ct, block_index);
        if type_width_bits == 8 {
            let mut types = vec![0_u8; block_variant_ct];
            file.read_exact(&mut types)
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            record_types.extend(types);
        } else {
            let mut packed_types = vec![0_u8; block_variant_ct.div_ceil(2)];
            file.read_exact(&mut packed_types)
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            for variant_in_block in 0..block_variant_ct {
                let byte = packed_types[variant_in_block / 2];
                let record_type = if variant_in_block % 2 == 0 {
                    byte & 0x0f
                } else {
                    byte >> 4
                };
                record_types.push(record_type);
            }
        }

        for _ in 0..block_variant_ct {
            let mut bytes = [0_u8; 4];
            file.read_exact(&mut bytes[..length_width])
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            record_lengths.push(u32::from_le_bytes(bytes));
        }
    }
    for record_type in &record_types {
        validate_supported_variable_record_type(path, *record_type)?;
    }
    let header_end = file.stream_position().map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    let mut record_offsets = Vec::with_capacity(variant_ct + 1);
    for (block_index, block_offset) in block_offsets.iter().enumerate() {
        let block_start = block_index * PGEN_VARIANT_BLOCK_SIZE;
        let block_end =
            (block_start + block_variant_count(variant_ct, block_index)).min(variant_ct);
        let mut offset = *block_offset;
        if record_offsets.len() == block_start {
            if block_index == 0 && offset != header_end {
                return Err(GenoioError::invalid_source(
                    path,
                    "pgen first variant-block offset does not match header length",
                ));
            }
            record_offsets.push(offset);
        } else if record_offsets
            .get(block_start)
            .is_none_or(|expected_offset| *expected_offset != offset)
        {
            return Err(GenoioError::invalid_source(
                path,
                "pgen variant-block offset does not match preceding record lengths",
            ));
        }
        for length in &record_lengths[block_start..block_end] {
            offset = offset.checked_add(u64::from(*length)).ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen record offset is out of range")
            })?;
            record_offsets.push(offset);
        }
    }
    if record_offsets.len() != variant_ct + 1 {
        return Err(GenoioError::invalid_source(
            path,
            "pgen variable-width header did not yield one offset per variant",
        ));
    }
    let actual_len = file
        .metadata()
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if record_offsets[variant_ct] > actual_len {
        return Err(GenoioError::invalid_source(
            path,
            "pgen variable-width records extend past end of file",
        ));
    }

    Ok((record_types, record_offsets))
}

fn read_variable_width_header_body_prefix(
    path: &Path,
    file: &mut File,
    variant_ct: usize,
    header_format: u8,
    prefix_variant_ct: usize,
) -> Result<(Vec<u8>, Vec<u64>)> {
    let type_length_format = header_format & 0x0f;
    let type_width_bits = match type_length_format {
        0..=3 => 4,
        4..=7 => 8,
        other => {
            return Err(GenoioError::invalid_source(
                path,
                format!("unsupported pgen variant-record type/length format {other}"),
            ));
        }
    };
    let length_width = usize::from((type_length_format & 0x03) + 1);
    let allele_count_format = (header_format >> 4) & 0x03;
    if allele_count_format != 0 {
        return Err(GenoioError::invalid_source(
            path,
            "unsupported pgen allele-count table; multiallelic PGEN decode is not implemented",
        ));
    }

    let block_ct = variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let mut block_offsets = Vec::with_capacity(block_ct);
    for _ in 0..block_ct {
        let mut bytes = [0_u8; 8];
        file.read_exact(&mut bytes)
            .map_err(|source| GenoioError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        block_offsets.push(u64::from_le_bytes(bytes));
    }

    let header_end = variable_width_header_end(path, variant_ct, type_width_bits, length_width)?;
    let prefix_block_ct = prefix_variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let mut record_types = Vec::with_capacity(prefix_variant_ct);
    let mut record_offsets = Vec::with_capacity(prefix_variant_ct.saturating_add(1));
    for (block_index, block_offset) in block_offsets
        .iter()
        .take(prefix_block_ct)
        .copied()
        .enumerate()
    {
        let block_variant_ct = block_variant_count(variant_ct, block_index);
        let block_start = block_index * PGEN_VARIANT_BLOCK_SIZE;
        let needed_in_block = prefix_variant_ct
            .saturating_sub(block_start)
            .min(block_variant_ct);
        // Type and length tables are block-grouped in the PGEN header. We
        // still have to skip through unneeded entries in the last touched
        // block so the file cursor reaches the matching length table.
        read_variable_record_type_prefix(
            path,
            file,
            type_width_bits,
            block_variant_ct,
            needed_in_block,
            &mut record_types,
        )?;
        if record_offsets.is_empty() {
            if block_index == 0 && block_offset != header_end {
                return Err(GenoioError::invalid_source(
                    path,
                    "pgen first variant-block offset does not match header length",
                ));
            }
            record_offsets.push(block_offset);
        } else if record_offsets
            .get(block_start)
            .is_none_or(|expected_offset| *expected_offset != block_offset)
        {
            return Err(GenoioError::invalid_source(
                path,
                "pgen variant-block offset does not match preceding record lengths",
            ));
        }
        let mut offset = block_offset;
        for _ in 0..needed_in_block {
            let mut bytes = [0_u8; 4];
            file.read_exact(&mut bytes[..length_width])
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            offset = offset
                .checked_add(u64::from(u32::from_le_bytes(bytes)))
                .ok_or_else(|| {
                    GenoioError::invalid_source(path, "pgen record offset is out of range")
                })?;
            record_offsets.push(offset);
        }
        let remaining_lengths = block_variant_ct - needed_in_block;
        skip_bytes(path, file, remaining_lengths * length_width)?;
    }
    // Only validate the prefix that may be decoded for this block. Unsupported
    // later records should not prevent first-block reads from succeeding.
    for record_type in &record_types {
        validate_supported_variable_record_type(path, *record_type)?;
    }
    if let Some(prefix_end) = record_offsets.last() {
        let actual_len = file
            .metadata()
            .map_err(|source| GenoioError::Io {
                path: path.to_path_buf(),
                source,
            })?
            .len();
        if *prefix_end > actual_len {
            return Err(GenoioError::invalid_source(
                path,
                "pgen variable-width records extend past end of file",
            ));
        }
    }
    Ok((record_types, record_offsets))
}

fn variable_width_header_end(
    path: &Path,
    variant_ct: usize,
    type_width_bits: usize,
    length_width: usize,
) -> Result<u64> {
    let block_ct = variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let block_offsets_len = block_ct.checked_mul(8).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen variable-width header is out of range")
    })?;
    let mut header_end = PGEN_HEADER_LEN
        .checked_add(u64::try_from(block_offsets_len).map_err(|_| {
            GenoioError::invalid_source(path, "pgen variable-width header is out of range")
        })?)
        .ok_or_else(|| {
            GenoioError::invalid_source(path, "pgen variable-width header is out of range")
        })?;
    for block_index in 0..block_ct {
        let block_variant_ct = block_variant_count(variant_ct, block_index);
        let type_table_len = if type_width_bits == 8 {
            block_variant_ct
        } else {
            block_variant_ct.div_ceil(2)
        };
        let length_table_len = block_variant_ct.checked_mul(length_width).ok_or_else(|| {
            GenoioError::invalid_source(path, "pgen variable-width header is out of range")
        })?;
        let table_len = type_table_len
            .checked_add(length_table_len)
            .ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen variable-width header is out of range")
            })?;
        header_end = header_end
            .checked_add(u64::try_from(table_len).map_err(|_| {
                GenoioError::invalid_source(path, "pgen variable-width header is out of range")
            })?)
            .ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen variable-width header is out of range")
            })?;
    }
    Ok(header_end)
}

fn read_variable_record_type_prefix(
    path: &Path,
    file: &mut File,
    type_width_bits: usize,
    block_variant_ct: usize,
    needed_in_block: usize,
    record_types: &mut Vec<u8>,
) -> Result<()> {
    if type_width_bits == 8 {
        let mut types = vec![0_u8; needed_in_block];
        file.read_exact(&mut types)
            .map_err(|source| GenoioError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        record_types.extend(types);
        skip_bytes(path, file, block_variant_ct - needed_in_block)?;
        return Ok(());
    }

    let packed_needed = needed_in_block.div_ceil(2);
    let mut packed_types = vec![0_u8; packed_needed];
    file.read_exact(&mut packed_types)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    for variant_in_block in 0..needed_in_block {
        let byte = packed_types[variant_in_block / 2];
        let record_type = if variant_in_block % 2 == 0 {
            byte & 0x0f
        } else {
            byte >> 4
        };
        record_types.push(record_type);
    }
    // Four-bit type tables pack two variants per byte, so skipping must use
    // packed byte counts rather than raw variant counts.
    skip_bytes(path, file, block_variant_ct.div_ceil(2) - packed_needed)?;
    Ok(())
}

fn skip_bytes(path: &Path, file: &mut File, len: usize) -> Result<()> {
    let offset = i64::try_from(len)
        .map_err(|_| GenoioError::invalid_source(path, "pgen skip is out of range"))?;
    file.seek(SeekFrom::Current(offset))
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(())
}

fn block_variant_count(variant_ct: usize, block_index: usize) -> usize {
    let block_start = block_index * PGEN_VARIANT_BLOCK_SIZE;
    (variant_ct - block_start).min(PGEN_VARIANT_BLOCK_SIZE)
}

fn validate_supported_variable_record_type(path: &Path, record_type: u8) -> Result<()> {
    if record_type & 0x08 != 0 {
        return Err(GenoioError::invalid_source(
            path,
            "unsupported pgen multiallelic hard-call patch set",
        ));
    }
    let dosage_bits = (record_type >> 5) & 0x03;
    if dosage_bits != 0 && record_type & 0x10 != 0 {
        return Err(GenoioError::invalid_source(
            path,
            "unsupported pgen hardcall-phase track with dosage",
        ));
    }
    if record_type & 0x80 != 0 && dosage_bits != 2 {
        return Err(GenoioError::invalid_source(
            path,
            "unsupported pgen phased-dosage track without full dosage track",
        ));
    }
    match record_type & 0x07 {
        0 | 1 | 2 | 3 | 4 | 6 | 7 => Ok(()),
        compression => Err(GenoioError::invalid_source(
            path,
            format!("unsupported pgen main-track compression type {compression}"),
        )),
    }
}

fn fixed_width_dosage_record_len(sample_ct: usize) -> usize {
    sample_ct.div_ceil(4) + sample_ct * 2
}

fn fixed_width_phased_dosage_record_len(sample_ct: usize) -> usize {
    sample_ct.div_ceil(4) + sample_ct * 4
}

fn fixed_width_record_len(header: &PgenHeader) -> usize {
    match header.layout {
        PgenLayout::FixedWidth => header.bytes_per_variant,
        PgenLayout::FixedWidthDosage => fixed_width_dosage_record_len(header.sample_ct),
        PgenLayout::FixedWidthPhasedDosage => {
            fixed_width_phased_dosage_record_len(header.sample_ct)
        }
        PgenLayout::VariableWidth => header.bytes_per_variant,
    }
}

fn validate_fixed_width_pgen_payload_len(
    path: &Path,
    file: &File,
    variant_ct: usize,
    bytes_per_record: usize,
) -> Result<()> {
    let payload_len = variant_ct
        .checked_mul(bytes_per_record)
        .ok_or_else(|| GenoioError::invalid_source(path, "pgen payload length is out of range"))?;
    let expected_len = PGEN_HEADER_LEN
        .checked_add(u64::try_from(payload_len).map_err(|_| {
            GenoioError::invalid_source(path, "pgen payload length is out of range")
        })?)
        .ok_or_else(|| GenoioError::invalid_source(path, "pgen payload length is out of range"))?;
    let actual_len = file
        .metadata()
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if actual_len != expected_len {
        return Err(GenoioError::invalid_source(
            path,
            format!("pgen payload length {actual_len} does not match fixed-width header"),
        ));
    }
    Ok(())
}

pub(super) fn validate_plink2_dimensions(
    path: &Path,
    header: &PgenHeader,
    sample_ct: usize,
    variant_ct: usize,
) -> Result<()> {
    validate_plink2_sample_count(path, header, sample_ct)?;
    if header.variant_ct != variant_ct {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "pgen variant count {} does not match pvar variant count {variant_ct}",
                header.variant_ct
            ),
        ));
    }
    Ok(())
}

pub(super) fn validate_plink2_sample_count(
    path: &Path,
    header: &PgenHeader,
    sample_ct: usize,
) -> Result<()> {
    if header.sample_ct != sample_ct {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "pgen sample count {} does not match psam sample count {sample_ct}",
                header.sample_ct
            ),
        ));
    }
    Ok(())
}

pub(super) fn open_pgen_payload(path: &Path) -> Result<File> {
    let mut file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    file.seek(SeekFrom::Start(PGEN_HEADER_LEN))
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(file)
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
        &mut decoder_state.missing,
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

pub(super) fn read_plink2_variant_packed(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    match header.layout {
        PgenLayout::FixedWidth
        | PgenLayout::FixedWidthDosage
        | PgenLayout::FixedWidthPhasedDosage => {
            read_fixed_width_variant_packed(path, file, header, variant_index, decoder_state)
        }
        PgenLayout::VariableWidth => {
            read_variable_width_variant_packed(path, file, header, variant_index, decoder_state)
        }
    }
}

fn read_fixed_width_variant_packed(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    seek_fixed_width_variant_record(path, file, header, variant_index)?;
    read_fixed_width_variant_packed_sequential(path, file, header, decoder_state)
}

pub(super) fn seek_fixed_width_variant_record(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
) -> Result<()> {
    let payload_offset = variant_index
        .checked_mul(fixed_width_record_len(header))
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
    Ok(())
}

pub(super) fn read_fixed_width_variant_packed_sequential(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    decoder_state
        .record
        .resize(fixed_width_record_len(header), 0);
    file.read_exact(&mut decoder_state.record)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    decoder_state.packed.load_pgen_payload(
        &decoder_state.record[..header.bytes_per_variant],
        header.sample_ct,
    );
    Ok(())
}

fn read_variable_width_variant_packed(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
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
    let compression = record_type & 0x07;
    match compression {
        0 => {
            if record.len() < header.bytes_per_variant {
                return Err(GenoioError::invalid_source(
                    path,
                    "pgen uncompressed record is shorter than expected",
                ));
            }
            decoder_state
                .packed
                .load_pgen_payload(&record[..header.bytes_per_variant], header.sample_ct);
        }
        1 => decode_one_bit_record(path, record, header.sample_ct, &mut decoder_state.packed)?,
        2 | 3 => {
            if !decoder_state.has_previous_non_ld {
                return Err(GenoioError::invalid_source(
                    path,
                    "pgen LD-compressed record appears before any non-LD record",
                ));
            }
            decode_ld_compressed_record(
                path,
                record,
                header.sample_ct,
                &decoder_state.previous_non_ld_packed,
                compression == 3,
                &mut decoder_state.packed,
            )?;
        }
        4 => decode_difflist_record(path, record, header.sample_ct, 0, &mut decoder_state.packed)?,
        6 => decode_difflist_record(path, record, header.sample_ct, 2, &mut decoder_state.packed)?,
        7 => decode_difflist_record(path, record, header.sample_ct, 3, &mut decoder_state.packed)?,
        other => {
            return Err(GenoioError::invalid_source(
                path,
                format!("unsupported pgen main-track compression type {other}"),
            ));
        }
    }
    if decoder_state.packed.sample_ct() != header.sample_ct {
        return Err(GenoioError::invalid_source(
            path,
            "pgen decoded category count does not match sample count",
        ));
    }
    if !matches!(compression, 2 | 3) {
        decoder_state
            .previous_non_ld_packed
            .copy_from(&decoder_state.packed);
        decoder_state.has_previous_non_ld = true;
    }
    Ok(())
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
        &mut decoder_state.missing,
    );
    overlay_fixed_width_dosages(
        path,
        &decoder_state.record[header.bytes_per_variant..],
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing,
    )
}

fn read_fixed_width_phased_dosage_variant_record(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<usize> {
    let record_len = fixed_width_phased_dosage_record_len(header.sample_ct);
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
    if decoder_state.packed.sample_ct() != header.sample_ct {
        return Err(GenoioError::invalid_source(
            path,
            "pgen decoded category count does not match sample count",
        ));
    }
    Ok(header.bytes_per_variant)
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
        &mut decoder_state.missing,
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
        &mut decoder_state.missing,
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
        &mut decoder_state.missing,
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
            missing: &mut decoder_state.missing,
        },
    )
}

pub(super) fn read_plink2_variant_haplotype_main_track(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<usize> {
    validate_variable_width_haplotype_layout(path, header)?;
    read_variable_width_record(path, file, header, variant_index, decoder_state)?;
    let cursor = decode_variable_width_main_track(
        path,
        decoder_state.record.as_slice(),
        header.record_types[variant_index],
        header.sample_ct,
        &decoder_state.previous_non_ld_packed,
        decoder_state.has_previous_non_ld,
        &mut decoder_state.packed,
    )?;
    update_variable_width_ld_state(header.record_types[variant_index], decoder_state);
    Ok(cursor)
}

pub(super) fn read_plink2_variant_haplotype_dosage_track(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<usize> {
    match header.layout {
        PgenLayout::FixedWidthPhasedDosage => read_fixed_width_phased_dosage_variant_record(
            path,
            file,
            header,
            variant_index,
            decoder_state,
        ),
        PgenLayout::VariableWidth => read_plink2_variant_haplotype_main_track(
            path,
            file,
            header,
            variant_index,
            decoder_state,
        ),
        PgenLayout::FixedWidth | PgenLayout::FixedWidthDosage => Err(GenoioError::unsupported(
            "plink2 haplotype dosage reads require explicit phased dosage records",
        )),
    }
}

fn validate_variable_width_haplotype_layout(_path: &Path, header: &PgenHeader) -> Result<()> {
    if matches!(header.layout, PgenLayout::VariableWidth) {
        Ok(())
    } else {
        Err(GenoioError::unsupported(
            "plink2 haplotype reads require variable-width explicit phased records",
        ))
    }
}

pub(super) fn decode_plink2_haplotype_hardcall_aux(
    path: &Path,
    header: &PgenHeader,
    variant_index: usize,
    cursor: usize,
    source_indices: &[usize],
    decoder_state: &PgenDecoderState,
    haplotype_state: &mut PgenHaplotypeDecodeState,
) -> Result<()> {
    let record_type = header.record_types[variant_index];
    if ((record_type >> 5) & 0x03) != 0 || record_type & 0x80 != 0 {
        return Err(GenoioError::unsupported(
            "pgen haplotype hardcall read does not accept dosage records",
        ));
    }
    if record_type & 0x10 == 0 {
        return Err(GenoioError::unsupported(
            "unphased pgen hardcall record retained in haplotype read",
        ));
    }
    decode_hardcall_phase_track(
        path,
        decoder_state.record.as_slice(),
        cursor,
        source_indices,
        &decoder_state.packed,
        haplotype_state,
    )
}

pub(super) fn decode_plink2_haplotype_dosage_aux(
    path: &Path,
    header: &PgenHeader,
    variant_index: usize,
    cursor: usize,
    source_indices: &[usize],
    decoder_state: &PgenDecoderState,
    haplotype_state: &mut PgenHaplotypeDecodeState,
) -> Result<()> {
    if matches!(header.layout, PgenLayout::FixedWidthPhasedDosage) {
        return decode_full_phased_dosage_tracks(
            path,
            decoder_state.record.as_slice(),
            cursor,
            header.sample_ct,
            source_indices,
            haplotype_state,
        );
    }
    let record_type = header.record_types[variant_index];
    let dosage_bits = (record_type >> 5) & 0x03;
    if record_type & 0x80 == 0 {
        return Err(GenoioError::unsupported(
            "pgen record does not contain explicit phased dosage values",
        ));
    }
    if dosage_bits != 2 {
        return Err(GenoioError::unsupported(
            "unsupported pgen phased dosage representation; only full dosage tracks are supported",
        ));
    }
    decode_full_phased_dosage_tracks(
        path,
        decoder_state.record.as_slice(),
        cursor,
        header.sample_ct,
        source_indices,
        haplotype_state,
    )
}

fn read_variable_width_record(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    if variant_index >= header.variant_ct || variant_index + 1 >= header.record_offsets.len() {
        return Err(GenoioError::invalid_source(
            path,
            "pgen variant index is outside variable-width record table",
        ));
    }
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
    Ok(())
}

fn update_variable_width_ld_state(record_type: u8, decoder_state: &mut PgenDecoderState) {
    if !matches!(record_type & 0x07, 2 | 3) {
        decoder_state
            .previous_non_ld_packed
            .copy_from(&decoder_state.packed);
        decoder_state.has_previous_non_ld = true;
    }
}

fn decode_hardcall_phase_track(
    path: &Path,
    record: &[u8],
    cursor: usize,
    source_indices: &[usize],
    packed: &PackedGenotypes,
    haplotype_state: &mut PgenHaplotypeDecodeState,
) -> Result<()> {
    let heterozygote_ct = (0..packed.sample_ct())
        .filter(|sample_index| packed.get(*sample_index) == 1)
        .count();
    ensure_record_bytes(path, record, cursor, 1)?;
    let phasepresent_stored = bit_is_set(&record[cursor..], 0);
    let phaseinfo_start_bit = if phasepresent_stored {
        let phasepresent_bits = 1 + heterozygote_ct;
        ensure_record_bytes(path, record, cursor, phasepresent_bits.div_ceil(8))?;
        cursor
            .checked_add(phasepresent_bits.div_ceil(8))
            .ok_or_else(|| GenoioError::invalid_source(path, "pgen phase offset is out of range"))?
            * 8
    } else {
        cursor
            .checked_mul(8)
            .and_then(|bit| bit.checked_add(1))
            .ok_or_else(|| GenoioError::invalid_source(path, "pgen phase offset is out of range"))?
    };
    ensure_record_bits(path, record, phaseinfo_start_bit, heterozygote_ct)?;

    haplotype_state.selected_haplotype_values.clear();
    haplotype_state.selected_haplotype_missing.clear();
    haplotype_state.selected_collapsed_values.clear();
    haplotype_state.selected_collapsed_missing.clear();
    haplotype_state
        .selected_haplotype_values
        .resize(source_indices.len() * 2, 0.0);
    haplotype_state
        .selected_haplotype_missing
        .resize(source_indices.len() * 2, false);

    let mut selected_cursor = SelectedSampleCursor::new(source_indices);
    let mut heterozygote_index = 0_usize;
    let mut phased_heterozygote_index = 0_usize;
    for sample_index in 0..packed.sample_ct() {
        let category = packed.get(sample_index);
        let selected_index = selected_cursor.selected_index_for(sample_index);
        match category {
            0 => {
                if let Some(selected_index) = selected_index {
                    set_selected_haplotype_pair(haplotype_state, selected_index, 0.0, 0.0, false);
                }
            }
            1 => {
                let phase_present = if phasepresent_stored {
                    bit_is_set_from_abs(record, cursor * 8 + 1 + heterozygote_index)
                } else {
                    true
                };
                heterozygote_index += 1;
                if !phase_present && selected_index.is_some() {
                    return Err(GenoioError::unsupported(
                        "unphased pgen heterozygous hardcall retained in haplotype read",
                    ));
                }
                if phase_present {
                    let swapped = bit_is_set_from_abs(
                        record,
                        phaseinfo_start_bit + phased_heterozygote_index,
                    );
                    phased_heterozygote_index += 1;
                    if let Some(selected_index) = selected_index {
                        if swapped {
                            set_selected_haplotype_pair(
                                haplotype_state,
                                selected_index,
                                1.0,
                                0.0,
                                false,
                            );
                        } else {
                            set_selected_haplotype_pair(
                                haplotype_state,
                                selected_index,
                                0.0,
                                1.0,
                                false,
                            );
                        }
                    }
                }
            }
            2 => {
                if let Some(selected_index) = selected_index {
                    set_selected_haplotype_pair(haplotype_state, selected_index, 1.0, 1.0, false);
                }
            }
            3 => {
                if let Some(selected_index) = selected_index {
                    set_selected_haplotype_pair(haplotype_state, selected_index, 0.0, 0.0, true);
                }
            }
            _ => unreachable!("two-bit hard-call code should be masked"),
        }
    }
    let end_bit = phaseinfo_start_bit + phased_heterozygote_index;
    if end_bit.div_ceil(8) != record.len() {
        return Err(GenoioError::invalid_source(
            path,
            "pgen phased hardcall record has trailing or missing bytes",
        ));
    }
    for source_index in source_indices {
        match packed.get(*source_index) {
            0 => push_collapsed_dosage(haplotype_state, 0.0, false),
            1 => push_collapsed_dosage(haplotype_state, 1.0, false),
            2 => push_collapsed_dosage(haplotype_state, 2.0, false),
            3 => push_collapsed_dosage(haplotype_state, 0.0, true),
            _ => unreachable!("two-bit hard-call code should be masked"),
        }
    }
    Ok(())
}

fn decode_full_phased_dosage_tracks(
    path: &Path,
    record: &[u8],
    cursor: usize,
    sample_ct: usize,
    source_indices: &[usize],
    haplotype_state: &mut PgenHaplotypeDecodeState,
) -> Result<()> {
    let dosage_bytes_len = sample_ct.checked_mul(2).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen dosage byte count is out of range")
    })?;
    ensure_record_bytes(path, record, cursor, dosage_bytes_len)?;
    let phase_cursor = cursor.checked_add(dosage_bytes_len).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen phased dosage offset is out of range")
    })?;
    ensure_record_bytes(path, record, phase_cursor, dosage_bytes_len)?;
    if phase_cursor + dosage_bytes_len != record.len() {
        return Err(GenoioError::invalid_source(
            path,
            "pgen phased dosage record has trailing or missing bytes",
        ));
    }

    haplotype_state.selected_haplotype_values.clear();
    haplotype_state.selected_haplotype_missing.clear();
    haplotype_state.selected_collapsed_values.clear();
    haplotype_state.selected_collapsed_missing.clear();
    for source_index in source_indices {
        let dosage_offset = cursor + source_index * 2;
        let phase_offset = phase_cursor + source_index * 2;
        let dosage = decode_phased_dosage_total(
            path,
            u16::from_le_bytes([record[dosage_offset], record[dosage_offset + 1]]),
        )?;
        let phase_delta = decode_phased_dosage_delta(
            path,
            i16::from_le_bytes([record[phase_offset], record[phase_offset + 1]]),
        )?;
        let (Some(total), Some(delta)) = (dosage, phase_delta) else {
            haplotype_state.selected_haplotype_values.extend([0.0, 0.0]);
            haplotype_state
                .selected_haplotype_missing
                .extend([true, true]);
            push_collapsed_dosage(haplotype_state, 0.0, true);
            continue;
        };
        let left = (total + delta) * 0.5;
        let right = (total - delta) * 0.5;
        validate_phased_dosage_haplotype_components(path, left, right)?;
        haplotype_state
            .selected_haplotype_values
            .extend([left, right]);
        haplotype_state
            .selected_haplotype_missing
            .extend([false, false]);
        push_collapsed_dosage(haplotype_state, total, false);
    }
    Ok(())
}

fn decode_phased_dosage_total(path: &Path, raw: u16) -> Result<Option<f32>> {
    if raw == u16::MAX {
        return Ok(None);
    }
    if raw > PGEN_MAX_DOSAGE_RAW {
        return Err(GenoioError::invalid_source(
            path,
            format!("pgen phased dosage total raw value {raw} exceeds 32768"),
        ));
    }
    Ok(Some(f32::from(raw) * PGEN_DOSAGE_SCALE))
}

fn decode_phased_dosage_delta(path: &Path, raw: i16) -> Result<Option<f32>> {
    if raw == i16::MIN {
        return Ok(None);
    }
    if !(PGEN_MIN_PHASE_RAW..=PGEN_MAX_PHASE_RAW).contains(&raw) {
        return Err(GenoioError::invalid_source(
            path,
            format!("pgen phased dosage phase raw value {raw} is outside [-16384, 16384]"),
        ));
    }
    Ok(Some(f32::from(raw) * PGEN_PHASE_SCALE))
}

fn validate_phased_dosage_haplotype_components(path: &Path, left: f32, right: f32) -> Result<()> {
    if (0.0..=1.0).contains(&left) && (0.0..=1.0).contains(&right) {
        return Ok(());
    }
    Err(GenoioError::invalid_source(
        path,
        format!(
            "pgen phased dosage haplotype components are outside [0, 1]: left={left}, right={right}"
        ),
    ))
}

fn set_selected_haplotype_pair(
    haplotype_state: &mut PgenHaplotypeDecodeState,
    selected_index: usize,
    left: f32,
    right: f32,
    missing: bool,
) {
    let offset = selected_index * 2;
    haplotype_state.selected_haplotype_values[offset] = left;
    haplotype_state.selected_haplotype_values[offset + 1] = right;
    haplotype_state.selected_haplotype_missing[offset] = missing;
    haplotype_state.selected_haplotype_missing[offset + 1] = missing;
}

fn push_collapsed_dosage(
    haplotype_state: &mut PgenHaplotypeDecodeState,
    value: f32,
    missing: bool,
) {
    haplotype_state.selected_collapsed_values.push(value);
    haplotype_state.selected_collapsed_missing.push(missing);
}

fn decode_variable_width_main_track(
    path: &Path,
    record: &[u8],
    record_type: u8,
    sample_ct: usize,
    previous_non_ld_packed: &PackedGenotypes,
    has_previous_non_ld: bool,
    packed: &mut PackedGenotypes,
) -> Result<usize> {
    let compression = record_type & 0x07;
    match compression {
        0 => {
            let bytes_per_variant = sample_ct.div_ceil(4);
            if record.len() < bytes_per_variant {
                return Err(GenoioError::invalid_source(
                    path,
                    "pgen uncompressed record is shorter than expected",
                ));
            }
            packed.load_pgen_payload(&record[..bytes_per_variant], sample_ct);
            Ok(bytes_per_variant)
        }
        1 => decode_one_bit_record_with_cursor(path, record, sample_ct, packed),
        2 | 3 => {
            if !has_previous_non_ld {
                return Err(GenoioError::invalid_source(
                    path,
                    "pgen LD-compressed record appears before any non-LD record",
                ));
            }
            let mut cursor = 0;
            let entries = decode_difflist(path, record, &mut cursor, sample_ct, true)?;
            if previous_non_ld_packed.sample_ct() != sample_ct {
                return Err(GenoioError::invalid_source(
                    path,
                    "pgen LD state length does not match sample count",
                ));
            }
            packed.copy_from(previous_non_ld_packed);
            for (sample_index, category) in entries {
                packed.set(sample_index, category);
            }
            if compression == 3 {
                packed.invert_0_2();
            }
            Ok(cursor)
        }
        4 => decode_difflist_record_with_cursor(path, record, sample_ct, 0, packed),
        6 => decode_difflist_record_with_cursor(path, record, sample_ct, 2, packed),
        7 => decode_difflist_record_with_cursor(path, record, sample_ct, 3, packed),
        other => Err(GenoioError::invalid_source(
            path,
            format!("unsupported pgen main-track compression type {other}"),
        )),
    }
}

fn decode_one_bit_record_with_cursor(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    packed: &mut PackedGenotypes,
) -> Result<usize> {
    let common_categories = *record.first().ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen 1-bit record is missing common-category byte")
    })?;
    let (low_category, high_category) = match common_categories {
        1 => (0, 1),
        2 => (0, 2),
        3 => (0, 3),
        5 => (1, 2),
        6 => (1, 3),
        9 => (2, 3),
        other => {
            return Err(GenoioError::invalid_source(
                path,
                format!("invalid pgen 1-bit common-category byte {other}"),
            ));
        }
    };
    let bitarray_len = sample_ct.div_ceil(8);
    if record.len() < 1 + bitarray_len {
        return Err(GenoioError::invalid_source(
            path,
            "pgen 1-bit record is shorter than expected",
        ));
    }
    let bitarray = &record[1..1 + bitarray_len];
    packed.resize(sample_ct);
    packed.clear_to(low_category);
    for sample_index in 0..sample_ct {
        if bit_is_set(bitarray, sample_index) {
            packed.set(sample_index, high_category);
        }
    }
    let mut cursor = 1 + bitarray_len;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    Ok(cursor)
}

fn decode_difflist_record_with_cursor(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    base_category: u8,
    packed: &mut PackedGenotypes,
) -> Result<usize> {
    packed.resize(sample_ct);
    packed.clear_to(base_category);
    let mut cursor = 0;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    Ok(cursor)
}

fn overlay_variable_width_dosages(
    path: &Path,
    record: &[u8],
    mut cursor: usize,
    dosage_bits: u8,
    sample_ct: usize,
    mut target: DosageOverlayTarget<'_>,
) -> Result<()> {
    let mut selected_samples = SelectedSampleCursor::new(target.source_indices);
    match dosage_bits {
        1 => overlay_difflist_dosages(
            path,
            record,
            &mut cursor,
            sample_ct,
            &mut selected_samples,
            &mut target,
        )?,
        2 => {
            let dosage_bytes_len = sample_ct.checked_mul(2).ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen dosage byte count is out of range")
            })?;
            ensure_record_bytes(path, record, cursor, dosage_bytes_len)?;
            for sample_index in 0..sample_ct {
                let byte_index = cursor + sample_index * 2;
                let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
                overlay_selected_pgen_dosage(sample_index, raw, &mut selected_samples, &mut target);
            }
        }
        3 => {
            let bitarray_len = sample_ct.div_ceil(8);
            ensure_record_bytes(path, record, cursor, bitarray_len)?;
            let bitarray = &record[cursor..cursor + bitarray_len];
            cursor += bitarray_len;
            let dosage_ct = (0..sample_ct)
                .filter(|sample_index| bit_is_set(bitarray, *sample_index))
                .count();
            let dosage_bytes_len = dosage_ct.checked_mul(2).ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen dosage byte count is out of range")
            })?;
            ensure_record_bytes(path, record, cursor, dosage_bytes_len)?;
            let mut dosage_index = 0;
            for sample_index in 0..sample_ct {
                if bit_is_set(bitarray, sample_index) {
                    let byte_index = cursor + dosage_index * 2;
                    let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
                    overlay_selected_pgen_dosage(
                        sample_index,
                        raw,
                        &mut selected_samples,
                        &mut target,
                    );
                    dosage_index += 1;
                }
            }
        }
        other => {
            return Err(GenoioError::invalid_source(
                path,
                format!("unsupported pgen dosage track type {other}"),
            ));
        }
    }
    Ok(())
}

fn overlay_difflist_dosages(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    selected_samples: &mut SelectedSampleCursor<'_>,
    target: &mut DosageOverlayTarget<'_>,
) -> Result<()> {
    let list_len = read_base128_varint(path, record, cursor)?;
    if list_len == 0 {
        return Ok(());
    }
    if list_len > sample_ct {
        return Err(GenoioError::invalid_source(
            path,
            "pgen difflist length exceeds sample count",
        ));
    }

    let group_ct = list_len.div_ceil(64);
    let sample_id_width = sample_id_width(sample_ct);
    let mut first_ids = Vec::with_capacity(group_ct);
    for _ in 0..group_ct {
        first_ids.push(read_fixed_width_sample_id(
            path,
            record,
            cursor,
            sample_id_width,
        )?);
    }
    ensure_record_bytes(path, record, *cursor, group_ct.saturating_sub(1))?;
    *cursor += group_ct.saturating_sub(1);

    let deltas_start = *cursor;
    let mut values_start = deltas_start;
    walk_difflist_ids(
        path,
        record,
        &mut values_start,
        sample_ct,
        list_len,
        &first_ids,
        |_, _| {},
    )?;

    let dosage_bytes_len = list_len.checked_mul(2).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen dosage byte count is out of range")
    })?;
    ensure_record_bytes(path, record, values_start, dosage_bytes_len)?;

    let mut ids_cursor = deltas_start;
    // The values follow the encoded sample IDs, so we first walk IDs to find
    // values_start, then walk them again while overlaying selected samples.
    walk_difflist_ids(
        path,
        record,
        &mut ids_cursor,
        sample_ct,
        list_len,
        &first_ids,
        |sample_index, dosage_index| {
            let byte_index = values_start + dosage_index * 2;
            let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
            overlay_selected_pgen_dosage(sample_index, raw, selected_samples, target);
        },
    )?;
    *cursor = values_start + dosage_bytes_len;
    Ok(())
}

fn walk_difflist_ids(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    list_len: usize,
    first_ids: &[usize],
    mut visit: impl FnMut(usize, usize),
) -> Result<()> {
    let mut previous_sample_id = None;
    let mut entry_index = 0;
    for (group_index, first_id) in first_ids.iter().copied().enumerate() {
        let group_len = (list_len - group_index * 64).min(64);
        let mut sample_id = first_id;
        validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
        visit(sample_id, entry_index);
        entry_index += 1;
        for _ in 1..group_len {
            let delta = read_base128_varint(path, record, cursor)?;
            sample_id = sample_id.checked_add(delta).ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen difflist sample id is out of range")
            })?;
            validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
            visit(sample_id, entry_index);
            entry_index += 1;
        }
    }
    Ok(())
}

fn overlay_selected_pgen_dosage(
    source_index: usize,
    raw: u16,
    selected_samples: &mut SelectedSampleCursor<'_>,
    target: &mut DosageOverlayTarget<'_>,
) {
    if let Some(selected_index) = selected_samples.selected_index_for(source_index) {
        apply_pgen_dosage(
            raw,
            &mut target.values[selected_index],
            &mut target.missing[selected_index],
        );
    }
}

fn overlay_fixed_width_dosages(
    path: &Path,
    dosage_bytes: &[u8],
    source_indices: &[usize],
    values: &mut [f32],
    missing: &mut [bool],
) -> Result<()> {
    for (selected_index, source_index) in source_indices.iter().copied().enumerate() {
        let byte_index = source_index.checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source(path, "pgen dosage offset is out of range")
        })?;
        ensure_record_bytes(path, dosage_bytes, byte_index, 2)?;
        let raw = u16::from_le_bytes([dosage_bytes[byte_index], dosage_bytes[byte_index + 1]]);
        apply_pgen_dosage(
            raw,
            &mut values[selected_index],
            &mut missing[selected_index],
        );
    }
    Ok(())
}

fn apply_pgen_dosage(raw: u16, value: &mut f32, is_missing: &mut bool) {
    if raw == u16::MAX {
        *value = 0.0;
        *is_missing = true;
        return;
    }
    *value = f32::from(raw) * (2.0 / 32768.0);
    *is_missing = false;
}

fn decode_one_bit_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    let common_categories = *record.first().ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen 1-bit record is missing common-category byte")
    })?;
    let (low_category, high_category) = match common_categories {
        1 => (0, 1),
        2 => (0, 2),
        3 => (0, 3),
        5 => (1, 2),
        6 => (1, 3),
        9 => (2, 3),
        other => {
            return Err(GenoioError::invalid_source(
                path,
                format!("invalid pgen 1-bit common-category byte {other}"),
            ));
        }
    };
    let bitarray_len = sample_ct.div_ceil(8);
    if record.len() < 1 + bitarray_len {
        return Err(GenoioError::invalid_source(
            path,
            "pgen 1-bit record is shorter than expected",
        ));
    }
    let bitarray = &record[1..1 + bitarray_len];
    packed.resize(sample_ct);
    packed.clear_to(low_category);
    for sample_index in 0..sample_ct {
        if bit_is_set(bitarray, sample_index) {
            packed.set(sample_index, high_category);
        }
    }
    let mut cursor = 1 + bitarray_len;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    Ok(())
}

fn decode_ld_compressed_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    previous_non_ld_packed: &PackedGenotypes,
    inverted: bool,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    if previous_non_ld_packed.sample_ct() != sample_ct {
        return Err(GenoioError::invalid_source(
            path,
            "pgen LD state length does not match sample count",
        ));
    }
    packed.copy_from(previous_non_ld_packed);
    let mut cursor = 0;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    if inverted {
        packed.invert_0_2();
    }
    Ok(())
}

fn decode_difflist_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    base_category: u8,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    packed.resize(sample_ct);
    packed.clear_to(base_category);
    let mut cursor = 0;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    Ok(())
}

fn decode_difflist(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    with_values: bool,
) -> Result<Vec<(usize, u8)>> {
    let list_len = read_base128_varint(path, record, cursor)?;
    if list_len == 0 {
        return Ok(Vec::new());
    }
    if list_len > sample_ct {
        return Err(GenoioError::invalid_source(
            path,
            "pgen difflist length exceeds sample count",
        ));
    }
    let group_ct = list_len.div_ceil(64);
    let sample_id_width = sample_id_width(sample_ct);
    let mut first_ids = Vec::with_capacity(group_ct);
    for _ in 0..group_ct {
        first_ids.push(read_fixed_width_sample_id(
            path,
            record,
            cursor,
            sample_id_width,
        )?);
    }
    ensure_record_bytes(path, record, *cursor, group_ct.saturating_sub(1))?;
    *cursor += group_ct.saturating_sub(1);

    let packed_values_start = *cursor;
    if with_values {
        ensure_record_bytes(path, record, *cursor, list_len.div_ceil(4))?;
        *cursor += list_len.div_ceil(4);
    }

    let mut entries = Vec::with_capacity(list_len);
    let mut previous_sample_id = None;
    for (group_index, first_id) in first_ids.into_iter().enumerate() {
        let group_len = (list_len - group_index * 64).min(64);
        let mut sample_id = first_id;
        validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
        entries.push((
            sample_id,
            packed_difflist_value(record, packed_values_start, entries.len(), with_values),
        ));
        for _ in 1..group_len {
            let delta = read_base128_varint(path, record, cursor)?;
            sample_id = sample_id.checked_add(delta).ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen difflist sample id is out of range")
            })?;
            validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
            entries.push((
                sample_id,
                packed_difflist_value(record, packed_values_start, entries.len(), with_values),
            ));
        }
    }
    Ok(entries)
}

fn read_base128_varint(path: &Path, record: &[u8], cursor: &mut usize) -> Result<usize> {
    let mut value = 0_usize;
    let mut shift = 0_u32;
    loop {
        ensure_record_bytes(path, record, *cursor, 1)?;
        let byte = record[*cursor];
        *cursor += 1;
        value |= usize::from(byte & 0x7f)
            .checked_shl(shift)
            .ok_or_else(|| GenoioError::invalid_source(path, "pgen varint is out of range"))?;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
        if shift >= usize::BITS {
            return Err(GenoioError::invalid_source(
                path,
                "pgen varint is out of range",
            ));
        }
    }
}

fn read_fixed_width_sample_id(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    width: usize,
) -> Result<usize> {
    ensure_record_bytes(path, record, *cursor, width)?;
    let mut value = 0_usize;
    for byte_index in 0..width {
        value |= usize::from(record[*cursor + byte_index]) << (8 * byte_index);
    }
    *cursor += width;
    Ok(value)
}

fn sample_id_width(sample_ct: usize) -> usize {
    if sample_ct <= 1 << 8 {
        1
    } else if sample_ct <= 1 << 16 {
        2
    } else if sample_ct <= 1 << 24 {
        3
    } else {
        4
    }
}

fn packed_difflist_value(record: &[u8], start: usize, index: usize, with_values: bool) -> u8 {
    if !with_values {
        return 0;
    }
    (record[start + index / 4] >> ((index % 4) * 2)) & 0b11
}

fn validate_difflist_sample_id(
    path: &Path,
    sample_id: usize,
    sample_ct: usize,
    previous_sample_id: &mut Option<usize>,
) -> Result<()> {
    if sample_id >= sample_ct {
        return Err(GenoioError::invalid_source(
            path,
            "pgen difflist sample id is outside sample count",
        ));
    }
    if previous_sample_id.is_some_and(|previous| sample_id <= previous) {
        return Err(GenoioError::invalid_source(
            path,
            "pgen difflist sample ids must be strictly increasing",
        ));
    }
    *previous_sample_id = Some(sample_id);
    Ok(())
}

fn ensure_record_bytes(path: &Path, record: &[u8], cursor: usize, len: usize) -> Result<()> {
    if cursor.checked_add(len).is_none_or(|end| end > record.len()) {
        return Err(GenoioError::invalid_source(
            path,
            "pgen record ended before expected data",
        ));
    }
    Ok(())
}

fn ensure_record_bits(path: &Path, record: &[u8], start_bit: usize, len: usize) -> Result<()> {
    let end_bit = start_bit
        .checked_add(len)
        .ok_or_else(|| GenoioError::invalid_source(path, "pgen bit range is out of range"))?;
    if end_bit > record.len() * 8 {
        return Err(GenoioError::invalid_source(
            path,
            "pgen record ended before expected bitarray data",
        ));
    }
    Ok(())
}

fn bit_is_set(bytes: &[u8], bit_index: usize) -> bool {
    bytes[bit_index / 8] & (1 << (bit_index % 8)) != 0
}

fn bit_is_set_from_abs(bytes: &[u8], bit_index: usize) -> bool {
    bytes[bit_index / 8] & (1 << (bit_index % 8)) != 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variable_record_helpers_write_packed_genotypes() {
        let path = Path::new("test.pgen");
        let mut packed = PackedGenotypes::default();

        decode_one_bit_record(path, &[2, 0b0000_1010, 0], 4, &mut packed)
            .expect("one-bit record should decode");
        assert_eq!(
            (0..4)
                .map(|sample_index| packed.get(sample_index))
                .collect::<Vec<_>>(),
            vec![0, 2, 0, 2]
        );

        decode_difflist_record(path, &[2, 1, 9, 2], 4, 0, &mut packed)
            .expect("difflist record should decode");
        assert_eq!(
            (0..4)
                .map(|sample_index| packed.get(sample_index))
                .collect::<Vec<_>>(),
            vec![0, 1, 0, 2]
        );

        let mut previous = PackedGenotypes::default();
        previous.resize(4);
        previous.clear_to(0);
        previous.set(1, 1);
        previous.set(2, 2);
        previous.set(3, 3);

        decode_ld_compressed_record(path, &[1, 2, 0], 4, &previous, true, &mut packed)
            .expect("LD-compressed record should decode");
        assert_eq!(
            (0..4)
                .map(|sample_index| packed.get(sample_index))
                .collect::<Vec<_>>(),
            vec![2, 1, 2, 3]
        );
    }

    #[test]
    fn variable_width_dosage_overlay_preserves_hardcall_inferred_values() {
        let path = Path::new("test.pgen");
        let mut record = vec![0b0000_0110];
        record.extend(100_u16.to_le_bytes());
        record.extend(200_u16.to_le_bytes());
        let mut values = vec![0.0, 2.0, 0.0];
        let mut missing = vec![false, false, true];

        overlay_variable_width_dosages(
            path,
            &record,
            0,
            3,
            4,
            DosageOverlayTarget {
                source_indices: &[0, 2, 3],
                values: &mut values,
                missing: &mut missing,
            },
        )
        .expect("dosage overlay should decode");

        assert_eq!(values, vec![0.0, f32::from(200_u16) * (2.0 / 32768.0), 0.0]);
        assert_eq!(missing, vec![false, false, true]);
    }

    #[test]
    fn variable_width_dosage_list_overlay_uses_source_order_without_dense_index() {
        let path = Path::new("test.pgen");
        let mut record = vec![3, 1, 3, 5];
        record.extend(100_u16.to_le_bytes());
        record.extend(200_u16.to_le_bytes());
        record.extend(300_u16.to_le_bytes());
        let mut values = vec![0.0, 1.0];
        let mut missing = vec![false, false];

        overlay_variable_width_dosages(
            path,
            &record,
            0,
            1,
            10,
            DosageOverlayTarget {
                source_indices: &[4, 9],
                values: &mut values,
                missing: &mut missing,
            },
        )
        .expect("dosage-list overlay should decode");

        assert_eq!(
            values,
            vec![
                f32::from(200_u16) * (2.0 / 32768.0),
                f32::from(300_u16) * (2.0 / 32768.0),
            ]
        );
        assert_eq!(missing, vec![false, false]);
    }
}
