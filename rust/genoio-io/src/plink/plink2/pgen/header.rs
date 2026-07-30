// pattern: Imperative Shell
//! PGEN header parsing and payload shape validation.
//!
//! Headers determine fixed-width versus variable-width layout, sample and
//! variant counts, per-variant record types, and block offsets. Unsupported PGEN
//! modes are rejected before callers allocate decode buffers.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;

use super::{
    PgenHeader, PgenLayout, PGEN_HEADER_LEN, PGEN_MAGIC, PGEN_MODE_FIXED_WIDTH_DOSAGE,
    PGEN_MODE_FIXED_WIDTH_HARDCALLS, PGEN_MODE_FIXED_WIDTH_PHASED_DOSAGE, PGEN_MODE_VARIABLE_WIDTH,
    PGEN_VARIANT_BLOCK_SIZE,
};

pub(in crate::plink::plink2) fn read_supported_pgen_header(path: &Path) -> Result<PgenHeader> {
    let mut file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    read_supported_pgen_header_inner(path, &mut file, None)
}

pub(in crate::plink::plink2) fn read_supported_pgen_header_prefix(
    path: &Path,
    requested_variant_ct: usize,
) -> Result<PgenHeader> {
    let mut file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    read_supported_pgen_header_inner(path, &mut file, Some(requested_variant_ct))
}

pub(in crate::plink::plink2) fn read_supported_pgen_header_from_file(
    path: &Path,
    file: &mut File,
) -> Result<PgenHeader> {
    read_supported_pgen_header_inner(path, file, None)
}

fn read_supported_pgen_header_inner(
    path: &Path,
    file: &mut File,
    requested_variant_ct: Option<usize>,
) -> Result<PgenHeader> {
    file.seek(SeekFrom::Start(0))
        .map_err(|source| GenoioError::Io {
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

    let counts = ParsedHeaderCounts::parse(path, &header)?;
    match header[2] {
        PGEN_MODE_FIXED_WIDTH_HARDCALLS => fixed_width_header(
            path,
            file,
            counts,
            header[11],
            FixedWidthHeaderSpec {
                layout: PgenLayout::FixedWidth,
                bytes_per_record: counts.bytes_per_variant,
                unsupported_flags_message: "unsupported pgen header flags; only fixed-width biallelic hardcalls without header extensions are supported",
            },
        ),
        PGEN_MODE_FIXED_WIDTH_DOSAGE => fixed_width_header(
            path,
            file,
            counts,
            header[11],
            FixedWidthHeaderSpec {
                layout: PgenLayout::FixedWidthDosage,
                bytes_per_record: fixed_width_dosage_record_len(counts.sample_ct),
                unsupported_flags_message: "unsupported pgen header flags; only fixed-width biallelic dosage without header extensions is supported",
            },
        ),
        PGEN_MODE_FIXED_WIDTH_PHASED_DOSAGE => fixed_width_header(
            path,
            file,
            counts,
            header[11],
            FixedWidthHeaderSpec {
                layout: PgenLayout::FixedWidthPhasedDosage,
                bytes_per_record: fixed_width_phased_dosage_record_len(counts.sample_ct),
                unsupported_flags_message: "unsupported pgen header flags; only fixed-width biallelic phased dosage without header extensions is supported",
            },
        ),
        PGEN_MODE_VARIABLE_WIDTH => {
            let prefix_variant_ct =
                requested_variant_ct.map(|requested| requested.min(counts.variant_ct));
            let (record_types, record_offsets) = match prefix_variant_ct {
                Some(prefix_variant_ct) => read_variable_width_header_body_prefix(
                    path,
                    file,
                    counts.variant_ct,
                    header[11],
                    prefix_variant_ct,
                )?,
                None => {
                    read_variable_width_header_body(path, file, counts.variant_ct, header[11])?
                }
            };
            Ok(PgenHeader {
                layout: PgenLayout::VariableWidth,
                variant_ct: counts.variant_ct,
                sample_ct: counts.sample_ct,
                bytes_per_variant: counts.bytes_per_variant,
                record_types,
                record_offsets,
            })
        }
        mode => Err(GenoioError::invalid_source(
            path,
            unsupported_mode_message(mode, requested_variant_ct.is_some()),
        )),
    }
}

fn fixed_width_header(
    path: &Path,
    file: &File,
    counts: ParsedHeaderCounts,
    header_flags: u8,
    spec: FixedWidthHeaderSpec,
) -> Result<PgenHeader> {
    if header_flags != 0 {
        return Err(GenoioError::invalid_source(
            path,
            spec.unsupported_flags_message,
        ));
    }
    validate_fixed_width_pgen_payload_len(path, file, counts.variant_ct, spec.bytes_per_record)?;
    Ok(PgenHeader {
        layout: spec.layout,
        variant_ct: counts.variant_ct,
        sample_ct: counts.sample_ct,
        bytes_per_variant: counts.bytes_per_variant,
        record_types: Vec::new(),
        record_offsets: Vec::new(),
    })
}

#[derive(Clone, Copy)]
struct ParsedHeaderCounts {
    variant_ct: usize,
    sample_ct: usize,
    bytes_per_variant: usize,
}

impl ParsedHeaderCounts {
    fn parse(path: &Path, header: &[u8; PGEN_HEADER_LEN as usize]) -> Result<Self> {
        let (variant_ct, sample_ct) = parse_pgen_header_counts(path, header)?;
        Ok(Self {
            variant_ct,
            sample_ct,
            bytes_per_variant: sample_ct.div_ceil(4),
        })
    }
}

struct FixedWidthHeaderSpec {
    layout: PgenLayout,
    bytes_per_record: usize,
    unsupported_flags_message: &'static str,
}

fn unsupported_mode_message(mode: u8, prefix: bool) -> String {
    if prefix {
        format!(
            "unsupported pgen mode 0x{mode:02x}; only fixed-width and variable-width biallelic hardcalls or unphased dosages are supported"
        )
    } else {
        format!(
            "unsupported pgen mode 0x{mode:02x}; only fixed-width and variable-width biallelic hardcalls are supported"
        )
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
    let (type_width_bits, length_width) = parse_variable_width_header_format(path, header_format)?;
    let block_ct = variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let block_offsets = read_block_offsets(path, file, block_ct)?;

    let mut record_types = Vec::with_capacity(variant_ct);
    let mut record_lengths = Vec::with_capacity(variant_ct);
    for block_index in 0..block_ct {
        let block_variant_ct = block_variant_count(variant_ct, block_index);
        read_variable_record_types(
            path,
            file,
            type_width_bits,
            block_variant_ct,
            &mut record_types,
        )?;

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
    for (variant_index, record_type) in record_types.iter().copied().enumerate() {
        validate_supported_variable_record_type_at(path, variant_index, record_type)?;
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
    validate_record_end(path, file, record_offsets[variant_ct])?;

    Ok((record_types, record_offsets))
}

fn read_variable_width_header_body_prefix(
    path: &Path,
    file: &mut File,
    variant_ct: usize,
    header_format: u8,
    prefix_variant_ct: usize,
) -> Result<(Vec<u8>, Vec<u64>)> {
    let (type_width_bits, length_width) = parse_variable_width_header_format(path, header_format)?;
    let block_ct = variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let block_offsets = read_block_offsets(path, file, block_ct)?;

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
    for (variant_index, record_type) in record_types.iter().copied().enumerate() {
        validate_supported_variable_record_type_at(path, variant_index, record_type)?;
    }
    if let Some(prefix_end) = record_offsets.last() {
        validate_record_end(path, file, *prefix_end)?;
    }
    Ok((record_types, record_offsets))
}

fn parse_variable_width_header_format(path: &Path, header_format: u8) -> Result<(usize, usize)> {
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
    Ok((type_width_bits, length_width))
}

fn read_block_offsets(path: &Path, file: &mut File, block_ct: usize) -> Result<Vec<u64>> {
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
    Ok(block_offsets)
}

fn read_variable_record_types(
    path: &Path,
    file: &mut File,
    type_width_bits: usize,
    block_variant_ct: usize,
    record_types: &mut Vec<u8>,
) -> Result<()> {
    if type_width_bits == 8 {
        let mut types = vec![0_u8; block_variant_ct];
        file.read_exact(&mut types)
            .map_err(|source| GenoioError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        record_types.extend(types);
        return Ok(());
    }

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
    Ok(())
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
    if record_type & 0x80 != 0 && dosage_bits == 0 {
        return Err(GenoioError::invalid_source(
            path,
            "pgen phased-dosage track is present without a dosage track",
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

fn validate_supported_variable_record_type_at(
    path: &Path,
    variant_index: usize,
    record_type: u8,
) -> Result<()> {
    validate_supported_variable_record_type(path, record_type)?;
    if variant_index.is_multiple_of(PGEN_VARIANT_BLOCK_SIZE) && matches!(record_type & 0x07, 2 | 3)
    {
        return Err(GenoioError::invalid_source(
            path,
            "pgen LD-compressed record appears before any non-LD record in its variant block; LD compression is forbidden for the first record of a variant block",
        ));
    }
    Ok(())
}

pub(super) fn fixed_width_dosage_record_len(sample_ct: usize) -> usize {
    sample_ct.div_ceil(4) + sample_ct * 2
}

pub(super) fn fixed_width_phased_dosage_record_len(sample_ct: usize) -> usize {
    sample_ct.div_ceil(4) + sample_ct * 4
}

pub(super) fn fixed_width_record_len(header: &PgenHeader) -> usize {
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

fn validate_record_end(path: &Path, file: &File, record_end: u64) -> Result<()> {
    let actual_len = file
        .metadata()
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if record_end > actual_len {
        return Err(GenoioError::invalid_source(
            path,
            "pgen variable-width records extend past end of file",
        ));
    }
    Ok(())
}

pub(in crate::plink::plink2) fn validate_plink2_dimensions(
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

pub(in crate::plink::plink2) fn validate_plink2_sample_count(
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

pub(in crate::plink::plink2) fn open_pgen_payload(path: &Path) -> Result<File> {
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
