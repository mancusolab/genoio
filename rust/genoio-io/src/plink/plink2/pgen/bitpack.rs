// pattern: Functional Core
//! Bit and variable-width integer helpers for PGEN record payloads.

use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;

pub(super) fn read_base128_varint(path: &Path, record: &[u8], cursor: &mut usize) -> Result<usize> {
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

pub(super) fn read_fixed_width_sample_id(
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

pub(super) fn sample_id_width(sample_ct: usize) -> usize {
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

pub(super) fn packed_difflist_value(
    record: &[u8],
    start: usize,
    index: usize,
    with_values: bool,
) -> u8 {
    if !with_values {
        return 0;
    }
    (record[start + index / 4] >> ((index % 4) * 2)) & 0b11
}

pub(super) fn validate_difflist_sample_id(
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

pub(super) fn ensure_record_bytes(
    path: &Path,
    record: &[u8],
    cursor: usize,
    len: usize,
) -> Result<()> {
    if cursor.checked_add(len).is_none_or(|end| end > record.len()) {
        return Err(GenoioError::invalid_source(
            path,
            "pgen record ended before expected data",
        ));
    }
    Ok(())
}

pub(super) fn ensure_record_bits(
    path: &Path,
    record: &[u8],
    start_bit: usize,
    len: usize,
) -> Result<()> {
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

pub(super) fn bit_is_set(bytes: &[u8], bit_index: usize) -> bool {
    bytes[bit_index / 8] & (1 << (bit_index % 8)) != 0
}
