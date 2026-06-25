// pattern: Functional Core
//! Bit and variable-width integer helpers for PGEN record payloads.
//!
//! PGEN variable-width records share base-128 varints, fixed-width sample IDs,
//! packed two-bit values, and bounds checks. Keeping them here makes malformed
//! record handling consistent across main, dosage, and haplotype tracks.

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

pub(super) fn packed_difflist_value(record: &[u8], start: usize, index: usize) -> u8 {
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

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    fn assert_error_contains(error: genoio_core::GenoioError, expected: &str) {
        let message = error.to_string();
        assert!(
            message.contains(expected),
            "expected error to contain {expected:?}, got {message:?}"
        );
    }

    #[test]
    fn base128_varint_decodes_multibyte_values_and_advances_cursor() {
        let path = Path::new("test.pgen");
        let record = [0xac, 0x02, 0xff];
        let mut cursor = 0;

        let value = read_base128_varint(path, &record, &mut cursor).unwrap();

        assert_eq!(value, 300);
        assert_eq!(cursor, 2);
    }

    #[test]
    fn base128_varint_rejects_truncated_continuation() {
        let path = Path::new("test.pgen");
        let mut cursor = 0;

        let error = read_base128_varint(path, &[0x80], &mut cursor)
            .expect_err("unterminated varint should fail");

        assert_error_contains(error, "record ended before expected data");
    }

    #[test]
    fn base128_varint_rejects_out_of_range_values() {
        let path = Path::new("test.pgen");
        let record = vec![0x80; (usize::BITS as usize / 7) + 1];
        let mut cursor = 0;

        let error = read_base128_varint(path, &record, &mut cursor)
            .expect_err("oversized varint should fail");

        assert_error_contains(error, "pgen varint is out of range");
    }

    #[test]
    fn fixed_width_sample_id_reads_little_endian_widths() {
        let path = Path::new("test.pgen");
        let mut cursor = 0;

        let value = read_fixed_width_sample_id(path, &[0x34, 0x12, 0xef], &mut cursor, 2).unwrap();

        assert_eq!(value, 0x1234);
        assert_eq!(cursor, 2);
    }

    #[test]
    fn record_bounds_helpers_reject_short_and_overflowing_ranges() {
        let path = Path::new("test.pgen");

        let short_bytes =
            ensure_record_bytes(path, &[1, 2], 1, 2).expect_err("short byte range should fail");
        assert_error_contains(short_bytes, "record ended before expected data");

        let overflowing_bytes = ensure_record_bytes(path, &[1, 2], usize::MAX, 1)
            .expect_err("overflowing byte range should fail");
        assert_error_contains(overflowing_bytes, "record ended before expected data");

        let short_bits =
            ensure_record_bits(path, &[0], 7, 2).expect_err("short bit range should fail");
        assert_error_contains(short_bits, "record ended before expected bitarray data");

        let overflowing_bits = ensure_record_bits(path, &[0], usize::MAX, 1)
            .expect_err("overflowing bit range should fail");
        assert_error_contains(overflowing_bits, "pgen bit range is out of range");
    }

    #[test]
    fn difflist_sample_id_validation_rejects_bad_order_and_bounds() {
        let path = Path::new("test.pgen");
        let mut previous = None;

        validate_difflist_sample_id(path, 2, 4, &mut previous).unwrap();

        let duplicate = validate_difflist_sample_id(path, 2, 4, &mut previous)
            .expect_err("duplicate sample id should fail");
        assert_error_contains(duplicate, "strictly increasing");

        let mut previous = Some(2);
        let out_of_order = validate_difflist_sample_id(path, 1, 4, &mut previous)
            .expect_err("out-of-order sample id should fail");
        assert_error_contains(out_of_order, "strictly increasing");

        let mut previous = None;
        let out_of_range = validate_difflist_sample_id(path, 4, 4, &mut previous)
            .expect_err("out-of-range sample id should fail");
        assert_error_contains(out_of_range, "outside sample count");
    }
}
