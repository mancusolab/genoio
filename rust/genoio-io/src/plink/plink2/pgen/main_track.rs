// pattern: Functional Core
//! PGEN hard-call main-track decompression.
//!
//! Main-track records provide the hard-call baseline for dense, sparse, dosage,
//! and haplotype reads. This module handles uncompressed, one-bit, LD-compressed,
//! and difflist encodings.

use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;
use crate::hardcall::PackedHardcalls as PackedGenotypes;

use super::bitpack::{
    bit_is_set, ensure_record_bytes, packed_difflist_value, read_base128_varint,
    read_fixed_width_sample_id, sample_id_width, validate_difflist_sample_id,
};

pub(super) fn decode_variable_width_main_track(
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
        2 | 3 => decode_ld_compressed_record_with_cursor(
            path,
            record,
            sample_ct,
            previous_non_ld_packed,
            compression == 3,
            has_previous_non_ld,
            packed,
        ),
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
    apply_difflist_entries(path, record, &mut cursor, sample_ct, packed)?;
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
    apply_difflist_entries(path, record, &mut cursor, sample_ct, packed)?;
    Ok(cursor)
}

pub(super) fn decode_one_bit_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    decode_one_bit_record_with_cursor(path, record, sample_ct, packed).map(|_| ())
}

pub(super) fn decode_ld_compressed_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    previous_non_ld_packed: &PackedGenotypes,
    inverted: bool,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    decode_ld_compressed_record_with_cursor(
        path,
        record,
        sample_ct,
        previous_non_ld_packed,
        inverted,
        true,
        packed,
    )
    .map(|_| ())
}

fn decode_ld_compressed_record_with_cursor(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    previous_non_ld_packed: &PackedGenotypes,
    inverted: bool,
    has_previous_non_ld: bool,
    packed: &mut PackedGenotypes,
) -> Result<usize> {
    if !has_previous_non_ld {
        return Err(GenoioError::invalid_source(
            path,
            "pgen LD-compressed record appears before any non-LD record",
        ));
    }
    if previous_non_ld_packed.sample_ct() != sample_ct {
        return Err(GenoioError::invalid_source(
            path,
            "pgen LD state length does not match sample count",
        ));
    }
    packed.copy_from(previous_non_ld_packed);
    let mut cursor = 0;
    apply_difflist_entries(path, record, &mut cursor, sample_ct, packed)?;
    if inverted {
        packed.invert_0_2();
    }
    Ok(cursor)
}

pub(super) fn decode_difflist_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    base_category: u8,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    decode_difflist_record_with_cursor(path, record, sample_ct, base_category, packed).map(|_| ())
}

fn apply_difflist_entries(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    visit_difflist_entries(path, record, cursor, sample_ct, |sample_index, category| {
        packed.set(sample_index, category);
    })
}

fn visit_difflist_entries(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    mut visit: impl FnMut(usize, u8),
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
    // PGEN stores first IDs for every group before packed values and deltas.
    // Keep only offsets into the record so decoding does not allocate per row.
    let first_ids_start = *cursor;
    let first_ids_len = group_ct.checked_mul(sample_id_width).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen difflist sample id range is out of range")
    })?;
    ensure_record_bytes(path, record, *cursor, first_ids_len)?;
    *cursor += first_ids_len;
    ensure_record_bytes(path, record, *cursor, group_ct.saturating_sub(1))?;
    *cursor += group_ct.saturating_sub(1);

    let packed_values_start = *cursor;
    ensure_record_bytes(path, record, *cursor, list_len.div_ceil(4))?;
    *cursor += list_len.div_ceil(4);

    let mut previous_sample_id = None;
    let mut entry_index = 0_usize;
    for group_index in 0..group_ct {
        let group_len = (list_len - group_index * 64).min(64);
        let mut first_id_cursor = first_ids_start + group_index * sample_id_width;
        let mut sample_id =
            read_fixed_width_sample_id(path, record, &mut first_id_cursor, sample_id_width)?;
        validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
        visit(
            sample_id,
            packed_difflist_value(record, packed_values_start, entry_index),
        );
        entry_index += 1;
        for _ in 1..group_len {
            let delta = read_base128_varint(path, record, cursor)?;
            sample_id = sample_id.checked_add(delta).ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen difflist sample id is out of range")
            })?;
            validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
            visit(
                sample_id,
                packed_difflist_value(record, packed_values_start, entry_index),
            );
            entry_index += 1;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    fn decoded_categories(packed: &PackedGenotypes) -> Vec<u8> {
        (0..packed.sample_ct())
            .map(|sample_index| packed.get(sample_index))
            .collect()
    }

    fn assert_error_contains(error: genoio_core::GenoioError, expected: &str) {
        let message = error.to_string();
        assert!(
            message.contains(expected),
            "expected error to contain {expected:?}, got {message:?}"
        );
    }

    #[test]
    fn variable_record_helpers_write_packed_genotypes() {
        let path = Path::new("test.pgen");
        let mut packed = PackedGenotypes::default();

        decode_one_bit_record(path, &[2, 0b0000_1010, 0], 4, &mut packed)
            .expect("one-bit record should decode");
        assert_eq!(decoded_categories(&packed), vec![0, 2, 0, 2]);

        decode_difflist_record(path, &[2, 1, 9, 2], 4, 0, &mut packed)
            .expect("difflist record should decode");
        assert_eq!(decoded_categories(&packed), vec![0, 1, 0, 2]);

        let mut previous = PackedGenotypes::default();
        previous.resize(4);
        previous.clear_to(0);
        previous.set(1, 1);
        previous.set(2, 2);
        previous.set(3, 3);

        decode_ld_compressed_record(path, &[1, 2, 0], 4, &previous, true, &mut packed)
            .expect("LD-compressed record should decode");
        assert_eq!(decoded_categories(&packed), vec![2, 1, 2, 3]);
    }

    #[test]
    fn difflist_record_decodes_four_packed_values() {
        let path = Path::new("test.pgen");
        let mut packed = PackedGenotypes::default();

        decode_difflist_record(path, &[4, 0, 0b11_10_01_00, 1, 1, 1], 4, 0, &mut packed)
            .expect("four-entry difflist record should decode");

        assert_eq!(decoded_categories(&packed), vec![0, 1, 2, 3]);
    }

    #[test]
    fn difflist_record_decodes_entries_across_group_boundary() {
        let path = Path::new("test.pgen");
        let mut record = vec![65, 0, 64, 0];
        record.extend(std::iter::repeat_n(0b01_01_01_01, 16));
        record.push(0b01);
        record.extend(std::iter::repeat_n(1, 63));
        let mut packed = PackedGenotypes::default();

        decode_difflist_record(path, &record, 130, 0, &mut packed)
            .expect("multi-group difflist record should decode");

        let categories = decoded_categories(&packed);
        assert!(categories[..65].iter().all(|category| *category == 1));
        assert!(categories[65..].iter().all(|category| *category == 0));
    }

    #[test]
    fn difflist_record_rejects_duplicate_and_out_of_range_sample_ids() {
        let path = Path::new("test.pgen");
        let mut packed = PackedGenotypes::default();

        let duplicate = decode_difflist_record(path, &[2, 1, 0, 0], 4, 0, &mut packed)
            .expect_err("duplicate difflist sample id should fail");
        assert_error_contains(duplicate, "strictly increasing");

        let out_of_range = decode_difflist_record(path, &[1, 4, 0], 4, 0, &mut packed)
            .expect_err("out-of-range difflist sample id should fail");
        assert_error_contains(out_of_range, "outside sample count");
    }

    #[test]
    fn difflist_record_rejects_non_increasing_group_boundary() {
        let path = Path::new("test.pgen");
        let mut record = vec![65, 0, 63, 0];
        record.extend(std::iter::repeat_n(0, 17));
        record.extend(std::iter::repeat_n(1, 63));
        let mut packed = PackedGenotypes::default();

        let error = decode_difflist_record(path, &record, 130, 0, &mut packed)
            .expect_err("non-increasing group boundary should fail");

        assert_error_contains(error, "strictly increasing");
    }

    #[test]
    fn one_bit_record_decodes_with_and_without_difflist_overlay() {
        let path = Path::new("test.pgen");
        let mut packed = PackedGenotypes::default();

        decode_one_bit_record(path, &[2, 0b0000_1010, 0], 4, &mut packed)
            .expect("one-bit record without overlay should decode");
        assert_eq!(decoded_categories(&packed), vec![0, 2, 0, 2]);

        decode_one_bit_record(path, &[2, 0b0000_1010, 1, 2, 1], 4, &mut packed)
            .expect("one-bit record with difflist overlay should decode");
        assert_eq!(decoded_categories(&packed), vec![0, 2, 1, 2]);
    }
}
