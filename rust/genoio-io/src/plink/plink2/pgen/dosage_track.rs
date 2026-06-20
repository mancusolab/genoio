// pattern: Functional Core
//! PGEN dosage-track overlay helpers.

use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;

use super::bitpack::{
    bit_is_set, ensure_record_bytes, read_base128_varint, read_fixed_width_sample_id,
    sample_id_width, validate_difflist_sample_id,
};
use super::SelectedSampleCursor;

pub(super) struct DosageOverlayTarget<'a> {
    pub(super) source_indices: &'a [usize],
    pub(super) values: &'a mut [f32],
    pub(super) missing: &'a mut [bool],
}

pub(super) fn overlay_variable_width_dosages(
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

pub(super) fn overlay_fixed_width_dosages(
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

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

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
