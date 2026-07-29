// pattern: Functional Core
//! PGEN dosage-track overlay helpers.
//!
//! Dosage tracks store sparse deviations from hard-call-derived values. These
//! helpers apply fixed-width and variable-width overlays only to selected
//! samples while preserving missingness bits.

use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;

use super::bitpack::{
    bit_is_set, ensure_record_bytes, read_base128_varint, read_fixed_width_sample_id,
    sample_id_width, validate_difflist_sample_id,
};
use super::{insert_sorted_unique_index, remove_sorted_index, SelectedSampleCursor};

pub(super) struct DosageOverlayTarget<'a> {
    pub(super) source_indices: &'a [usize],
    pub(super) values: &'a mut [f32],
    pub(super) missing_indices: &'a mut Vec<usize>,
}

pub(super) fn overlay_variable_width_dosages(
    path: &Path,
    record: &[u8],
    mut cursor: usize,
    dosage_bits: u8,
    sample_ct: usize,
    mut target: DosageOverlayTarget<'_>,
    mut dosage_source_indices: Option<&mut Vec<usize>>,
) -> Result<usize> {
    if let Some(indices) = dosage_source_indices.as_deref_mut() {
        indices.clear();
    }
    let mut selected_samples = SelectedSampleCursor::new(target.source_indices);
    match dosage_bits {
        1 => overlay_difflist_dosages(
            path,
            record,
            &mut cursor,
            sample_ct,
            &mut selected_samples,
            &mut target,
            dosage_source_indices.as_deref_mut(),
        )?,
        2 => {
            let dosage_bytes_len = sample_ct.checked_mul(2).ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen dosage byte count is out of range")
            })?;
            ensure_record_bytes(path, record, cursor, dosage_bytes_len)?;
            for sample_index in 0..sample_ct {
                if let Some(indices) = dosage_source_indices.as_deref_mut() {
                    indices.push(sample_index);
                }
                let byte_index = cursor + sample_index * 2;
                let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
                overlay_selected_pgen_dosage(
                    path,
                    sample_index,
                    raw,
                    true,
                    &mut selected_samples,
                    &mut target,
                )?;
            }
            cursor += dosage_bytes_len;
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
                    if let Some(indices) = dosage_source_indices.as_deref_mut() {
                        indices.push(sample_index);
                    }
                    let byte_index = cursor + dosage_index * 2;
                    let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
                    overlay_selected_pgen_dosage(
                        path,
                        sample_index,
                        raw,
                        false,
                        &mut selected_samples,
                        &mut target,
                    )?;
                    dosage_index += 1;
                }
            }
            cursor += dosage_bytes_len;
        }
        other => {
            return Err(GenoioError::invalid_source(
                path,
                format!("unsupported pgen dosage track type {other}"),
            ));
        }
    }
    Ok(cursor)
}

fn overlay_difflist_dosages(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    selected_samples: &mut SelectedSampleCursor<'_>,
    target: &mut DosageOverlayTarget<'_>,
    mut dosage_source_indices: Option<&mut Vec<usize>>,
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

    let layout = read_difflist_id_layout(path, record, cursor, sample_ct, list_len)?;
    let mut values_start = layout.deltas_start;
    walk_difflist_ids(
        path,
        record,
        &mut values_start,
        sample_ct,
        &layout,
        |_, _| Ok(()),
    )?;

    let dosage_bytes_len = list_len.checked_mul(2).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen dosage byte count is out of range")
    })?;
    ensure_record_bytes(path, record, values_start, dosage_bytes_len)?;

    let mut ids_cursor = layout.deltas_start;
    // The values follow the encoded sample IDs, so we first walk IDs to find
    // values_start, then walk them again while overlaying selected samples.
    walk_difflist_ids(
        path,
        record,
        &mut ids_cursor,
        sample_ct,
        &layout,
        |sample_index, dosage_index| {
            if let Some(indices) = dosage_source_indices.as_deref_mut() {
                indices.push(sample_index);
            }
            let byte_index = values_start + dosage_index * 2;
            let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
            overlay_selected_pgen_dosage(path, sample_index, raw, false, selected_samples, target)
        },
    )?;
    *cursor = values_start + dosage_bytes_len;
    Ok(())
}

struct DifflistIdLayout {
    first_ids_start: usize,
    deltas_start: usize,
    group_ct: usize,
    sample_id_width: usize,
    list_len: usize,
}

fn read_difflist_id_layout(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    list_len: usize,
) -> Result<DifflistIdLayout> {
    let group_ct = list_len.div_ceil(64);
    let sample_id_width = sample_id_width(sample_ct);
    let first_ids_start = *cursor;
    let first_ids_len = group_ct.checked_mul(sample_id_width).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen difflist sample id range is out of range")
    })?;
    ensure_record_bytes(path, record, *cursor, first_ids_len)?;
    *cursor += first_ids_len;
    ensure_record_bytes(path, record, *cursor, group_ct.saturating_sub(1))?;
    *cursor += group_ct.saturating_sub(1);

    Ok(DifflistIdLayout {
        first_ids_start,
        deltas_start: *cursor,
        group_ct,
        sample_id_width,
        list_len,
    })
}

fn walk_difflist_ids(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    layout: &DifflistIdLayout,
    mut visit: impl FnMut(usize, usize) -> Result<()>,
) -> Result<()> {
    let mut previous_sample_id = None;
    let mut entry_index = 0;
    for group_index in 0..layout.group_ct {
        let group_len = (layout.list_len - group_index * 64).min(64);
        let mut first_id_cursor = layout.first_ids_start + group_index * layout.sample_id_width;
        let mut sample_id =
            read_fixed_width_sample_id(path, record, &mut first_id_cursor, layout.sample_id_width)?;
        validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
        visit(sample_id, entry_index)?;
        entry_index += 1;
        for _ in 1..group_len {
            let delta = read_base128_varint(path, record, cursor)?;
            sample_id = sample_id.checked_add(delta).ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen difflist sample id is out of range")
            })?;
            validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
            visit(sample_id, entry_index)?;
            entry_index += 1;
        }
    }
    Ok(())
}

fn overlay_selected_pgen_dosage(
    path: &Path,
    source_index: usize,
    raw: u16,
    allow_missing_sentinel: bool,
    selected_samples: &mut SelectedSampleCursor<'_>,
    target: &mut DosageOverlayTarget<'_>,
) -> Result<()> {
    let dosage = decode_pgen_dosage(path, raw, allow_missing_sentinel)?;
    if let Some(selected_index) = selected_samples.selected_index_for(source_index) {
        apply_decoded_pgen_dosage(
            dosage,
            &mut target.values[selected_index],
            target.missing_indices,
            selected_index,
        );
    }
    Ok(())
}

pub(super) fn overlay_fixed_width_dosages(
    path: &Path,
    dosage_bytes: &[u8],
    source_indices: &[usize],
    values: &mut [f32],
    missing_indices: &mut Vec<usize>,
) -> Result<()> {
    missing_indices.clear();
    for (selected_index, source_index) in source_indices.iter().copied().enumerate() {
        let byte_index = source_index.checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source(path, "pgen dosage offset is out of range")
        })?;
        ensure_record_bytes(path, dosage_bytes, byte_index, 2)?;
        let raw = u16::from_le_bytes([dosage_bytes[byte_index], dosage_bytes[byte_index + 1]]);
        let dosage = decode_pgen_dosage(path, raw, true)?;
        apply_decoded_pgen_dosage(
            dosage,
            &mut values[selected_index],
            missing_indices,
            selected_index,
        );
    }
    Ok(())
}

fn decode_pgen_dosage(path: &Path, raw: u16, allow_missing_sentinel: bool) -> Result<Option<f32>> {
    if raw == u16::MAX && allow_missing_sentinel {
        return Ok(None);
    }
    if raw > super::PGEN_MAX_DOSAGE_RAW {
        let message = if raw == u16::MAX {
            "pgen sparse dosage entry uses the full-track missing sentinel".to_owned()
        } else {
            format!("pgen dosage raw value {raw} is reserved; expected 0..=32768")
        };
        return Err(GenoioError::invalid_source(path, message));
    }
    Ok(Some(f32::from(raw) * super::PGEN_DOSAGE_SCALE))
}

fn apply_decoded_pgen_dosage(
    dosage: Option<f32>,
    value: &mut f32,
    missing_indices: &mut Vec<usize>,
    selected_index: usize,
) {
    let Some(dosage) = dosage else {
        *value = 0.0;
        insert_sorted_unique_index(missing_indices, selected_index);
        return;
    };
    *value = dosage;
    remove_sorted_index(missing_indices, selected_index);
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    fn variable_record(dosage_bits: u8, raw: u16) -> Vec<u8> {
        let mut record = match dosage_bits {
            1 => vec![1, 0],
            2 => Vec::new(),
            3 => vec![1],
            other => panic!("unsupported test dosage type {other}"),
        };
        record.extend(raw.to_le_bytes());
        record
    }

    fn decode_variable_raw(dosage_bits: u8, raw: u16) -> Result<(f32, bool)> {
        let record = variable_record(dosage_bits, raw);
        let mut values = vec![0.0];
        let mut missing_indices = Vec::new();
        overlay_variable_width_dosages(
            Path::new("test.pgen"),
            &record,
            0,
            dosage_bits,
            1,
            DosageOverlayTarget {
                source_indices: &[0],
                values: &mut values,
                missing_indices: &mut missing_indices,
            },
            None,
        )?;
        Ok((values[0], !missing_indices.is_empty()))
    }

    #[test]
    fn variable_width_dosage_overlay_preserves_hardcall_inferred_values() {
        let path = Path::new("test.pgen");
        let mut record = vec![0b0000_0110];
        record.extend(100_u16.to_le_bytes());
        record.extend(200_u16.to_le_bytes());
        let mut values = vec![0.0, 2.0, 0.0];
        let mut missing_indices = vec![2];

        overlay_variable_width_dosages(
            path,
            &record,
            0,
            3,
            4,
            DosageOverlayTarget {
                source_indices: &[0, 2, 3],
                values: &mut values,
                missing_indices: &mut missing_indices,
            },
            None,
        )
        .expect("dosage overlay should decode");

        assert_eq!(values, vec![0.0, f32::from(200_u16) * (2.0 / 32768.0), 0.0]);
        assert_eq!(missing_indices, vec![2]);
    }

    #[test]
    fn variable_width_dosage_list_overlay_uses_source_order_without_dense_index() {
        let path = Path::new("test.pgen");
        let mut record = vec![3, 1, 3, 5];
        record.extend(100_u16.to_le_bytes());
        record.extend(200_u16.to_le_bytes());
        record.extend(300_u16.to_le_bytes());
        let mut values = vec![0.0, 1.0];
        let mut missing_indices = Vec::new();

        overlay_variable_width_dosages(
            path,
            &record,
            0,
            1,
            10,
            DosageOverlayTarget {
                source_indices: &[4, 9],
                values: &mut values,
                missing_indices: &mut missing_indices,
            },
            None,
        )
        .expect("dosage-list overlay should decode");

        assert_eq!(
            values,
            vec![
                f32::from(200_u16) * (2.0 / 32768.0),
                f32::from(300_u16) * (2.0 / 32768.0),
            ]
        );
        assert!(missing_indices.is_empty());
    }

    #[test]
    fn variable_width_dosage_list_overlay_crosses_group_boundary() {
        let path = Path::new("test.pgen");
        let mut record = vec![65, 0, 64, 0];
        record.extend(std::iter::repeat_n(1, 63));
        for value in 1_u16..=65 {
            record.extend(value.to_le_bytes());
        }
        let mut values = vec![0.0, 0.0, 9.0];
        let mut missing_indices = Vec::new();

        overlay_variable_width_dosages(
            path,
            &record,
            0,
            1,
            130,
            DosageOverlayTarget {
                source_indices: &[0, 64, 129],
                values: &mut values,
                missing_indices: &mut missing_indices,
            },
            None,
        )
        .expect("multi-group dosage-list overlay should decode");

        assert_eq!(
            values,
            vec![1.0 * (2.0 / 32768.0), 65.0 * (2.0 / 32768.0), 9.0,]
        );
        assert!(missing_indices.is_empty());
    }

    #[test]
    fn pbr_rust_plink2_003_raw_dosage_domain_is_representation_aware() {
        for dosage_bits in [1, 2, 3] {
            let (value, missing) =
                decode_variable_raw(dosage_bits, 32_768).expect("32768 must be valid");
            assert_eq!(value, 2.0);
            assert!(!missing);
        }

        let (_, missing) =
            decode_variable_raw(2, u16::MAX).expect("full-track sentinel must be valid");
        assert!(missing);

        for (dosage_bits, raw) in [
            (1, 32_769),
            (1, 65_534),
            (1, 65_535),
            (2, 32_769),
            (2, 65_534),
            (3, 32_769),
            (3, 65_534),
            (3, 65_535),
        ] {
            let error = decode_variable_raw(dosage_bits, raw)
                .expect_err("reserved or context-invalid raw dosage must fail");
            assert!(
                matches!(error, GenoioError::InvalidSource { .. }),
                "type {dosage_bits}, raw {raw}: {error}"
            );
        }

        for raw in [32_769_u16, 65_534_u16] {
            let mut values = vec![0.0];
            let mut missing_indices = Vec::new();
            let error = overlay_fixed_width_dosages(
                Path::new("test.pgen"),
                &raw.to_le_bytes(),
                &[0],
                &mut values,
                &mut missing_indices,
            )
            .expect_err("reserved fixed-width raw dosage must fail");
            assert!(matches!(error, GenoioError::InvalidSource { .. }));
        }
        let mut values = vec![0.0];
        let mut missing_indices = Vec::new();
        overlay_fixed_width_dosages(
            Path::new("test.pgen"),
            &u16::MAX.to_le_bytes(),
            &[0],
            &mut values,
            &mut missing_indices,
        )
        .expect("fixed-width full-track sentinel must be accepted");
        assert_eq!(missing_indices, vec![0]);
    }
}
