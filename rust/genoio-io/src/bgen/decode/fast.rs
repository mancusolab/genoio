// pattern: Functional Core

use std::path::Path;
use std::sync::OnceLock;

use genoio_core::GenoioError;

use crate::dosage_filter::DosageFilterCounts;
use crate::Result;

use super::generic::{
    decode_phased_a1_dosage, PackedProbabilityLayout, BIALLELIC_DIPLOID_STORED_PROBABILITY_COUNT,
};

pub(super) const UNPHASED_8_BIT_DEPTH: u8 = 8;
pub(super) const PHASED_16_BIT_DEPTH: u8 = 16;
const PHASED_16_BYTES_PER_SAMPLE: usize = 4;

/// Mutable column target for BGEN paths that can decode directly into the
/// caller's preallocated sample-major matrix.
pub(in crate::bgen) struct SampleMajorSlotMut<'a> {
    pub(in crate::bgen) values: &'a mut [f32],
    pub(in crate::bgen) row_width: usize,
    pub(in crate::bgen) variant_index: usize,
}

pub(super) fn decode_selected_unphased_8bit_a1_dosages(
    path: &Path,
    ploidies: &[u8],
    packed_probabilities: &[u8],
    packed_probability_layout: PackedProbabilityLayout,
    source_indices: &[usize],
    values: &mut Vec<f32>,
    missing_indices: &mut Vec<usize>,
) -> Result<()> {
    values.clear();
    missing_indices.clear();
    values.reserve(source_indices.len());
    debug_assert!(source_indices
        .windows(2)
        .all(|window| window[0] < window[1]));

    let mut selected_cursor = 0_usize;
    let mut probability_cursor = 0_usize;
    for (sample_index, &ploidy_byte) in ploidies.iter().enumerate() {
        if selected_cursor == source_indices.len() {
            break;
        }
        let is_selected = source_indices[selected_cursor] == sample_index;
        let is_missing = ploidy_byte & 0b1000_0000 != 0;
        if is_missing {
            if packed_probability_layout == PackedProbabilityLayout::AllSamples {
                probability_cursor = skip_missing_unphased_8bit_sample_probabilities(
                    path,
                    packed_probabilities,
                    probability_cursor,
                )?;
            }
            if is_selected {
                missing_indices.push(values.len());
                values.push(0.0);
                selected_cursor += 1;
            }
            continue;
        }

        let Some((&p_aa, &p_ab)) = packed_probabilities
            .get(probability_cursor)
            .zip(packed_probabilities.get(probability_cursor + 1))
        else {
            return Err(GenoioError::invalid_source(
                path,
                "bgen packed probability bytes are truncated",
            ));
        };
        probability_cursor += 2;

        if is_selected {
            values.push(unphased_8bit_a1_dosage(path, p_aa, p_ab)?);
            selected_cursor += 1;
        }
    }
    debug_assert_eq!(selected_cursor, source_indices.len());
    Ok(())
}

pub(super) fn decode_selected_called_unphased_8bit_a1_dosages(
    path: &Path,
    packed_probabilities: &[u8],
    source_indices: &[usize],
    values: &mut Vec<f32>,
    missing_indices: &mut Vec<usize>,
) -> Result<()> {
    values.clear();
    missing_indices.clear();
    values.reserve(source_indices.len());
    let lut = unphased_8bit_a1_dosage_lut();

    for &source_index in source_indices {
        let probability_cursor = source_index.checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen sample probability offset is out of range")
        })?;
        let Some((&p_aa, &p_ab)) = packed_probabilities
            .get(probability_cursor)
            .zip(packed_probabilities.get(probability_cursor + 1))
        else {
            return Err(GenoioError::invalid_source(
                path,
                "bgen packed probability bytes are truncated",
            ));
        };
        values.push(unphased_8bit_a1_dosage_from_lut(path, p_aa, p_ab, lut)?);
    }

    Ok(())
}

pub(super) fn decode_selected_called_unphased_8bit_a1_dosages_into_sample_major_slot(
    path: &Path,
    packed_probabilities: &[u8],
    source_indices: &[usize],
    slot: &mut SampleMajorSlotMut<'_>,
) -> Result<()> {
    validate_sample_major_slot_shape(source_indices.len(), slot)?;
    let lut = unphased_8bit_a1_dosage_lut();

    for (selected_cursor, &source_index) in source_indices.iter().enumerate() {
        let probability_cursor = source_index.checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen sample probability offset is out of range")
        })?;
        let Some((&p_aa, &p_ab)) = packed_probabilities
            .get(probability_cursor)
            .zip(packed_probabilities.get(probability_cursor + 1))
        else {
            return Err(GenoioError::invalid_source(
                path,
                "bgen packed probability bytes are truncated",
            ));
        };
        let target_index = selected_cursor * slot.row_width + slot.variant_index;
        slot.values[target_index] = unphased_8bit_a1_dosage_from_lut(path, p_aa, p_ab, lut)?;
    }

    Ok(())
}

pub(super) fn decode_selected_called_phased_16bit_a1_dosages(
    path: &Path,
    packed_probabilities: &[u8],
    source_indices: &[usize],
    values: &mut Vec<f32>,
    missing_indices: &mut Vec<usize>,
) -> Result<()> {
    values.clear();
    missing_indices.clear();
    values.reserve(source_indices.len());

    for &source_index in source_indices {
        let (first, second) = phased_16bit_raw_pair(path, packed_probabilities, source_index)?;
        values.push(decode_phased_16bit_a1_dosage(first, second));
    }

    Ok(())
}

pub(super) fn decode_selected_called_phased_16bit_a1_dosages_with_counts(
    path: &Path,
    packed_probabilities: &[u8],
    source_indices: &[usize],
    values: &mut Vec<f32>,
    missing_indices: &mut Vec<usize>,
) -> Result<DosageFilterCounts> {
    values.clear();
    missing_indices.clear();
    values.reserve(source_indices.len());
    let mut counts = DosageFilterCounts::default();

    for &source_index in source_indices {
        let (first, second) = phased_16bit_raw_pair(path, packed_probabilities, source_index)?;
        let dosage = decode_phased_16bit_a1_dosage(first, second);
        counts.record_called_dosage(dosage);
        values.push(dosage);
    }

    Ok(counts)
}

pub(super) fn decode_selected_called_phased_16bit_a1_dosages_into_sample_major_slot(
    path: &Path,
    packed_probabilities: &[u8],
    source_indices: &[usize],
    slot: &mut SampleMajorSlotMut<'_>,
) -> Result<()> {
    validate_sample_major_slot_shape(source_indices.len(), slot)?;

    for (selected_cursor, &source_index) in source_indices.iter().enumerate() {
        let (first, second) = phased_16bit_raw_pair(path, packed_probabilities, source_index)?;
        let target_index = selected_cursor * slot.row_width + slot.variant_index;
        slot.values[target_index] = decode_phased_16bit_a1_dosage(first, second);
    }

    Ok(())
}

pub(super) fn decode_phased_16bit_a1_dosage(first: u16, second: u16) -> f32 {
    decode_phased_a1_dosage(PHASED_16_BIT_DEPTH, u32::from(first), u32::from(second))
}

fn skip_missing_unphased_8bit_sample_probabilities(
    path: &Path,
    packed_probabilities: &[u8],
    probability_cursor: usize,
) -> Result<usize> {
    let probability_cursor = probability_cursor
        .checked_add(BIALLELIC_DIPLOID_STORED_PROBABILITY_COUNT)
        .ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen sample probability offset is out of range")
        })?;
    if probability_cursor > packed_probabilities.len() {
        return Err(GenoioError::invalid_source(
            path,
            "bgen packed probability bytes are truncated",
        ));
    }
    Ok(probability_cursor)
}

fn phased_16bit_raw_pair(
    path: &Path,
    packed_probabilities: &[u8],
    source_index: usize,
) -> Result<(u16, u16)> {
    let probability_cursor = source_index
        .checked_mul(PHASED_16_BYTES_PER_SAMPLE)
        .ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen sample probability offset is out of range")
        })?;
    let probability_end = probability_cursor
        .checked_add(PHASED_16_BYTES_PER_SAMPLE)
        .ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen sample probability offset is out of range")
        })?;
    let bytes = packed_probabilities
        .get(probability_cursor..probability_end)
        .ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen packed probability bytes are truncated")
        })?;
    Ok((
        u16::from_le_bytes([bytes[0], bytes[1]]),
        u16::from_le_bytes([bytes[2], bytes[3]]),
    ))
}

fn validate_sample_major_slot_shape(
    selected_sample_count: usize,
    slot: &SampleMajorSlotMut<'_>,
) -> Result<()> {
    let expected_len = selected_sample_count
        .checked_mul(slot.row_width)
        .ok_or_else(|| {
            GenoioError::internal_contract("sample-major dense matrix shape is out of range")
        })?;
    if slot.values.len() != expected_len {
        return Err(GenoioError::internal_contract(
            "sample-major dense buffer does not match declared shape",
        ));
    }
    if slot.variant_index >= slot.row_width {
        return Err(GenoioError::internal_contract(
            "sample-major variant index is outside row width",
        ));
    }
    Ok(())
}

fn unphased_8bit_a1_dosage(path: &Path, p_aa: u8, p_ab: u8) -> Result<f32> {
    unphased_8bit_a1_dosage_from_lut(path, p_aa, p_ab, unphased_8bit_a1_dosage_lut())
}

fn unphased_8bit_a1_dosage_from_lut(
    path: &Path,
    p_aa: u8,
    p_ab: u8,
    lut: &[f32; 65_536],
) -> Result<f32> {
    if u16::from(p_aa) + u16::from(p_ab) > 255 {
        return Err(GenoioError::invalid_source(
            path,
            "bgen malformed probability values produce invalid a1 dosage",
        ));
    }
    Ok(lut[usize::from(p_aa) << 8 | usize::from(p_ab)])
}

fn unphased_8bit_a1_dosage_lut() -> &'static [f32; 65_536] {
    static LUT: OnceLock<[f32; 65_536]> = OnceLock::new();
    LUT.get_or_init(|| {
        let mut values = [0.0_f32; 65_536];
        for p_aa in 0..=255_u16 {
            for p_ab in 0..=255_u16 {
                // Keep the same f32 operation order as the generic decoder so
                // existing exact-parity tests remain stable.
                let p_aa_f = f32::from(p_aa) / 255.0;
                let p_ab_f = f32::from(p_ab) / 255.0;
                let p_bb = 1.0 - p_aa_f - p_ab_f;
                values[usize::from(p_aa) << 8 | usize::from(p_ab)] =
                    (p_ab_f + 2.0 * p_bb).clamp(0.0, 2.0);
            }
        }
        values
    })
}
