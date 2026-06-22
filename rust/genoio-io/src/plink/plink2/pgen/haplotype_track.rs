// pattern: Mixed (unavoidable)
//! PGEN haplotype-track reads and auxiliary phase/dosage decode.
//!
//! The helpers decode phase-present bitmaps, phase values, and phased dosage
//! components, then map selected samples into haplotype-row output.
// Reason: Haplotype reads share the variable-width record buffer and LD state
// with main-track decode, so the IO cursor and auxiliary-track decode stay
// together to avoid copying record payloads.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;
use crate::hardcall::PackedHardcalls as PackedGenotypes;

use super::bitpack::{bit_is_set, ensure_record_bits, ensure_record_bytes};
use super::io::read_fixed_width_phased_dosage_variant_record;
use super::main_track::decode_variable_width_main_track;
use super::{
    PgenDecoderState, PgenHaplotypeDecodeState, PgenHeader, PgenLayout, SelectedSampleCursor,
    PGEN_DOSAGE_SCALE, PGEN_MAX_DOSAGE_RAW, PGEN_MAX_PHASE_RAW, PGEN_MIN_PHASE_RAW,
    PGEN_PHASE_SCALE,
};

pub(in crate::plink::plink2) fn read_plink2_variant_haplotype_main_track(
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

pub(in crate::plink::plink2) fn read_plink2_variant_haplotype_dosage_track(
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

pub(in crate::plink::plink2) fn decode_plink2_haplotype_hardcall_aux(
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

pub(in crate::plink::plink2) fn decode_plink2_haplotype_dosage_aux(
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
                    bit_is_set(record, cursor * 8 + 1 + heterozygote_index)
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
                    let swapped =
                        bit_is_set(record, phaseinfo_start_bit + phased_heterozygote_index);
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
