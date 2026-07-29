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
use super::dosage_track::{
    overlay_variable_width_dosages, DosageEntryScratch, DosageOverlayTarget,
};
use super::io::read_fixed_width_phased_dosage_variant_record;
use super::main_track::decode_variable_width_main_track;
use super::{
    insert_sorted_unique_index, PgenDecoderState, PgenHaplotypeDecodeState, PgenHeader, PgenLayout,
    SelectedSampleCursor, ValidatedPhasedDosage, PGEN_DOSAGE_SCALE, PGEN_MAX_DOSAGE_RAW,
    PGEN_MAX_PHASE_RAW, PGEN_MIN_PHASE_RAW, PGEN_PHASE_SCALE,
};

#[derive(Clone, Copy)]
struct HardcallPhaseLayout {
    phasepresent_start_bit: Option<usize>,
    phaseinfo_start_bit: usize,
    end_cursor: usize,
}

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
    if dosage_bits == 0 {
        return Err(GenoioError::unsupported(
            "pgen record does not contain dosage values",
        ));
    }

    let record = decoder_state.record.as_slice();
    let dosage_cursor = if record_type & 0x10 != 0 {
        decode_selected_hardcall_phase_prefix(
            path,
            record,
            cursor,
            source_indices,
            &decoder_state.packed,
            false,
            haplotype_state,
        )?
    } else {
        initialize_selected_hardcall_haplotypes(
            source_indices,
            &decoder_state.packed,
            None,
            record,
            haplotype_state,
        )?;
        cursor
    };
    let phase_cursor = overlay_variable_width_dosages(
        path,
        record,
        dosage_cursor,
        dosage_bits,
        header.sample_ct,
        DosageOverlayTarget {
            source_indices,
            values: &mut haplotype_state.selected_collapsed_values,
            missing_indices: &mut haplotype_state.selected_collapsed_missing_indices,
        },
        Some(DosageEntryScratch {
            source_indices: &mut haplotype_state.dosage_source_indices,
            totals: &mut haplotype_state.dosage_source_totals,
        }),
    )?;
    haplotype_state.selected_explicit_phase.clear();
    haplotype_state
        .selected_explicit_phase
        .resize(source_indices.len(), false);

    validate_variable_phased_dosage_record(
        path,
        record,
        record_type,
        phase_cursor,
        dosage_bits,
        (
            &haplotype_state.dosage_source_indices,
            &haplotype_state.dosage_source_totals,
        ),
        Some(&mut haplotype_state.validated_phased_dosages),
    )?;
    let mut selected_cursor = SelectedSampleCursor::new(source_indices);
    for dosage_index in 0..haplotype_state.dosage_source_indices.len() {
        let source_index = haplotype_state.dosage_source_indices[dosage_index];
        let Some(phased) = haplotype_state.validated_phased_dosages[dosage_index] else {
            continue;
        };
        if let Some(selected_index) = selected_cursor.selected_index_for(source_index) {
            apply_validated_explicit_phase(selected_index, phased, haplotype_state);
        }
    }
    finalize_implicit_haplotype_dosages(
        path,
        source_indices,
        &decoder_state.packed,
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
    let end_cursor = decode_selected_hardcall_phase_prefix(
        path,
        record,
        cursor,
        source_indices,
        packed,
        true,
        haplotype_state,
    )?;
    if end_cursor != record.len() {
        return Err(GenoioError::invalid_source(
            path,
            "pgen phased hardcall record has trailing or missing bytes",
        ));
    }
    Ok(())
}

pub(in crate::plink::plink2) fn skip_hardcall_phase_track(
    path: &Path,
    record: &[u8],
    cursor: usize,
    packed: &PackedGenotypes,
) -> Result<usize> {
    Ok(hardcall_phase_layout(path, record, cursor, packed)?.end_cursor)
}

fn hardcall_phase_layout(
    path: &Path,
    record: &[u8],
    cursor: usize,
    packed: &PackedGenotypes,
) -> Result<HardcallPhaseLayout> {
    let heterozygote_ct = (0..packed.sample_ct())
        .filter(|sample_index| packed.get(*sample_index) == 1)
        .count();
    ensure_record_bytes(path, record, cursor, 1)?;
    let phasepresent_stored = bit_is_set(&record[cursor..], 0);
    let (phasepresent_start_bit, phaseinfo_start_bit, phased_heterozygote_ct) =
        if phasepresent_stored {
            let phasepresent_bits = 1 + heterozygote_ct;
            let phasepresent_start_bit = cursor
                .checked_mul(8)
                .and_then(|bit| bit.checked_add(1))
                .ok_or_else(|| {
                GenoioError::invalid_source(path, "pgen phase offset is out of range")
            })?;
            ensure_record_bits(path, record, phasepresent_start_bit, heterozygote_ct)?;
            let phased_heterozygote_ct = (0..heterozygote_ct)
                .filter(|index| bit_is_set(record, phasepresent_start_bit + index))
                .count();
            let phaseinfo_start_bit = cursor
                .checked_add(phasepresent_bits.div_ceil(8))
                .and_then(|byte| byte.checked_mul(8))
                .ok_or_else(|| {
                    GenoioError::invalid_source(path, "pgen phase offset is out of range")
                })?;
            (
                Some(phasepresent_start_bit),
                phaseinfo_start_bit,
                phased_heterozygote_ct,
            )
        } else {
            let phaseinfo_start_bit = cursor
                .checked_mul(8)
                .and_then(|bit| bit.checked_add(1))
                .ok_or_else(|| {
                    GenoioError::invalid_source(path, "pgen phase offset is out of range")
                })?;
            (None, phaseinfo_start_bit, heterozygote_ct)
        };
    ensure_record_bits(path, record, phaseinfo_start_bit, phased_heterozygote_ct)?;
    let end_cursor = phaseinfo_start_bit
        .checked_add(phased_heterozygote_ct)
        .ok_or_else(|| GenoioError::invalid_source(path, "pgen phase offset is out of range"))?
        .div_ceil(8);
    Ok(HardcallPhaseLayout {
        phasepresent_start_bit,
        phaseinfo_start_bit,
        end_cursor,
    })
}

fn decode_selected_hardcall_phase_prefix(
    path: &Path,
    record: &[u8],
    cursor: usize,
    source_indices: &[usize],
    packed: &PackedGenotypes,
    require_selected_heterozygotes_phased: bool,
    haplotype_state: &mut PgenHaplotypeDecodeState,
) -> Result<usize> {
    let layout = hardcall_phase_layout(path, record, cursor, packed)?;
    initialize_selected_hardcall_haplotypes(
        source_indices,
        packed,
        Some(layout),
        record,
        haplotype_state,
    )?;
    if require_selected_heterozygotes_phased {
        for (selected_index, source_index) in source_indices.iter().copied().enumerate() {
            if packed.get(source_index) == 1
                && haplotype_state.selected_phase_swapped[selected_index].is_none()
            {
                return Err(GenoioError::unsupported(
                    "unphased pgen heterozygous hardcall retained in haplotype read",
                ));
            }
        }
    }
    Ok(layout.end_cursor)
}

fn initialize_selected_hardcall_haplotypes(
    source_indices: &[usize],
    packed: &PackedGenotypes,
    phase_layout: Option<HardcallPhaseLayout>,
    record: &[u8],
    haplotype_state: &mut PgenHaplotypeDecodeState,
) -> Result<()> {
    haplotype_state.selected_haplotype_values.clear();
    haplotype_state.selected_haplotype_missing_indices.clear();
    haplotype_state.selected_collapsed_values.clear();
    haplotype_state.selected_collapsed_missing_indices.clear();
    haplotype_state.selected_phase_swapped.clear();
    haplotype_state
        .selected_haplotype_values
        .resize(source_indices.len() * 2, 0.0);
    haplotype_state
        .selected_collapsed_values
        .resize(source_indices.len(), 0.0);
    haplotype_state
        .selected_phase_swapped
        .resize(source_indices.len(), None);

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
                    haplotype_state.selected_collapsed_values[selected_index] = 0.0;
                }
            }
            1 => {
                let phase_present = phase_layout.is_some_and(|layout| {
                    layout
                        .phasepresent_start_bit
                        .is_none_or(|start_bit| bit_is_set(record, start_bit + heterozygote_index))
                });
                heterozygote_index += 1;
                if phase_present {
                    let phaseinfo_start_bit = phase_layout
                        .ok_or_else(|| {
                            GenoioError::internal_contract(
                                "phased hardcall is missing phase-track layout",
                            )
                        })?
                        .phaseinfo_start_bit;
                    let swapped =
                        bit_is_set(record, phaseinfo_start_bit + phased_heterozygote_index);
                    phased_heterozygote_index += 1;
                    if let Some(selected_index) = selected_index {
                        haplotype_state.selected_phase_swapped[selected_index] = Some(swapped);
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
                if let Some(selected_index) = selected_index {
                    haplotype_state.selected_collapsed_values[selected_index] = 1.0;
                }
            }
            2 => {
                if let Some(selected_index) = selected_index {
                    set_selected_haplotype_pair(haplotype_state, selected_index, 1.0, 1.0, false);
                    haplotype_state.selected_collapsed_values[selected_index] = 2.0;
                }
            }
            3 => {
                if let Some(selected_index) = selected_index {
                    set_selected_haplotype_pair(haplotype_state, selected_index, 0.0, 0.0, true);
                    haplotype_state.selected_collapsed_values[selected_index] = 0.0;
                    insert_sorted_unique_index(
                        &mut haplotype_state.selected_collapsed_missing_indices,
                        selected_index,
                    );
                }
            }
            _ => unreachable!("two-bit hard-call code should be masked"),
        }
    }
    Ok(())
}

pub(in crate::plink::plink2) fn validate_variable_phased_dosage_record(
    path: &Path,
    record: &[u8],
    record_type: u8,
    phase_cursor: usize,
    dosage_bits: u8,
    dosage_entries: (&[usize], &[Option<f32>]),
    mut validated_phased_dosages: Option<&mut Vec<Option<ValidatedPhasedDosage>>>,
) -> Result<()> {
    let (dosage_source_indices, _) = dosage_entries;
    if let Some(validated) = validated_phased_dosages.as_deref_mut() {
        validated.clear();
        validated.resize(dosage_source_indices.len(), None);
    }
    let record_end = if record_type & 0x80 != 0 {
        validate_variable_explicit_phase_track(
            path,
            record,
            phase_cursor,
            dosage_bits,
            dosage_entries,
            validated_phased_dosages,
        )?
    } else {
        phase_cursor
    };
    if record_end != record.len() {
        return Err(GenoioError::invalid_source(
            path,
            "pgen phased dosage record has trailing or missing bytes",
        ));
    }
    Ok(())
}

fn validate_variable_explicit_phase_track(
    path: &Path,
    record: &[u8],
    cursor: usize,
    dosage_bits: u8,
    dosage_entries: (&[usize], &[Option<f32>]),
    mut validated_phased_dosages: Option<&mut Vec<Option<ValidatedPhasedDosage>>>,
) -> Result<usize> {
    let (dosage_source_indices, dosage_source_totals) = dosage_entries;
    if dosage_source_indices.len() != dosage_source_totals.len() {
        return Err(GenoioError::internal_contract(
            "pgen phased dosage totals and source indices are misaligned",
        ));
    }
    if dosage_bits == 2 {
        let phase_bytes_len = dosage_source_indices.len().checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source(path, "pgen phased dosage byte count is out of range")
        })?;
        ensure_record_bytes(path, record, cursor, phase_bytes_len)?;
        for dosage_index in 0..dosage_source_indices.len() {
            let byte_index = cursor + dosage_index * 2;
            let delta = decode_variable_phase_delta(
                path,
                i16::from_le_bytes([record[byte_index], record[byte_index + 1]]),
                true,
            )?;
            let total = dosage_source_totals[dosage_index];
            let phased = validate_phased_dosage_pair(path, total, delta)?;
            if let Some(validated) = validated_phased_dosages.as_deref_mut() {
                validated[dosage_index] = Some(phased);
            }
        }
        return Ok(cursor + phase_bytes_len);
    }

    let dosage_ct = dosage_source_indices.len();
    let phase_exists_len = dosage_ct.div_ceil(8);
    ensure_record_bytes(path, record, cursor, phase_exists_len)?;
    let phase_exists = &record[cursor..cursor + phase_exists_len];
    let phase_value_ct = (0..dosage_ct)
        .filter(|dosage_index| bit_is_set(phase_exists, *dosage_index))
        .count();
    let phase_values_cursor = cursor + phase_exists_len;
    let phase_bytes_len = phase_value_ct.checked_mul(2).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen phased dosage byte count is out of range")
    })?;
    ensure_record_bytes(path, record, phase_values_cursor, phase_bytes_len)?;

    let mut phase_value_index = 0_usize;
    for dosage_index in 0..dosage_ct {
        if !bit_is_set(phase_exists, dosage_index) {
            continue;
        }
        let byte_index = phase_values_cursor + phase_value_index * 2;
        let delta = decode_variable_phase_delta(
            path,
            i16::from_le_bytes([record[byte_index], record[byte_index + 1]]),
            false,
        )?;
        let total = dosage_source_totals[dosage_index];
        let phased = validate_phased_dosage_pair(path, total, delta)?;
        if let Some(validated) = validated_phased_dosages.as_deref_mut() {
            validated[dosage_index] = Some(phased);
        }
        phase_value_index += 1;
    }
    Ok(phase_values_cursor + phase_bytes_len)
}

fn decode_variable_phase_delta(path: &Path, raw: i16, allow_missing: bool) -> Result<Option<f32>> {
    if raw == i16::MIN {
        if allow_missing {
            return Ok(None);
        }
        return Err(GenoioError::invalid_source(
            path,
            "pgen sparse phased-dosage entry uses the full-track missing sentinel",
        ));
    }
    if !(PGEN_MIN_PHASE_RAW..=PGEN_MAX_PHASE_RAW).contains(&raw) {
        return Err(GenoioError::invalid_source(
            path,
            format!("pgen phased dosage phase raw value {raw} is outside [-16384, 16384]"),
        ));
    }
    Ok(Some(f32::from(raw) * PGEN_PHASE_SCALE))
}

fn validate_phased_dosage_pair(
    path: &Path,
    total: Option<f32>,
    delta: Option<f32>,
) -> Result<ValidatedPhasedDosage> {
    match (total, delta) {
        (None, None) => Ok(ValidatedPhasedDosage::Missing),
        (Some(total), Some(delta)) => {
            let left = (total + delta) * 0.5;
            let right = (total - delta) * 0.5;
            validate_phased_dosage_haplotype_components(path, left, right)?;
            Ok(ValidatedPhasedDosage::Present { total, left, right })
        }
        _ => Err(GenoioError::invalid_source(
            path,
            "pgen full dosage and phased-dosage missing sentinels are inconsistent",
        )),
    }
}

fn apply_validated_explicit_phase(
    selected_index: usize,
    phased: ValidatedPhasedDosage,
    haplotype_state: &mut PgenHaplotypeDecodeState,
) {
    match phased {
        ValidatedPhasedDosage::Missing => {
            set_selected_haplotype_pair(haplotype_state, selected_index, 0.0, 0.0, true);
        }
        ValidatedPhasedDosage::Present { left, right, .. } => {
            set_selected_haplotype_pair(haplotype_state, selected_index, left, right, false);
        }
    }
    haplotype_state.selected_explicit_phase[selected_index] = true;
}

fn finalize_implicit_haplotype_dosages(
    path: &Path,
    source_indices: &[usize],
    packed: &PackedGenotypes,
    haplotype_state: &mut PgenHaplotypeDecodeState,
) -> Result<()> {
    for (selected_index, source_index) in source_indices.iter().copied().enumerate() {
        if haplotype_state.selected_explicit_phase[selected_index] {
            continue;
        }
        if haplotype_state
            .selected_collapsed_missing_indices
            .binary_search(&selected_index)
            .is_ok()
        {
            set_selected_haplotype_pair(haplotype_state, selected_index, 0.0, 0.0, true);
            continue;
        }
        let has_stored_dosage = haplotype_state
            .dosage_source_indices
            .binary_search(&source_index)
            .is_ok();
        if has_stored_dosage {
            let Some(swapped) = haplotype_state.selected_phase_swapped[selected_index] else {
                return Err(GenoioError::unsupported(
                    "pgen phased dosage is unavailable: the dosage lacks explicit phase and cannot be implicitly phased from a phased heterozygous hardcall",
                ));
            };
            let total = haplotype_state.selected_collapsed_values[selected_index];
            let alt_haplotype = total.min(1.0);
            let ref_haplotype = (total - 1.0).max(0.0);
            let (left, right) = if swapped {
                (alt_haplotype, ref_haplotype)
            } else {
                (ref_haplotype, alt_haplotype)
            };
            validate_phased_dosage_haplotype_components(path, left, right)?;
            set_selected_haplotype_pair(haplotype_state, selected_index, left, right, false);
        } else if packed.get(source_index) == 1
            && haplotype_state.selected_phase_swapped[selected_index].is_none()
        {
            return Err(GenoioError::unsupported(
                "pgen phased dosage is unavailable for an unphased heterozygous hardcall",
            ));
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
    haplotype_state.dosage_source_totals.clear();
    haplotype_state.dosage_source_totals.reserve(sample_ct);
    for source_index in 0..sample_ct {
        let dosage_offset = cursor + source_index * 2;
        haplotype_state
            .dosage_source_totals
            .push(decode_phased_dosage_total(
                path,
                u16::from_le_bytes([record[dosage_offset], record[dosage_offset + 1]]),
            )?);
    }
    validate_full_phased_dosage_track(
        path,
        record,
        phase_cursor,
        sample_ct,
        &haplotype_state.dosage_source_totals,
        Some(&mut haplotype_state.validated_phased_dosages),
    )?;

    haplotype_state.selected_haplotype_values.clear();
    haplotype_state.selected_haplotype_missing_indices.clear();
    haplotype_state.selected_collapsed_values.clear();
    haplotype_state.selected_collapsed_missing_indices.clear();
    let mut selected_samples = SelectedSampleCursor::new(source_indices);
    for source_index in 0..sample_ct {
        let phased = haplotype_state.validated_phased_dosages[source_index].ok_or_else(|| {
            GenoioError::internal_contract("full pgen phased dosage validation is incomplete")
        })?;
        if selected_samples.selected_index_for(source_index).is_none() {
            continue;
        }
        match phased {
            ValidatedPhasedDosage::Missing => {
                let row_offset = haplotype_state.selected_haplotype_values.len();
                haplotype_state
                    .selected_haplotype_missing_indices
                    .extend([row_offset, row_offset + 1]);
                haplotype_state.selected_haplotype_values.extend([0.0, 0.0]);
                push_collapsed_dosage(haplotype_state, 0.0, true);
            }
            ValidatedPhasedDosage::Present { total, left, right } => {
                haplotype_state
                    .selected_haplotype_values
                    .extend([left, right]);
                push_collapsed_dosage(haplotype_state, total, false);
            }
        }
    }
    Ok(())
}

pub(in crate::plink::plink2) fn validate_full_phased_dosage_track(
    path: &Path,
    record: &[u8],
    phase_cursor: usize,
    sample_ct: usize,
    dosage_totals: &[Option<f32>],
    mut validated_phased_dosages: Option<&mut Vec<Option<ValidatedPhasedDosage>>>,
) -> Result<()> {
    if dosage_totals.len() != sample_ct {
        return Err(GenoioError::internal_contract(
            "pgen full phased dosage totals do not match the sample count",
        ));
    }
    let phase_bytes_len = sample_ct.checked_mul(2).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen phased dosage byte count is out of range")
    })?;
    ensure_record_bytes(path, record, phase_cursor, phase_bytes_len)?;
    let record_end = phase_cursor.checked_add(phase_bytes_len).ok_or_else(|| {
        GenoioError::invalid_source(path, "pgen phased dosage offset is out of range")
    })?;
    if record_end != record.len() {
        return Err(GenoioError::invalid_source(
            path,
            "pgen phased dosage record has trailing or missing bytes",
        ));
    }
    if let Some(validated) = validated_phased_dosages.as_deref_mut() {
        validated.clear();
        validated.resize(sample_ct, None);
    }
    for (source_index, total) in dosage_totals.iter().copied().enumerate() {
        let phase_offset = phase_cursor + source_index * 2;
        let phase_delta = decode_phased_dosage_delta(
            path,
            i16::from_le_bytes([record[phase_offset], record[phase_offset + 1]]),
        )?;
        let phased = validate_phased_dosage_pair(path, total, phase_delta)?;
        if let Some(validated) = validated_phased_dosages.as_deref_mut() {
            validated[source_index] = Some(phased);
        }
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
    if missing {
        insert_sorted_unique_index(
            &mut haplotype_state.selected_haplotype_missing_indices,
            offset,
        );
        insert_sorted_unique_index(
            &mut haplotype_state.selected_haplotype_missing_indices,
            offset + 1,
        );
    } else {
        super::remove_sorted_index(
            &mut haplotype_state.selected_haplotype_missing_indices,
            offset,
        );
        super::remove_sorted_index(
            &mut haplotype_state.selected_haplotype_missing_indices,
            offset + 1,
        );
    }
}

fn push_collapsed_dosage(
    haplotype_state: &mut PgenHaplotypeDecodeState,
    value: f32,
    missing: bool,
) {
    let index = haplotype_state.selected_collapsed_values.len();
    haplotype_state.selected_collapsed_values.push(value);
    if missing {
        haplotype_state
            .selected_collapsed_missing_indices
            .push(index);
    }
}
