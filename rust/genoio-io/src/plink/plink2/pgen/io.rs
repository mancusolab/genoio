// pattern: Imperative Shell
//! PGEN record seek and payload read helpers.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;

use super::header::{fixed_width_phased_dosage_record_len, fixed_width_record_len};
use super::main_track::{
    decode_difflist_record, decode_ld_compressed_record, decode_one_bit_record,
};
use super::{PgenDecoderState, PgenHeader, PgenLayout, PGEN_HEADER_LEN};

pub(in crate::plink::plink2) fn read_plink2_variant_packed(
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

pub(in crate::plink::plink2) fn seek_fixed_width_variant_record(
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

pub(in crate::plink::plink2) fn read_fixed_width_variant_packed_sequential(
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

pub(super) fn read_fixed_width_phased_dosage_variant_record(
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
