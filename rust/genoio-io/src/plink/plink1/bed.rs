// pattern: Imperative Shell
//! PLINK1 BED header, seek, and genotype payload decoding.
//!
//! BED stores one variant per packed byte slice in variant-major mode. This
//! module validates that layout and converts source records into packed
//! hard-call scratch owned by the caller.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::GenoioError;

use crate::error::Result;
use crate::hardcall::PackedHardcalls;

pub(super) fn open_bed_file(path: &Path) -> Result<File> {
    let mut file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = [0_u8; 3];
    file.read_exact(&mut header)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    validate_bed_header(path, &header)?;
    Ok(file)
}

fn validate_bed_header(path: &Path, header: &[u8; 3]) -> Result<()> {
    if header[0] != 0x6c || header[1] != 0x1b {
        return Err(GenoioError::invalid_source(path, "invalid bed magic bytes"));
    }
    if header[2] == 0x00 {
        return Err(GenoioError::invalid_source(
            path,
            "sample-major bed mode is not supported",
        ));
    }
    if header[2] != 0x01 {
        return Err(GenoioError::invalid_source(path, "invalid bed mode byte"));
    }
    Ok(())
}

pub(super) fn validate_bed_payload_len(
    path: &Path,
    file: &File,
    n_source_samples: usize,
    n_source_variants: usize,
    bytes_per_variant: usize,
) -> Result<()> {
    let expected_len = 3 + n_source_variants * bytes_per_variant;
    let actual_len = file
        .metadata()
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    let expected_len_u64 = u64::try_from(expected_len)
        .map_err(|_| GenoioError::invalid_source(path, "bed payload length is out of range"))?;
    if actual_len != expected_len_u64 {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "bed payload length {actual_len} does not match {n_source_samples} samples and {n_source_variants} variants"
            ),
        ));
    }
    Ok(())
}

pub(super) fn infer_bed_variant_count(
    path: &Path,
    file: &File,
    n_source_samples: usize,
    bytes_per_variant: usize,
) -> Result<usize> {
    let actual_len = file
        .metadata()
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if actual_len < 3 {
        return Err(GenoioError::invalid_source(
            path,
            "bed payload length is out of range",
        ));
    }
    let payload_len = actual_len - 3;
    let bytes_per_variant_u64 = u64::try_from(bytes_per_variant)
        .map_err(|_| GenoioError::invalid_source(path, "bed payload length is out of range"))?;
    if bytes_per_variant_u64 == 0 || payload_len % bytes_per_variant_u64 != 0 {
        return Err(GenoioError::invalid_source(
            path,
            format!("bed payload length {actual_len} does not match {n_source_samples} samples"),
        ));
    }
    usize::try_from(payload_len / bytes_per_variant_u64)
        .map_err(|_| GenoioError::invalid_source(path, "bed payload length is out of range"))
}

#[derive(Debug, Clone)]
pub(super) struct Plink1DecoderState {
    payload: Vec<u8>,
    pub(super) packed: PackedHardcalls,
    pub(super) values: Vec<f32>,
    pub(super) missing_indices: Vec<usize>,
}

impl Plink1DecoderState {
    pub(super) fn new(
        sample_ct: usize,
        bytes_per_variant: usize,
        selected_sample_ct: usize,
    ) -> Self {
        Self {
            payload: Vec::with_capacity(bytes_per_variant),
            packed: {
                let mut packed = PackedHardcalls::default();
                packed.resize(sample_ct);
                packed
            },
            values: Vec::with_capacity(selected_sample_ct),
            missing_indices: Vec::new(),
        }
    }
}

pub(super) fn read_plink1_variant_packed(
    path: &Path,
    file: &mut File,
    variant_index: usize,
    bytes_per_variant: usize,
    sample_ct: usize,
    decoder_state: &mut Plink1DecoderState,
) -> Result<()> {
    seek_plink1_variant(path, file, variant_index, bytes_per_variant)?;
    read_plink1_variant_packed_sequential(path, file, bytes_per_variant, sample_ct, decoder_state)
}

pub(super) fn seek_plink1_variant(
    path: &Path,
    file: &mut File,
    variant_index: usize,
    bytes_per_variant: usize,
) -> Result<()> {
    let offset = 3 + variant_index * bytes_per_variant;
    file.seek(SeekFrom::Start(u64::try_from(offset).map_err(|_| {
        GenoioError::invalid_source(path, "bed variant offset is out of range")
    })?))
    .map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

pub(super) fn read_plink1_variant_packed_sequential(
    path: &Path,
    file: &mut File,
    bytes_per_variant: usize,
    sample_ct: usize,
    decoder_state: &mut Plink1DecoderState,
) -> Result<()> {
    decoder_state.payload.resize(bytes_per_variant, 0);
    file.read_exact(&mut decoder_state.payload)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    decoder_state
        .packed
        .load_plink1_bed_payload(&decoder_state.payload, sample_ct);
    Ok(())
}
