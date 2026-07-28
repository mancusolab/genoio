// pattern: Imperative Shell

use std::io::Read;
use std::path::Path;

use flate2::read::ZlibDecoder;
use genoio_core::GenoioError;

use crate::Result;

use super::super::header::BgenCompression;
use super::super::io::{read_u32_le, skip_exact};

/// Reused scratch space for BGEN probability payload I/O.
///
/// Keeping compressed and decompressed buffers together makes call sites pass a
/// single owner while preserving allocation reuse across variants.
#[derive(Default)]
pub(in crate::bgen) struct ProbabilityPayloadBuffers {
    pub(in crate::bgen) payload: Vec<u8>,
    pub(in crate::bgen) compressed_payload: Vec<u8>,
}

pub(in crate::bgen) fn read_layout2_probability_payload_into(
    reader: &mut impl Read,
    path: &Path,
    compression: BgenCompression,
    buffers: &mut ProbabilityPayloadBuffers,
) -> Result<()> {
    let block_length = read_u32_le(reader, path)?;
    match compression {
        BgenCompression::None => {
            let payload_length = usize::try_from(block_length).map_err(|_| {
                GenoioError::invalid_source(
                    path,
                    "bgen uncompressed probability block is out of range",
                )
            })?;
            buffers.payload.clear();
            buffers.payload.resize(payload_length, 0);
            reader
                .read_exact(&mut buffers.payload)
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            Ok(())
        }
        BgenCompression::Zlib | BgenCompression::Zstd => {
            if block_length < 4 {
                return Err(GenoioError::invalid_source(
                    path,
                    "bgen compressed probability block length is smaller than decompressed length prefix",
                ));
            }
            let decompressed_block_length = read_u32_le(reader, path)?;
            let compressed_payload_length = usize::try_from(block_length - 4).map_err(|_| {
                GenoioError::invalid_source(
                    path,
                    "bgen compressed probability block is out of range",
                )
            })?;
            buffers.compressed_payload.clear();
            buffers
                .compressed_payload
                .resize(compressed_payload_length, 0);
            reader
                .read_exact(&mut buffers.compressed_payload)
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            decompress_probability_block_into(
                path,
                compression,
                &buffers.compressed_payload,
                decompressed_block_length,
                &mut buffers.payload,
            )
        }
        BgenCompression::Reserved => Err(GenoioError::invalid_source(
            path,
            "bgen compression value is reserved",
        )),
    }
}

/// Skip a Layout 2 probability payload by byte length only.
///
/// This is for metadata-only or otherwise discarded records. Retained matrix
/// records must call `read_layout2_probability_payload_into` so the decoded
/// probability contents are validated before use.
pub(in crate::bgen) fn skip_layout2_probability_payload_raw(
    reader: &mut impl Read,
    path: &Path,
    compression: BgenCompression,
) -> Result<()> {
    let block_length = read_u32_le(reader, path)?;
    match compression {
        BgenCompression::None | BgenCompression::Zlib | BgenCompression::Zstd => {
            skip_exact(reader, path, u64::from(block_length))
        }
        BgenCompression::Reserved => Err(GenoioError::invalid_source(
            path,
            "bgen compression value is reserved",
        )),
    }
}

fn decompress_probability_block_into(
    path: &Path,
    compression: BgenCompression,
    compressed_payload: &[u8],
    expected_decompressed_len: u32,
    decompressed: &mut Vec<u8>,
) -> Result<()> {
    let capacity = usize::try_from(expected_decompressed_len).map_err(|_| {
        GenoioError::invalid_source(
            path,
            "bgen decompressed probability block length is out of range",
        )
    })?;
    decompressed.clear();
    decompressed.reserve(capacity);
    match compression {
        BgenCompression::Zlib => {
            let mut decoder = ZlibDecoder::new(compressed_payload);
            decoder
                .read_to_end(decompressed)
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
        }
        BgenCompression::Zstd => {
            let mut decoder =
                zstd::stream::read::Decoder::new(compressed_payload).map_err(|source| {
                    GenoioError::Io {
                        path: path.to_path_buf(),
                        source,
                    }
                })?;
            decoder
                .read_to_end(decompressed)
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
        }
        BgenCompression::None | BgenCompression::Reserved => {
            return Err(GenoioError::invalid_source(
                path,
                "bgen compression value is not a compressed probability block",
            ));
        }
    }

    validate_decompressed_probability_block_len(
        path,
        decompressed.len(),
        expected_decompressed_len,
    )?;
    Ok(())
}

fn validate_decompressed_probability_block_len(
    path: &Path,
    actual_len: usize,
    expected_decompressed_len: u32,
) -> Result<()> {
    let expected_decompressed_len = usize::try_from(expected_decompressed_len).map_err(|_| {
        GenoioError::invalid_source(
            path,
            "bgen decompressed probability block length is out of range",
        )
    })?;
    if actual_len != expected_decompressed_len {
        return Err(GenoioError::invalid_source(
            path,
            "bgen decompressed probability block length does not match length prefix",
        ));
    }
    Ok(())
}
