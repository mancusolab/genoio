// pattern: Imperative Shell

use std::io::Read;
use std::path::Path;

use genoio_core::GenoioError;

use crate::Result;

pub(super) fn skip_exact(reader: &mut impl Read, path: &Path, mut len: u64) -> Result<()> {
    let mut buffer = [0_u8; 8192];
    while len > 0 {
        let chunk_len = buffer
            .len()
            .min(usize::try_from(len).unwrap_or(buffer.len()));
        reader
            .read_exact(&mut buffer[..chunk_len])
            .map_err(|source| GenoioError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        len -= chunk_len as u64;
    }
    Ok(())
}

pub(super) fn read_exact_vec(reader: &mut impl Read, path: &Path, len: usize) -> Result<Vec<u8>> {
    let mut bytes = vec![0_u8; len];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(bytes)
}

pub(super) fn read_len_prefixed_string_u16(
    reader: &mut impl Read,
    path: &Path,
    label: &str,
) -> Result<String> {
    let len = usize::from(read_u16_le(reader, path)?);
    read_utf8_string(reader, path, label, len)
}

pub(super) fn read_len_prefixed_string_u32(
    reader: &mut impl Read,
    path: &Path,
    label: &str,
) -> Result<String> {
    let len = usize::try_from(read_u32_le(reader, path)?).map_err(|_| {
        GenoioError::invalid_source(path, format!("bgen {label} length is out of range"))
    })?;
    read_utf8_string(reader, path, label, len)
}

fn read_utf8_string(
    reader: &mut impl Read,
    path: &Path,
    label: &str,
    len: usize,
) -> Result<String> {
    let mut bytes = vec![0_u8; len];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    String::from_utf8(bytes).map_err(|error| {
        GenoioError::invalid_source(path, format!("bgen {label} is not UTF-8: {error}"))
    })
}

pub(super) fn read_u16_le(reader: &mut impl Read, path: &Path) -> Result<u16> {
    let mut bytes = [0_u8; 2];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(u16::from_le_bytes(bytes))
}

pub(super) fn read_u32_le(reader: &mut impl Read, path: &Path) -> Result<u32> {
    let mut bytes = [0_u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(u32::from_le_bytes(bytes))
}

pub(super) fn read_u8(reader: &mut impl Read, path: &Path) -> Result<u8> {
    let mut byte = [0_u8; 1];
    reader
        .read_exact(&mut byte)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(byte[0])
}
