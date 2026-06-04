// pattern: Mixed (unavoidable)
// Reason: Format-local binary parsing is kept beside the filesystem entrypoint to match the
// existing reader module pattern in this crate.

use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use flate2::read::ZlibDecoder;
use genoio_core::{MetadataError, MetadataOutput, SampleRecord, SourceCapabilities, VariantRecord};

use crate::Result;

const BGEN_MAGIC: &[u8; 4] = b"bgen";
const ZERO_MAGIC: &[u8; 4] = &[0, 0, 0, 0];
const MIN_HEADER_LENGTH: u32 = 20;
const SAMPLE_IDENTIFIER_FLAG: u32 = 1 << 31;

pub fn read_bgen_metadata(bgen: &Path, sample: Option<&Path>) -> Result<MetadataOutput> {
    let mut reader = File::open(bgen).map_err(|source| MetadataError::Io {
        path: bgen.to_path_buf(),
        source,
    })?;
    let header = BgenHeader::read_from(&mut reader, bgen)?;
    header.validate(bgen)?;

    let samples = if header.flags.has_sample_identifiers {
        read_sample_identifier_block(&mut reader, bgen, header.sample_count)?
    } else if let Some(sample) = sample {
        read_companion_sample_file(sample, header.sample_count)?
    } else {
        return Err(MetadataError::parse(
            bgen,
            "bgen sample identifiers require embedded identifiers or a companion sample path",
        ));
    };

    reader
        .seek(SeekFrom::Start(u64::from(header.offset) + 4))
        .map_err(|source| MetadataError::Io {
            path: bgen.to_path_buf(),
            source,
        })?;
    let variants = read_layout2_variant_metadata(
        &mut reader,
        bgen,
        header.variant_count,
        header.sample_count,
        header.flags.compression,
    )?;

    Ok(MetadataOutput {
        samples,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}

struct BgenHeader {
    offset: u32,
    header_length: u32,
    variant_count: u32,
    sample_count: u32,
    flags: BgenFlags,
}

impl BgenHeader {
    fn read_from(reader: &mut impl Read, path: &Path) -> Result<Self> {
        let offset = read_u32_le(reader, path)?;
        let header_length = read_u32_le(reader, path)?;
        let variant_count = read_u32_le(reader, path)?;
        let sample_count = read_u32_le(reader, path)?;

        let mut magic = [0_u8; 4];
        reader
            .read_exact(&mut magic)
            .map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        if &magic != BGEN_MAGIC && &magic != ZERO_MAGIC {
            return Err(MetadataError::parse(path, "invalid bgen magic bytes"));
        }

        let free_data_length = header_length
            .checked_sub(MIN_HEADER_LENGTH)
            .ok_or_else(|| {
                MetadataError::parse(path, "bgen header length is smaller than 20 bytes")
            })?;
        skip_exact(reader, path, u64::from(free_data_length))?;
        let flags = BgenFlags::from_raw(read_u32_le(reader, path)?);
        Ok(Self {
            offset,
            header_length,
            variant_count,
            sample_count,
            flags,
        })
    }

    fn validate(&self, path: &Path) -> Result<()> {
        if self.header_length < MIN_HEADER_LENGTH {
            return Err(MetadataError::parse(
                path,
                "bgen header length is smaller than 20 bytes",
            ));
        }
        if self.header_length > self.offset {
            return Err(MetadataError::parse(
                path,
                "bgen header length exceeds variant data offset",
            ));
        }
        if self.flags.layout != BgenLayout::Layout2 {
            return Err(MetadataError::parse(
                path,
                "bgen metadata parsing requires layout 2",
            ));
        }
        if self.flags.compression == BgenCompression::Reserved {
            return Err(MetadataError::parse(
                path,
                "bgen compression value is reserved",
            ));
        }
        Ok(())
    }
}

struct BgenFlags {
    compression: BgenCompression,
    layout: BgenLayout,
    has_sample_identifiers: bool,
}

impl BgenFlags {
    fn from_raw(raw: u32) -> Self {
        Self {
            compression: BgenCompression::from_raw(raw & 0b11),
            layout: BgenLayout::from_raw((raw >> 2) & 0b1111),
            has_sample_identifiers: raw & SAMPLE_IDENTIFIER_FLAG != 0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BgenCompression {
    None,
    Zlib,
    Zstd,
    Reserved,
}

impl BgenCompression {
    fn from_raw(raw: u32) -> Self {
        match raw {
            0 => Self::None,
            1 => Self::Zlib,
            2 => Self::Zstd,
            _ => Self::Reserved,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BgenLayout {
    Layout1,
    Layout2,
    Reserved,
}

impl BgenLayout {
    fn from_raw(raw: u32) -> Self {
        match raw {
            1 => Self::Layout1,
            2 => Self::Layout2,
            _ => Self::Reserved,
        }
    }
}

fn read_sample_identifier_block(
    reader: &mut impl Read,
    path: &Path,
    expected_sample_count: u32,
) -> Result<Vec<SampleRecord>> {
    let block_length = read_u32_le(reader, path)?;
    if block_length < 8 {
        return Err(MetadataError::parse(
            path,
            "bgen sample identifiers block is too short",
        ));
    }

    let sample_count = read_u32_le(reader, path)?;
    if sample_count != expected_sample_count {
        return Err(MetadataError::parse(
            path,
            "bgen sample identifiers count does not match header sample count",
        ));
    }

    let mut records = Vec::with_capacity(usize::try_from(sample_count).map_err(|_| {
        MetadataError::parse(path, "bgen sample identifiers count is out of range")
    })?);
    let mut seen = HashSet::with_capacity(records.capacity());
    for _ in 0..sample_count {
        let sample_id = read_len_prefixed_string_u16(reader, path, "sample identifier")?;
        if sample_id.is_empty() {
            return Err(MetadataError::parse(
                path,
                "bgen sample identifier is empty",
            ));
        }
        if !seen.insert(sample_id.clone()) {
            return Err(MetadataError::parse(
                path,
                format!("bgen duplicate sample identifier: {sample_id}"),
            ));
        }
        records.push(sample_record(sample_id));
    }

    Ok(records)
}

fn read_companion_sample_file(
    path: &Path,
    expected_sample_count: u32,
) -> Result<Vec<SampleRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let capacity = usize::try_from(expected_sample_count)
        .map_err(|_| MetadataError::parse(path, "bgen sample count is out of range"))?;
    let mut records = Vec::with_capacity(capacity);
    let mut seen = HashSet::with_capacity(capacity);

    for line in contents.lines().skip(2) {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let sample_id = line.split_whitespace().next().ok_or_else(|| {
            MetadataError::parse(path, "bgen companion sample identifier is empty")
        })?;
        if sample_id.is_empty() {
            return Err(MetadataError::parse(
                path,
                "bgen companion sample identifier is empty",
            ));
        }
        if !seen.insert(sample_id.to_owned()) {
            return Err(MetadataError::parse(
                path,
                format!("bgen duplicate sample identifier: {sample_id}"),
            ));
        }
        records.push(sample_record(sample_id.to_owned()));
    }

    if records.len() != capacity {
        return Err(MetadataError::parse(
            path,
            format!(
                "bgen companion sample count does not match header sample count: expected {capacity}, found {}",
                records.len()
            ),
        ));
    }

    Ok(records)
}

fn sample_record(iid: String) -> SampleRecord {
    SampleRecord {
        fid: None,
        iid,
        father: None,
        mother: None,
        sex: None,
        phenotype: None,
        source_sample_index: None,
        haplotype_index: None,
    }
}

fn read_layout2_variant_metadata(
    reader: &mut impl Read,
    path: &Path,
    variant_count: u32,
    sample_count: u32,
    compression: BgenCompression,
) -> Result<Vec<VariantRecord>> {
    let mut variants = Vec::with_capacity(
        usize::try_from(variant_count)
            .map_err(|_| MetadataError::parse(path, "bgen variant count is out of range"))?,
    );

    for _ in 0..variant_count {
        variants.push(read_layout2_variant_identifying_data(reader, path)?);
        skip_layout2_probability_block(reader, path, sample_count, compression)?;
    }

    if variants.len() != usize::try_from(variant_count).unwrap_or(usize::MAX) {
        return Err(MetadataError::parse(
            path,
            "bgen parsed variant count does not match header variant count",
        ));
    }

    Ok(variants)
}

fn read_layout2_variant_identifying_data(
    reader: &mut impl Read,
    path: &Path,
) -> Result<VariantRecord> {
    let id = read_len_prefixed_string_u16(reader, path, "variant id")?;
    let rsid = read_len_prefixed_string_u16(reader, path, "variant rsid")?;
    let chrom = read_len_prefixed_string_u16(reader, path, "variant chromosome")?;
    let pos = read_u32_le(reader, path)?;
    let allele_count = read_u16_le(reader, path)?;
    if allele_count != 2 {
        return Err(MetadataError::parse(
            path,
            "unsupported bgen multiallelic variant metadata; only biallelic records are supported",
        ));
    }

    let mut alleles = Vec::with_capacity(usize::from(allele_count));
    for _ in 0..allele_count {
        alleles.push(read_len_prefixed_string_u32(
            reader,
            path,
            "variant allele",
        )?);
    }
    let a0 = alleles[0].clone();
    let a1 = alleles[1].clone();
    let id = if rsid.is_empty() { id } else { rsid };

    Ok(VariantRecord {
        chrom,
        pos,
        id,
        a0: a0.clone(),
        a1: a1.clone(),
        ref_allele: Some(a0.clone()),
        alt_allele: Some(a1.clone()),
        source_a0: a0,
        source_a1: a1,
        flipped: false,
        qual: None,
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
}

fn skip_layout2_probability_block(
    reader: &mut impl Read,
    path: &Path,
    expected_sample_count: u32,
    compression: BgenCompression,
) -> Result<()> {
    let block_length = read_u32_le(reader, path)?;
    match compression {
        BgenCompression::None => validate_and_skip_uncompressed_layout2_probability_block(
            reader,
            path,
            block_length,
            expected_sample_count,
        ),
        BgenCompression::Zlib | BgenCompression::Zstd => {
            if block_length < 4 {
                return Err(MetadataError::parse(
                    path,
                    "bgen compressed probability block length is smaller than decompressed length prefix",
                ));
            }
            let decompressed_block_length = read_u32_le(reader, path)?;
            let compressed_payload_length = usize::try_from(block_length - 4).map_err(|_| {
                MetadataError::parse(path, "bgen compressed probability block is out of range")
            })?;
            let compressed_payload = read_exact_vec(reader, path, compressed_payload_length)?;
            let decompressed_payload = decompress_probability_block(
                path,
                compression,
                &compressed_payload,
                decompressed_block_length,
            )?;
            validate_decompressed_layout2_probability_block(
                path,
                &decompressed_payload,
                expected_sample_count,
            )
        }
        BgenCompression::Reserved => Err(MetadataError::parse(
            path,
            "bgen compression value is reserved",
        )),
    }
}

fn validate_and_skip_uncompressed_layout2_probability_block(
    reader: &mut impl Read,
    path: &Path,
    block_length: u32,
    expected_sample_count: u32,
) -> Result<()> {
    let fixed_header_length = 10_u32
        .checked_add(expected_sample_count)
        .ok_or_else(|| MetadataError::parse(path, "bgen sample count is out of range"))?;
    if block_length < fixed_header_length {
        return Err(MetadataError::parse(
            path,
            "bgen uncompressed probability block is shorter than the layout 2 header",
        ));
    }

    let sample_count = read_u32_le(reader, path)?;
    if sample_count != expected_sample_count {
        return Err(MetadataError::parse(
            path,
            "bgen probability block sample count does not match header sample count",
        ));
    }

    let allele_count = read_u16_le(reader, path)?;
    if allele_count != 2 {
        return Err(MetadataError::parse(
            path,
            "unsupported bgen multiallelic probability block; only biallelic records are supported",
        ));
    }

    let min_ploidy = read_u8(reader, path)?;
    let max_ploidy = read_u8(reader, path)?;
    if min_ploidy != 2 || max_ploidy != 2 {
        return Err(MetadataError::parse(
            path,
            "unsupported bgen variable-ploidy probability block; only diploid records are supported",
        ));
    }

    for _ in 0..expected_sample_count {
        let ploidy = read_u8(reader, path)? & 0b0011_1111;
        if ploidy != 2 {
            return Err(MetadataError::parse(
                path,
                "unsupported bgen variable-ploidy probability block; only diploid records are supported",
            ));
        }
    }

    let phased = read_u8(reader, path)?;
    if phased != 0 {
        return Err(MetadataError::parse(
            path,
            "unsupported bgen phased probability block; only unphased records are supported",
        ));
    }

    let _bit_depth = read_u8(reader, path)?;
    let remaining = block_length - fixed_header_length;
    skip_exact(reader, path, u64::from(remaining))
}

fn decompress_probability_block(
    path: &Path,
    compression: BgenCompression,
    compressed_payload: &[u8],
    expected_decompressed_len: u32,
) -> Result<Vec<u8>> {
    let mut decompressed =
        Vec::with_capacity(usize::try_from(expected_decompressed_len).map_err(|_| {
            MetadataError::parse(
                path,
                "bgen decompressed probability block length is out of range",
            )
        })?);
    match compression {
        BgenCompression::Zlib => {
            let mut decoder = ZlibDecoder::new(compressed_payload);
            decoder
                .read_to_end(&mut decompressed)
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
        }
        BgenCompression::Zstd => {
            let mut decoder =
                zstd::stream::read::Decoder::new(compressed_payload).map_err(|source| {
                    MetadataError::Io {
                        path: path.to_path_buf(),
                        source,
                    }
                })?;
            decoder
                .read_to_end(&mut decompressed)
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
        }
        BgenCompression::None | BgenCompression::Reserved => {
            return Err(MetadataError::parse(
                path,
                "bgen compression value is not a compressed probability block",
            ));
        }
    }

    let expected_decompressed_len = usize::try_from(expected_decompressed_len).map_err(|_| {
        MetadataError::parse(
            path,
            "bgen decompressed probability block length is out of range",
        )
    })?;
    if decompressed.len() != expected_decompressed_len {
        return Err(MetadataError::parse(
            path,
            "bgen decompressed probability block length does not match length prefix",
        ));
    }
    Ok(decompressed)
}

fn validate_decompressed_layout2_probability_block(
    path: &Path,
    payload: &[u8],
    expected_sample_count: u32,
) -> Result<()> {
    let block_length = u32::try_from(payload.len()).map_err(|_| {
        MetadataError::parse(path, "bgen decompressed probability block is out of range")
    })?;
    validate_and_skip_uncompressed_layout2_probability_block(
        &mut &payload[..],
        path,
        block_length,
        expected_sample_count,
    )
}

fn skip_exact(reader: &mut impl Read, path: &Path, mut len: u64) -> Result<()> {
    let mut buffer = [0_u8; 8192];
    while len > 0 {
        let chunk_len = buffer
            .len()
            .min(usize::try_from(len).unwrap_or(buffer.len()));
        reader
            .read_exact(&mut buffer[..chunk_len])
            .map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        len -= u64::try_from(chunk_len).expect("skip chunk length should fit u64");
    }
    Ok(())
}

fn read_exact_vec(reader: &mut impl Read, path: &Path, len: usize) -> Result<Vec<u8>> {
    let mut bytes = vec![0_u8; len];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(bytes)
}

fn read_len_prefixed_string_u16(
    reader: &mut impl Read,
    path: &Path,
    label: &str,
) -> Result<String> {
    let len = usize::from(read_u16_le(reader, path)?);
    read_utf8_string(reader, path, label, len)
}

fn read_len_prefixed_string_u32(
    reader: &mut impl Read,
    path: &Path,
    label: &str,
) -> Result<String> {
    let len = usize::try_from(read_u32_le(reader, path)?)
        .map_err(|_| MetadataError::parse(path, format!("bgen {label} length is out of range")))?;
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
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    String::from_utf8(bytes)
        .map_err(|error| MetadataError::parse(path, format!("bgen {label} is not UTF-8: {error}")))
}

fn read_u16_le(reader: &mut impl Read, path: &Path) -> Result<u16> {
    let mut bytes = [0_u8; 2];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(u16::from_le_bytes(bytes))
}

fn read_u32_le(reader: &mut impl Read, path: &Path) -> Result<u32> {
    let mut bytes = [0_u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u8(reader: &mut impl Read, path: &Path) -> Result<u8> {
    let mut byte = [0_u8; 1];
    reader
        .read_exact(&mut byte)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(byte[0])
}
