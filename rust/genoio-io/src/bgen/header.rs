// pattern: Imperative Shell
//! Parse BGEN headers, sample identifiers, and variant metadata records.
//!
//! Header parsing validates the Layout 2 subset supported by this crate. Variant
//! metadata reads skip probability blocks so metadata-only calls avoid dosage
//! allocation and decompression work.

use std::collections::HashSet;
use std::fs;
use std::io::Read;
use std::path::Path;

use genoio_core::{GenoioError, SampleRecord, VariantMetadataArrowBuffers, VariantRecord};

use crate::Result;

use super::decode::skip_layout2_probability_payload_raw;
use super::io::{
    read_exact_vec, read_len_prefixed_string_u16, read_len_prefixed_string_u32,
    read_len_prefixed_utf8_u16_with, read_len_prefixed_utf8_u32_with, read_u16_le, read_u32_le,
    skip_exact, skip_len_prefixed_string_u16, skip_len_prefixed_string_u32,
};

const BGEN_MAGIC: &[u8; 4] = b"bgen";
const ZERO_MAGIC: &[u8; 4] = &[0, 0, 0, 0];
const MIN_HEADER_LENGTH: u32 = 20;
const SAMPLE_IDENTIFIER_FLAG: u32 = 1 << 31;

pub(super) fn read_bgen_samples(
    reader: &mut impl Read,
    bgen: &Path,
    sample: Option<&Path>,
    header: &BgenHeader,
) -> Result<Vec<SampleRecord>> {
    if header.flags.has_sample_identifiers {
        read_sample_identifier_block(reader, bgen, header.sample_count)
    } else if let Some(sample) = sample {
        read_companion_sample_file(sample, header.sample_count)
    } else {
        Err(GenoioError::invalid_source(
            bgen,
            "bgen sample identifiers require embedded identifiers or a companion sample path",
        ))
    }
}

pub(super) struct BgenHeader {
    pub(super) offset: u32,
    header_length: u32,
    pub(super) variant_count: u32,
    pub(super) sample_count: u32,
    pub(super) flags: BgenFlags,
}

impl BgenHeader {
    pub(super) fn read_from(reader: &mut impl Read, path: &Path) -> Result<Self> {
        let offset = read_u32_le(reader, path)?;
        let header_length = read_u32_le(reader, path)?;
        let variant_count = read_u32_le(reader, path)?;
        let sample_count = read_u32_le(reader, path)?;

        let mut magic = [0_u8; 4];
        reader
            .read_exact(&mut magic)
            .map_err(|source| GenoioError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        if &magic != BGEN_MAGIC && &magic != ZERO_MAGIC {
            return Err(GenoioError::invalid_source(
                path,
                "invalid bgen magic bytes",
            ));
        }

        let free_data_length = header_length
            .checked_sub(MIN_HEADER_LENGTH)
            .ok_or_else(|| {
                GenoioError::invalid_source(path, "bgen header length is smaller than 20 bytes")
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

    pub(super) fn validate(&self, path: &Path) -> Result<()> {
        if self.header_length < MIN_HEADER_LENGTH {
            return Err(GenoioError::invalid_source(
                path,
                "bgen header length is smaller than 20 bytes",
            ));
        }
        if self.header_length > self.offset {
            return Err(GenoioError::invalid_source(
                path,
                "bgen header length exceeds variant data offset",
            ));
        }
        if self.flags.layout != BgenLayout::Layout2 {
            return Err(GenoioError::unsupported(
                "bgen metadata parsing requires layout 2",
            ));
        }
        if self.flags.compression == BgenCompression::Reserved {
            return Err(GenoioError::unsupported(
                "bgen compression value is reserved",
            ));
        }
        Ok(())
    }
}

pub(super) struct BgenFlags {
    pub(super) compression: BgenCompression,
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
pub(super) enum BgenCompression {
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
        return Err(GenoioError::invalid_source(
            path,
            "bgen sample identifiers block is too short",
        ));
    }

    let sample_count = read_u32_le(reader, path)?;
    if sample_count != expected_sample_count {
        return Err(GenoioError::invalid_source(
            path,
            "bgen sample identifiers count does not match header sample count",
        ));
    }
    let body_len = usize::try_from(block_length - 8).map_err(|_| {
        GenoioError::invalid_source(path, "bgen sample identifiers block length is out of range")
    })?;
    let body = read_exact_vec(reader, path, body_len)?;
    let mut body_reader = body.as_slice();

    let mut records = Vec::with_capacity(usize::try_from(sample_count).map_err(|_| {
        GenoioError::invalid_source(path, "bgen sample identifiers count is out of range")
    })?);
    let mut seen = HashSet::with_capacity(records.capacity());
    for _ in 0..sample_count {
        let sample_id = read_len_prefixed_string_u16(&mut body_reader, path, "sample identifier")?;
        if sample_id.is_empty() {
            return Err(GenoioError::invalid_source(
                path,
                "bgen sample identifier is empty",
            ));
        }
        if !seen.insert(sample_id.clone()) {
            return Err(GenoioError::invalid_source(
                path,
                format!("bgen duplicate sample identifier: {sample_id}"),
            ));
        }
        records.push(sample_record(sample_id));
    }

    if !body_reader.is_empty() {
        return Err(GenoioError::invalid_source(
            path,
            "bgen sample identifiers block length does not match decoded sample identifiers",
        ));
    }

    Ok(records)
}

fn read_companion_sample_file(
    path: &Path,
    expected_sample_count: u32,
) -> Result<Vec<SampleRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let capacity = usize::try_from(expected_sample_count)
        .map_err(|_| GenoioError::invalid_source(path, "bgen sample count is out of range"))?;
    let mut records = Vec::with_capacity(capacity);
    let mut seen = HashSet::with_capacity(capacity);

    for line in contents.lines().skip(2) {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let sample_id = line.split_whitespace().next().ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen companion sample identifier is empty")
        })?;
        if sample_id.is_empty() {
            return Err(GenoioError::invalid_source(
                path,
                "bgen companion sample identifier is empty",
            ));
        }
        if !seen.insert(sample_id.to_owned()) {
            return Err(GenoioError::invalid_source(
                path,
                format!("bgen duplicate sample identifier: {sample_id}"),
            ));
        }
        records.push(sample_record(sample_id.to_owned()));
    }

    if records.len() != capacity {
        return Err(GenoioError::invalid_source(
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

pub(super) fn read_layout2_variant_metadata_arrow(
    reader: &mut impl Read,
    path: &Path,
    variant_count: u32,
    compression: BgenCompression,
) -> Result<VariantMetadataArrowBuffers> {
    let mut variants = VariantMetadataArrowBuffers::with_capacity(
        usize::try_from(variant_count)
            .map_err(|_| GenoioError::invalid_source(path, "bgen variant count is out of range"))?,
    );
    let mut string_scratch = Vec::new();

    for _ in 0..variant_count {
        read_layout2_variant_identifying_data_arrow(
            reader,
            path,
            &mut variants,
            &mut string_scratch,
        )?;
        skip_layout2_probability_payload_raw(reader, path, compression)?;
    }

    if variants.len() != usize::try_from(variant_count).unwrap_or(usize::MAX) {
        return Err(GenoioError::invalid_source(
            path,
            "bgen parsed variant count does not match header variant count",
        ));
    }

    Ok(variants)
}

pub(super) fn read_layout2_variant_identifying_data(
    reader: &mut impl Read,
    path: &Path,
) -> Result<VariantRecord> {
    let id = read_len_prefixed_string_u16(reader, path, "variant id")?;
    let rsid = read_len_prefixed_string_u16(reader, path, "variant rsid")?;
    let chrom = read_len_prefixed_string_u16(reader, path, "variant chromosome")?;
    let pos = read_u32_le(reader, path)?;
    let allele_count = read_u16_le(reader, path)?;
    if allele_count != 2 {
        return Err(GenoioError::unsupported(
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

fn read_layout2_variant_identifying_data_arrow(
    reader: &mut impl Read,
    path: &Path,
    variants: &mut VariantMetadataArrowBuffers,
    scratch: &mut Vec<u8>,
) -> Result<()> {
    let row_index = variants.len();
    read_len_prefixed_utf8_u16_with(reader, path, "variant id", scratch, |id| {
        variants.ids.append_value(id)
    })?;
    // Public BGEN IDs prefer rsid when present, but the on-disk order puts the
    // fallback ID first. Append the fallback, then replace the row only when an
    // rsid exists so the Arrow path does not materialize both strings.
    read_len_prefixed_utf8_u16_with(reader, path, "variant rsid", scratch, |rsid| {
        if !rsid.is_empty() {
            variants.ids.replace_value(row_index, rsid)?;
        }
        Ok(())
    })?;
    read_len_prefixed_utf8_u16_with(reader, path, "variant chromosome", scratch, |chrom| {
        variants.chroms.append_value(chrom)
    })?;
    let pos = i64::from(read_u32_le(reader, path)?);
    let allele_count = read_u16_le(reader, path)?;
    if allele_count != 2 {
        return Err(GenoioError::unsupported(
            "unsupported bgen multiallelic variant metadata; only biallelic records are supported",
        ));
    }

    read_len_prefixed_utf8_u32_with(reader, path, "variant allele", scratch, |a0| {
        variants.a0s.append_value(a0)?;
        variants.ref_alleles.push(Some(a0.to_owned()));
        variants.source_a0s.append_value(a0)
    })?;
    read_len_prefixed_utf8_u32_with(reader, path, "variant allele", scratch, |a1| {
        variants.a1s.append_value(a1)?;
        variants.alt_alleles.push(Some(a1.to_owned()));
        variants.source_a1s.append_value(a1)
    })?;

    variants.positions.push(pos);
    variants.flipped.push(false);
    variants.quals.push(None);
    variants.afs.push(None);
    variants.mafs.push(None);
    variants.macs.push(None);
    variants.missing_rates.push(None);
    variants.n_called.push(None);
    Ok(())
}

pub(super) fn skip_layout2_variant_identifying_data(
    reader: &mut impl Read,
    path: &Path,
) -> Result<()> {
    skip_len_prefixed_string_u16(reader, path)?;
    skip_len_prefixed_string_u16(reader, path)?;
    skip_len_prefixed_string_u16(reader, path)?;
    skip_exact(reader, path, 4)?;
    let allele_count = read_u16_le(reader, path)?;
    if allele_count != 2 {
        return Err(GenoioError::unsupported(
            "unsupported bgen multiallelic variant metadata; only biallelic records are supported",
        ));
    }
    for _ in 0..allele_count {
        skip_len_prefixed_string_u32(reader, path)?;
    }
    Ok(())
}
