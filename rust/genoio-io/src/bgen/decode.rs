// pattern: Imperative Shell

use std::io::Read;
use std::path::Path;

use flate2::read::ZlibDecoder;
use genoio_core::GenoioError;

use crate::Result;

use super::header::BgenCompression;
use super::io::{read_u16_le, read_u32_le, read_u8, skip_exact};

const DOSAGE_TOLERANCE: f32 = 1.0e-6;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BgenPhase {
    Unphased,
    Phased,
}

impl BgenPhase {
    fn from_raw(_path: &Path, raw: u8) -> Result<Self> {
        match raw {
            0 => Ok(Self::Unphased),
            1 => Ok(Self::Phased),
            _ => Err(GenoioError::unsupported(
                "unsupported bgen phased probability value; expected 0 or 1",
            )),
        }
    }
}

#[derive(Default)]
pub(super) struct DosageDecodeBuffers {
    pub(super) payload: Vec<u8>,
    pub(super) compressed_payload: Vec<u8>,
    pub(super) selected_values: Vec<f32>,
    pub(super) selected_missing: Vec<bool>,
}

#[derive(Default)]
pub(super) struct HaplotypeDecodeBuffers {
    pub(super) payload: Vec<u8>,
    pub(super) compressed_payload: Vec<u8>,
    pub(super) selected_haplotype_values: Vec<f32>,
    pub(super) selected_haplotype_missing: Vec<bool>,
    pub(super) selected_collapsed_values: Vec<f32>,
    pub(super) selected_collapsed_missing: Vec<bool>,
}

pub(super) fn decode_buffered_dosage_values(
    bgen: &Path,
    sample_count: u32,
    source_indices: &[usize],
    buffers: &mut DosageDecodeBuffers,
) -> Result<()> {
    let DosageDecodeBuffers {
        payload,
        selected_values,
        selected_missing,
        ..
    } = buffers;
    let decoded = DecodedDosageVariant::decode(bgen, payload, sample_count, 2)?;
    decoded.debug_assert_supported_subset();
    decoded.decode_selected_source_order(
        bgen,
        source_indices,
        selected_values,
        selected_missing,
    )?;
    Ok(())
}

pub(super) fn decode_buffered_haplotype_values(
    bgen: &Path,
    sample_count: u32,
    source_indices: &[usize],
    buffers: &mut HaplotypeDecodeBuffers,
) -> Result<()> {
    let HaplotypeDecodeBuffers {
        payload,
        selected_haplotype_values,
        selected_haplotype_missing,
        selected_collapsed_values,
        selected_collapsed_missing,
        ..
    } = buffers;
    let decoded = DecodedDosageVariant::decode(bgen, payload, sample_count, 2)?;
    decoded.debug_assert_supported_subset();
    decoded.decode_selected_phased_haplotypes_source_order(
        bgen,
        source_indices,
        selected_haplotype_values,
        selected_haplotype_missing,
        selected_collapsed_values,
        selected_collapsed_missing,
    )?;
    Ok(())
}

pub(super) fn skip_layout2_probability_block(
    reader: &mut impl Read,
    path: &Path,
    expected_sample_count: u32,
    variant_allele_count: u16,
    compression: BgenCompression,
) -> Result<()> {
    let payload = read_layout2_probability_payload(reader, path, compression)?;
    let decoded =
        DecodedDosageVariant::decode(path, &payload, expected_sample_count, variant_allele_count)?;
    decoded.debug_assert_supported_subset();
    let mut values = Vec::new();
    let mut missing = Vec::new();
    decoded.decode_source_order(path, &mut values, &mut missing)?;
    Ok(())
}

fn read_layout2_probability_payload(
    reader: &mut impl Read,
    path: &Path,
    compression: BgenCompression,
) -> Result<Vec<u8>> {
    let mut payload = Vec::new();
    let mut compressed_payload = Vec::new();
    read_layout2_probability_payload_into(
        reader,
        path,
        compression,
        &mut payload,
        &mut compressed_payload,
    )?;
    Ok(payload)
}

pub(super) fn read_layout2_probability_payload_into(
    reader: &mut impl Read,
    path: &Path,
    compression: BgenCompression,
    payload: &mut Vec<u8>,
    compressed_payload: &mut Vec<u8>,
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
            payload.clear();
            payload.resize(payload_length, 0);
            reader
                .read_exact(payload)
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
            compressed_payload.clear();
            compressed_payload.resize(compressed_payload_length, 0);
            reader
                .read_exact(compressed_payload)
                .map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            decompress_probability_block_into(
                path,
                compression,
                compressed_payload,
                decompressed_block_length,
                payload,
            )
        }
        BgenCompression::Reserved => Err(GenoioError::invalid_source(
            path,
            "bgen compression value is reserved",
        )),
    }
}

pub(super) fn skip_layout2_probability_payload(
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

#[derive(Debug, Clone)]
struct Layout2ProbabilityHeader {
    sample_count: u32,
    allele_count: u16,
    min_ploidy: u8,
    max_ploidy: u8,
    sample_ploidies: Vec<u8>,
    non_missing_sample_count: u32,
    phase: BgenPhase,
    bit_depth: u8,
    byte_len: usize,
}

impl Layout2ProbabilityHeader {
    fn decode(
        path: &Path,
        payload: &[u8],
        expected_sample_count: u32,
        variant_allele_count: u16,
    ) -> Result<Self> {
        let fixed_header_length = Self::fixed_header_length(path, expected_sample_count)?;
        if payload.len() < fixed_header_length {
            return Err(GenoioError::invalid_source(
                path,
                "bgen uncompressed probability block is shorter than the layout 2 header",
            ));
        }

        let reader = &mut &payload[..fixed_header_length];
        let sample_count = read_u32_le(reader, path)?;
        if sample_count != expected_sample_count {
            return Err(GenoioError::invalid_source(
                path,
                "bgen probability block sample count does not match header sample count",
            ));
        }

        let allele_count = read_u16_le(reader, path)?;
        if allele_count != variant_allele_count {
            return Err(GenoioError::invalid_source(
                path,
                "bgen probability block allele count does not match variant allele count",
            ));
        }
        if allele_count != 2 {
            return Err(GenoioError::unsupported(
                "unsupported bgen multiallelic probability block; only biallelic records are supported",
            ));
        }

        let min_ploidy = read_u8(reader, path)?;
        let max_ploidy = read_u8(reader, path)?;
        if min_ploidy != 2 || max_ploidy != 2 {
            return Err(GenoioError::unsupported(
                "unsupported bgen variable-ploidy probability block; only diploid records are supported",
            ));
        }

        let sample_count_usize = usize::try_from(expected_sample_count)
            .map_err(|_| GenoioError::invalid_source(path, "bgen sample count is out of range"))?;
        let mut sample_ploidies = Vec::with_capacity(sample_count_usize);
        let mut non_missing_sample_count = 0_u32;
        for _ in 0..expected_sample_count {
            let ploidy_byte = read_u8(reader, path)?;
            let is_missing = ploidy_byte & 0b1000_0000 != 0;
            let ploidy = ploidy_byte & 0b0011_1111;
            if !is_missing {
                if ploidy != 2 {
                    return Err(GenoioError::unsupported(
                        "unsupported bgen variable-ploidy probability block; only diploid records are supported",
                    ));
                }
                non_missing_sample_count =
                    non_missing_sample_count.checked_add(1).ok_or_else(|| {
                        GenoioError::invalid_source(
                            path,
                            "bgen non-missing sample count is out of range",
                        )
                    })?;
            }
            sample_ploidies.push(ploidy_byte);
        }

        let phase = BgenPhase::from_raw(path, read_u8(reader, path)?)?;

        let bit_depth = read_u8(reader, path)?;
        if !(1..=32).contains(&bit_depth) {
            return Err(GenoioError::unsupported(
                "unsupported bgen probability bit depth; expected 1..=32",
            ));
        }

        Ok(Self {
            sample_count,
            allele_count,
            min_ploidy,
            max_ploidy,
            sample_ploidies,
            non_missing_sample_count,
            phase,
            bit_depth,
            byte_len: fixed_header_length,
        })
    }

    fn fixed_header_length(path: &Path, sample_count: u32) -> Result<usize> {
        let sample_count = usize::try_from(sample_count)
            .map_err(|_| GenoioError::invalid_source(path, "bgen sample count is out of range"))?;
        10_usize.checked_add(sample_count).ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen probability header is out of range")
        })
    }

    fn required_packed_probability_bytes(&self, path: &Path) -> Result<usize> {
        let stored_probability_count = match self.phase {
            // For the supported biallelic diploid subset, unphased records
            // store P(A0/A0) and P(A0/A1); P(A1/A1) is inferred.
            BgenPhase::Unphased => 2,
            // Phased records store one non-last allele probability per
            // haplotype, so diploid biallelic records also store two values.
            BgenPhase::Phased => 2,
        };
        let bits = u64::from(self.non_missing_sample_count)
            .checked_mul(stored_probability_count)
            .and_then(|value| value.checked_mul(u64::from(self.bit_depth)))
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    path,
                    "bgen packed probability bit count is out of range",
                )
            })?;
        let bytes = bits.div_ceil(8);
        usize::try_from(bytes).map_err(|_| {
            GenoioError::invalid_source(path, "bgen packed probability bytes are out of range")
        })
    }
}

#[derive(Debug, Clone)]
struct DecodedDosageVariant<'a> {
    header: Layout2ProbabilityHeader,
    packed_probabilities: &'a [u8],
}

impl<'a> DecodedDosageVariant<'a> {
    fn decode(
        path: &Path,
        payload: &'a [u8],
        expected_sample_count: u32,
        variant_allele_count: u16,
    ) -> Result<Self> {
        let header = Layout2ProbabilityHeader::decode(
            path,
            payload,
            expected_sample_count,
            variant_allele_count,
        )?;
        let packed_probabilities = &payload[header.byte_len..];
        let required_packed_len = header.required_packed_probability_bytes(path)?;
        if packed_probabilities.len() < required_packed_len {
            return Err(GenoioError::invalid_source(
                path,
                "bgen probability block is truncated; packed probabilities are shorter than declared non-missing samples",
            ));
        }
        if packed_probabilities.len() > required_packed_len {
            return Err(GenoioError::invalid_source(
                path,
                "bgen probability block has trailing packed probability bytes",
            ));
        }

        Ok(Self {
            header,
            packed_probabilities,
        })
    }

    fn debug_assert_supported_subset(&self) {
        debug_assert_eq!(
            self.header.sample_count as usize,
            self.header.sample_ploidies.len()
        );
        debug_assert_eq!(self.header.allele_count, 2);
        debug_assert_eq!(self.header.min_ploidy, 2);
        debug_assert_eq!(self.header.max_ploidy, 2);
        debug_assert!((1..=32).contains(&self.header.bit_depth));
        debug_assert!(
            !self.packed_probabilities.is_empty() || self.header.non_missing_sample_count == 0
        );
    }

    fn decode_source_order(
        &self,
        path: &Path,
        values: &mut Vec<f32>,
        missing: &mut Vec<bool>,
    ) -> Result<()> {
        let sample_count = usize::try_from(self.header.sample_count)
            .map_err(|_| GenoioError::invalid_source(path, "bgen sample count is out of range"))?;
        values.clear();
        missing.clear();
        values.reserve(sample_count);
        missing.reserve(sample_count);

        let mut bit_reader = LittleEndianBitReader::new(self.packed_probabilities);
        for &ploidy_byte in &self.header.sample_ploidies {
            let is_missing = ploidy_byte & 0b1000_0000 != 0;
            if is_missing {
                values.push(0.0);
                missing.push(true);
                continue;
            }

            let first_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            let second_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            let value = match self.header.phase {
                BgenPhase::Unphased => {
                    decode_unphased_a1_dosage(path, self.header.bit_depth, first_raw, second_raw)?
                }
                BgenPhase::Phased => {
                    decode_phased_a1_dosage(self.header.bit_depth, first_raw, second_raw)
                }
            };
            values.push(value);
            missing.push(false);
        }

        Ok(())
    }

    fn decode_selected_source_order(
        &self,
        path: &Path,
        source_indices: &[usize],
        values: &mut Vec<f32>,
        missing: &mut Vec<bool>,
    ) -> Result<()> {
        values.clear();
        missing.clear();
        values.reserve(source_indices.len());
        missing.reserve(source_indices.len());

        let mut selected_cursor = 0_usize;
        let mut bit_reader = LittleEndianBitReader::new(self.packed_probabilities);
        for (sample_index, &ploidy_byte) in self.header.sample_ploidies.iter().enumerate() {
            let is_selected = source_indices
                .get(selected_cursor)
                .is_some_and(|&source_index| source_index == sample_index);
            let is_missing = ploidy_byte & 0b1000_0000 != 0;
            if is_missing {
                if is_selected {
                    values.push(0.0);
                    missing.push(true);
                    selected_cursor += 1;
                }
                continue;
            }

            let first_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            let second_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            if is_selected {
                let value = match self.header.phase {
                    BgenPhase::Unphased => decode_unphased_a1_dosage(
                        path,
                        self.header.bit_depth,
                        first_raw,
                        second_raw,
                    )?,
                    BgenPhase::Phased => {
                        decode_phased_a1_dosage(self.header.bit_depth, first_raw, second_raw)
                    }
                };
                values.push(value);
                missing.push(false);
                selected_cursor += 1;
            }
        }
        debug_assert_eq!(selected_cursor, source_indices.len());

        Ok(())
    }

    fn decode_selected_phased_haplotypes_source_order(
        &self,
        path: &Path,
        source_indices: &[usize],
        haplotype_values: &mut Vec<f32>,
        haplotype_missing: &mut Vec<bool>,
        collapsed_values: &mut Vec<f32>,
        collapsed_missing: &mut Vec<bool>,
    ) -> Result<()> {
        if self.header.phase != BgenPhase::Phased {
            return Err(GenoioError::unsupported(
                "unsupported bgen unphased probability block in retained haplotype dosage variant",
            ));
        }
        haplotype_values.clear();
        haplotype_missing.clear();
        collapsed_values.clear();
        collapsed_missing.clear();
        haplotype_values.reserve(source_indices.len() * 2);
        haplotype_missing.reserve(source_indices.len() * 2);
        collapsed_values.reserve(source_indices.len());
        collapsed_missing.reserve(source_indices.len());

        let mut selected_cursor = 0_usize;
        let mut bit_reader = LittleEndianBitReader::new(self.packed_probabilities);
        for (sample_index, &ploidy_byte) in self.header.sample_ploidies.iter().enumerate() {
            let is_selected = source_indices
                .get(selected_cursor)
                .is_some_and(|&source_index| source_index == sample_index);
            let is_missing = ploidy_byte & 0b1000_0000 != 0;
            if is_missing {
                if is_selected {
                    haplotype_values.extend_from_slice(&[0.0, 0.0]);
                    haplotype_missing.extend_from_slice(&[true, true]);
                    collapsed_values.push(0.0);
                    collapsed_missing.push(true);
                    selected_cursor += 1;
                }
                continue;
            }

            let first_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            let second_raw = bit_reader.read_u32(path, self.header.bit_depth)?;
            if is_selected {
                let first = decode_phased_a1_haplotype_dosage(self.header.bit_depth, first_raw);
                let second = decode_phased_a1_haplotype_dosage(self.header.bit_depth, second_raw);
                haplotype_values.extend_from_slice(&[first, second]);
                haplotype_missing.extend_from_slice(&[false, false]);
                collapsed_values.push((first + second).clamp(0.0, 2.0));
                collapsed_missing.push(false);
                selected_cursor += 1;
            }
        }
        debug_assert_eq!(selected_cursor, source_indices.len());

        Ok(())
    }
}

struct LittleEndianBitReader<'a> {
    bytes: &'a [u8],
    bit_offset: usize,
}

impl<'a> LittleEndianBitReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            bit_offset: 0,
        }
    }

    fn read_u32(&mut self, path: &Path, bit_count: u8) -> Result<u32> {
        debug_assert!((1..=32).contains(&bit_count));
        let available_bits = self.bytes.len().checked_mul(8).ok_or_else(|| {
            GenoioError::invalid_source(path, "bgen packed probability bit count is out of range")
        })?;
        let end_bit = self
            .bit_offset
            .checked_add(usize::from(bit_count))
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    path,
                    "bgen packed probability bit count is out of range",
                )
            })?;
        if end_bit > available_bits {
            return Err(GenoioError::invalid_source(
                path,
                "bgen packed probability bits are truncated",
            ));
        }

        if self.bit_offset.is_multiple_of(8) {
            let byte_offset = self.bit_offset / 8;
            let value = match bit_count {
                8 => u32::from(self.bytes[byte_offset]),
                16 => u32::from(u16::from_le_bytes([
                    self.bytes[byte_offset],
                    self.bytes[byte_offset + 1],
                ])),
                32 => u32::from_le_bytes([
                    self.bytes[byte_offset],
                    self.bytes[byte_offset + 1],
                    self.bytes[byte_offset + 2],
                    self.bytes[byte_offset + 3],
                ]),
                _ => {
                    let mut value = 0_u32;
                    for output_bit in 0..bit_count {
                        let source_bit = self.bit_offset + usize::from(output_bit);
                        let byte = self.bytes[source_bit / 8];
                        let bit = (byte >> (source_bit % 8)) & 1;
                        value |= u32::from(bit) << output_bit;
                    }
                    self.bit_offset = end_bit;
                    return Ok(value);
                }
            };
            self.bit_offset = end_bit;
            return Ok(value);
        }

        let mut value = 0_u32;
        for output_bit in 0..bit_count {
            let source_bit = self.bit_offset + usize::from(output_bit);
            let byte = self.bytes[source_bit / 8];
            let bit = (byte >> (source_bit % 8)) & 1;
            value |= u32::from(bit) << output_bit;
        }
        self.bit_offset = end_bit;
        Ok(value)
    }
}

fn decode_unphased_a1_dosage(
    path: &Path,
    bit_depth: u8,
    p_aa_raw: u32,
    p_ab_raw: u32,
) -> Result<f32> {
    let denominator = probability_denominator(bit_depth);
    let p_aa = p_aa_raw as f32 / denominator;
    let p_ab = p_ab_raw as f32 / denominator;
    let p_bb = 1.0 - p_aa - p_ab;
    let dosage = p_ab + 2.0 * p_bb;

    if p_bb < -DOSAGE_TOLERANCE || !(-DOSAGE_TOLERANCE..=2.0 + DOSAGE_TOLERANCE).contains(&dosage) {
        return Err(GenoioError::invalid_source(
            path,
            "bgen malformed probability values produce invalid a1 dosage",
        ));
    }

    Ok(dosage.clamp(0.0, 2.0))
}

fn decode_phased_a1_dosage(bit_depth: u8, p_hap0_a0_raw: u32, p_hap1_a0_raw: u32) -> f32 {
    let denominator = probability_denominator_f64(bit_depth);
    let p_hap0_a0 = p_hap0_a0_raw as f64 / denominator;
    let p_hap1_a0 = p_hap1_a0_raw as f64 / denominator;
    (2.0 - p_hap0_a0 - p_hap1_a0).clamp(0.0, 2.0) as f32
}

fn decode_phased_a1_haplotype_dosage(bit_depth: u8, p_hap_a0_raw: u32) -> f32 {
    let denominator = probability_denominator_f64(bit_depth);
    let p_hap_a0 = p_hap_a0_raw as f64 / denominator;
    (1.0 - p_hap_a0).clamp(0.0, 1.0) as f32
}

fn probability_denominator(bit_depth: u8) -> f32 {
    ((1_u64 << bit_depth) - 1) as f32
}

fn probability_denominator_f64(bit_depth: u8) -> f64 {
    ((1_u64 << bit_depth) - 1) as f64
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
            let decoded =
                zstd::stream::decode_all(compressed_payload).map_err(|source| GenoioError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            decompressed.extend_from_slice(&decoded);
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

#[cfg(test)]
mod tests {
    use super::*;

    fn layout2_payload(
        bit_depth: u8,
        ploidies: &[u8],
        calls: &[Option<(u32, u32)>],
        phased: u8,
    ) -> Vec<u8> {
        assert_eq!(ploidies.len(), calls.len());
        let mut payload = Vec::new();
        payload.extend_from_slice(
            &u32::try_from(calls.len())
                .expect("sample count should fit u32")
                .to_le_bytes(),
        );
        payload.extend_from_slice(&2_u16.to_le_bytes());
        payload.push(2);
        payload.push(2);
        payload.extend_from_slice(ploidies);
        payload.push(phased);
        payload.push(bit_depth);
        append_packed_probabilities(&mut payload, bit_depth, calls);
        payload
    }

    fn append_packed_probabilities(
        output: &mut Vec<u8>,
        bit_depth: u8,
        calls: &[Option<(u32, u32)>],
    ) {
        let mut current_byte = 0_u8;
        let mut bits_in_current_byte = 0_u8;
        for &(p_aa, p_ab) in calls.iter().flatten() {
            append_packed_probability_value(
                output,
                &mut current_byte,
                &mut bits_in_current_byte,
                bit_depth,
                p_aa,
            );
            append_packed_probability_value(
                output,
                &mut current_byte,
                &mut bits_in_current_byte,
                bit_depth,
                p_ab,
            );
        }
        if bits_in_current_byte > 0 {
            output.push(current_byte);
        }
    }

    fn append_packed_probability_value(
        output: &mut Vec<u8>,
        current_byte: &mut u8,
        bits_in_current_byte: &mut u8,
        bit_depth: u8,
        value: u32,
    ) {
        for bit_index in 0..bit_depth {
            let bit = ((value >> bit_index) & 1) as u8;
            *current_byte |= bit << *bits_in_current_byte;
            *bits_in_current_byte += 1;
            if *bits_in_current_byte == 8 {
                output.push(*current_byte);
                *current_byte = 0;
                *bits_in_current_byte = 0;
            }
        }
    }

    fn expected_dosage(bit_depth: u8, p_aa: u32, p_ab: u32) -> f32 {
        let denominator = ((1_u64 << bit_depth) - 1) as f32;
        let p_aa = p_aa as f32 / denominator;
        let p_ab = p_ab as f32 / denominator;
        p_ab + 2.0 * (1.0 - p_aa - p_ab)
    }

    #[test]
    fn decoded_dosage_variant_unpacks_little_endian_probabilities() {
        let payload = layout2_payload(
            3,
            &[2, 0b1000_0010, 2],
            &[Some((5, 2)), None, Some((1, 4))],
            0,
        );
        let decoded = DecodedDosageVariant::decode(Path::new("test.bgen"), &payload, 3, 2)
            .expect("payload should decode");

        let mut values = Vec::new();
        let mut missing = Vec::new();
        decoded
            .decode_source_order(Path::new("test.bgen"), &mut values, &mut missing)
            .expect("dosages should unpack");

        let expected = [expected_dosage(3, 5, 2), 0.0, expected_dosage(3, 1, 4)];
        assert_eq!(missing, vec![false, true, false]);
        for (observed, expected) in values.iter().zip(expected) {
            assert!((observed - expected).abs() <= f32::EPSILON);
        }
    }

    #[test]
    fn decoded_dosage_variant_collapses_phased_haplotype_probabilities() {
        let payload = layout2_payload(
            3,
            &[2, 0b1000_0010, 2],
            &[Some((7, 0)), None, Some((4, 2))],
            1,
        );
        let decoded = DecodedDosageVariant::decode(Path::new("test.bgen"), &payload, 3, 2)
            .expect("payload should decode");

        let mut values = Vec::new();
        let mut missing = Vec::new();
        decoded
            .decode_source_order(Path::new("test.bgen"), &mut values, &mut missing)
            .expect("phased dosages should unpack");

        let denominator = 7.0_f32;
        let expected = [1.0, 0.0, 2.0 - 4.0 / denominator - 2.0 / denominator];
        assert_eq!(missing, vec![false, true, false]);
        for (observed, expected) in values.iter().zip(expected) {
            assert!((observed - expected).abs() <= f32::EPSILON);
        }
    }

    #[test]
    fn decoded_dosage_variant_rejects_impossible_probability_sum() {
        let payload = layout2_payload(1, &[2], &[Some((1, 1))], 0);
        let decoded = DecodedDosageVariant::decode(Path::new("test.bgen"), &payload, 1, 2)
            .expect("payload header should decode");

        let mut values = Vec::new();
        let mut missing = Vec::new();
        let error = decoded
            .decode_source_order(Path::new("test.bgen"), &mut values, &mut missing)
            .expect_err("impossible probabilities should fail");

        assert!(error.to_string().contains("malformed probability"));
    }

    #[test]
    fn decoded_dosage_variant_rejects_trailing_packed_probability_bytes() {
        let mut payload = layout2_payload(8, &[2], &[Some((0, 0))], 0);
        payload.push(0);

        let error = DecodedDosageVariant::decode(Path::new("test.bgen"), &payload, 1, 2)
            .expect_err("trailing packed probability bytes should fail");

        assert!(error.to_string().contains("trailing"));
    }
}
