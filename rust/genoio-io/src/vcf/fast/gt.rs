//! Byte-level GT decoding for the compressed VCF fast path.
//!
//! This module intentionally handles only biallelic diploid hardcalls. Broader
//! FORMAT semantics stay on the htslib path, where correctness coverage is
//! better than a custom parser would be.

use std::path::Path;

use genoio_core::{variant_stats_from_counts, GenoioError, VariantStats};
use noodles_vcf as noodles;

use crate::error::Result;

pub(super) struct GtDecodeBuffers {
    values: Vec<f32>,
    missing: Vec<bool>,
    stats: Option<VariantStats>,
}

impl GtDecodeBuffers {
    /// Allocate per-record scratch once and reuse it for every retained record.
    pub(super) fn with_capacity(n_samples: usize) -> Self {
        Self {
            values: Vec::with_capacity(n_samples),
            missing: Vec::with_capacity(n_samples),
            stats: None,
        }
    }

    fn clear(&mut self) {
        self.values.clear();
        self.missing.clear();
        self.stats = None;
    }

    pub(super) fn values(&self) -> &[f32] {
        &self.values
    }

    pub(super) fn values_mut(&mut self) -> &mut [f32] {
        &mut self.values
    }

    pub(super) fn missing(&self) -> &[bool] {
        &self.missing
    }

    pub(super) fn stats(&self) -> Option<VariantStats> {
        self.stats
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum GtStatsMode {
    Skip,
    Compute,
}

impl GtStatsMode {
    pub(super) const fn from_needed(needed: bool) -> Self {
        if needed {
            Self::Compute
        } else {
            Self::Skip
        }
    }
}

pub(super) fn decode_gt_record(
    path: &Path,
    record: &noodles::Record,
    source_indices: &[usize],
    stats_mode: GtStatsMode,
    output: &mut GtDecodeBuffers,
) -> Result<()> {
    output.clear();
    let samples = record.samples();
    let sample_fields = samples.as_ref().as_bytes();
    let record_location = || {
        format!(
            "{}:{}",
            record.reference_sequence_name(),
            record_pos(record)
        )
    };
    let gt_index = gt_key_index(sample_fields).ok_or_else(|| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {} is missing FORMAT/GT", record_location()),
        )
    })?;

    match stats_mode {
        GtStatsMode::Compute => {
            let mut counts = GtCounts::default();
            scan_selected_gt(sample_fields, gt_index, source_indices, &mut |call| {
                counts.record(call);
                output.values.push(call.value);
                output.missing.push(call.is_missing);
            })
            .map_err(|reason| gt_error(path, &record_location(), reason))?;
            output.stats = Some(counts.variant_stats()?);
        }
        GtStatsMode::Skip => {
            scan_selected_gt(sample_fields, gt_index, source_indices, &mut |call| {
                output.values.push(call.value);
                output.missing.push(call.is_missing);
            })
            .map_err(|reason| gt_error(path, &record_location(), reason))?;
        }
    }
    Ok(())
}

fn gt_error(path: &Path, record_location: &str, reason: &str) -> GenoioError {
    GenoioError::invalid_source(
        path,
        format!("vcf record {record_location} has unsupported GT: {reason}"),
    )
}

fn record_pos(record: &noodles::Record) -> usize {
    record
        .variant_start()
        .and_then(|result| result.ok())
        .map(|pos| pos.get())
        .unwrap_or_default()
}

fn gt_key_index(sample_fields: &[u8]) -> Option<usize> {
    let format_end = sample_fields.iter().position(|&b| b == b'\t')?;
    sample_fields[..format_end]
        .split(|&b| b == b':')
        .position(|key| key == b"GT")
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct GtCall {
    value: f32,
    is_missing: bool,
    class: GtClass,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GtClass {
    HomRef,
    Het,
    HomAlt,
    Missing,
}

impl GtCall {
    const fn missing() -> Self {
        Self {
            value: 0.0,
            is_missing: true,
            class: GtClass::Missing,
        }
    }
}

#[derive(Default)]
struct GtCounts {
    hom_ref: u64,
    het: u64,
    hom_alt: u64,
    missing: u64,
}

impl GtCounts {
    fn record(&mut self, call: GtCall) {
        match call.class {
            GtClass::HomRef => self.hom_ref += 1,
            GtClass::Het => self.het += 1,
            GtClass::HomAlt => self.hom_alt += 1,
            GtClass::Missing => self.missing += 1,
        }
    }

    fn variant_stats(&self) -> Result<VariantStats> {
        variant_stats_from_counts(self.hom_ref, self.het, self.hom_alt, self.missing)
    }
}

fn scan_selected_gt(
    sample_fields: &[u8],
    gt_index: usize,
    source_indices: &[usize],
    emit: &mut impl FnMut(GtCall),
) -> std::result::Result<(), &'static str> {
    let Some(format_end) = sample_fields.iter().position(|&b| b == b'\t') else {
        return Err("record has FORMAT but no sample columns");
    };
    let mut selected_index = 0_usize;
    let mut sample_index = 0_usize;
    let mut field_start = format_end + 1;

    while selected_index < source_indices.len() {
        let target_index = source_indices[selected_index];
        if field_start > sample_fields.len() {
            return Err("selected sample index is outside the record");
        }
        let field_end = next_delimiter(sample_fields, field_start, b'\t');
        if sample_index == target_index {
            let sample = &sample_fields[field_start..field_end];
            emit(decode_gt_sample(sample, gt_index)?);
            selected_index += 1;
        }
        sample_index += 1;
        if field_end == sample_fields.len() {
            field_start = sample_fields.len() + 1;
        } else {
            field_start = field_end + 1;
        }
    }

    Ok(())
}

fn decode_gt_sample(sample: &[u8], gt_index: usize) -> std::result::Result<GtCall, &'static str> {
    let token = nth_colon_field(sample, gt_index).ok_or("sample is missing GT value")?;
    decode_gt_token(token)
}

fn nth_colon_field(sample: &[u8], index: usize) -> Option<&[u8]> {
    let mut field_start = 0_usize;
    for field_index in 0..=index {
        let field_end = next_delimiter(sample, field_start, b':');
        if field_index == index {
            return Some(&sample[field_start..field_end]);
        }
        if field_end == sample.len() {
            return None;
        }
        field_start = field_end + 1;
    }
    None
}

fn next_delimiter(buf: &[u8], start: usize, delimiter: u8) -> usize {
    buf[start..]
        .iter()
        .position(|&b| b == delimiter)
        .map_or(buf.len(), |offset| start + offset)
}

fn decode_gt_token(token: &[u8]) -> std::result::Result<GtCall, &'static str> {
    match token {
        b"0/0" | b"0|0" => Ok(GtCall {
            value: 0.0,
            is_missing: false,
            class: GtClass::HomRef,
        }),
        b"0/1" | b"1/0" | b"0|1" | b"1|0" => Ok(GtCall {
            value: 1.0,
            is_missing: false,
            class: GtClass::Het,
        }),
        b"1/1" | b"1|1" => Ok(GtCall {
            value: 2.0,
            is_missing: false,
            class: GtClass::HomAlt,
        }),
        b"./." | b".|." => Ok(GtCall::missing()),
        _ => Err("expected diploid biallelic hardcall"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_key_index_finds_gt_in_format_keys() {
        assert_eq!(gt_key_index(b"GT:DS\t0/1:0.4"), Some(0));
        assert_eq!(gt_key_index(b"DP:GT:GQ\t9:1/1:99"), Some(1));
        assert_eq!(gt_key_index(b"DP:GQ\t9:99"), None);
    }

    #[test]
    fn scan_selected_gt_decodes_source_order_subset() {
        let mut calls = Vec::new();
        scan_selected_gt(
            b"DP:GT\t9:0/0\t8:0/1\t7:1/1\t6:./.",
            1,
            &[1, 3],
            &mut |call| calls.push(call),
        )
        .expect("selected GTs should decode");

        assert_eq!(
            calls,
            vec![
                GtCall {
                    value: 1.0,
                    is_missing: false,
                    class: GtClass::Het,
                },
                GtCall::missing(),
            ]
        );
    }

    #[test]
    fn decode_gt_token_rejects_multiallelic_or_non_diploid_calls() {
        assert!(decode_gt_token(b"1/2").is_err());
        assert!(decode_gt_token(b"0").is_err());
        assert!(decode_gt_token(b"0/0/1").is_err());
    }
}
