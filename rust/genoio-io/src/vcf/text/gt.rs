//! Byte-level GT decoding for the text VCF backend.
//!
//! This module intentionally handles only biallelic diploid hardcalls. Broader
//! FORMAT semantics remain unsupported here until they are covered explicitly.

// pattern: Functional Core

use std::path::Path;

use genoio_core::{variant_stats_from_counts, GenoioError, VariantStats};
use noodles_vcf as noodles;

use crate::error::Result;

use super::format::{
    format_key_index, next_delimiter, nth_colon_field, scan_selected_format_tokens, FormatScanError,
};

pub(super) struct GtDecodeBuffers {
    values: Vec<f32>,
    missing: Vec<bool>,
    stats: Option<VariantStats>,
}

pub(super) struct HaplotypeSparseDecodeBuffers {
    a1_rows: Vec<usize>,
    n_rows: usize,
    has_missing: bool,
    stats: Option<VariantStats>,
}

pub(super) struct HaplotypeDenseDecodeBuffers {
    values: Vec<f32>,
    missing: Vec<bool>,
    stats: Option<VariantStats>,
}

impl HaplotypeDenseDecodeBuffers {
    /// Allocate per-record haplotype rows once and reuse them across records.
    pub(super) fn with_capacity(n_samples: usize) -> Self {
        Self {
            values: Vec::with_capacity(n_samples * 2),
            missing: Vec::with_capacity(n_samples * 2),
            stats: None,
        }
    }

    fn clear(&mut self) {
        self.values.clear();
        self.missing.clear();
        self.stats = None;
    }

    fn push_call(&mut self, call: HaplotypeCall) {
        for allele in call.alleles {
            // Missing alleles carry a zero placeholder plus the mask bit, which
            // keeps dense values rectangular without conflating missingness.
            self.values.push(f32::from(allele.unwrap_or(0)));
            self.missing.push(allele.is_none());
        }
    }

    pub(super) fn values(&self) -> &[f32] {
        &self.values
    }

    pub(super) fn missing(&self) -> &[bool] {
        &self.missing
    }
}

impl HaplotypeSparseDecodeBuffers {
    /// Allocate sparse haplotype scratch once and reuse it for every record.
    pub(super) fn with_capacity(n_samples: usize) -> Self {
        Self {
            a1_rows: Vec::with_capacity(n_samples * 2),
            n_rows: n_samples * 2,
            has_missing: false,
            stats: None,
        }
    }

    fn clear(&mut self) {
        self.a1_rows.clear();
        self.has_missing = false;
        self.stats = None;
    }

    fn push_call(&mut self, row_base: usize, call: HaplotypeCall) {
        for (offset, allele) in call.alleles.iter().enumerate() {
            match allele {
                Some(1) => self.a1_rows.push(row_base + offset),
                Some(_) => {}
                None => self.has_missing = true,
            }
        }
    }

    pub(super) fn a1_rows(&self) -> &[usize] {
        &self.a1_rows
    }

    pub(super) fn n_rows(&self) -> usize {
        self.n_rows
    }

    pub(super) fn has_missing(&self) -> bool {
        self.has_missing
    }
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
    let record_location = || format_record_location(record);
    let gt_index = gt_key_index(sample_fields).ok_or_else(|| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {} is missing FORMAT/GT", record_location()),
        )
    })?;

    match stats_mode {
        GtStatsMode::Compute => {
            let mut counts = GtCounts::default();
            scan_selected_gt_tokens(sample_fields, gt_index, source_indices, &mut |token| {
                let call = decode_gt_token(token)?;
                counts.record(call);
                output.values.push(call.value);
                output.missing.push(call.is_missing);
                Ok(())
            })
            .map_err(|reason| gt_error(path, &record_location(), reason))?;
            output.stats = Some(counts.variant_stats()?);
        }
        GtStatsMode::Skip => {
            scan_selected_gt_tokens(sample_fields, gt_index, source_indices, &mut |token| {
                let call = decode_gt_token(token)?;
                output.values.push(call.value);
                output.missing.push(call.is_missing);
                Ok(())
            })
            .map_err(|reason| gt_error(path, &record_location(), reason))?;
        }
    }
    Ok(())
}

pub(super) fn text_record_has_phased_genotype(record: &noodles::Record) -> bool {
    let samples = record.samples();
    let sample_fields = samples.as_ref().as_bytes();
    let Some(gt_index) = gt_key_index(sample_fields) else {
        return false;
    };
    let Some(format_end) = sample_fields.iter().position(|&b| b == b'\t') else {
        return false;
    };

    let mut field_start = format_end + 1;
    while field_start <= sample_fields.len() {
        let field_end = next_delimiter(sample_fields, field_start, b'\t');
        // Capability detection only needs evidence of phase separators. Avoid
        // decoding every allele during metadata scans.
        if nth_colon_field(&sample_fields[field_start..field_end], gt_index)
            .is_some_and(|gt| gt.contains(&b'|'))
        {
            return true;
        }
        if field_end == sample_fields.len() {
            break;
        }
        field_start = field_end + 1;
    }

    false
}

pub(super) fn decode_phased_gt_sparse_record(
    path: &Path,
    record: &noodles::Record,
    source_indices: &[usize],
    stats_mode: GtStatsMode,
    output: &mut HaplotypeSparseDecodeBuffers,
) -> Result<()> {
    output.clear();
    output.stats = decode_selected_phased_gt_record(
        path,
        record,
        source_indices,
        stats_mode,
        &mut |row_base, call| output.push_call(row_base, call),
    )?;
    Ok(())
}

pub(super) fn decode_phased_gt_dense_record(
    path: &Path,
    record: &noodles::Record,
    source_indices: &[usize],
    stats_mode: GtStatsMode,
    output: &mut HaplotypeDenseDecodeBuffers,
) -> Result<()> {
    output.clear();
    output.stats = decode_selected_phased_gt_record(
        path,
        record,
        source_indices,
        stats_mode,
        &mut |_row_base, call| output.push_call(call),
    )?;
    Ok(())
}

fn decode_selected_phased_gt_record(
    path: &Path,
    record: &noodles::Record,
    source_indices: &[usize],
    stats_mode: GtStatsMode,
    emit: &mut impl FnMut(usize, HaplotypeCall),
) -> Result<Option<VariantStats>> {
    let samples = record.samples();
    let sample_fields = samples.as_ref().as_bytes();
    let record_location = || format_record_location(record);
    let gt_index = gt_key_index(sample_fields).ok_or_else(|| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {} is missing FORMAT/GT", record_location()),
        )
    })?;

    match stats_mode {
        GtStatsMode::Compute => {
            let mut counts = GtCounts::default();
            // The callback row base lets dense and sparse callers share the
            // same selected-sample scan while preserving haplotype row order.
            let mut row_base = 0_usize;
            scan_selected_gt_tokens(sample_fields, gt_index, source_indices, &mut |token| {
                let call = decode_phased_gt_token(token)?;
                counts.record_class(call.genotype_class()?);
                emit(row_base, call);
                row_base += 2;
                Ok(())
            })
            .map_err(|reason| gt_error(path, &record_location(), reason))?;
            Ok(Some(counts.variant_stats()?))
        }
        GtStatsMode::Skip => {
            // Keep row accounting in the no-stats branch so output order is
            // identical regardless of filter shape.
            let mut row_base = 0_usize;
            scan_selected_gt_tokens(sample_fields, gt_index, source_indices, &mut |token| {
                emit(row_base, decode_phased_gt_token(token)?);
                row_base += 2;
                Ok(())
            })
            .map_err(|reason| gt_error(path, &record_location(), reason))?;
            Ok(None)
        }
    }
}

const UNPHASED_HAPLOTYPE_GT: &str =
    "contains an unphased GT separator in a retained haplotype variant";

fn gt_error(path: &Path, record_location: &str, reason: &str) -> GenoioError {
    if reason == UNPHASED_HAPLOTYPE_GT {
        return GenoioError::unsupported(format!(
            "vcf record {record_location} has unsupported GT: {reason}"
        ));
    }
    GenoioError::invalid_source(
        path,
        format!("vcf record {record_location} has unsupported GT: {reason}"),
    )
}

fn format_record_location(record: &noodles::Record) -> String {
    format!(
        "{}:{}",
        record.reference_sequence_name(),
        record_pos(record)
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
    format_key_index(sample_fields, b"GT")
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
        self.record_class(call.class);
    }

    fn record_class(&mut self, class: GtClass) {
        match class {
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

fn scan_selected_gt_tokens(
    sample_fields: &[u8],
    gt_index: usize,
    source_indices: &[usize],
    emit: &mut impl FnMut(&[u8]) -> std::result::Result<(), &'static str>,
) -> std::result::Result<(), &'static str> {
    scan_selected_format_tokens(
        sample_fields,
        gt_index,
        source_indices,
        "sample is missing GT value",
        emit,
    )
    .map_err(|error| match error {
        FormatScanError::Scan(reason) | FormatScanError::Emit(reason) => reason,
    })
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
        _ if is_non_diploid_gt_token(token) => Err("non-diploid GT"),
        _ if is_multiallelic_gt_token(token) => Err("multiallelic GT"),
        _ => Err("expected diploid biallelic hardcall"),
    }
}

fn is_non_diploid_gt_token(token: &[u8]) -> bool {
    let allele_count = token
        .split(|byte| matches!(byte, b'/' | b'|'))
        .filter(|allele| !allele.is_empty())
        .count();
    allele_count != 2
}

fn is_multiallelic_gt_token(token: &[u8]) -> bool {
    token
        .split(|byte| matches!(byte, b'/' | b'|'))
        .any(|allele| {
            std::str::from_utf8(allele)
                .ok()
                .and_then(|allele| allele.parse::<u32>().ok())
                .is_some_and(|index| index > 1)
        })
}

#[derive(Debug)]
struct HaplotypeCall {
    alleles: [Option<u8>; 2],
}

impl HaplotypeCall {
    fn genotype_class(&self) -> std::result::Result<GtClass, &'static str> {
        match self.alleles {
            [None, _] | [_, None] => Ok(GtClass::Missing),
            [Some(0), Some(0)] => Ok(GtClass::HomRef),
            [Some(1), Some(1)] => Ok(GtClass::HomAlt),
            [Some(0), Some(1)] | [Some(1), Some(0)] => Ok(GtClass::Het),
            _ => Err("expected diploid phased biallelic hardcall"),
        }
    }
}

fn decode_phased_gt_token(token: &[u8]) -> std::result::Result<HaplotypeCall, &'static str> {
    if token.len() != 3 {
        return Err("expected diploid phased biallelic hardcall");
    }
    if token[1] == b'/' {
        return Err(UNPHASED_HAPLOTYPE_GT);
    }
    if token[1] != b'|' {
        return Err("expected diploid phased biallelic hardcall");
    }

    let first = decode_phased_allele(token[0])?;
    let second = decode_phased_allele(token[2])?;
    Ok(HaplotypeCall {
        alleles: [first, second],
    })
}

fn decode_phased_allele(raw: u8) -> std::result::Result<Option<u8>, &'static str> {
    match raw {
        b'0' => Ok(Some(0)),
        b'1' => Ok(Some(1)),
        b'.' => Ok(None),
        _ => Err("expected diploid phased biallelic hardcall"),
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
        scan_selected_gt_tokens(
            b"DP:GT\t9:0/0\t8:0/1\t7:1/1\t6:./.",
            1,
            &[1, 3],
            &mut |token| {
                calls.push(decode_gt_token(token)?);
                Ok(())
            },
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

    #[test]
    fn decode_phased_gt_token_rejects_unphased_separator() {
        assert!(decode_phased_gt_token(b"0/1")
            .expect_err("unphased separator should fail")
            .contains("unphased"));
        let call = decode_phased_gt_token(b"1|0").expect("phased GT should decode");
        assert_eq!(call.alleles, [Some(1), Some(0)]);
        assert_eq!(
            call.genotype_class().expect("class should compute"),
            GtClass::Het
        );
    }
}
