//! Byte-level FORMAT/DS decoding for the text VCF backend.

use std::path::Path;

use genoio_core::GenoioError;
use noodles_vcf as noodles;

use crate::error::Result;

use super::format::{format_key_index, FormatScanError};

pub(super) struct DsDecodeBuffers {
    values: Vec<f32>,
    missing: Vec<bool>,
}

impl DsDecodeBuffers {
    /// Allocate per-record dosage scratch once and reuse it across records.
    pub(super) fn with_capacity(n_samples: usize) -> Self {
        Self {
            values: Vec::with_capacity(n_samples),
            missing: Vec::with_capacity(n_samples),
        }
    }

    fn clear(&mut self) {
        self.values.clear();
        self.missing.clear();
    }

    pub(super) fn values(&self) -> &[f32] {
        &self.values
    }

    pub(super) fn missing(&self) -> &[bool] {
        &self.missing
    }
}

pub(super) fn decode_ds_record(
    path: &Path,
    record: &noodles::Record,
    source_indices: &[usize],
    output: &mut DsDecodeBuffers,
) -> Result<()> {
    output.clear();
    let samples = record.samples();
    let sample_fields = samples.as_ref().as_bytes();
    let record_name = first_record_id(record);
    let ds_index = format_key_index(sample_fields, b"DS")
        .ok_or_else(|| GenoioError::unsupported("vcf dosage reads require FORMAT/DS values"))?;

    scan_selected_ds_tokens(
        path,
        &record_name,
        sample_fields,
        ds_index,
        source_indices,
        &mut |token| {
            let call = parse_ds_token(path, &record_name, token)?;
            output.values.push(call.value);
            output.missing.push(call.is_missing);
            Ok(())
        },
    )
}

fn scan_selected_ds_tokens(
    path: &Path,
    record_name: &str,
    sample_fields: &[u8],
    key_index: usize,
    source_indices: &[usize],
    emit: &mut impl FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    super::format::scan_selected_format_tokens(
        sample_fields,
        key_index,
        source_indices,
        "sample is missing DS value",
        emit,
    )
    .map_err(|error| match error {
        FormatScanError::Scan(reason) => GenoioError::invalid_source(
            path,
            format!("vcf record {record_name} has unsupported FORMAT/DS: {reason}"),
        ),
        FormatScanError::Emit(error) => error,
    })
}

fn parse_ds_token(path: &Path, record_name: &str, token: &[u8]) -> Result<DsCall> {
    if token == b"." {
        return Ok(DsCall {
            value: 0.0,
            is_missing: true,
        });
    }
    if token.contains(&b',') {
        return Err(GenoioError::invalid_source(
            path,
            format!("vcf record {record_name} has FORMAT/DS with multiple values for a sample; expected one"),
        ));
    }
    let raw = std::str::from_utf8(token).map_err(|error| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {record_name} has non-UTF8 FORMAT/DS value: {error}"),
        )
    })?;
    let value = raw.parse::<f32>().map_err(|error| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {record_name} has invalid FORMAT/DS value {raw}: {error}"),
        )
    })?;
    if !value.is_finite() || !(0.0..=2.0).contains(&value) {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf record {record_name} has invalid FORMAT/DS value {value}; expected finite value in [0, 2]"
            ),
        ));
    }
    Ok(DsCall {
        value,
        is_missing: false,
    })
}

fn first_record_id(record: &noodles::Record) -> String {
    let ids = record.ids();
    let id = ids.as_ref();
    if id.is_empty() {
        ".".to_string()
    } else {
        id.split(';').next().unwrap_or(".").to_string()
    }
}

struct DsCall {
    value: f32,
    is_missing: bool,
}
