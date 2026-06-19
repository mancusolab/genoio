//! Byte-level FORMAT/DS decoding for the compressed VCF fast path.

use std::path::Path;

use genoio_core::GenoioError;
use noodles_vcf as noodles;

use crate::error::Result;

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

    scan_selected_format_tokens(
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

fn format_key_index(sample_fields: &[u8], key: &[u8]) -> Option<usize> {
    let format_end = sample_fields.iter().position(|&b| b == b'\t')?;
    sample_fields[..format_end]
        .split(|&b| b == b':')
        .position(|candidate| candidate == key)
}

fn scan_selected_format_tokens(
    path: &Path,
    record_name: &str,
    sample_fields: &[u8],
    key_index: usize,
    source_indices: &[usize],
    emit: &mut impl FnMut(&[u8]) -> Result<()>,
) -> Result<()> {
    let Some(format_end) = sample_fields.iter().position(|&b| b == b'\t') else {
        return Err(GenoioError::invalid_source(
            path,
            format!("vcf record {record_name} has FORMAT but no sample columns"),
        ));
    };
    let mut selected_index = 0_usize;
    let mut sample_index = 0_usize;
    let mut field_start = format_end + 1;

    while selected_index < source_indices.len() {
        let target_index = source_indices[selected_index];
        if field_start > sample_fields.len() {
            return Err(GenoioError::invalid_source(
                path,
                format!("vcf record {record_name} is missing a selected sample column"),
            ));
        }
        let field_end = next_delimiter(sample_fields, field_start, b'\t');
        if sample_index == target_index {
            emit(
                nth_colon_field(&sample_fields[field_start..field_end], key_index).ok_or_else(
                    || {
                        GenoioError::invalid_source(
                            path,
                            format!("vcf record {record_name} sample is missing DS value"),
                        )
                    },
                )?,
            )?;
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

struct DsCall {
    value: f32,
    is_missing: bool,
}
