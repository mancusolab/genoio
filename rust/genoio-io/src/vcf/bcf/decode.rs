// pattern: Functional Core

use std::path::Path;

use genoio_core::{compute_dosage_variant_stats, GenoioError, VariantStats};
use noodles_bcf as bcf;
use noodles_vcf as noodles;
use noodles_vcf::variant::record::{
    samples::{keys::key, series::Value as NoodlesSampleValue},
    Ids as _,
};

use crate::error::Result;
use crate::hardcall::HardcallCounts;

/// Reusable dense BCF decode buffers for one retained variant.
///
/// Missing indices are sparse positions in `values` after sample selection.
pub(super) struct BcfDenseDecodeBuffers {
    pub(super) values: Vec<f32>,
    pub(super) missing_indices: Vec<usize>,
    pub(super) stats: Option<VariantStats>,
    pub(super) counts: Option<HardcallCounts>,
}

impl BcfDenseDecodeBuffers {
    pub(super) fn with_capacity(n_values: usize) -> Self {
        Self {
            values: Vec::with_capacity(n_values),
            missing_indices: Vec::new(),
            stats: None,
            counts: None,
        }
    }

    fn clear(&mut self) {
        self.values.clear();
        self.missing_indices.clear();
        self.stats = None;
        self.counts = None;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum BcfStatsMode {
    Skip,
    Counts,
    Compute,
}

pub(super) fn decode_gt_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    source_indices: &[usize],
    stats_mode: BcfStatsMode,
    decoded: &mut BcfDenseDecodeBuffers,
) -> Result<()> {
    decoded.clear();
    let samples = record.samples().map_err(|error| {
        GenoioError::invalid_source(path, format!("bcf samples error: {error}"))
    })?;
    let gt_series = samples
        .select(header, key::GENOTYPE)
        .ok_or_else(|| GenoioError::invalid_source(path, "bcf record is missing FORMAT/GT"))?
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf genotype series error: {error}"))
        })?;

    let mut counts = HardcallCounts::default();

    for source_index in source_indices {
        let call = decode_gt_call(path, header, record, &gt_series, *source_index)?;
        if !matches!(stats_mode, BcfStatsMode::Skip) {
            match call.class {
                BcfGtClass::HomRef => counts.record_hom_ref(),
                BcfGtClass::Het => counts.record_het(),
                BcfGtClass::HomAlt => counts.record_hom_alt(),
                BcfGtClass::Missing => counts.record_missing(),
            }
        }
        if call.is_missing() {
            decoded.missing_indices.push(decoded.values.len());
        }
        decoded.values.push(call.value);
    }

    if matches!(stats_mode, BcfStatsMode::Compute) {
        decoded.stats = Some(counts.variant_stats()?);
    }
    if matches!(stats_mode, BcfStatsMode::Counts) {
        decoded.counts = Some(counts);
    }
    Ok(())
}

pub(super) fn decode_ds_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    source_indices: &[usize],
    collect_stats: bool,
    decoded: &mut BcfDenseDecodeBuffers,
) -> Result<()> {
    decoded.clear();
    let samples = record.samples().map_err(|error| {
        GenoioError::invalid_source(path, format!("bcf samples error: {error}"))
    })?;
    let ds_series = samples
        .select(header, "DS")
        .ok_or_else(|| {
            GenoioError::unsupported("vcf dosage reads require FORMAT/DS values: missing DS")
        })?
        .map_err(|error| {
            GenoioError::unsupported(format!(
                "vcf dosage reads require FORMAT/DS values: {error}"
            ))
        })?;

    for source_index in source_indices {
        let value = ds_series
            .get(header, *source_index)
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    path,
                    format!(
                        "bcf record {} is missing a DS sample value",
                        record_id(record)
                    ),
                )
            })?
            .transpose()
            .map_err(|error| {
                GenoioError::unsupported(format!(
                    "vcf dosage reads require FORMAT/DS values: {error}"
                ))
            })?;

        let Some(value) = value else {
            decoded.missing_indices.push(decoded.values.len());
            decoded.values.push(0.0);
            continue;
        };
        let NoodlesSampleValue::Float(value) = value else {
            return Err(GenoioError::unsupported(
                "vcf dosage reads require scalar FORMAT/DS float values",
            ));
        };
        if !value.is_finite() || !(0.0..=2.0).contains(&value) {
            return Err(GenoioError::invalid_source(
                path,
                format!(
                    "vcf record {} has invalid FORMAT/DS value {value}; expected finite value in [0, 2]",
                    record_id(record)
                ),
            ));
        }
        decoded.values.push(value);
    }

    if collect_stats {
        decoded.stats = Some(compute_dosage_variant_stats(
            &decoded.values,
            &decoded.missing_indices,
        )?);
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BcfGtClass {
    HomRef,
    Het,
    HomAlt,
    Missing,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct BcfGtCall {
    value: f32,
    class: BcfGtClass,
}

impl BcfGtCall {
    fn is_missing(self) -> bool {
        self.class == BcfGtClass::Missing
    }
}

fn decode_gt_call(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    gt_series: &bcf::record::samples::Series<'_>,
    source_index: usize,
) -> Result<BcfGtCall> {
    let value = gt_series
        .get(header, source_index)
        .ok_or_else(|| {
            GenoioError::invalid_source(
                path,
                format!(
                    "bcf record {} is missing a GT sample value",
                    record_id(record)
                ),
            )
        })?
        .transpose()
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf genotype value error: {error}"))
        })?;

    let Some(NoodlesSampleValue::Genotype(genotype)) = value else {
        return Ok(BcfGtCall {
            value: 0.0,
            class: BcfGtClass::Missing,
        });
    };

    let mut alt_count = 0_u8;
    let mut allele_count = 0_usize;
    for result in genotype.iter() {
        let (allele, _) = result.map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf genotype allele error: {error}"))
        })?;
        if allele_count >= 2 {
            return Err(GenoioError::invalid_source(
                path,
                format!(
                    "vcf record {} has non-diploid GT with at least 3 alleles",
                    record_id(record)
                ),
            ));
        }
        allele_count += 1;
        let Some(allele) = allele else {
            return Ok(BcfGtCall {
                value: 0.0,
                class: BcfGtClass::Missing,
            });
        };
        match allele {
            0 => {}
            1 => alt_count += 1,
            other => {
                return Err(GenoioError::invalid_source(
                    path,
                    format!(
                        "vcf record {} has multiallelic GT allele index {other}",
                        record_id(record)
                    ),
                ));
            }
        }
    }
    if allele_count != 2 {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf record {} has non-diploid GT with {allele_count} alleles",
                record_id(record),
            ),
        ));
    }

    let class = match alt_count {
        0 => BcfGtClass::HomRef,
        1 => BcfGtClass::Het,
        2 => BcfGtClass::HomAlt,
        _ => unreachable!("two diploid GT alleles can only produce dosage 0, 1, or 2"),
    };
    Ok(BcfGtCall {
        value: f32::from(alt_count),
        class,
    })
}

fn record_id(record: &bcf::Record) -> String {
    record.ids().iter().next().unwrap_or(".").to_string()
}
