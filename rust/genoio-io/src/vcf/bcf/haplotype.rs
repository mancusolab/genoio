// pattern: Functional Core

use std::path::Path;

use genoio_core::GenoioError;
use noodles_bcf as bcf;
use noodles_vcf as noodles;
use noodles_vcf::variant::record::{
    samples::series::value::genotype::Phasing as NoodlesGenotypePhasing,
    samples::{keys::key, series::Value as NoodlesSampleValue},
    Ids as _,
};

use crate::error::Result;

/// Reusable dense phased BCF haplotype decode buffers for one retained variant.
///
/// Missing indices are haplotype-row positions in `values`.
pub(super) struct BcfHaplotypeDecodeBuffers {
    pub(super) values: Vec<f32>,
    pub(super) missing_indices: Vec<usize>,
}

impl BcfHaplotypeDecodeBuffers {
    pub(super) fn with_capacity(n_samples: usize) -> Self {
        Self {
            values: Vec::with_capacity(n_samples * 2),
            missing_indices: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.values.clear();
        self.missing_indices.clear();
    }
}

pub(super) fn decode_phased_haplotype_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    source_indices: &[usize],
    decoded: &mut BcfHaplotypeDecodeBuffers,
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

    for source_index in source_indices {
        let (sample_values, sample_missing) =
            decode_phased_haplotype_call(path, header, record, &gt_series, *source_index)?;
        let row_offset = decoded.values.len();
        decoded.values.extend(sample_values);
        for (offset, is_missing) in sample_missing.into_iter().enumerate() {
            if is_missing {
                decoded.missing_indices.push(row_offset + offset);
            }
        }
    }

    Ok(())
}

fn decode_phased_haplotype_call(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    gt_series: &bcf::record::samples::Series<'_>,
    source_index: usize,
) -> Result<([f32; 2], [bool; 2])> {
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
        return Ok(([0.0, 0.0], [true, true]));
    };
    let mut values = [0.0, 0.0];
    let mut missing = [false, false];
    let mut allele_count = 0_usize;
    for result in genotype.iter() {
        let (allele, phasing) = result.map_err(|error| {
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
        let allele_index = allele_count;
        allele_count += 1;
        if allele_index > 0 && phasing == NoodlesGenotypePhasing::Unphased {
            return Err(GenoioError::unsupported(format!(
                "vcf haplotype read record {} contains an unphased GT separator in a retained haplotype variant",
                record_id(record)
            )));
        }
        match allele {
            None => missing[allele_index] = true,
            Some(0) => {}
            Some(1) => values[allele_index] = 1.0,
            Some(other) => {
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

    Ok((values, missing))
}

fn record_id(record: &bcf::Record) -> String {
    record.ids().iter().next().unwrap_or(".").to_string()
}
