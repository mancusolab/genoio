//! Convert noodles lazy records into genoio metadata records.

use std::path::Path;

use genoio_core::{GenoioError, VariantRecord};
use noodles_vcf as noodles;

use crate::error::Result;

use super::super::finite_qual;

pub(super) fn variant_record_from_record(
    path: &Path,
    record: &noodles::Record,
) -> Result<VariantRecord> {
    let variant = metadata_variant_record_from_record(path, record)?;
    validate_biallelic_record(
        path,
        &variant.chrom,
        variant.pos,
        variant.alt_allele.as_deref(),
    )?;
    Ok(variant)
}

pub(super) fn metadata_variant_record_from_record(
    path: &Path,
    record: &noodles::Record,
) -> Result<VariantRecord> {
    let chrom = record.reference_sequence_name().to_string();
    if chrom.is_empty() {
        return Err(GenoioError::invalid_source(
            path,
            "vcf record is missing a chromosome id",
        ));
    }
    let pos = record
        .variant_start()
        .transpose()
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf position error: {error}")))?
        .ok_or_else(|| GenoioError::invalid_source(path, "vcf record position is missing"))?;
    let pos = u32::try_from(pos.get())
        .map_err(|_| GenoioError::invalid_source(path, "vcf record position is out of range"))?;

    let ref_allele = record.reference_bases();
    let alternate_bases = record.alternate_bases();
    let alt_allele = alternate_bases.as_ref();
    let ref_allele = ref_allele.to_string();
    let alt_allele = alt_allele.to_string();
    let first_alt = first_alt_allele(&alt_allele).ok_or_else(|| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {chrom}:{pos} has fewer than two alleles"),
        )
    })?;

    Ok(VariantRecord {
        chrom,
        pos,
        id: first_id(record),
        a0: ref_allele.clone(),
        a1: first_alt.clone(),
        ref_allele: Some(ref_allele.clone()),
        alt_allele: Some(alt_allele.clone()),
        source_a0: ref_allele,
        source_a1: first_alt,
        flipped: false,
        qual: record
            .quality_score()
            .transpose()
            .map_err(|error| GenoioError::invalid_source(path, format!("vcf qual error: {error}")))?
            .and_then(finite_qual),
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
}

fn first_id(record: &noodles::Record) -> String {
    let ids = record.ids();
    let id = ids.as_ref();
    if id.is_empty() {
        ".".to_string()
    } else {
        id.split(';').next().unwrap_or(".").to_string()
    }
}

fn first_alt_allele(alt: &str) -> Option<String> {
    if alt.is_empty() {
        None
    } else {
        Some(alt.split(',').next().unwrap_or(alt).to_string())
    }
}

fn validate_biallelic_record(path: &Path, chrom: &str, pos: u32, alt: Option<&str>) -> Result<()> {
    if alt.is_some_and(|alt| !alt.is_empty() && !alt.contains(',')) {
        return Ok(());
    }
    Err(GenoioError::invalid_source(
        path,
        format!("vcf dense reads require biallelic records; record {chrom}:{pos} is not biallelic"),
    ))
}
