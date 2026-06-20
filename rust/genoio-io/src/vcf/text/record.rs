//! Convert noodles lazy records into genoio metadata records.
//!
//! The conversion keeps REF/ALT orientation, finite QUAL values, and concrete
//! region checks in one place for text VCF and BCF callers.

// pattern: Functional Core

use std::path::Path;

use genoio_core::{GenoioError, RegionPredicate, VariantRecord};
use noodles_vcf as noodles;

use crate::error::Result;

use super::super::finite_qual;

/// Convert a noodles VCF/BCF record into genoio's metadata shape.
pub(in crate::vcf) fn variant_record_from_noodles_variant_record<R>(
    path: &Path,
    header: &noodles::Header,
    record: &R,
) -> Result<VariantRecord>
where
    R: noodles::variant::Record + ?Sized,
{
    let chrom = record
        .reference_sequence_name(header)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf chrom error: {error}")))?
        .to_string();
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

    let ref_allele = reference_bases_string(path, record)?;
    let alt_values = alternate_bases_strings(path, record)?;
    let first_alt = alt_values.first().cloned().ok_or_else(|| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {chrom}:{pos} has fewer than two alleles"),
        )
    })?;
    let alt_allele = alt_values.join(",");

    Ok(VariantRecord {
        chrom,
        pos,
        id: first_trait_record_id(record),
        a0: ref_allele.clone(),
        a1: first_alt.clone(),
        ref_allele: Some(ref_allele.clone()),
        alt_allele: Some(alt_allele),
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

pub(super) fn variant_in_region(variant: &VariantRecord, region: &RegionPredicate) -> bool {
    variant.chrom == region.chrom && variant.pos >= region.start && variant.pos <= region.end
}

pub(super) fn skip_variant_for_region(
    variant: &VariantRecord,
    region: Option<&RegionPredicate>,
) -> bool {
    region.is_some_and(|region| !variant_in_region(variant, region))
}

pub(super) fn variant_record_from_text_record(
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

fn first_trait_record_id<R>(record: &R) -> String
where
    R: noodles::variant::Record + ?Sized,
{
    record.ids().iter().next().unwrap_or(".").to_string()
}

fn reference_bases_string<R>(path: &Path, record: &R) -> Result<String>
where
    R: noodles::variant::Record + ?Sized,
{
    let reference_bases = record.reference_bases();
    let mut bases = Vec::with_capacity(reference_bases.len());
    for base in reference_bases.iter() {
        bases.push(base.map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf reference allele error: {error}"))
        })?);
    }
    String::from_utf8(bases).map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf reference allele is not UTF-8: {error}"))
    })
}

fn alternate_bases_strings<R>(path: &Path, record: &R) -> Result<Vec<String>>
where
    R: noodles::variant::Record + ?Sized,
{
    let alternate_bases = record.alternate_bases();
    let mut values = Vec::with_capacity(alternate_bases.len());
    for result in alternate_bases.iter() {
        values.push(result.map(str::to_string).map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf alternate allele error: {error}"))
        })?);
    }
    Ok(values)
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

pub(super) fn validate_biallelic_variant(path: &Path, variant: &VariantRecord) -> Result<()> {
    validate_biallelic_record(
        path,
        &variant.chrom,
        variant.pos,
        variant.alt_allele.as_deref(),
    )
}

fn validate_biallelic_record(path: &Path, chrom: &str, pos: u32, alt: Option<&str>) -> Result<()> {
    if alt.is_some_and(|alt| !alt.is_empty() && !alt.contains(',')) {
        return Ok(());
    }
    if alt.is_some_and(|alt| alt.contains(',')) {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf dense reads require biallelic records; record {chrom}:{pos} has multi-ALT alleles: multi-ALT records are not supported"
            ),
        ));
    }
    Err(GenoioError::invalid_source(
        path,
        format!("vcf dense reads require biallelic records; record {chrom}:{pos} is not biallelic"),
    ))
}
