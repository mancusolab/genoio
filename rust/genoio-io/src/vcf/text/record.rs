//! Convert noodles lazy records into genoio metadata records.
//!
//! The conversion keeps REF/ALT orientation, finite QUAL values, and concrete
//! region checks in one place for text VCF and BCF callers.

// pattern: Functional Core

use std::path::Path;

use genoio_core::{
    DenseDiagnostics, GenoioError, PartialFilterDecision, RegionPredicate, VariantFilter,
    VariantMetadataBuffers, VariantMetadataView,
};
use noodles_vcf as noodles;

use crate::error::Result;
use crate::retention::{MetadataRetentionAction, RetainedVariantState};

use super::super::finite_qual;

/// Borrowed metadata view over a noodles text VCF record.
///
/// Text VCF matrix reads use this to evaluate metadata predicates and append
/// returned variant rows without allocating a full [`VariantRecord`] first.
pub(super) struct TextVariantView<'a> {
    chrom: &'a str,
    pos: u32,
    ids: noodles::record::Ids<'a>,
    ref_allele: &'a str,
    alternate_bases: noodles::record::AlternateBases<'a>,
    qual: Option<f32>,
}

pub(super) struct PreparedTextCandidate<'a> {
    pub(super) variant: TextVariantView<'a>,
    pub(super) needs_genotype_decision: bool,
}

pub(super) enum TextCandidateAction<'a> {
    Skip,
    Stop,
    Decode(PreparedTextCandidate<'a>),
}

impl VariantMetadataView for TextVariantView<'_> {
    fn chrom(&self) -> &str {
        self.chrom
    }

    fn pos(&self) -> u32 {
        self.pos
    }

    fn id(&self) -> &str {
        first_id_str(self.ids.as_ref())
    }

    fn a0(&self) -> &str {
        self.ref_allele
    }

    fn a1(&self) -> &str {
        first_alt_allele_str(self.alternate_bases.as_ref()).unwrap_or("")
    }

    fn ref_allele(&self) -> Option<&str> {
        Some(self.ref_allele)
    }

    fn alt_allele(&self) -> Option<&str> {
        Some(self.alternate_bases.as_ref())
    }

    fn qual(&self) -> Option<f32> {
        self.qual
    }
}

pub(in crate::vcf) fn append_public_variant_metadata_from_noodles_variant_record<R>(
    path: &Path,
    header: &noodles::Header,
    record: &R,
    variants: &mut VariantMetadataBuffers,
) -> Result<()>
where
    R: noodles::variant::Record + ?Sized,
{
    let chrom = record
        .reference_sequence_name(header)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf chrom error: {error}")))?;
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
    let pos = i64::try_from(pos.get())
        .map_err(|_| GenoioError::invalid_source(path, "vcf record position is out of range"))?;

    let ref_allele = reference_bases_string(path, record)?;
    let first_alt = first_alternate_base_string(path, record)?.ok_or_else(|| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {chrom}:{pos} has fewer than two alleles"),
        )
    })?;
    let id = first_trait_record_id(record);

    variants.push(chrom, pos, &id, &ref_allele, &first_alt)
}

pub(super) fn variant_in_region_view<V: VariantMetadataView + ?Sized>(
    variant: &V,
    region: &RegionPredicate,
) -> bool {
    variant.chrom() == region.chrom && variant.pos() >= region.start && variant.pos() <= region.end
}

pub(super) fn skip_variant_for_region<V: VariantMetadataView + ?Sized>(
    variant: &V,
    region: Option<&RegionPredicate>,
) -> bool {
    region.is_some_and(|region| !variant_in_region_view(variant, region))
}

pub(super) fn text_variant_view_from_text_record<'a>(
    path: &Path,
    record: &'a noodles::Record,
) -> Result<TextVariantView<'a>> {
    let chrom = text_record_chrom(path, record)?;
    let pos = text_record_pos(path, record)?;
    let ref_allele = record.reference_bases();
    let alternate_bases = record.alternate_bases();
    let alt_allele = alternate_bases.as_ref();
    // Validate the first ALT while the path and coordinates are available, but
    // keep the full ALT field borrowed for metadata output.
    let _ = first_alt_allele_str(alt_allele).ok_or_else(|| {
        GenoioError::invalid_source(
            path,
            format!("vcf record {chrom}:{pos} has fewer than two alleles"),
        )
    })?;
    let qual = record
        .quality_score()
        .transpose()
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf qual error: {error}")))?
        .and_then(finite_qual);

    Ok(TextVariantView {
        chrom,
        pos,
        ids: record.ids(),
        ref_allele,
        alternate_bases,
        qual,
    })
}

/// Apply candidate-local metadata policy shared by stateless and persistent text reads.
pub(super) fn prepare_text_candidate<'a>(
    path: &Path,
    record: &'a noodles::Record,
    region: Option<&RegionPredicate>,
    variant_filter: Option<&VariantFilter>,
    retention: &mut RetainedVariantState,
    diagnostics: &mut DenseDiagnostics,
) -> Result<TextCandidateAction<'a>> {
    let variant = text_variant_view_from_text_record(path, record)?;
    if skip_variant_for_region(&variant, region) {
        return Ok(TextCandidateAction::Skip);
    }
    let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
        filter.partial_decision_view(&variant)
    });
    match retention.metadata_decision(partial_decision, diagnostics) {
        MetadataRetentionAction::Skip => return Ok(TextCandidateAction::Skip),
        MetadataRetentionAction::Stop => return Ok(TextCandidateAction::Stop),
        MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
    }
    validate_biallelic_variant(path, &variant)?;
    Ok(TextCandidateAction::Decode(PreparedTextCandidate {
        variant,
        needs_genotype_decision: matches!(partial_decision, PartialFilterDecision::NeedGenotypes),
    }))
}

pub(super) fn append_public_variant_metadata_from_text_record(
    path: &Path,
    record: &noodles::Record,
    variants: &mut VariantMetadataBuffers,
) -> Result<()> {
    let view = text_variant_view_from_text_record(path, record)?;
    append_public_variant_metadata_from_text_view(&view, variants)
}

pub(super) fn append_public_variant_metadata_from_text_view(
    variant: &TextVariantView<'_>,
    variants: &mut VariantMetadataBuffers,
) -> Result<()> {
    // Metadata-only text reads share the same borrowed append path as matrix
    // reads so REF/ALT/id handling stays consistent.
    variants.push_view(variant)
}

fn text_record_chrom<'a>(path: &Path, record: &'a noodles::Record) -> Result<&'a str> {
    let chrom = record.reference_sequence_name();
    if chrom.is_empty() {
        return Err(GenoioError::invalid_source(
            path,
            "vcf record is missing a chromosome id",
        ));
    }
    Ok(chrom)
}

fn text_record_pos(path: &Path, record: &noodles::Record) -> Result<u32> {
    let pos = record
        .variant_start()
        .transpose()
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf position error: {error}")))?
        .ok_or_else(|| GenoioError::invalid_source(path, "vcf record position is missing"))?;
    u32::try_from(pos.get())
        .map_err(|_| GenoioError::invalid_source(path, "vcf record position is out of range"))
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

fn first_alternate_base_string<R>(path: &Path, record: &R) -> Result<Option<String>>
where
    R: noodles::variant::Record + ?Sized,
{
    let alternate_bases = record.alternate_bases();
    let first = match alternate_bases.iter().next() {
        Some(result) => result.map(str::to_string).map(Some).map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf alternate allele error: {error}"))
        }),
        None => Ok(None),
    };
    first
}

fn first_id_str(id: &str) -> &str {
    if id.is_empty() {
        return ".";
    }
    match id.split(';').next() {
        Some(first) if !first.is_empty() => first,
        _ => ".",
    }
}

fn first_alt_allele_str(alt: &str) -> Option<&str> {
    if alt.is_empty() {
        return None;
    }
    match alt.split(',').next() {
        Some(first) if !first.is_empty() => Some(first),
        _ => None,
    }
}

pub(super) fn validate_biallelic_variant<V: VariantMetadataView + ?Sized>(
    path: &Path,
    variant: &V,
) -> Result<()> {
    validate_biallelic_record(path, variant.chrom(), variant.pos(), variant.alt_allele())
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

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use genoio_core::VariantFilter;

    use super::*;

    fn read_text_record(line: &str) -> noodles::Record {
        let mut reader = noodles::io::Reader::new(Cursor::new(line.as_bytes()));
        let mut record = noodles::Record::default();
        reader
            .read_record(&mut record)
            .expect("record should parse");
        record
    }

    #[test]
    fn text_variant_view_borrows_fields_for_filtering_and_metadata_append() {
        let record = read_text_record("chr1\t42\trs1;rs2\tA\tG,T\t30\tPASS\t.\tGT\t0/1\n");
        let view = text_variant_view_from_text_record(Path::new("fixture.vcf"), &record)
            .expect("view should borrow text record metadata");
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "and",
            "left": {"op": "predicate", "name": "chrom", "params": {"value": "chr1"}},
            "right": {"op": "predicate", "name": "qual", "params": {"min": 20.0}}
        }))
        .expect("filter should parse");
        let mut variants = VariantMetadataBuffers::with_capacity(1);

        append_public_variant_metadata_from_text_view(&view, &mut variants)
            .expect("borrowed view should append to metadata buffers");

        assert_eq!(filter.metadata_decision_view(&view), Some(true));
        assert!(skip_variant_for_region(
            &view,
            Some(&RegionPredicate {
                chrom: "chr2".to_string(),
                start: 1,
                end: 100,
            })
        ));
        assert_eq!(variants.positions, vec![42]);
        assert_eq!(variants.chroms.values, b"chr1");
        assert_eq!(variants.ids.values, b"rs1");
        assert_eq!(variants.a0s.values, b"A");
        assert_eq!(variants.a1s.values, b"G");
    }
}
