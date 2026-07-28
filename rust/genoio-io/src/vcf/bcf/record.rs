// pattern: Functional Core

use std::path::Path;

use genoio_core::{GenoioError, VariantMetadataBuffers, VariantMetadataView, VariantStats};
use noodles_bcf as bcf;
use noodles_vcf as noodles;
use noodles_vcf::variant::record::{AlternateBases as _, Ids as _};

use crate::error::Result;

use super::super::finite_qual;

/// Borrowed metadata view over a lazy BCF record.
///
/// BCF matrix loops use this for metadata filters, validation, and metadata append
/// without building an owned `VariantRecord` on the normal retained path.
pub(super) struct BcfVariantView<'a> {
    chrom: &'a str,
    pos: u32,
    ids: bcf::record::Ids<'a>,
    ref_allele: bcf::record::ReferenceBases<'a>,
    alternate_bases: bcf::record::AlternateBases<'a>,
    multi_alt_allele: Option<String>,
    qual: Option<f32>,
}

impl VariantMetadataView for BcfVariantView<'_> {
    fn chrom(&self) -> &str {
        self.chrom
    }

    fn pos(&self) -> u32 {
        self.pos
    }

    fn id(&self) -> &str {
        self.ids.iter().next().unwrap_or(".")
    }

    fn a0(&self) -> &str {
        std::str::from_utf8(self.ref_allele.as_ref()).unwrap_or("")
    }

    fn a1(&self) -> &str {
        match self.alternate_bases.iter().next().transpose() {
            Ok(Some(alt)) => alt,
            _ => "",
        }
    }

    fn ref_allele(&self) -> Option<&str> {
        Some(self.a0())
    }

    fn alt_allele(&self) -> Option<&str> {
        Some(
            self.multi_alt_allele
                .as_deref()
                .unwrap_or_else(|| self.a1()),
        )
    }

    fn qual(&self) -> Option<f32> {
        self.qual
    }
}

pub(super) fn bcf_variant_view_from_record<'a>(
    path: &Path,
    header: &'a noodles::Header,
    record: &'a bcf::Record,
) -> Result<BcfVariantView<'a>> {
    let chrom = <bcf::Record as noodles::variant::Record>::reference_sequence_name(record, header)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf chrom error: {error}")))?;
    if chrom.is_empty() {
        return Err(GenoioError::invalid_source(
            path,
            "bcf record is missing a chromosome id",
        ));
    }

    let pos = record
        .variant_start()
        .transpose()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf position error: {error}")))?
        .ok_or_else(|| GenoioError::invalid_source(path, "bcf record position is missing"))?;
    let pos = u32::try_from(pos.get())
        .map_err(|_| GenoioError::invalid_source(path, "bcf record position is out of range"))?;

    let ref_allele = record.reference_bases();
    std::str::from_utf8(ref_allele.as_ref()).map_err(|error| {
        GenoioError::invalid_source(path, format!("bcf reference allele is not UTF-8: {error}"))
    })?;
    let ids = record.ids();
    let alternate_bases = record.alternate_bases();
    let multi_alt_allele = bcf_multi_alt_allele_string(path, &alternate_bases)?;
    let qual = record
        .quality_score()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf qual error: {error}")))?
        .and_then(finite_qual);

    Ok(BcfVariantView {
        chrom,
        pos,
        ids,
        ref_allele,
        alternate_bases,
        multi_alt_allele,
        qual,
    })
}

pub(super) fn push_bcf_variant_row(
    variants: &mut Option<VariantMetadataBuffers>,
    variant: &BcfVariantView<'_>,
    stats: Option<VariantStats>,
    flipped: bool,
) -> Result<()> {
    let Some(variants) = variants.as_mut() else {
        return Ok(());
    };
    let row_index = variants.len();
    if flipped {
        variants.push_flipped_view(variant)?;
    } else {
        variants.push_view(variant)?;
    }
    if let Some(stats) = stats {
        variants.attach_stats(row_index, stats)?;
    }
    Ok(())
}

fn bcf_multi_alt_allele_string(
    path: &Path,
    alternate_bases: &bcf::record::AlternateBases<'_>,
) -> Result<Option<String>> {
    let mut iter = alternate_bases.iter();
    let first = iter
        .next()
        .transpose()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf alt error: {error}")))?
        .ok_or_else(|| GenoioError::invalid_source(path, "bcf record is missing ALT"))?;
    let Some(second) = iter.next() else {
        return Ok(None);
    };

    let mut joined = String::from(first);
    joined.push(',');
    joined.push_str(
        second.map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf alt error: {error}"))
        })?,
    );
    for alt in iter {
        joined.push(',');
        joined.push_str(alt.map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf alt error: {error}"))
        })?);
    }
    Ok(Some(joined))
}
