// pattern: Imperative Shell

use std::path::Path;

use genoio_core::{
    MetadataError, MetadataOutput, SampleRecord, SourceCapabilities, VariantRecord,
};
use rust_htslib::bcf::{record::GenotypeAllele, Read, Reader};

use crate::error::Result;

pub fn read_vcf_metadata(path: &Path) -> Result<MetadataOutput> {
    let mut reader = Reader::from_path(path)
        .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
    let header = reader.header().clone();
    let samples = header
        .samples()
        .iter()
        .map(|sample| SampleRecord {
            fid: None,
            iid: String::from_utf8_lossy(sample).into_owned(),
            father: None,
            mother: None,
            sex: None,
            phenotype: None,
        })
        .collect::<Vec<_>>();

    let mut variants = Vec::new();
    let mut has_phased_genotype_evidence = false;
    for record_result in reader.records() {
        let record = record_result
            .map_err(|error| MetadataError::parse(path, format!("vcf record error: {error}")))?;
        if !has_phased_genotype_evidence && record_has_phased_genotype(&record, samples.len())? {
            has_phased_genotype_evidence = true;
        }
        let rid = record
            .rid()
            .ok_or_else(|| MetadataError::parse(path, "vcf record is missing a chromosome id"))?;
        let chrom = String::from_utf8_lossy(
            header
                .rid2name(rid)
                .map_err(|error| MetadataError::parse(path, format!("vcf rid error: {error}")))?,
        )
        .into_owned();
        let pos = u32::try_from(record.pos() + 1)
            .map_err(|_| MetadataError::parse(path, "vcf record position is out of range"))?;
        let id = String::from_utf8_lossy(&record.id()).into_owned();
        let alleles = record.alleles();
        if alleles.len() < 2 {
            return Err(MetadataError::parse(
                path,
                format!("vcf record {chrom}:{pos} has fewer than two alleles"),
            ));
        }
        let ref_allele = String::from_utf8_lossy(alleles[0]).into_owned();
        let alt_values = alleles[1..]
            .iter()
            .map(|allele| String::from_utf8_lossy(allele).into_owned())
            .collect::<Vec<_>>();
        let first_alt = alt_values[0].clone();
        let alt_allele = alt_values.join(",");

        variants.push(VariantRecord {
            chrom,
            pos,
            id,
            a0: ref_allele.clone(),
            a1: first_alt.clone(),
            ref_allele: Some(ref_allele.clone()),
            alt_allele: Some(alt_allele),
            source_a0: ref_allele,
            source_a1: first_alt,
            flipped: false,
        });
    }

    let capabilities = if has_phased_genotype_evidence {
        SourceCapabilities::phased_genotypes()
    } else {
        SourceCapabilities::genotype_only()
    };

    Ok(MetadataOutput {
        samples,
        variants,
        capabilities,
    })
}

fn record_has_phased_genotype(
    record: &rust_htslib::bcf::Record,
    sample_count: usize,
) -> Result<bool> {
    let genotypes = match record.genotypes() {
        Ok(genotypes) => genotypes,
        Err(_) => return Ok(false),
    };
    Ok((0..sample_count).any(|sample_index| {
        genotypes.get(sample_index).iter().any(|allele| {
            matches!(
                allele,
                GenotypeAllele::Phased(_) | GenotypeAllele::PhasedMissing
            )
        })
    }))
}
