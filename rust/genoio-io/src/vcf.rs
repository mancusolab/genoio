// pattern: Imperative Shell

use std::path::Path;

use genoio_core::{
    attach_variant_stats, compute_variant_stats, select_samples_source_order,
    transpose_variant_major_to_sample_major, DenseGenotypeMatrix, MetadataError, MetadataOutput,
    SampleRecord, SourceCapabilities, VariantFilter, VariantRecord,
};
use rust_htslib::bcf::{record::GenotypeAllele, Read, Reader};

use crate::error::Result;

pub fn read_vcf_metadata(path: &Path) -> Result<MetadataOutput> {
    let mut reader = Reader::from_path(path)
        .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
    let header = reader.header().clone();
    let samples = sample_records_from_header(&header);

    let mut variants = Vec::new();
    let mut has_phased_genotype_evidence = false;
    for record_result in reader.records() {
        let record = record_result
            .map_err(|error| MetadataError::parse(path, format!("vcf record error: {error}")))?;
        if !has_phased_genotype_evidence && record_has_phased_genotype(&record, samples.len())? {
            has_phased_genotype_evidence = true;
        }
        variants.push(variant_record_from_record(path, &header, &record)?);
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

pub fn read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    reject_unindexed_compressed_region(path, variant_filter)?;
    let mut reader = Reader::from_path(path)
        .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
    for record_result in reader.records() {
        let record = record_result
            .map_err(|error| MetadataError::parse(path, format!("vcf record error: {error}")))?;
        let mut variant = variant_record_from_record(path, &header, &record)?;
        diagnostics.candidate_variants += 1;
        if variant_filter.and_then(|filter| filter.metadata_decision(&variant)) == Some(false) {
            diagnostics.dropped_metadata_variants += 1;
            continue;
        }

        validate_dense_biallelic_record(path, &record)?;
        let genotypes = record
            .genotypes()
            .map_err(|error| MetadataError::parse(path, format!("vcf genotype error: {error}")))?;
        let mut current_values = Vec::with_capacity(selection.source_indices.len());
        let mut current_missing = Vec::with_capacity(selection.source_indices.len());
        for source_index in &selection.source_indices {
            let (value, is_missing) =
                decode_diploid_gt(path, &record, &genotypes.get(*source_index))?;
            current_values.push(value);
            current_missing.push(is_missing);
        }

        let stats = if variant_filter.is_some_and(VariantFilter::requires_genotype_stats) {
            Some(compute_variant_stats(&current_values, &current_missing)?)
        } else {
            None
        };
        if variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref())) {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        variants.push(variant);
        variant_major_values.extend(current_values);
        variant_major_missing.extend(current_missing);
    }

    let n_samples = selection.samples.len();
    let n_variants = variants.len();
    diagnostics.retained_variants = n_variants;
    let values =
        transpose_variant_major_to_sample_major(&variant_major_values, n_samples, n_variants);
    let missing_mask =
        transpose_variant_major_to_sample_major(&variant_major_missing, n_samples, n_variants);

    DenseGenotypeMatrix::new(
        n_samples,
        n_variants,
        values,
        missing_mask,
        selection.samples,
        variants,
        diagnostics,
    )
}

fn validate_dense_biallelic_record(path: &Path, record: &rust_htslib::bcf::Record) -> Result<()> {
    let allele_count = record.alleles().len();
    if allele_count == 2 {
        return Ok(());
    }
    let record_id = record.id();
    let id = String::from_utf8_lossy(&record_id);
    let reason = if allele_count > 2 {
        "multi-ALT records are not supported"
    } else {
        "records with fewer than two alleles are not supported"
    };
    Err(MetadataError::parse(
        path,
        format!(
            "vcf dense reads require biallelic records; record {id} has {allele_count} alleles: {reason}"
        ),
    ))
}

fn sample_records_from_header(header: &rust_htslib::bcf::header::HeaderView) -> Vec<SampleRecord> {
    header
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
        .collect()
}

fn variant_record_from_record(
    path: &Path,
    header: &rust_htslib::bcf::header::HeaderView,
    record: &rust_htslib::bcf::Record,
) -> Result<VariantRecord> {
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

    Ok(VariantRecord {
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
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
}

fn reject_unindexed_compressed_region(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
) -> Result<()> {
    if !variant_filter.is_some_and(VariantFilter::has_region_predicate) || !is_compressed_vcf(path)
    {
        return Ok(());
    }
    let tbi = path.with_extension(format!(
        "{}.tbi",
        path.extension()
            .and_then(|extension| extension.to_str())
            .unwrap_or("")
    ));
    let csi = path.with_extension(format!(
        "{}.csi",
        path.extension()
            .and_then(|extension| extension.to_str())
            .unwrap_or("")
    ));
    if tbi.exists() || csi.exists() {
        return Ok(());
    }
    Err(MetadataError::parse(
        path,
        "region filter on compressed VCF requires an index",
    ))
}

fn is_compressed_vcf(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| matches!(extension, "gz" | "bgz"))
}

fn decode_diploid_gt(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    genotype: &rust_htslib::bcf::record::Genotype,
) -> Result<(f32, bool)> {
    if genotype.len() != 2 {
        return Err(MetadataError::parse(
            path,
            format!(
                "vcf record {} has non-diploid GT with {} alleles",
                String::from_utf8_lossy(&record.id()),
                genotype.len()
            ),
        ));
    }

    let mut dosage = 0.0;
    for allele in genotype.iter() {
        match allele {
            GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing => {
                return Ok((0.0, true));
            }
            GenotypeAllele::Unphased(index) | GenotypeAllele::Phased(index) => match index {
                0 => {}
                1 => dosage += 1.0,
                other => {
                    return Err(MetadataError::parse(
                        path,
                        format!(
                            "vcf record {} has multiallelic GT allele index {other}",
                            String::from_utf8_lossy(&record.id())
                        ),
                    ));
                }
            },
        }
    }

    Ok((dosage, false))
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
