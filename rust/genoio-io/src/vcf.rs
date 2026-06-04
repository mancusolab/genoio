// pattern: Imperative Shell

use std::path::{Path, PathBuf};

use genoio_core::{
    append_sparse_column, attach_variant_stats, compute_variant_stats, flip_values_to_minor_allele,
    reject_sparse_missing_values, select_samples_source_order,
    transpose_variant_major_to_sample_major, DenseGenotypeMatrix, MetadataError, MetadataOutput,
    PartialFilterDecision, RegionPredicate, SampleRecord, SourceCapabilities, SparseGenotypeMatrix,
    VariantFilter, VariantRecord, VariantWindow,
};
use rust_htslib::bcf::{
    record::{GenotypeAllele, Numeric},
    IndexedReader, Read, Reader,
};

use crate::error::Result;

/// Read VCF/BCF sample and variant metadata without returning genotypes.
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

/// Read retained VCF/BCF diploid genotypes as a dense sample-by-variant matrix.
pub fn read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_dense_windowed(path, requested_samples, variant_filter, None)
}

/// Read retained VCF/BCF diploid genotypes as dense values over an optional block window.
pub fn read_vcf_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<DenseGenotypeMatrix> {
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        let reader = Reader::from_path(path)
            .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
        return empty_vcf_dense(path, reader.header(), requested_samples);
    }

    // Region filters can be pushed into htslib only when the expression shape
    // is a concrete safe region and the compressed source has an index.
    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            return read_indexed_vcf_dense(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
            );
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    let mut reader = Reader::from_path(path)
        .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
    read_vcf_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        &mut reader,
    )
}

/// Read retained VCF/BCF FORMAT/DS values as dense sample-by-variant dosages.
pub fn read_vcf_dosage_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<DenseGenotypeMatrix> {
    if variant_filter.is_some_and(VariantFilter::requires_genotype_stats) {
        return Err(MetadataError::parse(
            path,
            "dosage-backed genotype reads do not support genotype-stat filters yet",
        ));
    }
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        let reader = Reader::from_path(path)
            .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
        return empty_vcf_dense(path, reader.header(), requested_samples);
    }

    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            return read_indexed_vcf_dosage_dense(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
            );
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    let mut reader = Reader::from_path(path)
        .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
    read_vcf_dosage_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        &mut reader,
    )
}

/// Read retained VCF/BCF diploid genotypes as sparse CSC.
pub fn read_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_vcf_sparse_windowed(path, requested_samples, variant_filter, None)
}

/// Read retained VCF/BCF diploid genotypes as sparse CSC over an optional block window.
pub fn read_vcf_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        let reader = Reader::from_path(path)
            .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
        return empty_vcf_sparse(path, reader.header(), requested_samples);
    }

    // Keep sparse and dense region behavior identical so both paths retain the
    // same variants and fail the same way for unindexed compressed inputs.
    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            return read_indexed_vcf_sparse(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
            );
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    let mut reader = Reader::from_path(path)
        .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
    read_vcf_sparse_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        &mut reader,
    )
}

/// Read phased VCF/BCF diploid genotypes as dense haplotype rows.
pub fn read_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_haplotypes_dense_windowed(path, requested_samples, variant_filter, None)
}

/// Read phased VCF/BCF diploid genotypes as dense haplotype rows over a block window.
pub fn read_vcf_haplotypes_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<DenseGenotypeMatrix> {
    reject_unindexed_compressed_region(path, variant_filter)?;
    let mut reader = Reader::from_path(path)
        .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
    read_vcf_haplotypes_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        &mut reader,
    )
}

/// Read phased VCF/BCF diploid genotypes as sparse haplotype rows.
pub fn read_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_vcf_haplotypes_sparse_windowed(path, requested_samples, variant_filter, None)
}

/// Read phased VCF/BCF diploid genotypes as sparse haplotype rows over a block window.
pub fn read_vcf_haplotypes_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    reject_unindexed_compressed_region(path, variant_filter)?;
    let mut reader = Reader::from_path(path)
        .map_err(|error| MetadataError::parse(path, format!("vcf reader error: {error}")))?;
    read_vcf_haplotypes_sparse_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        &mut reader,
    )
}

fn read_indexed_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
) -> Result<DenseGenotypeMatrix> {
    let mut reader = IndexedReader::from_path(path).map_err(|error| {
        MetadataError::parse(path, format!("indexed vcf reader error: {error}"))
    })?;
    let header = reader.header().clone();
    let rid = match header.name2rid(region.chrom.as_bytes()) {
        Ok(rid) => rid,
        Err(_) => return empty_vcf_dense(path, &header, requested_samples),
    };
    reader
        .fetch(
            rid,
            u64::from(region.start - 1),
            Some(u64::from(region.end - 1)),
        )
        .map_err(|error| MetadataError::parse(path, format!("vcf region fetch error: {error}")))?;

    read_vcf_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        &mut reader,
    )
}

fn read_indexed_vcf_dosage_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
) -> Result<DenseGenotypeMatrix> {
    let mut reader = IndexedReader::from_path(path).map_err(|error| {
        MetadataError::parse(path, format!("indexed vcf reader error: {error}"))
    })?;
    let header = reader.header().clone();
    let rid = match header.name2rid(region.chrom.as_bytes()) {
        Ok(rid) => rid,
        Err(_) => return empty_vcf_dense(path, &header, requested_samples),
    };
    reader
        .fetch(
            rid,
            u64::from(region.start - 1),
            Some(u64::from(region.end - 1)),
        )
        .map_err(|error| MetadataError::parse(path, format!("vcf region fetch error: {error}")))?;

    read_vcf_dosage_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        &mut reader,
    )
}

fn read_indexed_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
) -> Result<SparseGenotypeMatrix> {
    let mut reader = IndexedReader::from_path(path).map_err(|error| {
        MetadataError::parse(path, format!("indexed vcf reader error: {error}"))
    })?;
    let header = reader.header().clone();
    let rid = match header.name2rid(region.chrom.as_bytes()) {
        Ok(rid) => rid,
        Err(_) => return empty_vcf_sparse(path, &header, requested_samples),
    };
    reader
        .fetch(
            rid,
            u64::from(region.start - 1),
            Some(u64::from(region.end - 1)),
        )
        .map_err(|error| MetadataError::parse(path, format!("vcf region fetch error: {error}")))?;

    read_vcf_sparse_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        &mut reader,
    )
}

fn read_vcf_dense_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    reader: &mut R,
) -> Result<DenseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
    let mut retained_index = 0_usize;
    for record_result in reader.records() {
        // Check before pulling another record; otherwise block reads still pay
        // to scan one extra variant after the requested retained window.
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            break;
        }
        let record = record_result
            .map_err(|error| MetadataError::parse(path, format!("vcf record error: {error}")))?;
        let mut variant = variant_record_from_record(path, &header, &record)?;
        diagnostics.candidate_variants += 1;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                    break;
                }
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {}
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        validate_dense_biallelic_record(path, &record)?;
        let genotypes = record
            .format(b"GT")
            .integer()
            .map_err(|error| MetadataError::parse(path, format!("vcf genotype error: {error}")))?;
        let mut current_values = Vec::with_capacity(selection.source_indices.len());
        let mut current_missing = Vec::with_capacity(selection.source_indices.len());
        for source_index in &selection.source_indices {
            let (value, is_missing) =
                decode_raw_diploid_gt(path, &record, genotypes[*source_index])?;
            current_values.push(value);
            current_missing.push(is_missing);
        }

        let stats = if needs_genotype_decision {
            Some(compute_variant_stats(&current_values, &current_missing)?)
        } else {
            None
        };
        if needs_genotype_decision
            && variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref()))
        {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if needs_genotype_decision {
            if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                break;
            }
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                continue;
            }
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

fn read_vcf_dosage_dense_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    reader: &mut R,
) -> Result<DenseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
    let mut retained_index = 0_usize;
    for record_result in reader.records() {
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            break;
        }
        let record = record_result
            .map_err(|error| MetadataError::parse(path, format!("vcf record error: {error}")))?;
        let variant = variant_record_from_record(path, &header, &record)?;
        diagnostics.candidate_variants += 1;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {
                return Err(MetadataError::parse(
                    path,
                    "dosage-backed genotype reads do not support genotype-stat filters yet",
                ));
            }
        }

        validate_dense_biallelic_record(path, &record)?;
        let decoded = decode_ds_dosage_record(path, &record, &selection.source_indices)?;
        variants.push(variant);
        variant_major_values.extend(decoded.values);
        variant_major_missing.extend(decoded.missing);
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

fn read_vcf_haplotypes_dense_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    reader: &mut R,
) -> Result<DenseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
    let mut retained_index = 0_usize;
    for record_result in reader.records() {
        // Haplotype reads use the same retained-window semantics as genotype
        // reads, but each retained sample contributes two output rows.
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            break;
        }
        let record = record_result
            .map_err(|error| MetadataError::parse(path, format!("vcf record error: {error}")))?;
        let mut variant = variant_record_from_record(path, &header, &record)?;
        diagnostics.candidate_variants += 1;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                    break;
                }
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {}
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        validate_dense_biallelic_record(path, &record)?;
        let stats = if needs_genotype_decision {
            let decoded = decode_diploid_genotype_record(path, &record, &selection.source_indices)?;
            Some(compute_variant_stats(&decoded.values, &decoded.missing)?)
        } else {
            None
        };
        if needs_genotype_decision
            && variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref()))
        {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if needs_genotype_decision {
            if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                break;
            }
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                continue;
            }
        }
        let decoded = decode_phased_haplotype_record(path, &record, &selection.source_indices)?;
        variants.push(variant);
        variant_major_values.extend(decoded.haplotype_values);
        variant_major_missing.extend(decoded.haplotype_missing);
    }

    let samples = haplotype_sample_records(&selection.samples, &selection.source_indices);
    let n_samples = samples.len();
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
        samples,
        variants,
        diagnostics,
    )
}

fn read_vcf_haplotypes_sparse_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    reader: &mut R,
) -> Result<SparseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let samples = haplotype_sample_records(&selection.samples, &selection.source_indices);
    let n_samples = samples.len();
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();
    let mut retained_index = 0_usize;
    for record_result in reader.records() {
        // Stop before reading the next record so sparse blocks do not scan the
        // remainder of large VCF/BCF files after the block is filled.
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            break;
        }
        let record = record_result
            .map_err(|error| MetadataError::parse(path, format!("vcf record error: {error}")))?;
        let mut variant = variant_record_from_record(path, &header, &record)?;
        diagnostics.candidate_variants += 1;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {}
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        validate_dense_biallelic_record(path, &record)?;
        let stats = if needs_genotype_decision {
            let decoded = decode_diploid_genotype_record(path, &record, &selection.source_indices)?;
            Some(compute_variant_stats(&decoded.values, &decoded.missing)?)
        } else {
            None
        };
        if needs_genotype_decision
            && variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref()))
        {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if needs_genotype_decision {
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                continue;
            }
        }
        let decoded = decode_phased_haplotype_record(path, &record, &selection.source_indices)?;
        reject_sparse_missing_values(&decoded.haplotype_missing)?;
        append_sparse_column(
            &mut indptr,
            &mut indices,
            &mut data,
            &decoded.haplotype_values,
        );
        variants.push(variant);
    }

    let n_variants = variants.len();
    diagnostics.retained_variants = n_variants;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        samples,
        variants,
        diagnostics,
    )
}

fn read_vcf_sparse_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    reader: &mut R,
) -> Result<SparseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let n_samples = selection.samples.len();
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();
    let mut retained_index = 0_usize;
    for record_result in reader.records() {
        // Stop before reading the next record so sparse blocks do not scan the
        // remainder of large VCF/BCF files after the block is filled.
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            break;
        }
        let record = record_result
            .map_err(|error| MetadataError::parse(path, format!("vcf record error: {error}")))?;
        let mut variant = variant_record_from_record(path, &header, &record)?;
        diagnostics.candidate_variants += 1;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {}
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        validate_dense_biallelic_record(path, &record)?;
        let genotypes = record
            .format(b"GT")
            .integer()
            .map_err(|error| MetadataError::parse(path, format!("vcf genotype error: {error}")))?;
        let mut current_values = Vec::with_capacity(selection.source_indices.len());
        let mut current_missing = Vec::with_capacity(selection.source_indices.len());
        for source_index in &selection.source_indices {
            let (value, is_missing) =
                decode_raw_diploid_gt(path, &record, genotypes[*source_index])?;
            current_values.push(value);
            current_missing.push(is_missing);
        }

        let stats = if needs_genotype_decision {
            Some(compute_variant_stats(&current_values, &current_missing)?)
        } else {
            None
        };
        if needs_genotype_decision
            && variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref()))
        {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if needs_genotype_decision {
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                continue;
            }
        }
        reject_sparse_missing_values(&current_missing)?;
        flip_values_to_minor_allele(&mut current_values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &current_values);
        variants.push(variant);
    }

    let n_variants = variants.len();
    diagnostics.retained_variants = n_variants;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        selection.samples,
        variants,
        diagnostics,
    )
}

fn empty_vcf_dense(
    path: &Path,
    header: &rust_htslib::bcf::header::HeaderView,
    requested_samples: Option<&[String]>,
) -> Result<DenseGenotypeMatrix> {
    let all_samples = sample_records_from_header(header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    DenseGenotypeMatrix::new(
        selection.samples.len(),
        0,
        Vec::new(),
        Vec::new(),
        selection.samples,
        Vec::new(),
        diagnostics,
    )
}

fn empty_vcf_sparse(
    path: &Path,
    header: &rust_htslib::bcf::header::HeaderView,
    requested_samples: Option<&[String]>,
) -> Result<SparseGenotypeMatrix> {
    let all_samples = sample_records_from_header(header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    SparseGenotypeMatrix::new(
        selection.samples.len(),
        0,
        vec![0],
        Vec::new(),
        Vec::new(),
        selection.samples,
        Vec::new(),
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
            source_sample_index: None,
            haplotype_index: None,
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
        qual: finite_qual(record.qual()),
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
}

fn finite_qual(qual: f32) -> Option<f32> {
    qual.is_finite().then_some(qual)
}

fn reject_unindexed_compressed_region(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
) -> Result<()> {
    if variant_filter
        .and_then(VariantFilter::concrete_region_pushdown)
        .is_none()
        || !is_compressed_vcf(path)
    {
        return Ok(());
    }
    if has_vcf_index(path) {
        return Ok(());
    }
    Err(MetadataError::parse(
        path,
        "region filter on compressed VCF requires an index",
    ))
}

fn has_vcf_index(path: &Path) -> bool {
    companion_index_path(path, "tbi").exists() || companion_index_path(path, "csi").exists()
}

fn companion_index_path(path: &Path, index_extension: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}", path.to_string_lossy(), index_extension))
}

fn is_compressed_vcf(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| matches!(extension, "gz" | "bgz"))
}

fn decode_raw_diploid_gt(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    genotype: &[i32],
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
    for encoded in genotype {
        match decode_raw_gt_allele(*encoded) {
            RawGtAllele::Missing => return Ok((0.0, true)),
            RawGtAllele::Reference => {}
            RawGtAllele::Alternate => dosage += 1.0,
            RawGtAllele::Unsupported(other) => {
                return Err(MetadataError::parse(
                    path,
                    format!(
                        "vcf record {} has multiallelic GT allele index {other}",
                        String::from_utf8_lossy(&record.id())
                    ),
                ));
            }
        }
    }

    Ok((dosage, false))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RawGtAllele {
    Missing,
    Reference,
    Alternate,
    Unsupported(i32),
}

fn decode_raw_gt_allele(encoded: i32) -> RawGtAllele {
    match encoded {
        0 | 1 => RawGtAllele::Missing,
        value if value > 1 => match (value >> 1) - 1 {
            0 => RawGtAllele::Reference,
            1 => RawGtAllele::Alternate,
            allele => RawGtAllele::Unsupported(allele),
        },
        value => RawGtAllele::Unsupported(value),
    }
}

struct DecodedDiploidGenotypeRecord {
    values: Vec<f32>,
    missing: Vec<bool>,
}

fn decode_diploid_genotype_record(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    source_indices: &[usize],
) -> Result<DecodedDiploidGenotypeRecord> {
    let genotypes = record
        .format(b"GT")
        .integer()
        .map_err(|error| MetadataError::parse(path, format!("vcf genotype error: {error}")))?;
    let mut values = Vec::with_capacity(source_indices.len());
    let mut missing = Vec::with_capacity(source_indices.len());

    for source_index in source_indices {
        let (value, is_missing) = decode_raw_diploid_gt(path, record, genotypes[*source_index])?;
        values.push(value);
        missing.push(is_missing);
    }

    Ok(DecodedDiploidGenotypeRecord { values, missing })
}

fn decode_ds_dosage_record(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    source_indices: &[usize],
) -> Result<DecodedDiploidGenotypeRecord> {
    let dosages = record.format(b"DS").float().map_err(|error| {
        MetadataError::parse(
            path,
            format!("vcf dosage reads require FORMAT/DS values: {error}"),
        )
    })?;
    let mut values = Vec::with_capacity(source_indices.len());
    let mut missing = Vec::with_capacity(source_indices.len());

    for source_index in source_indices {
        let dosage = dosages[*source_index];
        if dosage.len() != 1 {
            return Err(MetadataError::parse(
                path,
                format!(
                    "vcf record {} has FORMAT/DS with {} values for a sample; expected one",
                    String::from_utf8_lossy(&record.id()),
                    dosage.len()
                ),
            ));
        }
        let value = dosage[0];
        if value.is_missing() {
            values.push(0.0);
            missing.push(true);
            continue;
        }
        if !value.is_finite() || !(0.0..=2.0).contains(&value) {
            return Err(MetadataError::parse(
                path,
                format!(
                    "vcf record {} has invalid FORMAT/DS value {value}; expected finite value in [0, 2]",
                    String::from_utf8_lossy(&record.id())
                ),
            ));
        }
        values.push(value);
        missing.push(false);
    }

    Ok(DecodedDiploidGenotypeRecord { values, missing })
}

struct PhasedHaplotypeRecord {
    haplotype_values: Vec<f32>,
    haplotype_missing: Vec<bool>,
}

fn decode_phased_haplotype_record(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    source_indices: &[usize],
) -> Result<PhasedHaplotypeRecord> {
    let genotypes = record
        .genotypes()
        .map_err(|error| MetadataError::parse(path, format!("vcf genotype error: {error}")))?;
    let mut haplotype_values = Vec::with_capacity(source_indices.len() * 2);
    let mut haplotype_missing = Vec::with_capacity(source_indices.len() * 2);

    for source_index in source_indices {
        let genotype = genotypes.get(*source_index);
        let (values, missing) = decode_phased_diploid_gt(path, record, &genotype)?;
        haplotype_values.extend(values);
        haplotype_missing.extend(missing);
    }

    Ok(PhasedHaplotypeRecord {
        haplotype_values,
        haplotype_missing,
    })
}

fn decode_phased_diploid_gt(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    genotype: &rust_htslib::bcf::record::Genotype,
) -> Result<([f32; 2], [bool; 2])> {
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

    let mut values = [0.0, 0.0];
    let mut missing = [false, false];
    for (allele_index, allele) in genotype.iter().enumerate() {
        if allele_index > 0
            && matches!(
                allele,
                GenotypeAllele::UnphasedMissing | GenotypeAllele::Unphased(_)
            )
        {
            return Err(MetadataError::parse(
                path,
                format!(
                    "vcf record {} contains an unphased GT separator in a retained haplotype variant",
                    String::from_utf8_lossy(&record.id())
                ),
            ));
        }
        match allele {
            GenotypeAllele::PhasedMissing | GenotypeAllele::UnphasedMissing => {
                missing[allele_index] = true;
            }
            GenotypeAllele::Phased(index) | GenotypeAllele::Unphased(index) => match index {
                0 => {}
                1 => values[allele_index] = 1.0,
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

    Ok((values, missing))
}

fn haplotype_sample_records(
    samples: &[SampleRecord],
    source_indices: &[usize],
) -> Vec<SampleRecord> {
    let mut haplotype_samples = Vec::with_capacity(samples.len() * 2);
    for (sample, source_index) in samples.iter().zip(source_indices) {
        for haplotype_index in 0..2 {
            let mut haplotype_sample = sample.clone();
            haplotype_sample.source_sample_index = Some(*source_index);
            haplotype_sample.haplotype_index = Some(haplotype_index);
            haplotype_samples.push(haplotype_sample);
        }
    }
    haplotype_samples
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
