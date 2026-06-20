//! Sparse CSC output for the text VCF backend.
//!
//! Sparse text reads decode selected GT or haplotype tokens and append retained
//! variants directly to CSC buffers. Missing retained calls are rejected because
//! the sparse matrix format stores values only.

// pattern: Mixed (unavoidable)
// Reason: This hot path combines lazy VCF record IO with direct CSC emission to
// avoid dense per-record staging.

use std::io::BufRead;
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, flip_values_to_minor_allele,
    flip_variant_metadata_to_minor_allele, reject_sparse_missing, reject_sparse_missing_values,
    should_flip_haplotype_to_minor_allele, DenseSampleSelection, GenoioError,
    PartialFilterDecision, RegionPredicate, SparseGenotypeMatrix, VariantFilter, VariantRecord,
    VariantWindow,
};
use noodles_vcf as noodles;

use crate::error::Result;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::super::haplotype_sample_records;
use super::gt::{
    decode_gt_record, decode_phased_gt_sparse_record, GtDecodeBuffers, GtStatsMode,
    HaplotypeSparseDecodeBuffers,
};
use super::record::{
    metadata_variant_record_from_record, skip_variant_for_region, validate_biallelic_variant,
};

pub(super) fn read_sparse_records<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    source_region: Option<&RegionPredicate>,
    selection: DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<SparseGenotypeMatrix> {
    let DenseSampleSelection {
        source_indices,
        samples,
        mut diagnostics,
    } = selection;
    let n_samples = samples.len();
    let variant_capacity = variant_window.map_or(0, |window| window.len);
    let mut indptr = Vec::with_capacity(variant_capacity.saturating_add(1));
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::with_capacity(variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = noodles::Record::default();
    let mut decoded = GtDecodeBuffers::with_capacity(source_indices.len());

    loop {
        if retention.window_is_satisfied() {
            break;
        }
        if reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("text VCF record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = metadata_variant_record_from_record(path, &record)?;
        if skip_variant_for_region(&variant, source_region) {
            continue;
        }
        let partial_decision = variant_filter
            .map(|filter| filter.partial_decision(&variant))
            .unwrap_or(PartialFilterDecision::Accept);
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }
        validate_biallelic_variant(path, &variant)?;

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        // Decode into reusable dense scratch first. CSC output still needs the
        // dense column briefly for missing-value rejection and minor-allele
        // flipping to preserve the public sparse contract.
        decode_gt_record(
            path,
            &record,
            &source_indices,
            GtStatsMode::from_needed(needs_genotype_decision),
            &mut decoded,
        )?;

        if needs_genotype_decision {
            let stats = decoded.stats();
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if let Some(stats) = stats {
                attach_variant_stats(&mut variant, stats);
            }
        }

        reject_sparse_missing_values(decoded.missing())?;
        // Genotype sparse columns store minor-allele dosage by convention.
        flip_values_to_minor_allele(decoded.values_mut(), &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, decoded.values());
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

pub(super) fn read_haplotype_sparse_records<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    source_region: Option<&RegionPredicate>,
    selection: DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<SparseGenotypeMatrix> {
    let DenseSampleSelection {
        source_indices,
        samples,
        mut diagnostics,
    } = selection;
    let n_samples = samples.len() * 2;
    let samples = haplotype_sample_records(&samples, &source_indices);
    let variant_capacity = variant_window.map_or(0, |window| window.len);
    let mut indptr = Vec::with_capacity(variant_capacity.saturating_add(1));
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::with_capacity(variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = noodles::Record::default();
    let mut decoded = HaplotypeSparseDecodeBuffers::with_capacity(source_indices.len());
    let mut stats_decoded = GtDecodeBuffers::with_capacity(source_indices.len());

    loop {
        if retention.window_is_satisfied() {
            break;
        }
        if reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("text VCF record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = metadata_variant_record_from_record(path, &record)?;
        if skip_variant_for_region(&variant, source_region) {
            continue;
        }
        let partial_decision = variant_filter
            .map(|filter| filter.partial_decision(&variant))
            .unwrap_or(PartialFilterDecision::Accept);
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }
        validate_biallelic_variant(path, &variant)?;

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        if needs_genotype_decision {
            // Apply genotype-stat filters before phased decoding so rejected
            // unphased records do not fail a haplotype sparse read.
            decode_gt_record(
                path,
                &record,
                &source_indices,
                GtStatsMode::Compute,
                &mut stats_decoded,
            )?;
            let stats = stats_decoded.stats();
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if let Some(stats) = stats {
                attach_variant_stats(&mut variant, stats);
            }
        }

        decode_phased_gt_sparse_record(
            path,
            &record,
            &source_indices,
            GtStatsMode::Skip,
            &mut decoded,
        )?;
        reject_sparse_missing(decoded.has_missing())?;
        append_haplotype_minor_sparse_column(
            &mut indptr,
            &mut indices,
            &mut data,
            &mut variant,
            &decoded,
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

fn append_haplotype_minor_sparse_column(
    indptr: &mut Vec<usize>,
    indices: &mut Vec<usize>,
    data: &mut Vec<f32>,
    variant: &mut VariantRecord,
    decoded: &HaplotypeSparseDecodeBuffers,
) {
    let a1_rows = decoded.a1_rows();
    if should_flip_haplotype_to_minor_allele(a1_rows.len(), decoded.n_rows()) {
        flip_variant_metadata_to_minor_allele(variant);
        // `a1_rows` is sorted because selected samples are scanned in source
        // order and each phased call contributes rows in haplotype order.
        append_haplotype_complement_column(indices, data, decoded.n_rows(), a1_rows);
    } else {
        append_haplotype_rows(indices, data, a1_rows);
    }
    indptr.push(indices.len());
}

fn append_haplotype_rows(indices: &mut Vec<usize>, data: &mut Vec<f32>, rows: &[usize]) {
    indices.extend_from_slice(rows);
    data.extend(std::iter::repeat_n(1.0, rows.len()));
}

fn append_haplotype_complement_column(
    indices: &mut Vec<usize>,
    data: &mut Vec<f32>,
    n_rows: usize,
    a1_rows: &[usize],
) {
    let mut next_a1 = 0_usize;
    for row in 0..n_rows {
        if next_a1 < a1_rows.len() && a1_rows[next_a1] == row {
            next_a1 += 1;
        } else {
            indices.push(row);
            data.push(1.0);
        }
    }
}
