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
    append_sparse_column, append_sparse_value, finish_sparse_column, reject_sparse_missing,
    should_flip_haplotype_to_minor_allele, DenseSampleSelection, GenoioError, RegionPredicate,
    SampleMetadataBuffers, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};
use noodles_vcf as noodles;

use crate::error::Result;
use crate::retention::{RetainedVariantState, RetentionAction};

use super::super::haplotype_sample_records;
use super::gt::{
    decode_gt_record, decode_phased_gt_sparse_record, GtDecodeBuffers, GtStatsMode,
    HaplotypeSparseDecodeBuffers,
};
use super::record::{prepare_text_candidate, TextCandidateAction};
use super::{
    dense_output_variant_capacity, VariantMetadataSink, VariantMetadataSinkKind, VcfMetadataReturn,
};

#[expect(
    clippy::too_many_arguments,
    reason = "sparse VCF loop receives prevalidated output mode, selection, and reader state"
)]
pub(super) fn read_sparse_records_with_metadata<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    source_region: Option<&RegionPredicate>,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    selection: DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<SparseGenotypeMatrix> {
    let DenseSampleSelection {
        source_indices,
        samples,
        mut diagnostics,
    } = selection;
    let n_samples = samples.len();
    let variant_capacity = dense_output_variant_capacity(variant_window);
    let mut indptr = Vec::with_capacity(variant_capacity.saturating_add(1));
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = VariantMetadataSink::new(variant_sink_kind, variant_capacity);
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

        let prepared = match prepare_text_candidate(
            path,
            &record,
            source_region,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            TextCandidateAction::Skip => continue,
            TextCandidateAction::Stop => break,
            TextCandidateAction::Decode(prepared) => prepared,
        };
        let variant = prepared.variant;
        let needs_genotype_decision = prepared.needs_genotype_decision;
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

        let mut stats_to_attach = None;
        if needs_genotype_decision {
            let stats = decoded.stats();
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate_view(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            stats_to_attach = stats;
        }

        reject_sparse_missing(!decoded.missing_indices().is_empty())?;
        // Genotype sparse columns store minor-allele dosage by convention.
        let flipped = flip_values_to_minor_allele(decoded.values_mut());
        append_sparse_column(&mut indptr, &mut indices, &mut data, decoded.values())?;
        variants.push_view_with_optional_stats_and_orientation(
            &variant,
            stats_to_attach,
            flipped,
        )?;
    }

    let n_variants = indptr.len().saturating_sub(1);
    diagnostics.retained_variants = n_variants;
    let samples =
        SampleMetadataBuffers::optional_from_records(&samples, metadata_return.samples, false)?;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        samples,
        variants.into_output()?,
        diagnostics,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "sparse VCF loop receives prevalidated output mode, selection, and reader state"
)]
pub(super) fn read_haplotype_sparse_records_with_metadata<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    source_region: Option<&RegionPredicate>,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    selection: DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<SparseGenotypeMatrix> {
    let DenseSampleSelection {
        source_indices,
        samples,
        mut diagnostics,
    } = selection;
    let n_samples = samples.len() * 2;
    let haplotype_samples = haplotype_sample_records(&samples, &source_indices);
    let output_samples = SampleMetadataBuffers::optional_from_records(
        &haplotype_samples,
        metadata_return.samples,
        true,
    )?;
    let variant_capacity = dense_output_variant_capacity(variant_window);
    let mut indptr = Vec::with_capacity(variant_capacity.saturating_add(1));
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = VariantMetadataSink::new(variant_sink_kind, variant_capacity);
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

        let prepared = match prepare_text_candidate(
            path,
            &record,
            source_region,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            TextCandidateAction::Skip => continue,
            TextCandidateAction::Stop => break,
            TextCandidateAction::Decode(prepared) => prepared,
        };
        let variant = prepared.variant;
        let needs_genotype_decision = prepared.needs_genotype_decision;
        let mut stats_to_attach = None;
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
                variant_filter.is_none_or(|filter| filter.evaluate_view(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            stats_to_attach = stats;
        }

        decode_phased_gt_sparse_record(
            path,
            &record,
            &source_indices,
            GtStatsMode::Skip,
            &mut decoded,
        )?;
        reject_sparse_missing(decoded.has_missing())?;
        let flipped =
            append_haplotype_minor_sparse_column(&mut indptr, &mut indices, &mut data, &decoded)?;
        variants.push_view_with_optional_stats_and_orientation(
            &variant,
            stats_to_attach,
            flipped,
        )?;
    }

    let n_variants = indptr.len().saturating_sub(1);
    diagnostics.retained_variants = n_variants;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        output_samples,
        variants.into_output()?,
        diagnostics,
    )
}

/// Append a haplotype sparse column and report whether metadata must be flipped.
pub(super) fn append_haplotype_minor_sparse_column(
    indptr: &mut Vec<i32>,
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
    decoded: &HaplotypeSparseDecodeBuffers,
) -> Result<bool> {
    let a1_rows = decoded.a1_rows();
    let flipped = should_flip_haplotype_to_minor_allele(a1_rows.len(), decoded.n_rows());
    if flipped {
        // `a1_rows` is sorted because selected samples are scanned in source
        // order and each phased call contributes rows in haplotype order.
        append_haplotype_complement_column(indices, data, decoded.n_rows(), a1_rows)?;
    } else {
        append_haplotype_rows(indices, data, a1_rows)?;
    }
    finish_sparse_column(indptr, data.len())?;
    Ok(flipped)
}

/// Flip genotype dosages in-place when allele 1 is the major allele.
pub(super) fn flip_values_to_minor_allele(values: &mut [f32]) -> bool {
    let a1_count = values.iter().sum::<f32>();
    let a0_count = 2.0 * values.len() as f32 - a1_count;
    if a1_count <= a0_count {
        return false;
    }
    for value in values {
        *value = 2.0 - *value;
    }
    true
}

fn append_haplotype_rows(
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
    rows: &[usize],
) -> Result<()> {
    for &row in rows {
        append_sparse_value(indices, data, row, 1.0)?;
    }
    Ok(())
}

fn append_haplotype_complement_column(
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
    n_rows: usize,
    a1_rows: &[usize],
) -> Result<()> {
    let mut next_a1 = 0_usize;
    for row in 0..n_rows {
        if next_a1 < a1_rows.len() && a1_rows[next_a1] == row {
            next_a1 += 1;
        } else {
            append_sparse_value(indices, data, row, 1.0)?;
        }
    }
    Ok(())
}
