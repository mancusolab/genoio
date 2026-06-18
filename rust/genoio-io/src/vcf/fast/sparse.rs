//! Sparse CSC output for the compressed VCF fast path.

use std::io::BufRead;
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, flip_values_to_minor_allele,
    reject_sparse_missing_values, DenseSampleSelection, GenoioError, PartialFilterDecision,
    SparseGenotypeMatrix, VariantFilter, VariantWindow,
};
use noodles_vcf as noodles;

use crate::error::Result;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::gt::{decode_gt_record, GtDecodeBuffers, GtStatsMode};
use super::record::variant_record_from_record;

pub(super) fn read_sparse_records<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
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
            GenoioError::invalid_source(path, format!("vcf fast record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = variant_record_from_record(path, &record)?;
        let partial_decision = variant_filter
            .map(|filter| filter.partial_decision(&variant))
            .unwrap_or(PartialFilterDecision::Accept);
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

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
