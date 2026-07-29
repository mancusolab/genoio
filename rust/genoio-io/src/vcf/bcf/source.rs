//! Lazy BCF readers built on noodles-bcf records.
//!
//! BCF exposes typed sample values, so this module deliberately does not reuse
//! the text VCF byte scanners. The record loop mirrors the public VCF/BCF
//! semantics while keeping one lazy `bcf::Record` buffer alive across variants.

// pattern: Mixed
// Reason: BCF setup, lazy record iteration, and decode routing share ownership
// of the same reusable record buffer.

use std::fs::File;
use std::path::Path;

use genoio_core::{
    append_sparse_column, reject_sparse_missing, select_samples_source_order, DenseGenotypeMatrix,
    DenseLayout, DenseMissingPolicy, DenseSampleSelection, GenoioError, MetadataOutput,
    SampleMetadataBuffers, SourceCapabilities, SparseGenotypeMatrix, VariantFilter,
    VariantMetadataBuffers, VariantMetadataView, VariantStats, VariantWindow,
};
use noodles_bcf as bcf;
use noodles_vcf as noodles;

use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::hardcall::evaluate_hardcall_counts_filter;
use crate::matrix::apply_dense_missing_policy_to_variant;
use crate::retention::{RetainedVariantState, RetentionAction};

use super::super::text::append_public_variant_metadata_from_noodles_variant_record;
use super::super::{
    haplotype_sample_records, sample_records_from_noodles_header,
    variant_record_has_phased_genotype,
};
use super::decode::{decode_ds_record, decode_gt_record, BcfDenseDecodeBuffers, BcfStatsMode};
use super::haplotype::{decode_phased_haplotype_record, BcfHaplotypeDecodeBuffers};
use super::record::{prepare_bcf_candidate, push_bcf_variant_row, BcfCandidateAction};

const BCF_METADATA_INITIAL_VARIANT_CAPACITY: usize = 4096;

pub(in crate::vcf) fn read_metadata(path: &Path) -> Result<MetadataOutput> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    let samples = sample_records_from_noodles_header(&header);

    let mut variants = VariantMetadataBuffers::with_capacity(BCF_METADATA_INITIAL_VARIANT_CAPACITY);
    let mut has_phased_genotype_evidence = false;
    // Reuse noodles' lazy BCF record buffer and append public columns directly.
    let mut record = bcf::Record::default();
    loop {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        if !has_phased_genotype_evidence
            && variant_record_has_phased_genotype(path, &header, &record)?
        {
            has_phased_genotype_evidence = true;
        }
        append_public_variant_metadata_from_noodles_variant_record(
            path,
            &header,
            &record,
            &mut variants,
        )?;
    }

    let capabilities = if has_phased_genotype_evidence {
        SourceCapabilities::phased_genotypes()
    } else {
        SourceCapabilities::genotype_only()
    };

    Ok(MetadataOutput {
        samples: SampleMetadataBuffers::from_records(&samples, false)?,
        variants,
        capabilities,
    })
}

pub(super) struct BcfInput {
    pub(super) reader: bcf::io::Reader<noodles_bgzf::io::Reader<File>>,
    pub(super) header: noodles::Header,
    pub(super) selection: DenseSampleSelection,
}

pub(super) fn open_bcf_input(
    path: &Path,
    requested_samples: Option<&[String]>,
) -> Result<BcfInput> {
    open_bcf_input_with_hooks(path, requested_samples, || {}, || {})
}

pub(super) fn open_bcf_input_with_hooks(
    path: &Path,
    requested_samples: Option<&[String]>,
    on_source_open: impl FnOnce(),
    on_header_parse: impl FnOnce(),
) -> Result<BcfInput> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    on_source_open();
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    on_header_parse();
    let source_samples = sample_records_from_noodles_header(&header);
    let selection = select_samples_source_order(&source_samples, requested_samples, path)?;
    Ok(BcfInput {
        reader,
        header,
        selection,
    })
}

pub(in crate::vcf) fn read_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    let BcfInput {
        mut reader,
        header,
        selection,
    } = open_bcf_input(path, requested_samples)?;
    let DenseSampleSelection {
        source_indices,
        samples,
        mut diagnostics,
    } = selection;

    let n_rows = samples.len();
    let mut indptr = vec![0_i32];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = return_variants.then(|| VariantMetadataBuffers::with_capacity(0));
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = bcf::Record::default();
    let mut decoded = BcfDenseDecodeBuffers::with_capacity(source_indices.len());

    while !retention.window_is_satisfied() {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        let prepared = match prepare_bcf_candidate(
            path,
            &header,
            &record,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            BcfCandidateAction::Skip => continue,
            BcfCandidateAction::Stop => break,
            BcfCandidateAction::Decode(prepared) => prepared,
        };
        let variant = prepared.variant;
        let needs_genotype_decision = prepared.needs_genotype_decision;
        decode_gt_record(
            path,
            &header,
            &record,
            &source_indices,
            BcfStatsMode::from_needed(needs_genotype_decision),
            &mut decoded,
        )?;
        if needs_genotype_decision {
            match retention.genotype_decision(
                variant_filter
                    .is_none_or(|filter| filter.evaluate_view(&variant, decoded.stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }

        reject_sparse_missing(!decoded.missing_indices.is_empty())?;
        let flipped = flip_values_to_minor_allele(decoded.values.as_mut_slice());
        append_sparse_column(&mut indptr, &mut indices, &mut data, &decoded.values)?;
        push_bcf_variant_row(&mut variants, &variant, decoded.stats, flipped)?;
    }

    let n_cols = indptr.len() - 1;
    diagnostics.retained_variants = n_cols;
    let samples = SampleMetadataBuffers::optional_from_records(&samples, return_samples, false)?;
    SparseGenotypeMatrix::new(
        n_rows,
        n_cols,
        indptr,
        indices,
        data,
        samples,
        variants,
        diagnostics,
    )
}

pub(in crate::vcf) fn read_haplotypes_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    let BcfInput {
        mut reader,
        header,
        selection,
    } = open_bcf_input(path, requested_samples)?;
    let DenseSampleSelection {
        source_indices,
        samples: selected_samples,
        mut diagnostics,
    } = selection;

    let mut variants = return_variants.then(|| VariantMetadataBuffers::with_capacity(0));
    let mut variant_major_values = Vec::new();
    let mut n_variants = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = bcf::Record::default();
    let mut decoded_gt = BcfDenseDecodeBuffers::with_capacity(source_indices.len());
    let mut decoded = BcfHaplotypeDecodeBuffers::with_capacity(source_indices.len());

    while !retention.window_is_satisfied() {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        let prepared = match prepare_bcf_candidate(
            path,
            &header,
            &record,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            BcfCandidateAction::Skip => continue,
            BcfCandidateAction::Stop => break,
            BcfCandidateAction::Decode(prepared) => prepared,
        };
        let variant = prepared.variant;
        let needs_genotype_decision = prepared.needs_genotype_decision;
        let filter_result = if needs_genotype_decision {
            let filter = variant_filter.ok_or_else(|| {
                GenoioError::internal_contract("genotype decision requires a variant filter")
            })?;
            decode_gt_record(
                path,
                &header,
                &record,
                &source_indices,
                if return_variants {
                    BcfStatsMode::Compute
                } else {
                    BcfStatsMode::Counts
                },
                &mut decoded_gt,
            )?;
            Some(evaluate_bcf_gt_filter(
                &decoded_gt,
                filter,
                &variant,
                return_variants,
                "haplotype",
            )?)
        } else {
            None
        };
        if let Some((retain_variant, _)) = filter_result {
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }

        let stats_to_attach = filter_result.and_then(|(_, stats)| stats);
        push_bcf_variant_row(&mut variants, &variant, stats_to_attach, false)?;

        decode_phased_haplotype_record(path, &header, &record, &source_indices, &mut decoded)?;
        apply_dense_missing_policy_to_variant(
            &mut decoded.values,
            &decoded.missing_indices,
            missing_policy,
        )?;
        n_variants += 1;
        variant_major_values.extend_from_slice(&decoded.values);
    }

    let haplotype_samples = haplotype_sample_records(&selected_samples, &source_indices);
    let samples =
        SampleMetadataBuffers::optional_from_records(&haplotype_samples, return_samples, true)?;
    let n_samples = selected_samples.len() * 2;
    diagnostics.retained_variants = n_variants;
    DenseGenotypeMatrix::new_with_layout(
        n_samples,
        n_variants,
        variant_major_values,
        DenseLayout::VariantMajor,
        samples,
        variants,
        diagnostics,
    )
}

pub(in crate::vcf) fn read_haplotypes_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    let BcfInput {
        mut reader,
        header,
        selection,
    } = open_bcf_input(path, requested_samples)?;
    let DenseSampleSelection {
        source_indices,
        samples: selected_samples,
        mut diagnostics,
    } = selection;

    let haplotype_samples = haplotype_sample_records(&selected_samples, &source_indices);
    let output_samples =
        SampleMetadataBuffers::optional_from_records(&haplotype_samples, return_samples, true)?;
    let n_rows = selected_samples.len() * 2;
    let mut indptr = vec![0_i32];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = return_variants.then(|| VariantMetadataBuffers::with_capacity(0));
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = bcf::Record::default();
    let mut decoded = BcfHaplotypeDecodeBuffers::with_capacity(source_indices.len());
    let mut stats_decoded = BcfDenseDecodeBuffers::with_capacity(source_indices.len());

    while !retention.window_is_satisfied() {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        let prepared = match prepare_bcf_candidate(
            path,
            &header,
            &record,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            BcfCandidateAction::Skip => continue,
            BcfCandidateAction::Stop => break,
            BcfCandidateAction::Decode(prepared) => prepared,
        };
        let variant = prepared.variant;
        let needs_genotype_decision = prepared.needs_genotype_decision;
        let stats = if needs_genotype_decision {
            decode_gt_record(
                path,
                &header,
                &record,
                &source_indices,
                BcfStatsMode::Compute,
                &mut stats_decoded,
            )?;
            Some(
                stats_decoded
                    .stats
                    .ok_or_else(|| GenoioError::internal_contract("bcf GT stats missing"))?,
            )
        } else {
            None
        };
        if needs_genotype_decision {
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate_view(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }

        decode_phased_haplotype_record(path, &header, &record, &source_indices, &mut decoded)?;
        reject_sparse_missing(!decoded.missing_indices.is_empty())?;
        let flipped = flip_haplotype_values_to_minor_allele(decoded.values.as_mut_slice());
        append_sparse_column(&mut indptr, &mut indices, &mut data, &decoded.values)?;
        push_bcf_variant_row(&mut variants, &variant, stats, flipped)?;
    }

    let n_cols = indptr.len() - 1;
    diagnostics.retained_variants = n_cols;
    SparseGenotypeMatrix::new(
        n_rows,
        n_cols,
        indptr,
        indices,
        data,
        output_samples,
        variants,
        diagnostics,
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DenseField {
    Gt,
    Ds,
}

pub(in crate::vcf) fn read_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    read_dense_windowed_with_field(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        return_samples,
        return_variants,
        DenseField::Gt,
    )
}

pub(in crate::vcf) fn read_dosage_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    read_dense_windowed_with_field(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        return_samples,
        return_variants,
        DenseField::Ds,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "shared BCF dense loop mirrors read options plus field and metadata return choices"
)]
fn read_dense_windowed_with_field(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    field: DenseField,
) -> Result<DenseGenotypeMatrix> {
    let BcfInput {
        mut reader,
        header,
        selection,
    } = open_bcf_input(path, requested_samples)?;
    let DenseSampleSelection {
        source_indices,
        samples,
        mut diagnostics,
    } = selection;

    let mut variants = return_variants.then(|| VariantMetadataBuffers::with_capacity(0));
    let mut variant_major_values = Vec::new();
    let mut n_variants = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = bcf::Record::default();
    let mut decoded = BcfDenseDecodeBuffers::with_capacity(source_indices.len());

    while !retention.window_is_satisfied() {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        let prepared = match prepare_bcf_candidate(
            path,
            &header,
            &record,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            BcfCandidateAction::Skip => continue,
            BcfCandidateAction::Stop => break,
            BcfCandidateAction::Decode(prepared) => prepared,
        };
        let variant = prepared.variant;
        let needs_genotype_decision = prepared.needs_genotype_decision;
        match field {
            DenseField::Gt => decode_gt_record(
                path,
                &header,
                &record,
                &source_indices,
                match (needs_genotype_decision, return_variants) {
                    (true, false) => BcfStatsMode::Counts,
                    (true, true) => BcfStatsMode::Compute,
                    (false, _) => BcfStatsMode::Skip,
                },
                &mut decoded,
            )?,
            DenseField::Ds => {
                decode_ds_record(path, &header, &record, &source_indices, false, &mut decoded)?
            }
        };

        let mut stats_to_attach = None;
        if needs_genotype_decision {
            let (retain_variant, stats) = match field {
                DenseField::Gt => {
                    let filter = variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?;
                    evaluate_bcf_gt_filter(&decoded, filter, &variant, return_variants, "GT")?
                }
                DenseField::Ds => {
                    let filter = variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?;
                    evaluate_dosage_filter(
                        &decoded.values,
                        &decoded.missing_indices,
                        filter,
                        &variant,
                        return_variants,
                    )?
                }
            };
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            stats_to_attach = stats;
        }
        push_bcf_variant_row(&mut variants, &variant, stats_to_attach, false)?;

        n_variants += 1;
        apply_dense_missing_policy_to_variant(
            &mut decoded.values,
            &decoded.missing_indices,
            missing_policy,
        )?;
        variant_major_values.extend_from_slice(&decoded.values);
    }

    let n_samples = samples.len();
    diagnostics.retained_variants = n_variants;
    let samples = SampleMetadataBuffers::optional_from_records(&samples, return_samples, false)?;
    DenseGenotypeMatrix::new_with_layout(
        n_samples,
        n_variants,
        variant_major_values,
        DenseLayout::VariantMajor,
        samples,
        variants,
        diagnostics,
    )
}

pub(super) fn evaluate_bcf_gt_filter<V: VariantMetadataView + ?Sized>(
    decoded: &BcfDenseDecodeBuffers,
    filter: &VariantFilter,
    variant: &V,
    require_stats: bool,
    context: &str,
) -> Result<(bool, Option<VariantStats>)> {
    if !require_stats {
        let counts = decoded.counts.ok_or_else(|| {
            GenoioError::internal_contract(format!("bcf {context} filter fast path missing counts"))
        })?;
        return evaluate_hardcall_counts_filter(
            counts,
            filter,
            filter.genotype_filter_plan(),
            Some(variant),
            false,
        );
    }

    let stats = decoded.stats.ok_or_else(|| {
        GenoioError::internal_contract(format!("bcf {context} filter missing stats"))
    })?;
    Ok((filter.evaluate_view(variant, Some(&stats)), Some(stats)))
}

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

pub(super) fn flip_haplotype_values_to_minor_allele(values: &mut [f32]) -> bool {
    let a1_count = values.iter().sum::<f32>();
    let a0_count = values.len() as f32 - a1_count;
    if a1_count <= a0_count {
        return false;
    }
    for value in values {
        *value = 1.0 - *value;
    }
    true
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::super::super::read_vcf_dense_windowed_with_threads as read_vcf_dense_windowed_with_threads_facade;
    use super::*;
    use genoio_core::{
        DenseGenotypeMatrix, DenseLayout, SparseGenotypeMatrix, StringColumnBuffers,
        VariantMetadataBuffers,
    };
    use noodles_core::Position;
    use noodles_vcf::{
        self as noodles_vcf,
        header::record::value::{
            map::{
                format::{Number, Type},
                Contig, Format,
            },
            Map,
        },
        variant::{
            io::Write as _,
            record::samples::keys::key,
            record_buf::{samples::sample::Value, samples::Keys, AlternateBases, Ids, Samples},
        },
    };

    fn assert_values_with_nan(actual: &[f32], expected: &[f32]) {
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected) {
            if expected.is_nan() {
                assert!(actual.is_nan());
            } else {
                assert_eq!(actual, expected);
            }
        }
    }

    fn variants(output: &Option<VariantMetadataBuffers>) -> &VariantMetadataBuffers {
        output
            .as_ref()
            .expect("variant metadata buffers should be returned")
    }

    fn string_at(column: &StringColumnBuffers, index: usize) -> &str {
        let start = column.offsets[index] as usize;
        let end = column.offsets[index + 1] as usize;
        std::str::from_utf8(&column.values[start..end]).expect("string column should be UTF-8")
    }

    fn variant_ids(variants: &VariantMetadataBuffers) -> Vec<&str> {
        (0..variants.len())
            .map(|index| string_at(&variants.ids, index))
            .collect()
    }

    fn variant_id(variants: &VariantMetadataBuffers, index: usize) -> &str {
        string_at(&variants.ids, index)
    }

    fn variant_a0(variants: &VariantMetadataBuffers, index: usize) -> &str {
        string_at(&variants.a0s, index)
    }

    fn variant_a1(variants: &VariantMetadataBuffers, index: usize) -> &str {
        string_at(&variants.a1s, index)
    }

    fn read_vcf_dense_windowed_with_threads_for_test(
        path: &Path,
        requested_samples: Option<&[String]>,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        matrix_only: bool,
        threads: Option<usize>,
    ) -> Result<DenseGenotypeMatrix> {
        read_vcf_dense_windowed_with_threads_facade(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            DenseMissingPolicy::Nan,
            !matrix_only,
            !matrix_only,
            threads,
        )
    }

    fn read_dense_windowed_for_test(
        path: &Path,
        requested_samples: Option<&[String]>,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        read_dense_windowed_with_field(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            missing_policy,
            !matrix_only,
            !matrix_only,
            DenseField::Gt,
        )
    }

    fn read_dosage_dense_windowed_for_test(
        path: &Path,
        requested_samples: Option<&[String]>,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        super::read_dosage_dense_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            missing_policy,
            !matrix_only,
            !matrix_only,
        )
    }

    fn read_sparse_windowed_for_test(
        path: &Path,
        requested_samples: Option<&[String]>,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
    ) -> Result<SparseGenotypeMatrix> {
        super::read_sparse_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            true,
            true,
        )
    }

    fn read_haplotypes_dense_windowed_for_test(
        path: &Path,
        requested_samples: Option<&[String]>,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        super::read_haplotypes_dense_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            missing_policy,
            !matrix_only,
            !matrix_only,
        )
    }

    fn read_haplotypes_sparse_windowed_for_test(
        path: &Path,
        requested_samples: Option<&[String]>,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
    ) -> Result<SparseGenotypeMatrix> {
        super::read_haplotypes_sparse_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            true,
            true,
        )
    }

    fn write_test_bcf(path: &Path) {
        let file = fs::File::create(path).expect("test BCF should be created");
        let mut writer = noodles_bcf::io::Writer::new(file);
        let header = noodles_vcf::Header::builder()
            .add_contig("1", Map::<Contig>::new())
            .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
            .add_sample_name("s1")
            .add_sample_name("s2")
            .build();

        writer
            .write_header(&header)
            .expect("test BCF header should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs1", 10, "A", &["G"], ["0/0", "0/1"]),
            )
            .expect("first BCF record should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs2", 20, "C", &["T"], ["1/1", "./."]),
            )
            .expect("second BCF record should be written");
    }

    fn write_test_bcf_with_ds(path: &Path) {
        let file = fs::File::create(path).expect("test BCF should be created");
        let mut writer = noodles_bcf::io::Writer::new(file);
        let ds_format = Map::<Format>::builder()
            .set_number(Number::Count(1))
            .set_type(Type::Float)
            .set_description("Expected alternate allele dosage")
            .build()
            .expect("DS format should build");
        let header = noodles_vcf::Header::builder()
            .add_contig("1", Map::<Contig>::new())
            .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
            .add_format("DS", ds_format)
            .add_sample_name("s1")
            .add_sample_name("s2")
            .build();

        writer
            .write_header(&header)
            .expect("test BCF header should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_ds_test_record(
                    "rs1",
                    10,
                    "A",
                    &["G"],
                    [("0/0", Some(0.2)), ("0/1", Some(1.4))],
                ),
            )
            .expect("first BCF record should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_ds_test_record("rs2", 20, "C", &["T"], [("1/1", Some(2.0)), ("./.", None)]),
            )
            .expect("second BCF record should be written");
    }

    fn write_test_bcf_phased(path: &Path) {
        let file = fs::File::create(path).expect("test BCF should be created");
        let mut writer = noodles_bcf::io::Writer::new(file);
        let header = noodles_vcf::Header::builder()
            .add_contig("1", Map::<Contig>::new())
            .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
            .add_sample_name("s1")
            .add_sample_name("s2")
            .build();

        writer
            .write_header(&header)
            .expect("test BCF header should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs1", 10, "A", &["G"], ["0|1", "1|1"]),
            )
            .expect("first phased BCF record should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs2", 20, "C", &["T"], ["1|0", "0|0"]),
            )
            .expect("second phased BCF record should be written");
    }

    fn write_test_bcf_mixed_phase_for_stat_filter(path: &Path) {
        let file = fs::File::create(path).expect("test BCF should be created");
        let mut writer = noodles_bcf::io::Writer::new(file);
        let header = noodles_vcf::Header::builder()
            .add_contig("1", Map::<Contig>::new())
            .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
            .add_sample_name("s1")
            .add_sample_name("s2")
            .build();

        writer
            .write_header(&header)
            .expect("test BCF header should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs_phased", 10, "A", &["G"], ["0|1", "1|0"]),
            )
            .expect("first mixed-phase BCF record should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs_unphased_monomorphic", 20, "C", &["T"], ["0/0", "0/0"]),
            )
            .expect("second mixed-phase BCF record should be written");
    }

    fn bcf_test_record(
        id: &str,
        pos: usize,
        reference_bases: &str,
        alternate_bases: &[&str],
        genotypes: [&str; 2],
    ) -> noodles_vcf::variant::RecordBuf {
        let ids: Ids = [id.to_string()].into_iter().collect();
        let keys: Keys = [String::from(key::GENOTYPE)].into_iter().collect();
        let samples = Samples::new(
            keys,
            genotypes
                .into_iter()
                .map(|gt| vec![Some(Value::from(gt))])
                .collect(),
        );

        noodles_vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("1")
            .set_variant_start(Position::try_from(pos).expect("position should be valid"))
            .set_ids(ids)
            .set_reference_bases(reference_bases)
            .set_alternate_bases(AlternateBases::from(
                alternate_bases
                    .iter()
                    .copied()
                    .map(str::to_string)
                    .collect::<Vec<_>>(),
            ))
            .set_samples(samples)
            .build()
    }

    fn bcf_ds_test_record(
        id: &str,
        pos: usize,
        reference_bases: &str,
        alternate_bases: &[&str],
        calls: [(&str, Option<f32>); 2],
    ) -> noodles_vcf::variant::RecordBuf {
        let ids: Ids = [id.to_string()].into_iter().collect();
        let keys: Keys = [String::from(key::GENOTYPE), "DS".to_string()]
            .into_iter()
            .collect();
        let samples = Samples::new(
            keys,
            calls
                .into_iter()
                .map(|(gt, ds)| vec![Some(Value::from(gt)), ds.map(Value::from)])
                .collect(),
        );

        noodles_vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("1")
            .set_variant_start(Position::try_from(pos).expect("position should be valid"))
            .set_ids(ids)
            .set_reference_bases(reference_bases)
            .set_alternate_bases(AlternateBases::from(
                alternate_bases
                    .iter()
                    .copied()
                    .map(str::to_string)
                    .collect::<Vec<_>>(),
            ))
            .set_samples(samples)
            .build()
    }

    #[test]
    fn bcf_dense_gt_reads_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let dense = read_dense_windowed_for_test(
            file.path(),
            None,
            None,
            None,
            DenseMissingPolicy::Nan,
            false,
        )
        .expect("BCF dense GT should decode");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 2);
        assert_values_with_nan(&dense.values, &[0.0, 1.0, 2.0, f32::NAN]);
        assert_eq!(dense.layout, DenseLayout::VariantMajor);
        assert_eq!(variant_ids(variants(&dense.variants)), vec!["rs1", "rs2"]);
    }

    #[test]
    fn bcf_dense_gt_rejects_threaded_reads() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let error = read_vcf_dense_windowed_with_threads_for_test(
            file.path(),
            None,
            None,
            None,
            false,
            Some(2),
        )
        .expect_err("BCF should reject explicit thread count");

        assert!(error
            .to_string()
            .contains("threaded BCF reads are not supported"));
    }

    #[test]
    fn bcf_dense_gt_rejects_non_diploid_calls_without_collecting_alleles() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        {
            let out = fs::File::create(file.path()).expect("test BCF should be created");
            let mut writer = noodles_bcf::io::Writer::new(out);
            let header = noodles_vcf::Header::builder()
                .add_contig("1", Map::<Contig>::new())
                .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
                .add_sample_name("s1")
                .add_sample_name("s2")
                .build();
            writer
                .write_header(&header)
                .expect("test BCF header should be written");
            writer
                .write_variant_record(
                    &header,
                    &bcf_test_record("rs_triploid", 10, "A", &["G"], ["0/1/1", "0/0"]),
                )
                .expect("triploid BCF record should be written");
        }

        let error = read_dense_windowed_for_test(
            file.path(),
            None,
            None,
            None,
            DenseMissingPolicy::Nan,
            true,
        )
        .expect_err("non-diploid BCF GT should fail");

        assert!(error.to_string().contains("non-diploid GT"), "{error}");
    }

    #[test]
    fn bcf_dense_gt_rejects_multiallelic_allele_indexes_without_collecting_alleles() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        {
            let out = fs::File::create(file.path()).expect("test BCF should be created");
            let mut writer = noodles_bcf::io::Writer::new(out);
            let header = noodles_vcf::Header::builder()
                .add_contig("1", Map::<Contig>::new())
                .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
                .add_sample_name("s1")
                .add_sample_name("s2")
                .build();
            writer
                .write_header(&header)
                .expect("test BCF header should be written");
            writer
                .write_variant_record(
                    &header,
                    &bcf_test_record("rs_gt2", 10, "A", &["G"], ["0/2", "0/0"]),
                )
                .expect("multiallelic GT-index BCF record should be written");
        }

        let error = read_dense_windowed_for_test(
            file.path(),
            None,
            None,
            None,
            DenseMissingPolicy::Nan,
            true,
        )
        .expect_err("multiallelic BCF GT allele index should fail");

        assert!(
            error.to_string().contains("multiallelic GT allele index"),
            "{error}"
        );
    }

    #[test]
    fn bcf_sparse_gt_rejects_missing_calls_after_streaming_decode() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let error = read_sparse_windowed_for_test(
            file.path(),
            None,
            None,
            Some(VariantWindow { start: 1, len: 1 }),
        )
        .expect_err("sparse BCF GT should reject retained missing calls");

        assert!(error.to_string().contains("sparse missing values"));
    }

    #[test]
    fn bcf_dense_gt_applies_retained_windows() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let dense = read_dense_windowed_for_test(
            file.path(),
            None,
            None,
            Some(VariantWindow { start: 1, len: 1 }),
            DenseMissingPolicy::Nan,
            false,
        )
        .expect("BCF dense GT window should decode");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 1);
        assert_values_with_nan(&dense.values, &[2.0, f32::NAN]);
        assert_eq!(variant_id(variants(&dense.variants), 0), "rs2");
    }

    #[test]
    fn bcf_dense_gt_filters_stats_after_sample_selection() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());
        let samples = vec!["s2".to_string()];
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "mac",
            "params": {"min": 1}
        }))
        .expect("filter should parse");

        let dense = read_dense_windowed_for_test(
            file.path(),
            Some(&samples),
            Some(&filter),
            None,
            DenseMissingPolicy::Raise,
            false,
        )
        .expect("BCF dense GT filter should decode");

        assert_eq!(dense.n_samples, 1);
        assert_eq!(dense.n_variants, 1);
        assert_eq!(dense.values, vec![1.0]);
        let samples = dense.samples.as_ref().expect("sample metadata");
        assert_eq!(
            samples.iter().next().expect("first sample").iid.as_str(),
            "s2"
        );
        let variant_metadata = variants(&dense.variants);
        assert_eq!(variant_id(variant_metadata, 0), "rs1");
    }

    #[test]
    fn bcf_dense_gt_matrix_only_filters_from_hardcall_counts() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "mac",
            "params": {"min": 1}
        }))
        .expect("filter should parse");

        let dense = read_dense_windowed_for_test(
            file.path(),
            None,
            Some(&filter),
            None,
            DenseMissingPolicy::Raise,
            true,
        )
        .expect("matrix-only BCF dense GT should filter");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 1);
        assert!(dense.samples.is_none());
        assert!(dense.variants.is_none());
        assert_eq!(dense.values, vec![0.0, 1.0]);
    }

    #[test]
    fn bcf_dense_gt_can_return_samples_without_variant_buffers() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let dense = super::read_dense_windowed(
            file.path(),
            None,
            None,
            None,
            DenseMissingPolicy::Nan,
            true,
            false,
        )
        .expect("BCF dense GT output should decode without variant buffers");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 2);
        assert_eq!(dense.layout, DenseLayout::VariantMajor);
        let samples = dense.samples.as_ref().expect("sample metadata");
        assert_eq!(
            samples.iter().next().expect("first sample").iid.as_str(),
            "s1"
        );
        assert!(dense.variants.is_none());
        assert_values_with_nan(&dense.values, &[0.0, 1.0, 2.0, f32::NAN]);
    }

    #[test]
    fn bcf_dense_ds_reads_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_with_ds(file.path());

        let dense = read_dosage_dense_windowed_for_test(
            file.path(),
            None,
            None,
            None,
            DenseMissingPolicy::Nan,
            false,
        )
        .expect("BCF dense DS should decode");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 2);
        assert_values_with_nan(&dense.values, &[0.2, 1.4, 2.0, f32::NAN]);
        assert_eq!(dense.layout, DenseLayout::VariantMajor);
        assert_eq!(variant_ids(variants(&dense.variants)), vec!["rs1", "rs2"]);
    }

    #[test]
    fn bcf_dense_ds_matrix_only_filters_from_dosage_stats() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_with_ds(file.path());
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "mac",
            "params": {"max": 1}
        }))
        .expect("filter should parse");

        let dense = read_dosage_dense_windowed_for_test(
            file.path(),
            None,
            Some(&filter),
            None,
            DenseMissingPolicy::Nan,
            true,
        )
        .expect("matrix-only BCF dense DS should filter");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 1);
        assert!(dense.samples.is_none());
        assert!(dense.variants.is_none());
        assert_values_with_nan(&dense.values, &[2.0, f32::NAN]);
    }

    #[test]
    fn bcf_sparse_gt_reads_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "mac",
            "params": {"min": 1}
        }))
        .expect("filter should parse");

        let sparse = read_sparse_windowed_for_test(file.path(), None, Some(&filter), None)
            .expect("BCF sparse GT should decode");

        assert_eq!(sparse.n_rows, 2);
        assert_eq!(sparse.n_cols, 1);
        assert_eq!(sparse.indptr, vec![0, 1]);
        assert_eq!(sparse.indices, vec![1]);
        assert_eq!(sparse.data, vec![1.0]);
        assert_eq!(variant_id(variants(&sparse.variants), 0), "rs1");
    }

    #[test]
    fn bcf_sparse_gt_can_omit_samples_and_variant_buffers() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "mac",
            "params": {"min": 1}
        }))
        .expect("filter should parse");

        let sparse =
            super::read_sparse_windowed(file.path(), None, Some(&filter), None, false, false)
                .expect("BCF sparse GT output should decode without metadata buffers");

        assert_eq!(sparse.n_rows, 2);
        assert_eq!(sparse.n_cols, 1);
        assert!(sparse.samples.is_none());
        assert!(sparse.variants.is_none());
        assert_eq!(sparse.indptr, vec![0, 1]);
        assert_eq!(sparse.indices, vec![1]);
        assert_eq!(sparse.data, vec![1.0]);
    }

    #[test]
    fn bcf_dense_haplotypes_read_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_phased(file.path());

        let dense = read_haplotypes_dense_windowed_for_test(
            file.path(),
            None,
            None,
            None,
            DenseMissingPolicy::Raise,
            false,
        )
        .expect("BCF dense haplotypes should decode");

        assert_eq!(dense.n_samples, 4);
        assert_eq!(dense.n_variants, 2);
        assert_eq!(dense.values, vec![0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]);
        assert_eq!(dense.layout, DenseLayout::VariantMajor);
        let samples = dense.samples.as_ref().expect("sample metadata");
        let mut samples = samples.iter();
        let first = samples.next().expect("first haplotype sample");
        let second = samples.next().expect("second haplotype sample");
        assert_eq!(first.iid.as_str(), "s1");
        assert_eq!(first.haplotype_index, Some(0));
        assert_eq!(second.haplotype_index, Some(1));
    }

    #[test]
    fn bcf_dense_haplotypes_matrix_only_filter_drops_unphased_before_decode() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_mixed_phase_for_stat_filter(file.path());
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "maf",
            "params": {"min": 0.1}
        }))
        .expect("filter should parse");

        let dense = read_haplotypes_dense_windowed_for_test(
            file.path(),
            None,
            Some(&filter),
            None,
            DenseMissingPolicy::Raise,
            true,
        )
        .expect("matrix-only BCF haplotypes should prefilter before phased decode");

        assert_eq!(dense.n_samples, 4);
        assert_eq!(dense.n_variants, 1);
        assert!(dense.samples.is_none());
        assert!(dense.variants.is_none());
        assert_eq!(dense.values, vec![0.0, 1.0, 1.0, 0.0]);
    }

    #[test]
    fn bcf_sparse_haplotypes_read_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_phased(file.path());

        let sparse = read_haplotypes_sparse_windowed_for_test(file.path(), None, None, None)
            .expect("BCF sparse haplotypes should decode");

        assert_eq!(sparse.n_rows, 4);
        assert_eq!(sparse.n_cols, 2);
        assert_eq!(sparse.indptr, vec![0, 1, 2]);
        assert_eq!(sparse.indices, vec![0, 0]);
        assert_eq!(sparse.data, vec![1.0, 1.0]);
        let variant_metadata = variants(&sparse.variants);
        assert_eq!(variant_a0(variant_metadata, 0), "G");
        assert_eq!(variant_a1(variant_metadata, 0), "A");
        assert_eq!(variant_a0(variant_metadata, 1), "C");
        assert_eq!(variant_a1(variant_metadata, 1), "T");
    }
}
