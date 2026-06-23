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
    append_sparse_column, attach_variant_stats, compute_dosage_variant_stats,
    flip_haplotype_values_to_minor_allele, flip_values_to_minor_allele,
    reject_sparse_missing_values, select_samples_source_order, DenseGenotypeMatrix,
    DenseMissingPolicy, DenseSampleSelection, GenoioError, MetadataOutput, PartialFilterDecision,
    SourceCapabilities, SparseGenotypeMatrix, VariantFilter, VariantRecord, VariantStats,
    VariantWindow,
};
use noodles_bcf as bcf;
use noodles_vcf as noodles;
use noodles_vcf::variant::record::{
    samples::series::value::genotype::Phasing as NoodlesGenotypePhasing,
    samples::{keys::key, series::Value as NoodlesSampleValue},
    AlternateBases as _, Ids as _,
};

use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::hardcall::{evaluate_hardcall_counts_filter, HardcallCounts};
use crate::matrix::{
    apply_dense_missing_policy_to_variant, finish_variant_major_dense_matrix,
    missing_indices_from_mask, VariantMajorDenseParts,
};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::text::variant_record_from_noodles_variant_record;
use super::{
    haplotype_sample_records, sample_records_from_noodles_header,
    variant_record_has_phased_genotype,
};

pub(super) fn read_metadata(path: &Path) -> Result<MetadataOutput> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    let samples = sample_records_from_noodles_header(&header);

    let mut variants = Vec::new();
    let mut has_phased_genotype_evidence = false;
    // Reuse noodles' lazy BCF record buffer so metadata scans do not allocate a
    // full RecordBuf for each variant before genotype decoding exists.
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
        variants.push(variant_record_from_noodles_variant_record(
            path, &header, &record,
        )?);
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

pub(super) fn read_dense_windowed(
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
        matrix_only,
        DenseField::Gt,
    )
}

pub(super) fn read_dosage_dense_windowed(
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
        matrix_only,
        DenseField::Ds,
    )
}

struct BcfInput {
    reader: bcf::io::Reader<noodles_bgzf::io::Reader<File>>,
    header: noodles::Header,
    selection: DenseSampleSelection,
}

fn open_bcf_input(path: &Path, requested_samples: Option<&[String]>) -> Result<BcfInput> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    let source_samples = sample_records_from_noodles_header(&header);
    let selection = select_samples_source_order(&source_samples, requested_samples, path)?;
    Ok(BcfInput {
        reader,
        header,
        selection,
    })
}

pub(super) fn read_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
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
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = bcf::Record::default();

    while !retention.window_is_satisfied() {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        let mut variant = variant_record_from_noodles_variant_record(path, &header, &record)?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        validate_biallelic_variant(path, &variant)?;
        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        let decoded = decode_gt_record(
            path,
            &header,
            &record,
            &source_indices,
            BcfStatsMode::from_needed(needs_genotype_decision),
        )?;
        if needs_genotype_decision {
            match retention.genotype_decision(
                variant_filter
                    .is_none_or(|filter| filter.evaluate(&variant, decoded.stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }

        if let Some(stats) = decoded.stats {
            attach_variant_stats(&mut variant, stats);
        }
        reject_sparse_missing_values(&decoded.missing)?;
        let mut values = decoded.values;
        flip_values_to_minor_allele(&mut values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &values);
        variants.push(variant);
    }

    let n_cols = variants.len();
    diagnostics.retained_variants = n_cols;
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

pub(super) fn read_haplotypes_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    matrix_only: bool,
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

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut missing_indices = Vec::new();
    let mut n_variants = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = bcf::Record::default();

    while !retention.window_is_satisfied() {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        let variant = if !matrix_only || variant_filter.is_some() {
            Some(variant_record_from_noodles_variant_record(
                path, &header, &record,
            )?)
        } else {
            None
        };
        let partial_decision = match (variant_filter, variant.as_ref()) {
            (Some(filter), Some(variant)) => filter.partial_decision(variant),
            (Some(_), None) => {
                return Err(GenoioError::internal_contract(
                    "bcf filter requires variant metadata",
                ));
            }
            (None, _) => PartialFilterDecision::Accept,
        };
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        if let Some(variant) = variant.as_ref() {
            validate_biallelic_variant(path, variant)?;
        } else {
            validate_biallelic_lazy_record(path, &header, &record)?;
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        let filter_result = if needs_genotype_decision {
            let variant = variant.as_ref().ok_or_else(|| {
                GenoioError::internal_contract("bcf filter requires variant metadata")
            })?;
            let filter = variant_filter.ok_or_else(|| {
                GenoioError::internal_contract("genotype decision requires a variant filter")
            })?;
            let decoded_gt = decode_gt_record(
                path,
                &header,
                &record,
                &source_indices,
                if matrix_only {
                    BcfStatsMode::Counts
                } else {
                    BcfStatsMode::Compute
                },
            )?;
            Some(evaluate_bcf_gt_filter(
                &decoded_gt,
                filter,
                variant,
                matrix_only,
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

        if !matrix_only {
            let mut variant = variant.ok_or_else(|| {
                GenoioError::internal_contract("bcf metadata output requires variant metadata")
            })?;
            if let Some((_, Some(stats))) = filter_result {
                attach_variant_stats(&mut variant, stats);
            }
            variants.push(variant);
        }

        let mut decoded = decode_phased_haplotype_record(path, &header, &record, &source_indices)?;
        missing_indices_from_mask(&decoded.missing, &mut missing_indices);
        apply_dense_missing_policy_to_variant(
            &mut decoded.values,
            &missing_indices,
            missing_policy,
        )?;
        n_variants += 1;
        variant_major_values.extend(decoded.values);
    }

    let samples = if matrix_only {
        Vec::new()
    } else {
        haplotype_sample_records(&selected_samples, &source_indices)
    };
    let n_samples = selected_samples.len() * 2;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

pub(super) fn read_haplotypes_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
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

    let samples = haplotype_sample_records(&selected_samples, &source_indices);
    let n_rows = samples.len();
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = bcf::Record::default();

    while !retention.window_is_satisfied() {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        let mut variant = variant_record_from_noodles_variant_record(path, &header, &record)?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        validate_biallelic_variant(path, &variant)?;
        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        let stats = if needs_genotype_decision {
            Some(
                decode_gt_record(
                    path,
                    &header,
                    &record,
                    &source_indices,
                    BcfStatsMode::Compute,
                )?
                .stats
                .ok_or_else(|| GenoioError::internal_contract("bcf GT stats missing"))?,
            )
        } else {
            None
        };
        if needs_genotype_decision {
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }

        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        let decoded = decode_phased_haplotype_record(path, &header, &record, &source_indices)?;
        reject_sparse_missing_values(&decoded.missing)?;
        let mut values = decoded.values;
        flip_haplotype_values_to_minor_allele(&mut values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &values);
        variants.push(variant);
    }

    let n_cols = variants.len();
    diagnostics.retained_variants = n_cols;
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DenseField {
    Gt,
    Ds,
}

fn read_dense_windowed_with_field(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    matrix_only: bool,
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

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut missing_indices = Vec::new();
    let mut n_variants = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = bcf::Record::default();

    while !retention.window_is_satisfied() {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        let variant = if !matrix_only || variant_filter.is_some() {
            Some(variant_record_from_noodles_variant_record(
                path, &header, &record,
            )?)
        } else {
            None
        };
        let partial_decision = match (variant_filter, variant.as_ref()) {
            (Some(filter), Some(variant)) => filter.partial_decision(variant),
            (Some(_), None) => {
                return Err(GenoioError::internal_contract(
                    "bcf filter requires variant metadata",
                ));
            }
            (None, _) => PartialFilterDecision::Accept,
        };
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        if let Some(variant) = variant.as_ref() {
            validate_biallelic_variant(path, variant)?;
        } else {
            validate_biallelic_lazy_record(path, &header, &record)?;
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        let mut decoded = match field {
            DenseField::Gt => decode_gt_record(
                path,
                &header,
                &record,
                &source_indices,
                match (needs_genotype_decision, matrix_only) {
                    (true, true) => BcfStatsMode::Counts,
                    (true, false) => BcfStatsMode::Compute,
                    (false, _) => BcfStatsMode::Skip,
                },
            )?,
            DenseField::Ds => decode_ds_record(path, &header, &record, &source_indices, false)?,
        };

        if needs_genotype_decision {
            let variant_ref = variant.as_ref().ok_or_else(|| {
                GenoioError::internal_contract("bcf filter requires variant metadata")
            })?;
            let (retain_variant, stats) = match field {
                DenseField::Gt => {
                    let filter = variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?;
                    evaluate_bcf_gt_filter(&decoded, filter, variant_ref, matrix_only, "GT")?
                }
                DenseField::Ds => {
                    let filter = variant_filter.ok_or_else(|| {
                        GenoioError::internal_contract(
                            "genotype decision requires a variant filter",
                        )
                    })?;
                    evaluate_dosage_filter(
                        &decoded.values,
                        &decoded.missing,
                        filter,
                        variant_ref,
                        !matrix_only,
                    )?
                }
            };
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if !matrix_only {
                let mut variant = variant.ok_or_else(|| {
                    GenoioError::internal_contract("bcf metadata output requires variant metadata")
                })?;
                if let Some(stats) = stats {
                    attach_variant_stats(&mut variant, stats);
                }
                variants.push(variant);
            }
        } else if !matrix_only {
            let variant = variant.ok_or_else(|| {
                GenoioError::internal_contract("bcf metadata output requires variant metadata")
            })?;
            variants.push(variant);
        }

        n_variants += 1;
        missing_indices_from_mask(&decoded.missing, &mut missing_indices);
        apply_dense_missing_policy_to_variant(
            &mut decoded.values,
            &missing_indices,
            missing_policy,
        )?;
        variant_major_values.extend(decoded.values);
    }

    let n_samples = samples.len();
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

fn evaluate_bcf_gt_filter(
    decoded: &DecodedBcfDenseValues,
    filter: &VariantFilter,
    variant: &VariantRecord,
    matrix_only: bool,
    context: &str,
) -> Result<(bool, Option<VariantStats>)> {
    if matrix_only {
        let counts = decoded.counts.ok_or_else(|| {
            GenoioError::internal_contract(format!(
                "matrix-only bcf {context} filter missing counts"
            ))
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
    Ok((filter.evaluate(variant, Some(&stats)), Some(stats)))
}

fn validate_biallelic_lazy_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
) -> Result<()> {
    if record.alternate_bases().len() == 1 {
        return Ok(());
    }

    let variant = variant_record_from_noodles_variant_record(path, header, record)?;
    validate_biallelic_variant(path, &variant)
}

fn validate_biallelic_variant(path: &Path, variant: &VariantRecord) -> Result<()> {
    if variant
        .alt_allele
        .as_deref()
        .is_some_and(|alt| !alt.is_empty() && !alt.contains(','))
    {
        return Ok(());
    }

    if variant
        .alt_allele
        .as_deref()
        .is_some_and(|alt| alt.contains(','))
    {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf dense reads require biallelic records; record {}:{} has multi-ALT alleles: multi-ALT records are not supported",
                variant.chrom, variant.pos
            ),
        ));
    }

    Err(GenoioError::invalid_source(
        path,
        format!(
            "vcf dense reads require biallelic records; record {}:{} is not biallelic",
            variant.chrom, variant.pos
        ),
    ))
}

struct DecodedBcfDenseValues {
    values: Vec<f32>,
    missing: Vec<bool>,
    stats: Option<VariantStats>,
    counts: Option<HardcallCounts>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BcfStatsMode {
    Skip,
    Counts,
    Compute,
}

impl BcfStatsMode {
    const fn from_needed(needed: bool) -> Self {
        if needed {
            Self::Compute
        } else {
            Self::Skip
        }
    }
}

fn decode_gt_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    source_indices: &[usize],
    stats_mode: BcfStatsMode,
) -> Result<DecodedBcfDenseValues> {
    let samples = record.samples().map_err(|error| {
        GenoioError::invalid_source(path, format!("bcf samples error: {error}"))
    })?;
    let gt_series = samples
        .select(header, key::GENOTYPE)
        .ok_or_else(|| GenoioError::invalid_source(path, "bcf record is missing FORMAT/GT"))?
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf genotype series error: {error}"))
        })?;

    let mut values = Vec::with_capacity(source_indices.len());
    let mut missing = Vec::with_capacity(source_indices.len());
    let mut counts = HardcallCounts::default();

    for source_index in source_indices {
        let call = decode_gt_call(path, header, record, &gt_series, *source_index)?;
        if !matches!(stats_mode, BcfStatsMode::Skip) {
            match call.class {
                BcfGtClass::HomRef => counts.record_hom_ref(),
                BcfGtClass::Het => counts.record_het(),
                BcfGtClass::HomAlt => counts.record_hom_alt(),
                BcfGtClass::Missing => counts.record_missing(),
            }
        }
        values.push(call.value);
        missing.push(call.is_missing());
    }

    let stats = if matches!(stats_mode, BcfStatsMode::Compute) {
        Some(counts.variant_stats()?)
    } else {
        None
    };
    let counts = if matches!(stats_mode, BcfStatsMode::Counts) {
        Some(counts)
    } else {
        None
    };
    Ok(DecodedBcfDenseValues {
        values,
        missing,
        stats,
        counts,
    })
}

fn decode_ds_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    source_indices: &[usize],
    collect_stats: bool,
) -> Result<DecodedBcfDenseValues> {
    let samples = record.samples().map_err(|error| {
        GenoioError::invalid_source(path, format!("bcf samples error: {error}"))
    })?;
    let ds_series = samples
        .select(header, "DS")
        .ok_or_else(|| {
            GenoioError::unsupported("vcf dosage reads require FORMAT/DS values: missing DS")
        })?
        .map_err(|error| {
            GenoioError::unsupported(format!(
                "vcf dosage reads require FORMAT/DS values: {error}"
            ))
        })?;

    let mut values = Vec::with_capacity(source_indices.len());
    let mut missing = Vec::with_capacity(source_indices.len());
    for source_index in source_indices {
        let value = ds_series
            .get(header, *source_index)
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    path,
                    format!(
                        "bcf record {} is missing a DS sample value",
                        record_id(record)
                    ),
                )
            })?
            .transpose()
            .map_err(|error| {
                GenoioError::unsupported(format!(
                    "vcf dosage reads require FORMAT/DS values: {error}"
                ))
            })?;

        let Some(value) = value else {
            values.push(0.0);
            missing.push(true);
            continue;
        };
        let NoodlesSampleValue::Float(value) = value else {
            return Err(GenoioError::unsupported(
                "vcf dosage reads require scalar FORMAT/DS float values",
            ));
        };
        if !value.is_finite() || !(0.0..=2.0).contains(&value) {
            return Err(GenoioError::invalid_source(
                path,
                format!(
                    "vcf record {} has invalid FORMAT/DS value {value}; expected finite value in [0, 2]",
                    record_id(record)
                ),
            ));
        }
        values.push(value);
        missing.push(false);
    }

    let stats = if collect_stats {
        Some(compute_dosage_variant_stats(&values, &missing)?)
    } else {
        None
    };
    Ok(DecodedBcfDenseValues {
        values,
        missing,
        stats,
        counts: None,
    })
}

struct DecodedBcfHaplotypes {
    values: Vec<f32>,
    missing: Vec<bool>,
}

fn decode_phased_haplotype_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    source_indices: &[usize],
) -> Result<DecodedBcfHaplotypes> {
    let samples = record.samples().map_err(|error| {
        GenoioError::invalid_source(path, format!("bcf samples error: {error}"))
    })?;
    let gt_series = samples
        .select(header, key::GENOTYPE)
        .ok_or_else(|| GenoioError::invalid_source(path, "bcf record is missing FORMAT/GT"))?
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf genotype series error: {error}"))
        })?;

    let mut values = Vec::with_capacity(source_indices.len() * 2);
    let mut missing = Vec::with_capacity(source_indices.len() * 2);
    for source_index in source_indices {
        let (sample_values, sample_missing) =
            decode_phased_haplotype_call(path, header, record, &gt_series, *source_index)?;
        values.extend(sample_values);
        missing.extend(sample_missing);
    }

    Ok(DecodedBcfHaplotypes { values, missing })
}

fn decode_phased_haplotype_call(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    gt_series: &bcf::record::samples::Series<'_>,
    source_index: usize,
) -> Result<([f32; 2], [bool; 2])> {
    let value = gt_series
        .get(header, source_index)
        .ok_or_else(|| {
            GenoioError::invalid_source(
                path,
                format!(
                    "bcf record {} is missing a GT sample value",
                    record_id(record)
                ),
            )
        })?
        .transpose()
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf genotype value error: {error}"))
        })?;

    let Some(NoodlesSampleValue::Genotype(genotype)) = value else {
        return Ok(([0.0, 0.0], [true, true]));
    };
    let alleles = genotype
        .iter()
        .map(|result| {
            result.map_err(|error| {
                GenoioError::invalid_source(path, format!("bcf genotype allele error: {error}"))
            })
        })
        .collect::<Result<Vec<_>>>()?;
    if alleles.len() != 2 {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf record {} has non-diploid GT with {} alleles",
                record_id(record),
                alleles.len()
            ),
        ));
    }

    let mut values = [0.0, 0.0];
    let mut missing = [false, false];
    for (allele_index, (allele, phasing)) in alleles.into_iter().enumerate() {
        if allele_index > 0 && phasing == NoodlesGenotypePhasing::Unphased {
            return Err(GenoioError::unsupported(format!(
                "vcf haplotype read record {} contains an unphased GT separator in a retained haplotype variant",
                record_id(record)
            )));
        }
        match allele {
            None => missing[allele_index] = true,
            Some(0) => {}
            Some(1) => values[allele_index] = 1.0,
            Some(other) => {
                return Err(GenoioError::invalid_source(
                    path,
                    format!(
                        "vcf record {} has multiallelic GT allele index {other}",
                        record_id(record)
                    ),
                ));
            }
        }
    }

    Ok((values, missing))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BcfGtClass {
    HomRef,
    Het,
    HomAlt,
    Missing,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct BcfGtCall {
    value: f32,
    class: BcfGtClass,
}

impl BcfGtCall {
    fn is_missing(self) -> bool {
        self.class == BcfGtClass::Missing
    }
}

fn decode_gt_call(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    gt_series: &bcf::record::samples::Series<'_>,
    source_index: usize,
) -> Result<BcfGtCall> {
    let value = gt_series
        .get(header, source_index)
        .ok_or_else(|| {
            GenoioError::invalid_source(
                path,
                format!(
                    "bcf record {} is missing a GT sample value",
                    record_id(record)
                ),
            )
        })?
        .transpose()
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf genotype value error: {error}"))
        })?;

    let Some(NoodlesSampleValue::Genotype(genotype)) = value else {
        return Ok(BcfGtCall {
            value: 0.0,
            class: BcfGtClass::Missing,
        });
    };

    let alleles = genotype
        .iter()
        .map(|result| {
            result.map(|(position, _)| position).map_err(|error| {
                GenoioError::invalid_source(path, format!("bcf genotype allele error: {error}"))
            })
        })
        .collect::<Result<Vec<_>>>()?;
    if alleles.len() != 2 {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf record {} has non-diploid GT with {} alleles",
                record_id(record),
                alleles.len()
            ),
        ));
    }

    let mut alt_count = 0_u8;
    for allele in alleles {
        let Some(allele) = allele else {
            return Ok(BcfGtCall {
                value: 0.0,
                class: BcfGtClass::Missing,
            });
        };
        match allele {
            0 => {}
            1 => alt_count += 1,
            other => {
                return Err(GenoioError::invalid_source(
                    path,
                    format!(
                        "vcf record {} has multiallelic GT allele index {other}",
                        record_id(record)
                    ),
                ));
            }
        }
    }

    let class = match alt_count {
        0 => BcfGtClass::HomRef,
        1 => BcfGtClass::Het,
        2 => BcfGtClass::HomAlt,
        _ => unreachable!("two diploid GT alleles can only produce dosage 0, 1, or 2"),
    };
    Ok(BcfGtCall {
        value: f32::from(alt_count),
        class,
    })
}

fn record_id(record: &bcf::Record) -> String {
    record.ids().iter().next().unwrap_or(".").to_string()
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::super::read_vcf_dense_windowed_with_threads;
    use super::*;
    use genoio_core::DenseLayout;
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

        let dense = read_dense_windowed(
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
        assert_eq!(
            dense
                .variants
                .iter()
                .map(|variant| variant.id.as_str())
                .collect::<Vec<_>>(),
            vec!["rs1", "rs2"]
        );
    }

    #[test]
    fn bcf_dense_gt_rejects_threaded_reads() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let error =
            read_vcf_dense_windowed_with_threads(file.path(), None, None, None, false, Some(2))
                .expect_err("BCF should reject explicit thread count");

        assert!(error
            .to_string()
            .contains("threaded BCF reads are not supported"));
    }

    #[test]
    fn bcf_dense_gt_applies_retained_windows() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let dense = read_dense_windowed(
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
        assert_eq!(dense.variants[0].id, "rs2");
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

        let dense = read_dense_windowed(
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
        assert_eq!(dense.samples[0].iid, "s2");
        assert_eq!(dense.variants[0].id, "rs1");
        assert_eq!(dense.variants[0].mac, Some(1));
        assert_eq!(dense.variants[0].n_called, Some(1));
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

        let dense = read_dense_windowed(
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
        assert!(dense.samples.is_empty());
        assert!(dense.variants.is_empty());
        assert_eq!(dense.values, vec![0.0, 1.0]);
    }

    #[test]
    fn bcf_dense_ds_reads_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_with_ds(file.path());

        let dense = read_dosage_dense_windowed(
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
        assert_eq!(
            dense
                .variants
                .iter()
                .map(|variant| variant.id.as_str())
                .collect::<Vec<_>>(),
            vec!["rs1", "rs2"]
        );
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

        let dense = read_dosage_dense_windowed(
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
        assert!(dense.samples.is_empty());
        assert!(dense.variants.is_empty());
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

        let sparse = read_sparse_windowed(file.path(), None, Some(&filter), None)
            .expect("BCF sparse GT should decode");

        assert_eq!(sparse.n_rows, 2);
        assert_eq!(sparse.n_cols, 1);
        assert_eq!(sparse.indptr, vec![0, 1]);
        assert_eq!(sparse.indices, vec![1]);
        assert_eq!(sparse.data, vec![1.0]);
        assert_eq!(sparse.variants[0].id, "rs1");
    }

    #[test]
    fn bcf_dense_haplotypes_read_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_phased(file.path());

        let dense = read_haplotypes_dense_windowed(
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
        assert_eq!(dense.samples[0].iid, "s1");
        assert_eq!(dense.samples[0].haplotype_index, Some(0));
        assert_eq!(dense.samples[1].haplotype_index, Some(1));
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

        let dense = read_haplotypes_dense_windowed(
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
        assert!(dense.samples.is_empty());
        assert!(dense.variants.is_empty());
        assert_eq!(dense.values, vec![0.0, 1.0, 1.0, 0.0]);
    }

    #[test]
    fn bcf_sparse_haplotypes_read_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_phased(file.path());

        let sparse = read_haplotypes_sparse_windowed(file.path(), None, None, None)
            .expect("BCF sparse haplotypes should decode");

        assert_eq!(sparse.n_rows, 4);
        assert_eq!(sparse.n_cols, 2);
        assert_eq!(sparse.indptr, vec![0, 1, 2]);
        assert_eq!(sparse.indices, vec![0, 0]);
        assert_eq!(sparse.data, vec![1.0, 1.0]);
        assert!(sparse.variants[0].flipped);
        assert!(!sparse.variants[1].flipped);
    }
}
