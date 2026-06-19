//! Lazy BCF readers built on noodles-bcf records.
//!
//! BCF exposes typed sample values, so this module deliberately does not reuse
//! the text VCF byte scanners. The record loop mirrors the htslib path while
//! keeping one lazy `bcf::Record` buffer alive across variants.

// pattern: Mixed
// Reason: BCF setup, lazy record iteration, and decode routing share ownership
// of the same reusable record buffer.

use std::fs::File;
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, compute_dosage_variant_stats,
    flip_haplotype_values_to_minor_allele, flip_values_to_minor_allele,
    reject_sparse_missing_values, select_samples_source_order, variant_stats_from_counts,
    DenseGenotypeMatrix, GenoioError, PartialFilterDecision, SparseGenotypeMatrix, VariantFilter,
    VariantRecord, VariantStats, VariantWindow,
};
use noodles_bcf as bcf;
use noodles_vcf as noodles;
use noodles_vcf::variant::record::{
    samples::series::value::genotype::Phasing as NoodlesGenotypePhasing,
    samples::{keys::key, series::Value as NoodlesSampleValue},
    AlternateBases as _, Ids as _,
};

use crate::error::Result;
use crate::matrix::{finish_variant_major_dense_matrix, VariantMajorDenseParts};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::fast::metadata_variant_record_from_variant_record;
use super::{haplotype_sample_records, sample_records_from_noodles_header};

pub(super) fn read_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_dense_windowed_with_field(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
        DenseField::Gt,
    )
}

pub(super) fn read_dosage_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_dense_windowed_with_field(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
        DenseField::Ds,
    )
}

pub(super) fn read_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    let source_samples = sample_records_from_noodles_header(&header);
    let selection = select_samples_source_order(&source_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let n_rows = selection.samples.len();
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

        let mut variant = metadata_variant_record_from_variant_record(path, &header, &record)?;
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
            &selection.source_indices,
            needs_genotype_decision,
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
        selection.samples,
        variants,
        diagnostics,
    )
}

pub(super) fn read_haplotypes_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    let source_samples = sample_records_from_noodles_header(&header);
    let selection = select_samples_source_order(&source_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
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
            Some(metadata_variant_record_from_variant_record(
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
        let stats = if needs_genotype_decision {
            Some(
                decode_gt_record(path, &header, &record, &selection.source_indices, true)?
                    .stats
                    .ok_or_else(|| GenoioError::internal_contract("bcf GT stats missing"))?,
            )
        } else {
            None
        };
        if needs_genotype_decision {
            let variant = variant.as_ref().ok_or_else(|| {
                GenoioError::internal_contract("bcf filter requires variant metadata")
            })?;
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
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

        let decoded =
            decode_phased_haplotype_record(path, &header, &record, &selection.source_indices)?;
        n_variants += 1;
        variant_major_values.extend(decoded.values);
        variant_major_missing.extend(decoded.missing);
    }

    let samples = if matrix_only {
        Vec::new()
    } else {
        haplotype_sample_records(&selection.samples, &selection.source_indices)
    };
    let n_samples = selection.samples.len() * 2;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
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
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    let source_samples = sample_records_from_noodles_header(&header);
    let selection = select_samples_source_order(&source_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let samples = haplotype_sample_records(&selection.samples, &selection.source_indices);
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

        let mut variant = metadata_variant_record_from_variant_record(path, &header, &record)?;
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
                decode_gt_record(path, &header, &record, &selection.source_indices, true)?
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
        let decoded =
            decode_phased_haplotype_record(path, &header, &record, &selection.source_indices)?;
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
    matrix_only: bool,
    field: DenseField,
) -> Result<DenseGenotypeMatrix> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    let source_samples = sample_records_from_noodles_header(&header);
    let selection = select_samples_source_order(&source_samples, requested_samples, path)?;
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
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
            Some(metadata_variant_record_from_variant_record(
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
        let decoded = match field {
            DenseField::Gt => decode_gt_record(
                path,
                &header,
                &record,
                &selection.source_indices,
                needs_genotype_decision,
            )?,
            DenseField::Ds => decode_ds_record(
                path,
                &header,
                &record,
                &selection.source_indices,
                needs_genotype_decision,
            )?,
        };

        if needs_genotype_decision {
            let variant = variant.as_ref().ok_or_else(|| {
                GenoioError::internal_contract("bcf filter requires variant metadata")
            })?;
            match retention.genotype_decision(
                variant_filter
                    .is_none_or(|filter| filter.evaluate(variant, decoded.stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }

        if !matrix_only {
            let mut variant = variant.ok_or_else(|| {
                GenoioError::internal_contract("bcf metadata output requires variant metadata")
            })?;
            if let Some(stats) = decoded.stats {
                attach_variant_stats(&mut variant, stats);
            }
            variants.push(variant);
        }

        n_variants += 1;
        variant_major_values.extend(decoded.values);
        variant_major_missing.extend(decoded.missing);
    }

    let n_samples = selection.samples.len();
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

fn validate_biallelic_lazy_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
) -> Result<()> {
    if record.alternate_bases().len() == 1 {
        return Ok(());
    }

    let variant = metadata_variant_record_from_variant_record(path, header, record)?;
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
}

fn decode_gt_record(
    path: &Path,
    header: &noodles::Header,
    record: &bcf::Record,
    source_indices: &[usize],
    collect_stats: bool,
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
    let mut hom_ref_count = 0_u64;
    let mut het_count = 0_u64;
    let mut hom_alt_count = 0_u64;
    let mut missing_count = 0_u64;

    for source_index in source_indices {
        let call = decode_gt_call(path, header, record, &gt_series, *source_index)?;
        if collect_stats {
            match call.class {
                BcfGtClass::HomRef => hom_ref_count += 1,
                BcfGtClass::Het => het_count += 1,
                BcfGtClass::HomAlt => hom_alt_count += 1,
                BcfGtClass::Missing => missing_count += 1,
            }
        }
        values.push(call.value);
        missing.push(call.is_missing());
    }

    let stats = if collect_stats {
        Some(variant_stats_from_counts(
            hom_ref_count,
            het_count,
            hom_alt_count,
            missing_count,
        )?)
    } else {
        None
    };
    Ok(DecodedBcfDenseValues {
        values,
        missing,
        stats,
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
