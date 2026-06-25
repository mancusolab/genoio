// pattern: Imperative Shell
//! VCF and BCF reader facade.
//!
//! The facade routes `.bcf` paths to the typed BCF backend and all other paths
//! to the text VCF backend. Public functions preserve one contract for dense,
//! sparse, dosage, haplotype, metadata, threaded, and windowed reads.

use std::path::Path;

use genoio_core::{
    DenseGenotypeMatrix, DenseMissingPolicy, GenoioError, MetadataOutput, SampleRecord,
    SparseGenotypeMatrix, VariantFilter, VariantWindow,
};
use noodles_vcf as noodles;
use noodles_vcf::variant::record::samples::{
    keys::key,
    series::{value::genotype::Phasing as NoodlesGenotypePhasing, Value as NoodlesSampleValue},
};

use self::policy::read_text_vcf_with_optional_index;
use crate::error::Result;

mod bcf;
mod policy;
mod text;

/// Read VCF/BCF sample metadata and public variant metadata as columnar buffers.
pub fn read_vcf_public_metadata(path: &Path) -> Result<MetadataOutput> {
    if is_bcf_path(path) {
        return bcf::read_metadata(path);
    }
    text::read_vcf_public_metadata(path)
}

pub fn read_vcf_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_dense_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        return_samples,
        return_variants,
        None,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "Python metadata bridge passes sample and variant return choices explicitly"
)]
pub fn read_vcf_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return bcf::read_dense_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            missing_policy,
            return_samples,
            return_variants,
        );
    }

    read_text_vcf_with_optional_index(
        path,
        variant_filter,
        || {
            text::empty_vcf_dense(
                path,
                requested_samples,
                return_samples,
                return_variants,
                threads,
            )
        },
        |region| {
            text::read_vcf_dense_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                region,
                missing_policy,
                return_samples,
                return_variants,
                threads,
            )
        },
        || {
            text::read_vcf_dense(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                missing_policy,
                return_samples,
                return_variants,
                threads,
            )
        },
    )
}

pub fn read_vcf_dosage_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_dosage_dense_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        return_samples,
        return_variants,
        None,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "Python metadata bridge passes sample and variant return choices explicitly"
)]
pub fn read_vcf_dosage_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return bcf::read_dosage_dense_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            missing_policy,
            return_samples,
            return_variants,
        );
    }

    read_text_vcf_with_optional_index(
        path,
        variant_filter,
        || {
            text::empty_vcf_dense(
                path,
                requested_samples,
                return_samples,
                return_variants,
                threads,
            )
        },
        |region| {
            text::read_vcf_dosage_dense_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                region,
                missing_policy,
                return_samples,
                return_variants,
                threads,
            )
        },
        || {
            text::read_vcf_dosage_dense(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                missing_policy,
                return_samples,
                return_variants,
                threads,
            )
        },
    )
}

pub fn read_vcf_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    read_vcf_sparse_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        return_samples,
        return_variants,
        None,
    )
}

pub fn read_vcf_sparse_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return bcf::read_sparse_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            return_samples,
            return_variants,
        );
    }

    read_text_vcf_with_optional_index(
        path,
        variant_filter,
        || {
            text::empty_vcf_sparse(
                path,
                requested_samples,
                return_samples,
                return_variants,
                threads,
            )
        },
        |region| {
            text::read_vcf_sparse_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                region,
                return_samples,
                return_variants,
                threads,
            )
        },
        || {
            text::read_vcf_sparse(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                return_samples,
                return_variants,
                threads,
            )
        },
    )
}

pub fn read_vcf_haplotypes_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_haplotypes_dense_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        return_samples,
        return_variants,
        None,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "Python metadata bridge passes sample and variant return choices explicitly"
)]
pub fn read_vcf_haplotypes_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return bcf::read_haplotypes_dense_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            missing_policy,
            return_samples,
            return_variants,
        );
    }

    read_text_vcf_with_optional_index(
        path,
        variant_filter,
        || {
            text::empty_vcf_haplotypes_dense(
                path,
                requested_samples,
                return_samples,
                return_variants,
                threads,
            )
        },
        |region| {
            text::read_vcf_haplotypes_dense_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                region,
                missing_policy,
                return_samples,
                return_variants,
                threads,
            )
        },
        || {
            text::read_vcf_haplotypes_dense(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                missing_policy,
                return_samples,
                return_variants,
                threads,
            )
        },
    )
}

pub fn read_vcf_haplotypes_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    read_vcf_haplotypes_sparse_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        return_samples,
        return_variants,
        None,
    )
}

pub fn read_vcf_haplotypes_sparse_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return bcf::read_haplotypes_sparse_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            return_samples,
            return_variants,
        );
    }

    read_text_vcf_with_optional_index(
        path,
        variant_filter,
        || {
            text::empty_vcf_haplotypes_sparse(
                path,
                requested_samples,
                return_samples,
                return_variants,
                threads,
            )
        },
        |region| {
            text::read_vcf_haplotypes_sparse_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                region,
                return_samples,
                return_variants,
                threads,
            )
        },
        || {
            text::read_vcf_haplotypes_sparse(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                return_samples,
                return_variants,
                threads,
            )
        },
    )
}

fn validate_threaded_read_support(path: &Path, threads: Option<usize>) -> Result<()> {
    let Some(threads) = threads else {
        return Ok(());
    };
    if threads == 0 {
        return Err(GenoioError::invalid_source(
            path,
            "vcf thread count must be greater than zero",
        ));
    }
    if is_bcf_path(path) {
        return Err(GenoioError::invalid_source(
            path,
            "threaded BCF reads are not supported by the noodles backend",
        ));
    }
    if !is_compressed_vcf(path) {
        return Err(GenoioError::invalid_source(
            path,
            "threaded reads are only supported for compressed VCF sources",
        ));
    }
    Ok(())
}

fn sample_records_from_noodles_header(header: &noodles::Header) -> Vec<SampleRecord> {
    header
        .sample_names()
        .iter()
        .map(|sample| plain_sample_record(sample.to_string()))
        .collect()
}

fn plain_sample_record(iid: String) -> SampleRecord {
    SampleRecord {
        fid: None,
        iid,
        father: None,
        mother: None,
        sex: None,
        phenotype: None,
        source_sample_index: None,
        haplotype_index: None,
    }
}

fn variant_record_has_phased_genotype<R>(
    path: &Path,
    header: &noodles::Header,
    record: &R,
) -> Result<bool>
where
    R: noodles::variant::Record + ?Sized,
{
    let samples = record.samples().map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf samples error: {error}"))
    })?;
    let Some(gt_series_result) = samples.select(header, key::GENOTYPE) else {
        return Ok(false);
    };
    let gt_series = gt_series_result.map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf genotype series error: {error}"))
    })?;

    for value_result in gt_series.iter(header) {
        let Some(NoodlesSampleValue::Genotype(genotype)) = value_result.map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf genotype value error: {error}"))
        })?
        else {
            continue;
        };

        for allele_result in genotype.iter() {
            let (_, phasing) = allele_result.map_err(|error| {
                GenoioError::invalid_source(path, format!("vcf genotype allele error: {error}"))
            })?;
            if phasing == NoodlesGenotypePhasing::Phased {
                return Ok(true);
            }
        }
    }

    Ok(false)
}

fn finite_qual(qual: f32) -> Option<f32> {
    qual.is_finite().then_some(qual)
}

pub(super) fn is_compressed_vcf(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| matches!(extension, "gz" | "bgz"))
}

pub(super) fn is_bcf_path(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bcf"))
}

fn haplotype_sample_records(
    samples: &[SampleRecord],
    source_indices: &[usize],
) -> Vec<SampleRecord> {
    let mut haplotype_samples = Vec::with_capacity(samples.len() * 2);
    for (sample, source_index) in samples.iter().zip(source_indices) {
        let original_source_index = sample.source_sample_index.unwrap_or(*source_index);
        for haplotype_index in 0..2 {
            let mut haplotype_sample = sample.clone();
            haplotype_sample.source_sample_index = Some(original_source_index);
            haplotype_sample.haplotype_index = Some(haplotype_index);
            haplotype_samples.push(haplotype_sample);
        }
    }
    haplotype_samples
}
