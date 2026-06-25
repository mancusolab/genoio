// pattern: Imperative Shell
#![allow(dead_code, unused_imports)]

use std::path::Path;

use genoio_core::{
    DenseGenotypeMatrixArrowVariants, DenseMissingPolicy, SparseGenotypeMatrixArrowVariants,
    VariantFilter, VariantWindow,
};

pub(crate) use super::arrow::{
    dense_missing_sample_major_arrow, dense_values_sample_major_arrow, sparse_values_dense_arrow,
    string_at, strings, variant_a0, variant_a1, variant_chrom, variant_id, variant_ids, variants,
};
pub(crate) use ::genoio_io::Result;

pub(crate) fn read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    read_vcf_dense_windowed(path, requested_samples, variant_filter, None, false)
}

pub(crate) fn read_vcf_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_dense_windowed_with_arrow_variants(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_vcf_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_dense_windowed_with_threads_and_arrow_variants(
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

pub(crate) fn read_vcf_dosage_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_dosage_dense_windowed_with_arrow_variants(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_vcf_dosage_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_dosage_dense_windowed_with_threads_and_arrow_variants(
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

pub(crate) fn read_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    read_vcf_sparse_windowed(path, requested_samples, variant_filter, None)
}

pub(crate) fn read_vcf_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_sparse_windowed_with_arrow_variants(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )
}

pub(crate) fn read_vcf_sparse_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_sparse_windowed_with_threads_and_arrow_variants(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
        threads,
    )
}

pub(crate) fn read_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    read_vcf_haplotypes_dense_windowed(path, requested_samples, variant_filter, None, false)
}

pub(crate) fn read_vcf_haplotypes_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_haplotypes_dense_windowed_with_arrow_variants(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_vcf_haplotypes_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_haplotypes_dense_windowed_with_threads_and_arrow_variants(
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

pub(crate) fn read_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    read_vcf_haplotypes_sparse_windowed(path, requested_samples, variant_filter, None)
}

pub(crate) fn read_vcf_haplotypes_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_haplotypes_sparse_windowed_with_arrow_variants(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )
}

pub(crate) fn read_vcf_haplotypes_sparse_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    ::genoio_io::read_vcf_haplotypes_sparse_windowed_with_threads_and_arrow_variants(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
        threads,
    )
}
