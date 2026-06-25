// pattern: Imperative Shell
#![allow(dead_code)]

use std::path::Path;

use genoio_core::{
    DenseGenotypeMatrixArrowVariants, DenseLayout, DenseMissingPolicy,
    SparseGenotypeMatrixArrowVariants, StringColumnBuffers, VariantFilter,
    VariantMetadataArrowBuffers, VariantWindow,
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

pub(crate) fn dense_values_sample_major_arrow(
    matrix: &DenseGenotypeMatrixArrowVariants,
) -> Vec<f32> {
    dense_buffer_sample_major(
        &matrix.values,
        matrix.n_samples,
        matrix.n_variants,
        matrix.layout,
    )
}

pub(crate) fn dense_missing_sample_major_arrow(
    matrix: &DenseGenotypeMatrixArrowVariants,
) -> Vec<bool> {
    let missing = matrix
        .values
        .iter()
        .map(|value| value.is_nan())
        .collect::<Vec<_>>();
    dense_buffer_sample_major(&missing, matrix.n_samples, matrix.n_variants, matrix.layout)
}

pub(crate) fn sparse_values_dense_arrow(matrix: &SparseGenotypeMatrixArrowVariants) -> Vec<f32> {
    let mut dense = vec![0.0; matrix.n_rows * matrix.n_cols];
    for col in 0..matrix.n_cols {
        for offset in matrix.indptr[col]..matrix.indptr[col + 1] {
            let row = matrix.indices[offset];
            dense[row * matrix.n_cols + col] = matrix.data[offset];
        }
    }
    dense
}

pub(crate) fn variants(
    output: &Option<VariantMetadataArrowBuffers>,
) -> &VariantMetadataArrowBuffers {
    output
        .as_ref()
        .expect("variant Arrow buffers should be returned")
}

pub(crate) fn variant_ids(variants: &VariantMetadataArrowBuffers) -> Vec<&str> {
    strings(&variants.ids)
}

pub(crate) fn variant_id(variants: &VariantMetadataArrowBuffers, index: usize) -> &str {
    string_at(&variants.ids, index)
}

pub(crate) fn variant_chrom(variants: &VariantMetadataArrowBuffers, index: usize) -> &str {
    string_at(&variants.chroms, index)
}

pub(crate) fn variant_a0(variants: &VariantMetadataArrowBuffers, index: usize) -> &str {
    string_at(&variants.a0s, index)
}

pub(crate) fn variant_a1(variants: &VariantMetadataArrowBuffers, index: usize) -> &str {
    string_at(&variants.a1s, index)
}

pub(crate) fn string_at(column: &StringColumnBuffers, index: usize) -> &str {
    let start = column.offsets[index] as usize;
    let end = column.offsets[index + 1] as usize;
    std::str::from_utf8(&column.values[start..end]).expect("Arrow string column should be UTF-8")
}

pub(crate) fn strings(column: &StringColumnBuffers) -> Vec<&str> {
    (0..column.len())
        .map(|index| string_at(column, index))
        .collect()
}

fn dense_buffer_sample_major<T: Copy>(
    values: &[T],
    n_samples: usize,
    n_variants: usize,
    layout: DenseLayout,
) -> Vec<T> {
    match layout {
        DenseLayout::SampleMajor => values.to_vec(),
        DenseLayout::VariantMajor => {
            let mut sample_major = Vec::with_capacity(values.len());
            for sample_index in 0..n_samples {
                for variant_index in 0..n_variants {
                    sample_major.push(values[variant_index * n_samples + sample_index]);
                }
            }
            sample_major
        }
    }
}
