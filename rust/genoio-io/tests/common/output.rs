// pattern: Functional Core
#![allow(dead_code)]

use genoio_core::{
    DenseGenotypeMatrix, DenseLayout, SparseGenotypeMatrix, StringColumnBuffers,
    VariantMetadataBuffers,
};

pub(crate) fn dense_values_sample_major_output(matrix: &DenseGenotypeMatrix) -> Vec<f32> {
    dense_buffer_sample_major(
        &matrix.values,
        matrix.n_samples,
        matrix.n_variants,
        matrix.layout,
    )
}

pub(crate) fn dense_missing_sample_major_output(matrix: &DenseGenotypeMatrix) -> Vec<bool> {
    let missing = matrix
        .values
        .iter()
        .map(|value| value.is_nan())
        .collect::<Vec<_>>();
    dense_buffer_sample_major(&missing, matrix.n_samples, matrix.n_variants, matrix.layout)
}

pub(crate) fn sparse_values_dense_output(matrix: &SparseGenotypeMatrix) -> Vec<f32> {
    let mut dense = vec![0.0; matrix.n_rows * matrix.n_cols];
    for col in 0..matrix.n_cols {
        let start = usize::try_from(matrix.indptr[col]).expect("sparse pointer is nonnegative");
        let end = usize::try_from(matrix.indptr[col + 1]).expect("sparse pointer is nonnegative");
        for offset in start..end {
            let row = usize::try_from(matrix.indices[offset]).expect("sparse row is nonnegative");
            dense[row * matrix.n_cols + col] = matrix.data[offset];
        }
    }
    dense
}

pub(crate) fn variants(output: &Option<VariantMetadataBuffers>) -> &VariantMetadataBuffers {
    output
        .as_ref()
        .expect("variant metadata buffers should be returned")
}

pub(crate) fn variant_ids(variants: &VariantMetadataBuffers) -> Vec<&str> {
    strings(&variants.ids)
}

pub(crate) fn variant_id(variants: &VariantMetadataBuffers, index: usize) -> &str {
    string_at(&variants.ids, index)
}

pub(crate) fn variant_chrom(variants: &VariantMetadataBuffers, index: usize) -> &str {
    string_at(&variants.chroms, index)
}

pub(crate) fn variant_a0(variants: &VariantMetadataBuffers, index: usize) -> &str {
    string_at(&variants.a0s, index)
}

pub(crate) fn variant_a1(variants: &VariantMetadataBuffers, index: usize) -> &str {
    string_at(&variants.a1s, index)
}

pub(crate) fn string_at(column: &StringColumnBuffers, index: usize) -> &str {
    let start = column.offsets[index] as usize;
    let end = column.offsets[index + 1] as usize;
    std::str::from_utf8(&column.values[start..end]).expect("string column should be UTF-8")
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
