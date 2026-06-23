// pattern: Functional Core

use genoio_core::{DenseGenotypeMatrix, DenseLayout};

pub fn dense_values_sample_major(matrix: &DenseGenotypeMatrix) -> Vec<f32> {
    dense_buffer_sample_major(
        &matrix.values,
        matrix.n_samples,
        matrix.n_variants,
        matrix.layout,
    )
}

pub fn dense_missing_sample_major(matrix: &DenseGenotypeMatrix) -> Vec<bool> {
    dense_buffer_sample_major(
        &matrix.missing_mask,
        matrix.n_samples,
        matrix.n_variants,
        matrix.layout,
    )
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
            // Integration tests assert the public matrix contract. Production
            // readers may keep variant-major buffers internally to avoid a copy.
            for sample_index in 0..n_samples {
                for variant_index in 0..n_variants {
                    sample_major.push(values[variant_index * n_samples + sample_index]);
                }
            }
            sample_major
        }
    }
}
