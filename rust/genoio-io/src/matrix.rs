// pattern: Functional Core
//! Matrix assembly helpers shared by format readers.
//!
//! Format modules decode variants in the order that is cheapest for their
//! source files. These helpers validate dense buffers and preserve explicit
//! missing-value behavior without coupling readers to row metadata assembly.

use genoio_core::{DenseMissingPolicy, GenoioError};

use crate::error::Result;

pub(crate) fn shrink_sample_major_width<T: Copy>(
    values: &mut Vec<T>,
    n_samples: usize,
    old_width: usize,
    new_width: usize,
) {
    debug_assert!(new_width <= old_width);
    if old_width == new_width {
        return;
    }
    for sample_index in 1..n_samples {
        let source_start = sample_index * old_width;
        let target_start = sample_index * new_width;
        values.copy_within(source_start..source_start + new_width, target_start);
    }
    values.truncate(n_samples * new_width);
}

pub(crate) fn apply_dense_missing_policy_to_variant(
    values: &mut [f32],
    missing_indices: &[usize],
    policy: DenseMissingPolicy,
) -> Result<()> {
    if missing_indices.is_empty() {
        return Ok(());
    }
    validate_missing_indices(values.len(), missing_indices)?;
    match policy {
        DenseMissingPolicy::Raise => Err(GenoioError::missing_data(
            "missing genotype calls are present in retained data",
        )),
        DenseMissingPolicy::Nan => {
            for &index in missing_indices {
                values[index] = f32::NAN;
            }
            Ok(())
        }
        DenseMissingPolicy::Impute => {
            let called_count = values.len() - missing_indices.len();
            if called_count == 0 {
                return Err(GenoioError::missing_data(
                    "cannot impute all-missing variant",
                ));
            }
            let called_sum = sum_called_values(values, missing_indices);
            let mean = called_sum / called_count as f32;
            for &index in missing_indices {
                values[index] = mean;
            }
            Ok(())
        }
    }
}

fn sum_called_values(values: &[f32], missing_indices: &[usize]) -> f32 {
    let mut called_sum = 0.0_f32;
    let mut missing_cursor = 0_usize;
    for (index, &value) in values.iter().enumerate() {
        if missing_indices
            .get(missing_cursor)
            .is_some_and(|&missing_index| missing_index == index)
        {
            missing_cursor += 1;
        } else {
            called_sum += value;
        }
    }
    called_sum
}

fn validate_missing_indices(values_len: usize, missing_indices: &[usize]) -> Result<()> {
    let mut previous = None;
    for &index in missing_indices {
        if index >= values_len {
            return Err(GenoioError::internal_contract(
                "dense missing index is outside variant values",
            ));
        }
        if previous.is_some_and(|previous| index <= previous) {
            return Err(GenoioError::internal_contract(
                "dense missing indices must be sorted and unique",
            ));
        }
        previous = Some(index);
    }
    Ok(())
}

/// Write a variant-major vector into one column of preallocated sample-major
/// dense buffers.
pub(crate) fn write_sample_major_variant_slot(
    values: &mut [f32],
    n_samples: usize,
    row_width: usize,
    variant_index: usize,
    variant_values: &[f32],
) -> Result<()> {
    if variant_values.len() != n_samples {
        return Err(GenoioError::internal_contract(format!(
            "variant value count {} does not match sample count {n_samples}",
            variant_values.len(),
        )));
    }
    let expected_len = n_samples.checked_mul(row_width).ok_or_else(|| {
        GenoioError::internal_contract("sample-major dense matrix shape is out of range")
    })?;
    if values.len() != expected_len {
        return Err(GenoioError::internal_contract(
            "sample-major dense buffer does not match declared shape",
        ));
    }
    if variant_index >= row_width {
        return Err(GenoioError::internal_contract(
            "sample-major variant index is outside row width",
        ));
    }

    // Readers accept variants one at a time; this fills the corresponding
    // column in each preallocated sample row without a later full transpose.
    for (sample_index, &value) in variant_values.iter().enumerate().take(n_samples) {
        let target_index = sample_index * row_width + variant_index;
        values[target_index] = value;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn apply_dense_missing_policy_to_variant_imputes_missing_indices() {
        let mut values = vec![2.0, f32::NAN, 6.0];

        apply_dense_missing_policy_to_variant(&mut values, &[1], DenseMissingPolicy::Impute)
            .expect("single missing value should be imputed");

        assert_eq!(values, vec![2.0, 4.0, 6.0]);
    }

    #[test]
    fn apply_dense_missing_policy_to_variant_rejects_unsorted_indices() {
        let mut values = vec![0.0, 1.0, 2.0];
        let error =
            apply_dense_missing_policy_to_variant(&mut values, &[2, 1], DenseMissingPolicy::Nan)
                .expect_err("unsorted indices should fail");

        assert!(error
            .to_string()
            .contains("dense missing indices must be sorted and unique"));
    }

    #[test]
    fn write_sample_major_variant_slot_writes_one_column_per_sample() {
        let mut values = vec![0.0; 6];

        write_sample_major_variant_slot(&mut values, 2, 3, 1, &[4.0, 5.0])
            .expect("column write should succeed");

        assert_eq!(values, vec![0.0, 4.0, 0.0, 0.0, 5.0, 0.0]);
    }
}
