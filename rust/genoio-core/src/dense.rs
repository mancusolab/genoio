// pattern: Functional Core

use std::collections::HashSet;
use std::path::Path;

use crate::{MetadataError, SampleRecord, VariantRecord};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DenseDiagnostics {
    pub requested_samples: usize,
    pub retained_samples: usize,
    pub missing_samples: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DenseGenotypeMatrix {
    pub n_samples: usize,
    pub n_variants: usize,
    pub values: Vec<f32>,
    pub missing_mask: Vec<bool>,
    pub samples: Vec<SampleRecord>,
    pub variants: Vec<VariantRecord>,
    pub diagnostics: DenseDiagnostics,
}

impl DenseGenotypeMatrix {
    pub fn new(
        n_samples: usize,
        n_variants: usize,
        values: Vec<f32>,
        missing_mask: Vec<bool>,
        samples: Vec<SampleRecord>,
        variants: Vec<VariantRecord>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, MetadataError> {
        let expected_len = n_samples
            .checked_mul(n_variants)
            .ok_or_else(|| MetadataError::parse("<dense>", "dense matrix shape is out of range"))?;
        if values.len() != expected_len {
            return Err(MetadataError::parse(
                "<dense>",
                format!(
                    "dense values length {} does not match shape {n_samples} x {n_variants}",
                    values.len()
                ),
            ));
        }
        if missing_mask.len() != values.len() {
            return Err(MetadataError::parse(
                "<dense>",
                format!(
                    "dense missing mask length {} does not match values length {}",
                    missing_mask.len(),
                    values.len()
                ),
            ));
        }
        if samples.len() != n_samples {
            return Err(MetadataError::parse(
                "<dense>",
                format!(
                    "sample metadata length {} does not match n_samples {n_samples}",
                    samples.len()
                ),
            ));
        }
        if variants.len() != n_variants {
            return Err(MetadataError::parse(
                "<dense>",
                format!(
                    "variant metadata length {} does not match n_variants {n_variants}",
                    variants.len()
                ),
            ));
        }

        Ok(Self {
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples,
            variants,
            diagnostics,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseSampleSelection {
    pub source_indices: Vec<usize>,
    pub samples: Vec<SampleRecord>,
    pub diagnostics: DenseDiagnostics,
}

pub fn select_samples_source_order(
    source_samples: &[SampleRecord],
    requested: Option<&[String]>,
    error_path: &Path,
) -> Result<DenseSampleSelection, MetadataError> {
    let Some(requested) = requested else {
        return Ok(DenseSampleSelection {
            source_indices: (0..source_samples.len()).collect(),
            samples: source_samples.to_vec(),
            diagnostics: DenseDiagnostics {
                requested_samples: source_samples.len(),
                retained_samples: source_samples.len(),
                missing_samples: 0,
            },
        });
    };

    let requested_ids = requested.iter().collect::<HashSet<_>>();
    if requested_ids.len() != requested.len() {
        return Err(MetadataError::parse(
            error_path,
            "sample keep list must not contain duplicate sample IDs",
        ));
    }

    let mut source_indices = Vec::new();
    let mut samples = Vec::new();
    for (index, sample) in source_samples.iter().enumerate() {
        if requested_ids.contains(&sample.iid) {
            source_indices.push(index);
            samples.push(sample.clone());
        }
    }

    let missing_samples = requested.len() - samples.len();
    if missing_samples > 0 {
        return Err(MetadataError::parse(
            error_path,
            format!(
                "missing requested sample(s): requested={} retained={} missing={missing_samples}",
                requested.len(),
                samples.len()
            ),
        ));
    }

    Ok(DenseSampleSelection {
        source_indices,
        samples,
        diagnostics: DenseDiagnostics {
            requested_samples: requested.len(),
            retained_samples: requested.len(),
            missing_samples: 0,
        },
    })
}

pub fn transpose_variant_major_to_sample_major<T: Copy>(
    values: &[T],
    n_samples: usize,
    n_variants: usize,
) -> Vec<T> {
    let mut transposed = Vec::with_capacity(values.len());
    for sample_index in 0..n_samples {
        for variant_index in 0..n_variants {
            transposed.push(values[variant_index * n_samples + sample_index]);
        }
    }
    transposed
}
