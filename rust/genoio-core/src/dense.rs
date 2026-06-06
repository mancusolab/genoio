// pattern: Functional Core

use std::collections::HashSet;
use std::path::Path;

use crate::{GenoioError, SampleRecord, VariantRecord};

/// Counts describing source selection and variant filtering.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DenseDiagnostics {
    pub requested_samples: usize,
    pub retained_samples: usize,
    pub missing_samples: usize,
    pub candidate_variants: usize,
    pub retained_variants: usize,
    pub dropped_metadata_variants: usize,
    pub dropped_genotype_variants: usize,
}

/// Dense genotype matrix in sample-by-variant order.
///
/// `values` and `missing_mask` are both flat sample-major buffers with length
/// `n_samples * n_variants`.
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
    /// Build a dense matrix after validating shape and metadata lengths.
    pub fn new(
        n_samples: usize,
        n_variants: usize,
        values: Vec<f32>,
        missing_mask: Vec<bool>,
        samples: Vec<SampleRecord>,
        variants: Vec<VariantRecord>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        let expected_len = n_samples.checked_mul(n_variants).ok_or_else(|| {
            GenoioError::invalid_source("<dense>", "dense matrix shape is out of range")
        })?;
        if values.len() != expected_len {
            return Err(GenoioError::invalid_source(
                "<dense>",
                format!(
                    "dense values length {} does not match shape {n_samples} x {n_variants}",
                    values.len()
                ),
            ));
        }
        if missing_mask.len() != values.len() {
            return Err(GenoioError::invalid_source(
                "<dense>",
                format!(
                    "dense missing mask length {} does not match values length {}",
                    missing_mask.len(),
                    values.len()
                ),
            ));
        }
        if samples.len() != n_samples {
            return Err(GenoioError::invalid_source(
                "<dense>",
                format!(
                    "sample metadata length {} does not match n_samples {n_samples}",
                    samples.len()
                ),
            ));
        }
        if variants.len() != n_variants {
            return Err(GenoioError::invalid_source(
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

    /// Build a dense matrix when callers intentionally omitted metadata.
    pub fn new_matrix_only(
        n_samples: usize,
        n_variants: usize,
        values: Vec<f32>,
        missing_mask: Vec<bool>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        let expected_len = n_samples.checked_mul(n_variants).ok_or_else(|| {
            GenoioError::invalid_source("<dense>", "dense matrix shape is out of range")
        })?;
        if values.len() != expected_len {
            return Err(GenoioError::invalid_source(
                "<dense>",
                format!(
                    "dense values length {} does not match shape {n_samples} x {n_variants}",
                    values.len()
                ),
            ));
        }
        if missing_mask.len() != values.len() {
            return Err(GenoioError::invalid_source(
                "<dense>",
                format!(
                    "dense missing mask length {} does not match values length {}",
                    missing_mask.len(),
                    values.len()
                ),
            ));
        }

        Ok(Self {
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples: Vec::new(),
            variants: Vec::new(),
            diagnostics,
        })
    }
}

/// Result of applying an optional sample keep list to source sample metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseSampleSelection {
    pub source_indices: Vec<usize>,
    pub samples: Vec<SampleRecord>,
    pub diagnostics: DenseDiagnostics,
}

/// Select requested samples while preserving source row order.
pub fn select_samples_source_order(
    source_samples: &[SampleRecord],
    requested: Option<&[String]>,
    error_path: &Path,
) -> Result<DenseSampleSelection, GenoioError> {
    let Some(requested) = requested else {
        return Ok(DenseSampleSelection {
            source_indices: (0..source_samples.len()).collect(),
            samples: source_samples.to_vec(),
            diagnostics: DenseDiagnostics {
                requested_samples: source_samples.len(),
                retained_samples: source_samples.len(),
                missing_samples: 0,
                ..DenseDiagnostics::default()
            },
        });
    };

    let requested_ids = requested.iter().collect::<HashSet<_>>();
    if requested_ids.len() != requested.len() {
        return Err(GenoioError::invalid_source(
            error_path,
            "sample keep list must not contain duplicate sample IDs",
        ));
    }

    let mut source_indices = Vec::new();
    let mut samples = Vec::new();
    // Keep-list order is only used for membership. Matrix rows remain in
    // source order so VCF/PLINK reads agree and metadata aligns with values.
    for (index, sample) in source_samples.iter().enumerate() {
        if requested_ids.contains(&sample.iid) {
            source_indices.push(index);
            samples.push(sample.clone());
        }
    }

    let missing_samples = requested.len() - samples.len();
    if missing_samples > 0 {
        return Err(GenoioError::sample_filter(
            requested.len(),
            samples.len(),
            missing_samples,
        ));
    }

    Ok(DenseSampleSelection {
        source_indices,
        samples,
        diagnostics: DenseDiagnostics {
            requested_samples: requested.len(),
            retained_samples: requested.len(),
            missing_samples: 0,
            ..DenseDiagnostics::default()
        },
    })
}

/// Transpose a flat variant-major buffer into sample-major order.
pub fn transpose_variant_major_to_sample_major<T: Copy>(
    values: &[T],
    n_samples: usize,
    n_variants: usize,
) -> Vec<T> {
    // Readers append one variant at a time because source formats are
    // variant-major. Python callers expect sample rows, so transpose once at
    // the core boundary instead of reshaping incorrectly downstream.
    let mut transposed = Vec::with_capacity(values.len());
    for sample_index in 0..n_samples {
        for variant_index in 0..n_variants {
            transposed.push(values[variant_index * n_samples + sample_index]);
        }
    }
    transposed
}
