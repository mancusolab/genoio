// pattern: Functional Core

use std::collections::HashSet;
use std::path::Path;

use crate::{GenoioError, SampleRecord, VariantMetadataArrowBuffers, VariantRecord};

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

/// Flat buffer layout used by a dense genotype matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DenseLayout {
    /// Flat values are stored in public sample-by-variant order.
    SampleMajor,
    /// Flat values are stored as one complete variant after another.
    VariantMajor,
}

/// Dense missing-call policy applied while readers build matrix values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DenseMissingPolicy {
    /// Reject retained missing calls before returning the matrix.
    Raise,
    /// Store retained missing calls as `NaN` in the value buffer.
    Nan,
    /// Replace retained missing calls with the retained variant mean.
    Impute,
}

/// Dense genotype matrix with layout-tagged flat buffers.
///
/// `values` has length `n_samples * n_variants`. Consumers that read flat
/// buffers directly must inspect `layout`; Python assembly converts either
/// layout into the public sample-by-variant array shape.
#[derive(Debug, Clone, PartialEq)]
pub struct DenseGenotypeMatrix {
    pub n_samples: usize,
    pub n_variants: usize,
    pub values: Vec<f32>,
    pub layout: DenseLayout,
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
        samples: Vec<SampleRecord>,
        variants: Vec<VariantRecord>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        Self::new_with_layout(
            n_samples,
            n_variants,
            values,
            DenseLayout::SampleMajor,
            samples,
            variants,
            diagnostics,
        )
    }

    /// Build a dense matrix with an explicit flat-buffer layout.
    pub fn new_with_layout(
        n_samples: usize,
        n_variants: usize,
        values: Vec<f32>,
        layout: DenseLayout,
        samples: Vec<SampleRecord>,
        variants: Vec<VariantRecord>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        validate_dense_values(n_samples, n_variants, values.len())?;
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
            layout,
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
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        Self::new_matrix_only_with_layout(
            n_samples,
            n_variants,
            values,
            DenseLayout::SampleMajor,
            diagnostics,
        )
    }

    /// Build a dense matrix with explicit layout when callers omitted metadata.
    pub fn new_matrix_only_with_layout(
        n_samples: usize,
        n_variants: usize,
        values: Vec<f32>,
        layout: DenseLayout,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        validate_dense_values(n_samples, n_variants, values.len())?;

        Ok(Self {
            n_samples,
            n_variants,
            values,
            layout,
            samples: Vec::new(),
            variants: Vec::new(),
            diagnostics,
        })
    }
}

/// Dense genotype matrix with public variants staged in Arrow-compatible buffers.
#[derive(Debug, Clone, PartialEq)]
pub struct DenseGenotypeMatrixArrowVariants {
    pub n_samples: usize,
    pub n_variants: usize,
    pub values: Vec<f32>,
    pub layout: DenseLayout,
    pub samples: Vec<SampleRecord>,
    pub variants: Option<VariantMetadataArrowBuffers>,
    pub diagnostics: DenseDiagnostics,
}

impl DenseGenotypeMatrixArrowVariants {
    /// Convert a legacy dense matrix result into the columnar variant payload used by PyO3.
    pub fn from_matrix(
        matrix: DenseGenotypeMatrix,
        return_variants: bool,
    ) -> Result<Self, GenoioError> {
        let variants = if return_variants {
            Some(VariantMetadataArrowBuffers::from_records(&matrix.variants)?)
        } else {
            None
        };
        Self::new_with_layout(
            matrix.n_samples,
            matrix.n_variants,
            matrix.values,
            matrix.layout,
            matrix.samples,
            variants,
            matrix.diagnostics,
        )
    }

    /// Convert Arrow-buffered variant metadata back to the legacy row matrix shape.
    pub fn into_matrix(self) -> Result<DenseGenotypeMatrix, GenoioError> {
        let Self {
            n_samples,
            n_variants,
            values,
            layout,
            samples,
            variants,
            diagnostics,
        } = self;
        match variants {
            Some(variants) => DenseGenotypeMatrix::new_with_layout(
                n_samples,
                n_variants,
                values,
                layout,
                samples,
                variants.into_records()?,
                diagnostics,
            ),
            None if samples.is_empty() => DenseGenotypeMatrix::new_matrix_only_with_layout(
                n_samples,
                n_variants,
                values,
                layout,
                diagnostics,
            ),
            None => Err(GenoioError::internal_contract(
                "dense Arrow matrix cannot convert to row metadata without variant buffers",
            )),
        }
    }

    /// Build a dense matrix with optional sample metadata and optional Arrow variant metadata.
    pub fn new_with_layout(
        n_samples: usize,
        n_variants: usize,
        values: Vec<f32>,
        layout: DenseLayout,
        samples: Vec<SampleRecord>,
        variants: Option<VariantMetadataArrowBuffers>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        validate_dense_values(n_samples, n_variants, values.len())?;
        if !samples.is_empty() && samples.len() != n_samples {
            return Err(GenoioError::invalid_source(
                "<dense>",
                format!(
                    "sample metadata length {} does not match n_samples {n_samples}",
                    samples.len()
                ),
            ));
        }
        if let Some(variants) = variants.as_ref() {
            if variants.len() != n_variants {
                return Err(GenoioError::invalid_source(
                    "<dense>",
                    format!(
                        "variant metadata length {} does not match n_variants {n_variants}",
                        variants.len()
                    ),
                ));
            }
        }

        Ok(Self {
            n_samples,
            n_variants,
            values,
            layout,
            samples,
            variants,
            diagnostics,
        })
    }
}

fn validate_dense_values(
    n_samples: usize,
    n_variants: usize,
    values_len: usize,
) -> Result<(), GenoioError> {
    let expected_len = n_samples.checked_mul(n_variants).ok_or_else(|| {
        GenoioError::invalid_source("<dense>", "dense matrix shape is out of range")
    })?;
    if values_len != expected_len {
        return Err(GenoioError::invalid_source(
            "<dense>",
            format!(
                "dense values length {values_len} does not match shape {n_samples} x {n_variants}",
            ),
        ));
    }
    Ok(())
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
    // Keep this helper available for callers that truly need a physically
    // sample-major buffer. The Python bridge usually keeps variant-major
    // buffers and exposes the public shape with NumPy strides instead.
    let mut transposed = Vec::with_capacity(values.len());
    for sample_index in 0..n_samples {
        for variant_index in 0..n_variants {
            transposed.push(values[variant_index * n_samples + sample_index]);
        }
    }
    transposed
}
