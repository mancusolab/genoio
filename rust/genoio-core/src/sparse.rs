// pattern: Functional Core

use crate::{DenseDiagnostics, DenseGenotypeMatrix, MetadataError, SampleRecord, VariantRecord};

#[derive(Debug, Clone, PartialEq)]
pub struct SparseGenotypeMatrix {
    pub n_rows: usize,
    pub n_cols: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<usize>,
    pub data: Vec<f32>,
    pub samples: Vec<SampleRecord>,
    pub variants: Vec<VariantRecord>,
    pub diagnostics: DenseDiagnostics,
}

impl SparseGenotypeMatrix {
    pub fn new(
        n_rows: usize,
        n_cols: usize,
        indptr: Vec<usize>,
        indices: Vec<usize>,
        data: Vec<f32>,
        samples: Vec<SampleRecord>,
        variants: Vec<VariantRecord>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, MetadataError> {
        validate_csc_contract(n_rows, n_cols, &indptr, &indices, &data)?;
        if samples.len() != n_rows {
            return Err(MetadataError::parse(
                "<sparse>",
                format!("sample metadata length {} does not match n_rows {n_rows}", samples.len()),
            ));
        }
        if variants.len() != n_cols {
            return Err(MetadataError::parse(
                "<sparse>",
                format!(
                    "variant metadata length {} does not match n_cols {n_cols}",
                    variants.len()
                ),
            ));
        }

        Ok(Self {
            n_rows,
            n_cols,
            indptr,
            indices,
            data,
            samples,
            variants,
            diagnostics,
        })
    }
}

pub fn sparse_from_dense_minor_flipped(
    dense: DenseGenotypeMatrix,
) -> Result<SparseGenotypeMatrix, MetadataError> {
    reject_sparse_missing(&dense)?;

    let n_rows = dense.n_samples;
    let n_cols = dense.n_variants;
    let mut indptr = Vec::with_capacity(n_cols + 1);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = dense.variants;
    indptr.push(0);

    for col in 0..n_cols {
        let mut a1_count = 0.0_f32;
        for row in 0..n_rows {
            a1_count += dense.values[row * n_cols + col];
        }
        let a0_count = 2.0 * n_rows as f32 - a1_count;
        let flip = a1_count > a0_count;
        if flip {
            flip_variant(&mut variants[col]);
        }

        for row in 0..n_rows {
            let source_value = dense.values[row * n_cols + col];
            let value = if flip { 2.0 - source_value } else { source_value };
            if value != 0.0 {
                indices.push(row);
                data.push(value);
            }
        }
        indptr.push(indices.len());
    }

    SparseGenotypeMatrix::new(
        n_rows,
        n_cols,
        indptr,
        indices,
        data,
        dense.samples,
        variants,
        dense.diagnostics,
    )
}

fn validate_csc_contract(
    n_rows: usize,
    n_cols: usize,
    indptr: &[usize],
    indices: &[usize],
    data: &[f32],
) -> Result<(), MetadataError> {
    if indptr.is_empty() {
        return Err(MetadataError::parse("<sparse>", "sparse indptr must be nonempty"));
    }
    if indptr.len() != n_cols + 1 {
        return Err(MetadataError::parse(
            "<sparse>",
            format!(
                "sparse indptr length {} does not match n_cols + 1 ({})",
                indptr.len(),
                n_cols + 1
            ),
        ));
    }
    if indptr[0] != 0 {
        return Err(MetadataError::parse(
            "<sparse>",
            "sparse first pointer must be zero",
        ));
    }
    for pair in indptr.windows(2) {
        if pair[0] > pair[1] {
            return Err(MetadataError::parse(
                "<sparse>",
                "sparse pointers must be nondecreasing",
            ));
        }
    }
    let terminal_pointer = *indptr.last().expect("indptr is nonempty");
    if terminal_pointer != indices.len() || terminal_pointer != data.len() {
        return Err(MetadataError::parse(
            "<sparse>",
            format!(
                "sparse terminal pointer {terminal_pointer} must equal indices length {} and data length {}",
                indices.len(),
                data.len()
            ),
        ));
    }
    if let Some(row_index) = indices.iter().find(|row_index| **row_index >= n_rows) {
        return Err(MetadataError::parse(
            "<sparse>",
            format!("sparse row index {row_index} is outside n_rows {n_rows}"),
        ));
    }
    Ok(())
}

fn reject_sparse_missing(dense: &DenseGenotypeMatrix) -> Result<(), MetadataError> {
    if dense.missing_mask.iter().any(|missing| *missing) {
        return Err(MetadataError::parse(
            "<sparse>",
            "sparse missing values are not stored in this release",
        ));
    }
    Ok(())
}

fn flip_variant(variant: &mut VariantRecord) {
    std::mem::swap(&mut variant.a0, &mut variant.a1);
    variant.flipped = true;
    if let Some(af) = variant.af {
        variant.af = Some(1.0 - af);
    }
}
