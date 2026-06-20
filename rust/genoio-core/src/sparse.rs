// pattern: Functional Core

use crate::{DenseDiagnostics, GenoioError, SampleRecord, VariantRecord};

/// Sparse genotype matrix stored as CSC arrays.
///
/// Columns are variants and rows are samples. Python can expose the same data
/// as CSC or convert to CSR after crossing the FFI boundary.
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
    /// Build a sparse matrix after validating CSC and metadata contracts.
    #[expect(
        clippy::too_many_arguments,
        reason = "constructor mirrors validated CSC matrix fields"
    )]
    pub fn new(
        n_rows: usize,
        n_cols: usize,
        indptr: Vec<usize>,
        indices: Vec<usize>,
        data: Vec<f32>,
        samples: Vec<SampleRecord>,
        variants: Vec<VariantRecord>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        validate_csc_contract(n_rows, n_cols, &indptr, &indices, &data)?;
        if samples.len() != n_rows {
            return Err(GenoioError::invalid_source(
                "<sparse>",
                format!(
                    "sample metadata length {} does not match n_rows {n_rows}",
                    samples.len()
                ),
            ));
        }
        if variants.len() != n_cols {
            return Err(GenoioError::invalid_source(
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

/// Reject retained missing calls before constructing sparse matrices.
pub fn reject_sparse_missing_values(missing: &[bool]) -> Result<(), GenoioError> {
    reject_sparse_missing(missing.iter().any(|value| *value))
}

/// Reject a precomputed missing-call flag before sparse matrix construction.
pub fn reject_sparse_missing(has_missing: bool) -> Result<(), GenoioError> {
    if !has_missing {
        return Ok(());
    }
    Err(GenoioError::missing_data(
        "sparse missing values are not stored in this release",
    ))
}

/// Flip a biallelic 0/1/2 column so stored nonzeros represent the minor allele.
pub fn flip_values_to_minor_allele(values: &mut [f32], variant: &mut VariantRecord) {
    let a1_count = values.iter().sum::<f32>();
    let a0_count = 2.0 * values.len() as f32 - a1_count;
    if a1_count <= a0_count {
        return;
    }
    for value in values {
        *value = 2.0 - *value;
    }
    mark_variant_flipped(variant);
}

/// Flip a biallelic haplotype indicator column to encode the minor allele.
pub fn flip_haplotype_values_to_minor_allele(values: &mut [f32], variant: &mut VariantRecord) {
    let a1_count = values.iter().sum::<f32>();
    let a0_count = values.len() as f32 - a1_count;
    if a1_count <= a0_count {
        return;
    }
    for value in values {
        *value = 1.0 - *value;
    }
    mark_variant_flipped(variant);
}

/// Return true when allele 1 is the major allele among called haplotypes.
pub fn should_flip_haplotype_to_minor_allele(a1_count: usize, n_haplotypes: usize) -> bool {
    a1_count > n_haplotypes.saturating_sub(a1_count)
}

/// Flip only variant metadata after a backend emitted complement allele rows.
pub fn flip_variant_metadata_to_minor_allele(variant: &mut VariantRecord) {
    mark_variant_flipped(variant);
}

/// Append one dense variant column to CSC buffers, skipping zero entries.
pub fn append_sparse_column(
    indptr: &mut Vec<usize>,
    indices: &mut Vec<usize>,
    data: &mut Vec<f32>,
    values: &[f32],
) {
    for (row, value) in values.iter().enumerate() {
        if *value != 0.0 {
            indices.push(row);
            data.push(*value);
        }
    }
    indptr.push(indices.len());
}

fn validate_csc_contract(
    n_rows: usize,
    n_cols: usize,
    indptr: &[usize],
    indices: &[usize],
    data: &[f32],
) -> Result<(), GenoioError> {
    if indptr.is_empty() {
        return Err(GenoioError::invalid_source(
            "<sparse>",
            "sparse indptr must be nonempty",
        ));
    }
    if indptr.len() != n_cols + 1 {
        return Err(GenoioError::invalid_source(
            "<sparse>",
            format!(
                "sparse indptr length {} does not match n_cols + 1 ({})",
                indptr.len(),
                n_cols + 1
            ),
        ));
    }
    if indptr[0] != 0 {
        return Err(GenoioError::invalid_source(
            "<sparse>",
            "sparse first pointer must be zero",
        ));
    }
    for pair in indptr.windows(2) {
        if pair[0] > pair[1] {
            return Err(GenoioError::invalid_source(
                "<sparse>",
                "sparse pointers must be nondecreasing",
            ));
        }
    }
    let terminal_pointer = indptr[indptr.len() - 1];
    if terminal_pointer != indices.len() || terminal_pointer != data.len() {
        return Err(GenoioError::invalid_source(
            "<sparse>",
            format!(
                "sparse terminal pointer {terminal_pointer} must equal indices length {} and data length {}",
                indices.len(),
                data.len()
            ),
        ));
    }
    if let Some(row_index) = indices.iter().find(|row_index| **row_index >= n_rows) {
        return Err(GenoioError::invalid_source(
            "<sparse>",
            format!("sparse row index {row_index} is outside n_rows {n_rows}"),
        ));
    }
    Ok(())
}

fn mark_variant_flipped(variant: &mut VariantRecord) {
    std::mem::swap(&mut variant.a0, &mut variant.a1);
    variant.flipped = true;
    if let Some(af) = variant.af {
        variant.af = Some(1.0 - af);
    }
}
