// pattern: Functional Core

use crate::{
    DenseDiagnostics, GenoioError, SampleMetadataBuffers, VariantMetadataBuffers, VariantRecord,
};

/// Sparse genotype matrix stored as SciPy-compatible CSC arrays with optional metadata.
///
/// Columns are variants and rows are samples. Python can expose the same data
/// as CSC or convert to CSR after crossing the FFI boundary. `indptr` and
/// `indices` are `i32` because SciPy uses 32-bit sparse indices by default and
/// the Python bridge can transfer these buffers without a widening copy.
#[derive(Debug, Clone, PartialEq)]
pub struct SparseGenotypeMatrix {
    pub n_rows: usize,
    pub n_cols: usize,
    pub indptr: Vec<i32>,
    pub indices: Vec<i32>,
    pub data: Vec<f32>,
    /// Requested sample metadata; `None` means it was omitted.
    pub samples: Option<SampleMetadataBuffers>,
    /// Requested variant metadata; `None` means it was omitted.
    pub variants: Option<VariantMetadataBuffers>,
    pub diagnostics: DenseDiagnostics,
}

impl SparseGenotypeMatrix {
    /// Build a sparse matrix after validating CSC and optional metadata contracts.
    #[expect(
        clippy::too_many_arguments,
        reason = "constructor mirrors validated CSC matrix fields"
    )]
    pub fn new(
        n_rows: usize,
        n_cols: usize,
        indptr: Vec<i32>,
        indices: Vec<i32>,
        data: Vec<f32>,
        samples: Option<SampleMetadataBuffers>,
        variants: Option<VariantMetadataBuffers>,
        diagnostics: DenseDiagnostics,
    ) -> Result<Self, GenoioError> {
        validate_csc_contract(n_rows, n_cols, &indptr, &indices, &data)?;
        if let Some(samples) = samples.as_ref() {
            if samples.len() != n_rows {
                return Err(GenoioError::invalid_source(
                    "<sparse>",
                    format!(
                        "sample metadata length {} does not match n_rows {n_rows}",
                        samples.len()
                    ),
                ));
            }
        }
        if let Some(variants) = variants.as_ref() {
            if variants.len() != n_cols {
                return Err(GenoioError::invalid_source(
                    "<sparse>",
                    format!(
                        "variant metadata length {} does not match n_cols {n_cols}",
                        variants.len()
                    ),
                ));
            }
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
///
/// This is the common path for genotype sparse reads: parsers may decode one
/// retained variant into dense scratch for missing-value policy and allele
/// flipping, then this helper emits only the nonzero rows into CSC storage.
pub fn append_sparse_column(
    indptr: &mut Vec<i32>,
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
    values: &[f32],
) -> Result<(), GenoioError> {
    for (row, value) in values.iter().enumerate() {
        if *value != 0.0 {
            append_sparse_value(indices, data, row, *value)?;
        }
    }
    finish_sparse_column(indptr, data.len())
}

/// Append one nonzero sparse row/value pair after checking the row and nnz range.
///
/// Haplotype paths often already know the nonzero rows, so they call this
/// directly instead of materializing a dense temporary column. Keeping the
/// conversion here makes every sparse emitter use the same `i32` overflow
/// checks before data reaches the PyO3 boundary.
pub fn append_sparse_value(
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
    row: usize,
    value: f32,
) -> Result<(), GenoioError> {
    let row = i32::try_from(row).map_err(|_| sparse_i32_range_error("sparse row index", row))?;
    let next_nnz = data
        .len()
        .checked_add(1)
        .ok_or_else(|| sparse_i32_range_error("sparse nnz", usize::MAX))?;
    if next_nnz > SPARSE_INDEX_MAX {
        return Err(sparse_i32_range_error("sparse nnz", next_nnz));
    }
    indices.push(row);
    data.push(value);
    Ok(())
}

/// Finish a CSC column after checking the cumulative nonzero count range.
///
/// `nnz` becomes the next `indptr` value, so it must fit in the same `i32`
/// index payload as row indices.
pub fn finish_sparse_column(indptr: &mut Vec<i32>, nnz: usize) -> Result<(), GenoioError> {
    let nnz = i32::try_from(nnz).map_err(|_| sparse_i32_range_error("sparse nnz", nnz))?;
    indptr.push(nnz);
    Ok(())
}

fn validate_csc_contract(
    n_rows: usize,
    n_cols: usize,
    indptr: &[i32],
    indices: &[i32],
    data: &[f32],
) -> Result<(), GenoioError> {
    validate_sparse_i32_dimensions(n_rows, n_cols)?;
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
    if indptr.iter().any(|pointer| *pointer < 0) {
        return Err(GenoioError::invalid_source(
            "<sparse>",
            "sparse pointers must be nonnegative",
        ));
    }
    if indices.iter().any(|row_index| *row_index < 0) {
        return Err(GenoioError::invalid_source(
            "<sparse>",
            "sparse row indices must be nonnegative",
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
    let terminal_pointer = usize::try_from(indptr[indptr.len() - 1]).map_err(|_| {
        GenoioError::invalid_source("<sparse>", "sparse terminal pointer must be nonnegative")
    })?;
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
    if let Some(row_index) = indices
        .iter()
        .filter_map(|row_index| {
            usize::try_from(*row_index)
                .ok()
                .map(|row| (*row_index, row))
        })
        .find(|(_, row)| *row >= n_rows)
        .map(|(row_index, _)| row_index)
    {
        return Err(GenoioError::invalid_source(
            "<sparse>",
            format!("sparse row index {row_index} is outside n_rows {n_rows}"),
        ));
    }
    Ok(())
}

const SPARSE_INDEX_MAX: usize = i32::MAX as usize;

// Validate dimensions before inspecting `indptr` so huge shapes fail with the
// sparse-index contract error instead of overflowing `n_cols + 1` checks later.
fn validate_sparse_i32_dimensions(n_rows: usize, n_cols: usize) -> Result<(), GenoioError> {
    if n_rows > SPARSE_INDEX_MAX {
        return Err(sparse_i32_range_error("sparse n_rows", n_rows));
    }
    let indptr_len = n_cols
        .checked_add(1)
        .ok_or_else(|| sparse_i32_range_error("sparse n_cols + 1", usize::MAX))?;
    if indptr_len > SPARSE_INDEX_MAX {
        return Err(sparse_i32_range_error("sparse n_cols + 1", indptr_len));
    }
    Ok(())
}

fn sparse_i32_range_error(label: &str, value: usize) -> GenoioError {
    GenoioError::invalid_source(
        "<sparse>",
        format!("{label} {value} exceeds sparse int32 index range {SPARSE_INDEX_MAX}"),
    )
}

fn mark_variant_flipped(variant: &mut VariantRecord) {
    std::mem::swap(&mut variant.a0, &mut variant.a1);
    variant.flipped = true;
    if let Some(af) = variant.af {
        variant.af = Some(1.0 - af);
    }
}
