// pattern: Imperative Shell
//! Haplotype PLINK2 read orchestration.
//!
//! Haplotype readers use PGEN phase or phased-dosage auxiliary tracks and emit
//! two rows per selected sample. Sparse hard-call output rejects retained
//! missing haplotypes because sparse CSC has no missing-value channel.

use std::path::Path;

use genoio_core::{
    append_sparse_value, finish_sparse_column, reject_sparse_missing, DenseGenotypeMatrix,
    DenseMissingPolicy, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};

use crate::error::Result;

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors haplotype dense read options plus metadata return choices"
)]
pub fn read_plink2_haplotypes_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    let output = super::session::read_windowed(
        pgen,
        pvar,
        psam,
        crate::blocks::BlockReadOptions {
            matrix_kind: crate::blocks::MatrixKind::Haplotype,
            sparse: false,
            requested_samples: requested_samples.map(<[String]>::to_vec),
            variant_filter: variant_filter.cloned(),
            dosage_source: crate::blocks::DosageSource::Hardcall,
            missing_policy,
            return_samples,
            return_variants,
        },
        variant_window,
    )?;
    let crate::blocks::BlockOutput::Dense(output) = output else {
        return Err(genoio_core::GenoioError::internal_contract(
            "PLINK2 dense hardcall haplotype session returned sparse output",
        ));
    };
    Ok(output)
}

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors haplotype dosage read options plus metadata return choices"
)]
pub fn read_plink2_haplotypes_dosage_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    let output = super::session::read_windowed(
        pgen,
        pvar,
        psam,
        crate::blocks::BlockReadOptions {
            matrix_kind: crate::blocks::MatrixKind::Haplotype,
            sparse: false,
            requested_samples: requested_samples.map(<[String]>::to_vec),
            variant_filter: variant_filter.cloned(),
            dosage_source: crate::blocks::DosageSource::Dosage,
            missing_policy,
            return_samples,
            return_variants,
        },
        variant_window,
    )?;
    let crate::blocks::BlockOutput::Dense(output) = output else {
        return Err(genoio_core::GenoioError::internal_contract(
            "PLINK2 dense dosage haplotype session returned sparse output",
        ));
    };
    Ok(output)
}

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors haplotype sparse read options plus metadata return choices"
)]
pub fn read_plink2_haplotypes_sparse_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    let output = super::session::read_windowed(
        pgen,
        pvar,
        psam,
        crate::blocks::BlockReadOptions {
            matrix_kind: crate::blocks::MatrixKind::Haplotype,
            sparse: true,
            requested_samples: requested_samples.map(<[String]>::to_vec),
            variant_filter: variant_filter.cloned(),
            dosage_source: crate::blocks::DosageSource::Hardcall,
            missing_policy: DenseMissingPolicy::Raise,
            return_samples,
            return_variants,
        },
        variant_window,
    )?;
    let crate::blocks::BlockOutput::Sparse(output) = output else {
        return Err(genoio_core::GenoioError::internal_contract(
            "PLINK2 sparse hardcall haplotype session returned dense output",
        ));
    };
    Ok(output)
}

#[inline]
pub(super) fn append_haplotype_sparse_column(
    indptr: &mut Vec<i32>,
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
    values: &[f32],
    missing_indices: &[usize],
) -> Result<()> {
    reject_sparse_missing(!missing_indices.is_empty())?;
    // Haplotype decode already returns row-major values, so append nonzeros
    // directly while reusing the core sparse index range checks.
    for (row, &value) in values.iter().enumerate() {
        if value != 0.0 {
            append_sparse_value(indices, data, row, value)?;
        }
    }
    finish_sparse_column(indptr, data.len())
}
