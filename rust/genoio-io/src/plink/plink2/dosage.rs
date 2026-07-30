// pattern: Imperative Shell
//! Dense dosage PLINK2 read orchestration.
//!
//! Dosage reads start from hard-call main-track values and apply PGEN dosage
//! overlays before filter evaluation. Output is staged variant-major because
//! dosage overlays naturally decode one retained variant at a time.

use std::path::Path;

use genoio_core::{DenseGenotypeMatrix, DenseMissingPolicy, VariantFilter, VariantWindow};

use crate::error::Result;

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors dense dosage read options plus metadata return choices"
)]
pub fn read_plink2_dosage_dense_windowed(
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
            matrix_kind: crate::blocks::MatrixKind::Genotype,
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
            "PLINK2 dense dosage session returned sparse output",
        ));
    };
    Ok(output)
}
