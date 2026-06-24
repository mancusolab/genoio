// pattern: Imperative Shell
//! BGEN reader facade and metadata entry point.
//!
//! The facade exposes the public BGEN API while submodules handle header
//! parsing, index lookup, probability decoding, genotype-filter stats, and dense
//! output assembly.

use std::path::Path;

use crate::Result;
use genoio_core::{
    DenseGenotypeMatrixArrowVariants, DenseMissingPolicy, MetadataArrowOutput, MetadataOutput,
    SourceCapabilities, VariantFilter, VariantWindow,
};

use crate::matrix::dense_matrix_to_arrow_variants;

mod decode;
mod dense;
mod filter;
mod haplotype;
mod header;
mod index;
mod io;
mod session;

pub use dense::{
    read_bgen_dosage_dense, read_bgen_dosage_dense_windowed,
    read_bgen_dosage_dense_windowed_with_missing_policy,
};
pub use haplotype::{
    read_bgen_haplotypes_dosage_dense_windowed,
    read_bgen_haplotypes_dosage_dense_windowed_with_missing_policy,
};

use session::BgenReadSession;

/// Read BGEN sample and variant metadata without returning dosages.
pub fn read_bgen_metadata(bgen: &Path, sample: Option<&Path>) -> Result<MetadataOutput> {
    let mut session = BgenReadSession::open(bgen)?;
    let samples = session.read_samples(sample)?;
    session.seek_to_variants()?;
    let variants = session.read_all_variant_metadata()?;

    Ok(MetadataOutput {
        samples,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}

/// Read BGEN metadata with variant metadata staged as Arrow-compatible buffers.
pub fn read_bgen_metadata_arrow(bgen: &Path, sample: Option<&Path>) -> Result<MetadataArrowOutput> {
    read_bgen_metadata(bgen, sample).and_then(MetadataArrowOutput::from_metadata)
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors dense dosage read options plus metadata return choices"
)]
pub fn read_bgen_dosage_dense_windowed_with_arrow_variants(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let matrix_only = !return_samples && !return_variants;
    read_bgen_dosage_dense_windowed_with_missing_policy(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        matrix_only,
    )
    .and_then(|matrix| dense_matrix_to_arrow_variants(matrix, return_samples, return_variants))
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors haplotype dosage read options plus metadata return choices"
)]
pub fn read_bgen_haplotypes_dosage_dense_windowed_with_arrow_variants(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let matrix_only = !return_samples && !return_variants;
    read_bgen_haplotypes_dosage_dense_windowed_with_missing_policy(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        matrix_only,
    )
    .and_then(|matrix| dense_matrix_to_arrow_variants(matrix, return_samples, return_variants))
}
