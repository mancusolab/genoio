// pattern: Imperative Shell
//! BGEN reader facade and metadata entry point.
//!
//! The facade exposes the public BGEN API while submodules handle header
//! parsing, index lookup, probability decoding, genotype-filter stats, and dense
//! output assembly.

use std::path::Path;

use crate::Result;
use genoio_core::{MetadataOutput, SourceCapabilities};

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
