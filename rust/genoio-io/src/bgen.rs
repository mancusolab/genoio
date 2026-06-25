// pattern: Imperative Shell
//! BGEN reader facade and metadata entry point.
//!
//! The facade exposes the public BGEN API while submodules handle header
//! parsing, index lookup, probability decoding, genotype-filter stats, and dense
//! output assembly.

use std::path::Path;

use crate::Result;
use genoio_core::{MetadataArrowOutput, SourceCapabilities};

mod decode;
mod dense;
mod filter;
mod haplotype;
mod header;
mod index;
mod io;
mod session;

pub use dense::read_bgen_dosage_dense_windowed_with_arrow_variants;
pub use haplotype::read_bgen_haplotypes_dosage_dense_windowed_with_arrow_variants;

use session::BgenReadSession;

/// Read BGEN metadata with variant metadata staged as Arrow-compatible buffers.
pub fn read_bgen_metadata_arrow(bgen: &Path, sample: Option<&Path>) -> Result<MetadataArrowOutput> {
    let mut session = BgenReadSession::open(bgen)?;
    let samples = session.read_samples(sample)?;
    session.seek_to_variants()?;
    let variants = session.read_all_variant_metadata_arrow()?;

    Ok(MetadataArrowOutput {
        samples,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}
