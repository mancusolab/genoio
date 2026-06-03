// pattern: Functional Core

//! Shared genotype metadata, filter, and matrix contracts.
//!
//! This crate is intentionally free of file IO and Python bindings. Reader
//! crates populate these validated containers, and the PyO3 layer converts
//! them into NumPy, SciPy, and Polars objects.

pub mod capabilities;
pub mod dense;
pub mod error;
pub mod filter;
pub mod metadata;
pub mod sparse;

pub use capabilities::SourceCapabilities;
pub use dense::{
    select_samples_source_order, transpose_variant_major_to_sample_major, DenseDiagnostics,
    DenseGenotypeMatrix, DenseSampleSelection,
};
pub use error::MetadataError;
pub use filter::{
    attach_variant_stats, compute_variant_stats, PartialFilterDecision, RegionPredicate,
    VariantFilter, VariantStats, VariantWindow,
};
pub use metadata::{MetadataOutput, SampleRecord, VariantRecord};
pub use sparse::{
    append_sparse_column, flip_values_to_minor_allele, reject_sparse_missing_values,
    SparseGenotypeMatrix,
};

/// Python package and Rust workspace package name.
pub const PACKAGE_NAME: &str = "genoio";
/// Cargo package version compiled into the extension.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// Build profile marker exposed for diagnostics.
pub const BUILD_PROFILE: &str = if cfg!(debug_assertions) {
    "debug"
} else {
    "release"
};
