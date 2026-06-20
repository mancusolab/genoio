// pattern: Functional Core

#![cfg_attr(
    not(test),
    deny(clippy::expect_used, clippy::panic, clippy::unwrap_used)
)]

//! Shared genotype metadata, filter, and matrix contracts.
//!
//! This crate is intentionally free of file IO and Python bindings. Reader
//! crates populate these validated containers, and the PyO3 layer converts
//! them into NumPy, SciPy, and Polars objects.
//!
//! The public Rust items here are workspace contracts between backend crates.
//! They are documented for maintainers, but they are not yet a standalone
//! stable Rust API. Compatibility promises are made at the Python package
//! boundary unless a future release explicitly publishes Rust API stability.

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
pub use error::GenoioError;
pub use filter::{
    attach_variant_stats, compute_dosage_variant_stats, compute_variant_stats,
    is_dosage_polymorphic, variant_stats_from_counts, GenotypeFilterConjunction,
    GenotypeFilterPlan, PartialFilterDecision, RegionPredicate, VariantFilter, VariantStats,
    VariantWindow,
};
pub use metadata::{MetadataOutput, SampleRecord, VariantRecord};
pub use sparse::{
    append_sparse_column, flip_haplotype_values_to_minor_allele, flip_values_to_minor_allele,
    flip_variant_metadata_to_minor_allele, reject_sparse_missing, reject_sparse_missing_values,
    should_flip_haplotype_to_minor_allele, SparseGenotypeMatrix,
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
