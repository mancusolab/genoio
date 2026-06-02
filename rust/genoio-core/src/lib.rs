// pattern: Functional Core

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
    attach_variant_stats, compute_variant_stats, RegionPredicate, VariantFilter, VariantStats,
};
pub use metadata::{MetadataOutput, SampleRecord, VariantRecord};
pub use sparse::{sparse_from_dense_minor_flipped, SparseGenotypeMatrix};

pub const PACKAGE_NAME: &str = "genoio";
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const BUILD_PROFILE: &str = if cfg!(debug_assertions) {
    "debug"
} else {
    "release"
};
