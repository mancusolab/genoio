// pattern: Functional Core

pub mod capabilities;
pub mod dense;
pub mod error;
pub mod filter;
pub mod metadata;

pub use capabilities::SourceCapabilities;
pub use dense::{
    select_samples_source_order, transpose_variant_major_to_sample_major, DenseDiagnostics,
    DenseGenotypeMatrix, DenseSampleSelection,
};
pub use error::MetadataError;
pub use filter::{attach_variant_stats, compute_variant_stats, VariantFilter, VariantStats};
pub use metadata::{MetadataOutput, SampleRecord, VariantRecord};

pub const PACKAGE_NAME: &str = "genoio";
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const BUILD_PROFILE: &str = if cfg!(debug_assertions) {
    "debug"
} else {
    "release"
};
