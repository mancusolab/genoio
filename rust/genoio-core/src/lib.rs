// pattern: Functional Core

pub mod capabilities;
pub mod error;
pub mod metadata;

pub use capabilities::SourceCapabilities;
pub use error::MetadataError;
pub use metadata::{MetadataOutput, SampleRecord, VariantRecord};

pub const PACKAGE_NAME: &str = "genoio";
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const BUILD_PROFILE: &str = if cfg!(debug_assertions) {
    "debug"
} else {
    "release"
};
