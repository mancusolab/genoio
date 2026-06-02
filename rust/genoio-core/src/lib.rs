// pattern: Functional Core

pub const PACKAGE_NAME: &str = "genoio";
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const BUILD_PROFILE: &str = if cfg!(debug_assertions) {
    "debug"
} else {
    "release"
};

pub mod source {}
pub mod contracts {}
