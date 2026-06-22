// pattern: Functional Core
//! Crate-local result alias for fallible IO and decode paths.
//!
//! `genoio-io` uses `genoio_core::GenoioError` directly so format readers share
//! one error surface with matrix validation and filter evaluation.

pub type Result<T> = std::result::Result<T, genoio_core::GenoioError>;
