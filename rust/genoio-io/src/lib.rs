// pattern: Imperative Shell

#![cfg_attr(
    not(test),
    deny(clippy::expect_used, clippy::panic, clippy::unwrap_used)
)]

//! File-format readers for `genoio-core` matrix and metadata contracts.
//!
//! Public functions in this crate perform filesystem IO and format parsing,
//! then return validated core structs. Python bindings call these functions
//! through `genoio-py`.

mod bgen;
mod error;
mod hardcall;
mod matrix;
mod plink1;
mod plink2;
mod plink_common;
mod retention;
mod vcf;

pub use bgen::{
    read_bgen_dosage_dense, read_bgen_dosage_dense_windowed,
    read_bgen_haplotypes_dosage_dense_windowed, read_bgen_metadata,
};
pub use error::Result;
pub use plink1::{
    read_plink1_dense, read_plink1_dense_windowed, read_plink1_metadata, read_plink1_sparse,
    read_plink1_sparse_windowed,
};
pub use plink2::{
    read_plink2_dense, read_plink2_dense_windowed, read_plink2_dosage_dense_windowed,
    read_plink2_haplotypes_dense_windowed, read_plink2_haplotypes_dosage_dense_windowed,
    read_plink2_haplotypes_sparse, read_plink2_haplotypes_sparse_windowed, read_plink2_metadata,
    read_plink2_sparse, read_plink2_sparse_windowed,
};
pub use vcf::{
    read_vcf_dense, read_vcf_dense_windowed, read_vcf_dense_windowed_with_threads,
    read_vcf_dosage_dense_windowed, read_vcf_dosage_dense_windowed_with_threads,
    read_vcf_haplotypes_dense, read_vcf_haplotypes_dense_windowed,
    read_vcf_haplotypes_dense_windowed_with_threads, read_vcf_haplotypes_sparse,
    read_vcf_haplotypes_sparse_windowed, read_vcf_haplotypes_sparse_windowed_with_threads,
    read_vcf_metadata, read_vcf_sparse, read_vcf_sparse_windowed,
    read_vcf_sparse_windowed_with_threads,
};

/// Return the compiled Rust IO backend name for diagnostics.
pub fn backend_name() -> &'static str {
    "genoio-io"
}
