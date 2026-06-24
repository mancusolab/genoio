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
//!
//! Naming grammar inside parser modules:
//!
//! - `read_*` performs source IO or advances a reader.
//! - `parse_*` turns text fields into typed records.
//! - `decode_*` turns encoded genotype or dosage payloads into values.
//! - `validate_*` checks source-data invariants.
//! - `ensure_*` checks caller preconditions before entering a backend.
//! - `skip_*` advances over data without materializing it.
//! - `Session` owns an open reader; `Context` is read-loop setup; `State` and
//!   `Buffers` are mutable decode scratch; `Parts` are final matrix assembly
//!   inputs.

mod bgen;
mod dosage_filter;
mod error;
mod hardcall;
mod matrix;
mod plink;
mod retention;
mod vcf;

pub use bgen::{
    read_bgen_dosage_dense, read_bgen_dosage_dense_windowed,
    read_bgen_dosage_dense_windowed_with_missing_policy,
    read_bgen_haplotypes_dosage_dense_windowed,
    read_bgen_haplotypes_dosage_dense_windowed_with_missing_policy, read_bgen_metadata,
};
pub use error::Result;
pub use plink::{
    read_plink1_dense, read_plink1_dense_windowed, read_plink1_dense_windowed_with_missing_policy,
    read_plink1_metadata, read_plink1_sparse, read_plink1_sparse_windowed, read_plink2_dense,
    read_plink2_dense_windowed, read_plink2_dense_windowed_with_missing_policy,
    read_plink2_dosage_dense_windowed, read_plink2_dosage_dense_windowed_with_missing_policy,
    read_plink2_haplotypes_dense_windowed,
    read_plink2_haplotypes_dense_windowed_with_missing_policy,
    read_plink2_haplotypes_dosage_dense_windowed,
    read_plink2_haplotypes_dosage_dense_windowed_with_missing_policy,
    read_plink2_haplotypes_sparse, read_plink2_haplotypes_sparse_windowed, read_plink2_metadata,
    read_plink2_sparse, read_plink2_sparse_windowed,
};
pub use vcf::{
    read_vcf_dense, read_vcf_dense_windowed, read_vcf_dense_windowed_with_missing_policy,
    read_vcf_dense_windowed_with_threads, read_vcf_dense_windowed_with_threads_and_missing_policy,
    read_vcf_dosage_dense_windowed, read_vcf_dosage_dense_windowed_with_missing_policy,
    read_vcf_dosage_dense_windowed_with_threads,
    read_vcf_dosage_dense_windowed_with_threads_and_missing_policy, read_vcf_haplotypes_dense,
    read_vcf_haplotypes_dense_windowed, read_vcf_haplotypes_dense_windowed_with_missing_policy,
    read_vcf_haplotypes_dense_windowed_with_threads,
    read_vcf_haplotypes_dense_windowed_with_threads_and_missing_policy, read_vcf_haplotypes_sparse,
    read_vcf_haplotypes_sparse_windowed, read_vcf_haplotypes_sparse_windowed_with_threads,
    read_vcf_metadata, read_vcf_sparse, read_vcf_sparse_windowed,
    read_vcf_sparse_windowed_with_threads,
};

/// Return the compiled Rust IO backend name for diagnostics.
pub fn backend_name() -> &'static str {
    "genoio-io"
}
