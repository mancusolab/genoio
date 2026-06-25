// pattern: Imperative Shell

//! PLINK-family readers.
//!
//! This facade keeps the crate-level API stable while the PLINK1 and PLINK2
//! implementations live under one format-family namespace.

mod common;
mod plink1;
mod plink2;

pub use plink1::{
    read_plink1_dense_windowed_with_arrow_variants, read_plink1_metadata_arrow,
    read_plink1_sparse_windowed_with_arrow_variants,
};
pub use plink2::{
    read_plink2_dense_windowed_with_arrow_variants,
    read_plink2_dosage_dense_windowed_with_arrow_variants,
    read_plink2_haplotypes_dense_windowed_with_arrow_variants,
    read_plink2_haplotypes_dosage_dense_windowed_with_arrow_variants,
    read_plink2_haplotypes_sparse_windowed_with_arrow_variants, read_plink2_metadata_arrow,
    read_plink2_sparse_windowed_with_arrow_variants,
};
