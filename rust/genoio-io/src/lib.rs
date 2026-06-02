// pattern: Imperative Shell

mod error;
mod plink1;
mod vcf;

pub use error::Result;
pub use plink1::{
    read_plink1_dense, read_plink1_dense_windowed, read_plink1_metadata, read_plink1_sparse,
    read_plink1_sparse_windowed,
};
pub use vcf::{
    read_vcf_dense, read_vcf_dense_windowed, read_vcf_haplotypes_dense,
    read_vcf_haplotypes_dense_windowed, read_vcf_haplotypes_sparse,
    read_vcf_haplotypes_sparse_windowed, read_vcf_metadata, read_vcf_sparse,
    read_vcf_sparse_windowed,
};

pub fn backend_name() -> &'static str {
    "genoio-io"
}
