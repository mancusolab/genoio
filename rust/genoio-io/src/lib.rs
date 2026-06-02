// pattern: Imperative Shell

mod error;
mod plink1;
mod vcf;

pub use error::Result;
pub use plink1::{read_plink1_dense, read_plink1_metadata, read_plink1_sparse};
pub use vcf::{read_vcf_dense, read_vcf_metadata, read_vcf_sparse};

pub fn backend_name() -> &'static str {
    "genoio-io"
}
