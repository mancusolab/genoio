// pattern: Functional Core

use crate::SourceCapabilities;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleRecord {
    pub fid: Option<String>,
    pub iid: String,
    pub father: Option<String>,
    pub mother: Option<String>,
    pub sex: Option<String>,
    pub phenotype: Option<String>,
    pub source_sample_index: Option<usize>,
    pub haplotype_index: Option<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct VariantRecord {
    pub chrom: String,
    pub pos: u32,
    pub id: String,
    pub a0: String,
    pub a1: String,
    pub ref_allele: Option<String>,
    pub alt_allele: Option<String>,
    pub source_a0: String,
    pub source_a1: String,
    pub flipped: bool,
    pub af: Option<f32>,
    pub maf: Option<f32>,
    pub mac: Option<u32>,
    pub missing_rate: Option<f32>,
    pub n_called: Option<u32>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetadataOutput {
    pub samples: Vec<SampleRecord>,
    pub variants: Vec<VariantRecord>,
    pub capabilities: SourceCapabilities,
}
