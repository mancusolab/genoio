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
}

#[derive(Debug, Clone, PartialEq, Eq)]
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
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MetadataOutput {
    pub samples: Vec<SampleRecord>,
    pub variants: Vec<VariantRecord>,
    pub capabilities: SourceCapabilities,
}
