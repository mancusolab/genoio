// pattern: Functional Core

/// Representation support advertised by a resolved genotype source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceCapabilities {
    pub supports_geno: bool,
    pub supports_haplo: bool,
    pub phased: bool,
}

impl SourceCapabilities {
    /// Capabilities for unphased genotype-only formats.
    pub const fn genotype_only() -> Self {
        Self {
            supports_geno: true,
            supports_haplo: false,
            phased: false,
        }
    }

    /// Capabilities for VCF sources with phased diploid genotype evidence.
    pub const fn phased_genotypes() -> Self {
        Self {
            supports_geno: true,
            supports_haplo: true,
            phased: true,
        }
    }
}
