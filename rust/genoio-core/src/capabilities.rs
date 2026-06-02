// pattern: Functional Core

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceCapabilities {
    pub supports_geno: bool,
    pub supports_haplo: bool,
    pub phased: bool,
}

impl SourceCapabilities {
    pub const fn genotype_only() -> Self {
        Self {
            supports_geno: true,
            supports_haplo: false,
            phased: false,
        }
    }

    pub const fn phased_genotypes() -> Self {
        Self {
            supports_geno: true,
            supports_haplo: true,
            phased: true,
        }
    }
}
