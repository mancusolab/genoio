// pattern: Functional Core

use crate::VariantRecord;

/// Metadata-only filter decision before genotype values are decoded.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartialFilterDecision {
    /// Metadata proves the variant passes the full expression.
    Accept,
    /// Metadata proves the variant fails the full expression.
    Reject,
    /// Genotype values and derived statistics are needed to decide.
    NeedGenotypes,
}

/// Concrete 1-based inclusive genomic region suitable for reader pushdown.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionPredicate {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

/// Per-variant statistics computed from called diploid genotype values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VariantStats {
    pub af: Option<f64>,
    pub maf: Option<f64>,
    pub mac: Option<f64>,
    pub missing_rate: f64,
    pub n_called: u32,
    pub polymorphic: bool,
}

/// Borrowed variant metadata contract for filter and validation hot paths.
///
/// Format readers can implement this view over parser-owned buffers, avoiding a
/// temporary [`VariantRecord`] when predicates only need borrowed fields. The
/// default methods describe metadata that many source formats do not provide or
/// only attach after genotype-stat filters retain a row.
pub trait VariantMetadataView {
    /// Source contig or chromosome label.
    fn chrom(&self) -> &str;
    /// 1-based source coordinate.
    fn pos(&self) -> u32;
    /// Public variant identifier after format-specific normalization.
    fn id(&self) -> &str;
    /// Public allele 0, potentially flipped for minor-allele sparse outputs.
    fn a0(&self) -> &str;
    /// Public allele 1, potentially flipped for minor-allele sparse outputs.
    fn a1(&self) -> &str;
    /// Original REF allele when the source format provides REF/ALT orientation.
    fn ref_allele(&self) -> Option<&str>;
    /// Original ALT allele string when available; comma-separated ALT marks multiallelic records.
    fn alt_allele(&self) -> Option<&str>;

    /// Source allele 0 before public allele flipping.
    fn source_a0(&self) -> &str {
        self.a0()
    }

    /// Source allele 1 before public allele flipping.
    fn source_a1(&self) -> &str {
        self.a1()
    }

    /// True when public `a0`/`a1` have been swapped relative to source alleles.
    fn flipped(&self) -> bool {
        false
    }

    /// Source quality score when the metadata format exposes one.
    fn qual(&self) -> Option<f32> {
        None
    }

    /// Attached allele frequency for retained genotype-stat-filtered variants.
    fn af(&self) -> Option<f32> {
        None
    }

    /// Attached minor allele frequency for retained genotype-stat-filtered variants.
    fn maf(&self) -> Option<f32> {
        None
    }

    /// Attached integer minor allele count for retained hard-call-compatible variants.
    fn mac(&self) -> Option<u32> {
        None
    }

    /// Attached missing-call rate for retained genotype-stat-filtered variants.
    fn missing_rate(&self) -> Option<f32> {
        None
    }

    /// Attached called genotype count for retained genotype-stat-filtered variants.
    fn n_called(&self) -> Option<u32> {
        None
    }
}

impl VariantMetadataView for VariantRecord {
    fn chrom(&self) -> &str {
        &self.chrom
    }

    fn pos(&self) -> u32 {
        self.pos
    }

    fn id(&self) -> &str {
        &self.id
    }

    fn a0(&self) -> &str {
        &self.a0
    }

    fn a1(&self) -> &str {
        &self.a1
    }

    fn ref_allele(&self) -> Option<&str> {
        self.ref_allele.as_deref()
    }

    fn alt_allele(&self) -> Option<&str> {
        self.alt_allele.as_deref()
    }

    fn source_a0(&self) -> &str {
        &self.source_a0
    }

    fn source_a1(&self) -> &str {
        &self.source_a1
    }

    fn flipped(&self) -> bool {
        self.flipped
    }

    fn qual(&self) -> Option<f32> {
        self.qual
    }

    fn af(&self) -> Option<f32> {
        self.af
    }

    fn maf(&self) -> Option<f32> {
        self.maf
    }

    fn mac(&self) -> Option<u32> {
        self.mac
    }

    fn missing_rate(&self) -> Option<f32> {
        self.missing_rate
    }

    fn n_called(&self) -> Option<u32> {
        self.n_called
    }
}

/// Retained-variant window for block reads.
///
/// `start` and `len` are expressed after filters have retained variants, not
/// necessarily in raw source-row coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantWindow {
    pub start: usize,
    pub len: usize,
}

impl VariantWindow {
    /// Return true when `retained_index` belongs to this window.
    pub fn contains(self, retained_index: usize) -> bool {
        retained_index >= self.start && retained_index < self.start.saturating_add(self.len)
    }

    /// Return true when no later retained variant can belong to this window.
    pub fn is_past(self, retained_index: usize) -> bool {
        retained_index >= self.start.saturating_add(self.len)
    }
}
