// pattern: Functional Core

use super::model::VariantStats;

/// Compiled genotype-dependent portion of a variant filter.
///
/// Backends use this to select a format-specific predicate kernel when the
/// filter shape is simple enough to avoid constructing full `VariantStats`.
/// The plan intentionally ignores metadata predicates because readers have
/// already handled them with `VariantFilter::partial_decision`. A non-generic
/// plan is therefore only a valid complete decision after partial evaluation has
/// returned `NeedGenotypes` for that source variant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GenotypeFilterPlan {
    /// Use the generic `VariantStats` evaluation path.
    Generic,
    /// Retain variants with nonzero minor allele count.
    Polymorphic,
    /// Retain variants with minor allele count within a closed range.
    MacRange { min: Option<u32>, max: Option<u32> },
    /// Retain variants with minor allele frequency within a closed range.
    MafRange { min: Option<f32>, max: Option<f32> },
    /// Retain variants with missing rate no greater than `max`.
    MissingRateMax { max: f32 },
    /// Retain variants satisfying a conjunction of simple genotype predicates.
    Conjunction(GenotypeFilterConjunction),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenotypeFilterConjunction {
    /// True when a nonzero reference and nonzero alternate allele are required.
    pub polymorphic: bool,
    /// Inclusive lower bound for minor allele count.
    pub mac_min: Option<u32>,
    /// Inclusive upper bound for minor allele count.
    pub mac_max: Option<u32>,
    /// Inclusive lower bound for minor allele frequency.
    pub maf_min: Option<f32>,
    /// Inclusive upper bound for minor allele frequency.
    pub maf_max: Option<f32>,
    /// Inclusive upper bound for missing-call rate.
    pub missing_rate_max: Option<f32>,
}

impl GenotypeFilterPlan {
    pub fn evaluate_stats(self, stats: &VariantStats) -> Option<bool> {
        match self {
            Self::Generic => None,
            Self::Polymorphic => Some(stats.polymorphic),
            Self::MacRange { min, max } => Some(mac_in_range(stats, min, max)),
            Self::MafRange { min, max } => Some(maf_in_range(stats, min, max)),
            Self::MissingRateMax { max } => Some(stats.missing_rate <= f64::from(max)),
            Self::Conjunction(plan) => Some(plan.evaluate_stats(stats)),
        }
    }
}

impl GenotypeFilterConjunction {
    pub(super) fn empty() -> Self {
        Self {
            polymorphic: false,
            mac_min: None,
            mac_max: None,
            maf_min: None,
            maf_max: None,
            missing_rate_max: None,
        }
    }

    pub(super) fn is_empty(self) -> bool {
        !self.polymorphic && !self.has_mac() && !self.has_maf() && !self.has_missing_rate()
    }

    fn has_mac(self) -> bool {
        self.mac_min.is_some() || self.mac_max.is_some()
    }

    fn has_maf(self) -> bool {
        self.maf_min.is_some() || self.maf_max.is_some()
    }

    fn has_missing_rate(self) -> bool {
        self.missing_rate_max.is_some()
    }

    pub(super) fn into_plan(self) -> GenotypeFilterPlan {
        if self.polymorphic && !self.has_mac() && !self.has_maf() && !self.has_missing_rate() {
            return GenotypeFilterPlan::Polymorphic;
        }
        if !self.polymorphic && self.has_mac() && !self.has_maf() && !self.has_missing_rate() {
            return GenotypeFilterPlan::MacRange {
                min: self.mac_min,
                max: self.mac_max,
            };
        }
        if !self.polymorphic && !self.has_mac() && self.has_maf() && !self.has_missing_rate() {
            return GenotypeFilterPlan::MafRange {
                min: self.maf_min,
                max: self.maf_max,
            };
        }
        if !self.polymorphic && !self.has_mac() && !self.has_maf() {
            if let Some(max) = self.missing_rate_max {
                return GenotypeFilterPlan::MissingRateMax { max };
            }
        }
        GenotypeFilterPlan::Conjunction(self)
    }

    fn evaluate_stats(self, stats: &VariantStats) -> bool {
        if self.polymorphic && !stats.polymorphic {
            return false;
        }
        if self.has_mac() && !mac_in_range(stats, self.mac_min, self.mac_max) {
            return false;
        }
        if self.has_maf() && !maf_in_range(stats, self.maf_min, self.maf_max) {
            return false;
        }
        if self
            .missing_rate_max
            .is_some_and(|max| stats.missing_rate > f64::from(max))
        {
            return false;
        }
        true
    }
}

fn mac_in_range(stats: &VariantStats, min: Option<u32>, max: Option<u32>) -> bool {
    stats.mac.is_some_and(|mac| {
        min.is_none_or(|threshold| mac >= f64::from(threshold))
            && max.is_none_or(|threshold| mac <= f64::from(threshold))
    })
}

fn maf_in_range(stats: &VariantStats, min: Option<f32>, max: Option<f32>) -> bool {
    stats.maf.is_some_and(|maf| {
        min.is_none_or(|threshold| maf >= f64::from(threshold))
            && max.is_none_or(|threshold| maf <= f64::from(threshold))
    })
}
