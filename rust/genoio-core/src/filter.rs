// pattern: Functional Core

//! Validated variant-filter expressions and genotype-statistics helpers.
//!
//! The implementation is split by responsibility behind this facade so callers
//! can continue importing the established `genoio_core::filter` API.

mod expression;
mod ir;
mod model;
mod plan;
mod stats;

pub use expression::VariantFilter;
pub use model::{
    PartialFilterDecision, RegionPredicate, VariantMetadataView, VariantStats, VariantWindow,
};
pub use plan::{GenotypeFilterConjunction, GenotypeFilterPlan};
pub use stats::{
    attach_variant_stats, compute_dosage_variant_stats, is_dosage_polymorphic,
    variant_stats_from_counts,
};
