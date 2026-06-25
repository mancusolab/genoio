// pattern: Functional Core

use std::collections::BTreeSet;

use serde::Deserialize;
use serde_json::Value;

use crate::{GenoioError, VariantRecord};

/// Serializable variant filter evaluated by Rust readers.
///
/// Python constructs a JSON-compatible filter IR. This type validates that IR,
/// tracks which predicates need genotype statistics, and exposes safe metadata
/// decisions for early record skipping.
///
/// Reader contract:
///
/// 1. Call [`VariantFilter::partial_decision`] after source metadata has been
///    parsed and before genotype values are decoded.
/// 2. Treat `Accept` and `Reject` as final decisions; they must not decode
///    genotypes just to re-evaluate the same filter.
/// 3. For `NeedGenotypes`, evaluate the genotype-dependent portion with the
///    most native backend representation available. [`GenotypeFilterPlan`]
///    describes the optimized shapes; `Generic` means the reader must construct
///    complete [`VariantStats`] and call [`VariantFilter::evaluate`].
///
/// Matrix-only reads may return only the retain/drop decision. Reads that return
/// variant metadata should attach complete stats for retained genotype-filtered
/// variants when the public output contract exposes those fields.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantFilter {
    expr: Expr,
}

/// Compiled genotype-dependent portion of a variant filter.
///
/// Backends use this to select a format-specific predicate kernel when the
/// filter shape is simple enough to avoid constructing full `VariantStats`.
/// The plan intentionally ignores metadata predicates because readers have
/// already handled them with [`VariantFilter::partial_decision`]. A non-generic
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
    fn empty() -> Self {
        Self {
            polymorphic: false,
            mac_min: None,
            mac_max: None,
            maf_min: None,
            maf_max: None,
            missing_rate_max: None,
        }
    }

    fn is_empty(self) -> bool {
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

    fn into_plan(self) -> GenotypeFilterPlan {
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

#[derive(Debug, Clone, PartialEq)]
enum Expr {
    AlwaysTrue,
    AlwaysFalse,
    Predicate(Predicate),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Not(Box<Expr>),
}

#[derive(Debug, Clone, PartialEq)]
enum Predicate {
    Chrom(String),
    Region { chrom: String, start: u32, end: u32 },
    IdIn(BTreeSet<String>),
    Snp,
    Biallelic,
    Qual { min: Option<f32>, max: Option<f32> },
    Maf { min: Option<f32>, max: Option<f32> },
    Mac { min: Option<u32>, max: Option<u32> },
    MissingRate { max: f32 },
    Polymorphic,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum RawExpr {
    Predicate {
        name: String,
        #[serde(default)]
        params: Value,
    },
    And {
        left: Box<RawExpr>,
        right: Box<RawExpr>,
    },
    Or {
        left: Box<RawExpr>,
        right: Box<RawExpr>,
    },
    Not {
        expr: Box<RawExpr>,
    },
}

impl VariantFilter {
    /// Parse, validate, and simplify the JSON-compatible filter IR from Python.
    ///
    /// Simplification is semantic-preserving: contradictory expressions become
    /// `AlwaysFalse`, repeated thresholds are tightened, and concrete region
    /// intersections are reduced before any reader sees the filter.
    pub fn from_json_value(value: Value) -> Result<Self, GenoioError> {
        let raw: RawExpr = serde_json::from_value(value).map_err(|error| {
            GenoioError::invalid_source("<filter>", format!("invalid filter IR: {error}"))
        })?;
        Ok(Self {
            expr: Expr::from_raw(raw)?.simplify(),
        })
    }

    /// Evaluate predicates that can be decided from variant metadata alone.
    ///
    /// Returns `None` when any part of the expression needs genotype-derived
    /// statistics such as MAF or missing rate.
    pub fn metadata_decision(&self, variant: &VariantRecord) -> Option<bool> {
        self.metadata_decision_view(variant)
    }

    /// Evaluate metadata-only predicates against borrowed variant fields.
    ///
    /// This is the hot-path entry point for parser-owned metadata. It keeps the
    /// old `VariantRecord` wrapper available while letting readers avoid owned
    /// strings when metadata alone can decide the filter.
    pub fn metadata_decision_view<V: VariantMetadataView + ?Sized>(
        &self,
        variant: &V,
    ) -> Option<bool> {
        self.expr.metadata_decision(variant)
    }

    /// Partially evaluate the filter using metadata available before GT decode.
    ///
    /// This preserves boolean semantics while letting readers avoid genotype
    /// statistics whenever metadata alone proves an accept or reject decision.
    pub fn partial_decision(&self, variant: &VariantRecord) -> PartialFilterDecision {
        self.partial_decision_view(variant)
    }

    /// Partially evaluate the filter using borrowed metadata fields.
    ///
    /// `NeedGenotypes` is the only result that should force a reader to decode
    /// genotypes solely for filter evaluation.
    pub fn partial_decision_view<V: VariantMetadataView + ?Sized>(
        &self,
        variant: &V,
    ) -> PartialFilterDecision {
        self.expr.partial_decision(variant)
    }

    /// Evaluate the complete filter against metadata and optional statistics.
    pub fn evaluate(&self, variant: &VariantRecord, stats: Option<&VariantStats>) -> bool {
        self.evaluate_view(variant, stats)
    }

    /// Evaluate the complete filter against borrowed metadata and optional statistics.
    ///
    /// Genotype-dependent predicates evaluate from `stats`; metadata predicates
    /// continue to read borrowed fields from `variant`.
    pub fn evaluate_view<V: VariantMetadataView + ?Sized>(
        &self,
        variant: &V,
        stats: Option<&VariantStats>,
    ) -> bool {
        self.expr.evaluate(variant, stats)
    }

    /// Return true when any predicate needs genotype statistics.
    pub fn requires_genotype_stats(&self) -> bool {
        self.expr.requires_genotype_stats()
    }

    /// Return true when the expression has no variant-metadata predicates.
    pub fn is_genotype_stats_only(&self) -> bool {
        self.expr.is_genotype_stats_only()
    }

    /// Evaluate a metadata-free expression from genotype statistics alone.
    ///
    /// Returns `None` when the expression contains any predicate requiring
    /// `VariantRecord` fields such as chromosome, position, ID, alleles, or
    /// quality.
    pub fn evaluate_genotype_stats(&self, stats: &VariantStats) -> Option<bool> {
        self.is_genotype_stats_only()
            .then(|| self.expr.evaluate_genotype_stats(stats))
    }

    /// Return a compiled plan for genotype-dependent filter evaluation.
    pub fn genotype_filter_plan(&self) -> GenotypeFilterPlan {
        self.expr.genotype_filter_plan()
    }

    /// Return true when the expression contains any region predicate.
    pub fn has_region_predicate(&self) -> bool {
        self.expr.has_region_predicate()
    }

    /// Return a region that can be safely pushed into an indexed reader.
    ///
    /// Only bare region predicates and conjunctions with a region are safe.
    /// Disjunctions and negation require full expression evaluation.
    pub fn concrete_region_pushdown(&self) -> Option<RegionPredicate> {
        self.expr.concrete_region_pushdown()
    }

    /// Return true when simplification proves the filter cannot retain variants.
    ///
    /// Readers use this as an early-empty fast path. It is not an error state:
    /// callers may intentionally ask for a chromosome, region, or ID set that
    /// has no overlap with a source.
    pub fn is_always_false(&self) -> bool {
        matches!(self.expr, Expr::AlwaysFalse)
    }
}

impl Expr {
    fn from_raw(raw: RawExpr) -> Result<Self, GenoioError> {
        match raw {
            RawExpr::Predicate { name, params } => {
                Ok(Self::Predicate(Predicate::from_raw(&name, params)?))
            }
            RawExpr::And { left, right } => Ok(Self::And(
                Box::new(Self::from_raw(*left)?),
                Box::new(Self::from_raw(*right)?),
            )),
            RawExpr::Or { left, right } => Ok(Self::Or(
                Box::new(Self::from_raw(*left)?),
                Box::new(Self::from_raw(*right)?),
            )),
            RawExpr::Not { expr } => Ok(Self::Not(Box::new(Self::from_raw(*expr)?))),
        }
    }

    fn metadata_decision<V: VariantMetadataView + ?Sized>(&self, variant: &V) -> Option<bool> {
        match self {
            Self::AlwaysTrue => Some(true),
            Self::AlwaysFalse => Some(false),
            Self::Predicate(predicate) => predicate.metadata_decision(variant),
            Self::And(left, right) => match left.metadata_decision(variant) {
                Some(false) => Some(false),
                Some(true) => right.metadata_decision(variant),
                None => match right.metadata_decision(variant) {
                    Some(false) => Some(false),
                    _ => None,
                },
            },
            Self::Or(left, right) => match left.metadata_decision(variant) {
                Some(true) => Some(true),
                Some(false) => right.metadata_decision(variant),
                None => match right.metadata_decision(variant) {
                    Some(true) => Some(true),
                    _ => None,
                },
            },
            Self::Not(expr) => expr.metadata_decision(variant).map(|decision| !decision),
        }
    }

    fn partial_decision<V: VariantMetadataView + ?Sized>(
        &self,
        variant: &V,
    ) -> PartialFilterDecision {
        match self.metadata_decision(variant) {
            Some(true) => PartialFilterDecision::Accept,
            Some(false) => PartialFilterDecision::Reject,
            None => PartialFilterDecision::NeedGenotypes,
        }
    }

    fn evaluate<V: VariantMetadataView + ?Sized>(
        &self,
        variant: &V,
        stats: Option<&VariantStats>,
    ) -> bool {
        match self {
            Self::AlwaysTrue => true,
            Self::AlwaysFalse => false,
            Self::Predicate(predicate) => predicate.evaluate(variant, stats),
            Self::And(left, right) => {
                left.evaluate(variant, stats) && right.evaluate(variant, stats)
            }
            Self::Or(left, right) => {
                left.evaluate(variant, stats) || right.evaluate(variant, stats)
            }
            Self::Not(expr) => !expr.evaluate(variant, stats),
        }
    }

    fn requires_genotype_stats(&self) -> bool {
        match self {
            Self::AlwaysTrue | Self::AlwaysFalse => false,
            Self::Predicate(predicate) => predicate.requires_genotype_stats(),
            Self::And(left, right) | Self::Or(left, right) => {
                left.requires_genotype_stats() || right.requires_genotype_stats()
            }
            Self::Not(expr) => expr.requires_genotype_stats(),
        }
    }

    fn is_genotype_stats_only(&self) -> bool {
        match self {
            Self::AlwaysTrue | Self::AlwaysFalse => true,
            Self::Predicate(predicate) => predicate.is_genotype_stats_only(),
            Self::And(left, right) | Self::Or(left, right) => {
                left.is_genotype_stats_only() && right.is_genotype_stats_only()
            }
            Self::Not(expr) => expr.is_genotype_stats_only(),
        }
    }

    fn evaluate_genotype_stats(&self, stats: &VariantStats) -> bool {
        match self {
            Self::AlwaysTrue => true,
            Self::AlwaysFalse => false,
            Self::Predicate(predicate) => predicate.evaluate_genotype_stats(stats),
            Self::And(left, right) => {
                left.evaluate_genotype_stats(stats) && right.evaluate_genotype_stats(stats)
            }
            Self::Or(left, right) => {
                left.evaluate_genotype_stats(stats) || right.evaluate_genotype_stats(stats)
            }
            Self::Not(expr) => !expr.evaluate_genotype_stats(stats),
        }
    }

    fn genotype_filter_plan(&self) -> GenotypeFilterPlan {
        let mut conjunction = GenotypeFilterConjunction::empty();
        if self.collect_genotype_conjunction(&mut conjunction) && !conjunction.is_empty() {
            conjunction.into_plan()
        } else {
            GenotypeFilterPlan::Generic
        }
    }

    fn collect_genotype_conjunction(&self, plan: &mut GenotypeFilterConjunction) -> bool {
        match self {
            Self::AlwaysTrue | Self::AlwaysFalse => true,
            Self::Predicate(predicate) => {
                predicate.collect_genotype_conjunction(plan);
                true
            }
            Self::And(left, right) => {
                left.collect_genotype_conjunction(plan) && right.collect_genotype_conjunction(plan)
            }
            Self::Or(_, _) | Self::Not(_) => false,
        }
    }

    fn has_region_predicate(&self) -> bool {
        match self {
            Self::AlwaysTrue | Self::AlwaysFalse => false,
            Self::Predicate(Predicate::Region { .. }) => true,
            Self::Predicate(_) => false,
            Self::And(left, right) | Self::Or(left, right) => {
                left.has_region_predicate() || right.has_region_predicate()
            }
            Self::Not(expr) => expr.has_region_predicate(),
        }
    }

    fn concrete_region_pushdown(&self) -> Option<RegionPredicate> {
        match self {
            Self::AlwaysTrue | Self::AlwaysFalse => None,
            Self::Predicate(Predicate::Region { chrom, start, end }) => Some(RegionPredicate {
                chrom: chrom.clone(),
                start: *start,
                end: *end,
            }),
            Self::And(left, right) => left
                .concrete_region_pushdown()
                .or_else(|| right.concrete_region_pushdown()),
            Self::Predicate(_) | Self::Or(_, _) | Self::Not(_) => None,
        }
    }

    fn simplify(self) -> Self {
        match self {
            Self::And(left, right) => simplify_and(*left, *right),
            Self::Or(left, right) => simplify_or(*left, *right),
            Self::Not(expr) => match expr.simplify() {
                Self::AlwaysTrue => Self::AlwaysFalse,
                Self::AlwaysFalse => Self::AlwaysTrue,
                Self::Not(inner) => *inner,
                simplified => Self::Not(Box::new(simplified)),
            },
            other => other,
        }
    }
}

fn simplify_and(left: Expr, right: Expr) -> Expr {
    let mut terms = Vec::new();
    flatten_and(left.simplify(), &mut terms);
    flatten_and(right.simplify(), &mut terms);

    let mut simplified = Vec::<Expr>::new();
    for term in terms {
        match term {
            Expr::AlwaysFalse => return Expr::AlwaysFalse,
            Expr::AlwaysTrue => {}
            Expr::Predicate(predicate) => {
                if !combine_predicate_term(&mut simplified, predicate, Predicate::and_combine) {
                    return Expr::AlwaysFalse;
                }
            }
            other => simplified.push(other),
        }
    }
    simplified.sort_by_key(and_term_cost);
    rebuild_conjunction(simplified)
}

fn simplify_or(left: Expr, right: Expr) -> Expr {
    let mut terms = Vec::new();
    flatten_or(left.simplify(), &mut terms);
    flatten_or(right.simplify(), &mut terms);

    let mut simplified = Vec::<Expr>::new();
    for term in terms {
        match term {
            Expr::AlwaysTrue => return Expr::AlwaysTrue,
            Expr::AlwaysFalse => {}
            Expr::Predicate(predicate) => {
                combine_predicate_term(&mut simplified, predicate, |existing, current| {
                    existing
                        .or_combine(current)
                        .map_or(PredicateCombine::Unchanged, PredicateCombine::Combined)
                });
            }
            other => simplified.push(other),
        }
    }
    rebuild_disjunction(simplified)
}

fn combine_predicate_term(
    terms: &mut Vec<Expr>,
    predicate: Predicate,
    combine: impl Fn(&Predicate, &Predicate) -> PredicateCombine,
) -> bool {
    // Combining two predicates can expose a new simplification with an earlier
    // term, so restart after each successful merge until the predicate settles.
    let mut pending = Some(predicate);
    let mut index = 0;
    while let Some(current) = pending.take() {
        if index >= terms.len() {
            pending = Some(current);
            break;
        }
        if let Expr::Predicate(existing) = &terms[index] {
            match combine(existing, &current) {
                PredicateCombine::Unchanged => {
                    pending = Some(current);
                    index += 1;
                }
                PredicateCombine::Combined(combined) => {
                    terms.remove(index);
                    pending = Some(combined);
                    index = 0;
                }
                PredicateCombine::AlwaysFalse => return false,
            }
        } else {
            pending = Some(current);
            index += 1;
        }
    }
    if let Some(predicate) = pending {
        terms.push(Expr::Predicate(predicate));
    }
    true
}

fn flatten_and(expr: Expr, out: &mut Vec<Expr>) {
    match expr {
        Expr::And(left, right) => {
            flatten_and(*left, out);
            flatten_and(*right, out);
        }
        other => out.push(other),
    }
}

fn flatten_or(expr: Expr, out: &mut Vec<Expr>) {
    match expr {
        Expr::Or(left, right) => {
            flatten_or(*left, out);
            flatten_or(*right, out);
        }
        other => out.push(other),
    }
}

fn rebuild_conjunction(mut terms: Vec<Expr>) -> Expr {
    if terms.is_empty() {
        return Expr::AlwaysTrue;
    }
    let mut expr = terms.remove(0);
    for term in terms {
        expr = Expr::And(Box::new(expr), Box::new(term));
    }
    expr
}

fn rebuild_disjunction(mut terms: Vec<Expr>) -> Expr {
    if terms.is_empty() {
        return Expr::AlwaysFalse;
    }
    let mut expr = terms.remove(0);
    for term in terms {
        expr = Expr::Or(Box::new(expr), Box::new(term));
    }
    expr
}

fn and_term_cost(expr: &Expr) -> u8 {
    match expr {
        Expr::AlwaysTrue | Expr::AlwaysFalse => 0,
        Expr::Predicate(Predicate::Chrom(_)) | Expr::Predicate(Predicate::Region { .. }) => 1,
        Expr::Predicate(Predicate::IdIn(_)) => 2,
        Expr::Predicate(Predicate::Qual { .. })
        | Expr::Predicate(Predicate::Snp)
        | Expr::Predicate(Predicate::Biallelic) => 3,
        Expr::Predicate(predicate) if predicate.requires_genotype_stats() => 10,
        Expr::Predicate(_) => 5,
        Expr::Not(expr) if !expr.requires_genotype_stats() => 6,
        Expr::And(_, _) | Expr::Or(_, _) | Expr::Not(_) => 8,
    }
}

enum PredicateCombine {
    Unchanged,
    Combined(Predicate),
    AlwaysFalse,
}

impl Predicate {
    fn from_raw(name: &str, params: Value) -> Result<Self, GenoioError> {
        match name {
            "chrom" => Ok(Self::Chrom(required_string(&params, "value")?)),
            "region" => {
                let region = required_string(&params, "value")?;
                let (chrom, start, end) = parse_region(&region)?;
                Ok(Self::Region { chrom, start, end })
            }
            "id_in" => Ok(Self::IdIn(required_string_set(&params, "values")?)),
            "snp" => {
                expect_no_params(&params)?;
                Ok(Self::Snp)
            }
            "biallelic" => {
                expect_no_params(&params)?;
                Ok(Self::Biallelic)
            }
            "qual" => {
                let min = optional_nonnegative_f32(&params, "min")?;
                let max = optional_nonnegative_f32(&params, "max")?;
                validate_range("qual", min, max)?;
                Ok(Self::Qual { min, max })
            }
            "maf" => {
                let min = optional_rate(&params, "min")?;
                let max = optional_rate(&params, "max")?;
                validate_range("maf", min, max)?;
                Ok(Self::Maf { min, max })
            }
            "mac" => {
                let min = optional_u32(&params, "min")?;
                let max = optional_u32(&params, "max")?;
                validate_range("mac", min, max)?;
                Ok(Self::Mac { min, max })
            }
            "missing_rate" => Ok(Self::MissingRate {
                max: required_rate(&params, "max")?,
            }),
            "polymorphic" => {
                expect_no_params(&params)?;
                Ok(Self::Polymorphic)
            }
            other => Err(GenoioError::invalid_source(
                "<filter>",
                format!("unknown predicate name: {other}"),
            )),
        }
    }

    fn metadata_decision<V: VariantMetadataView + ?Sized>(&self, variant: &V) -> Option<bool> {
        match self {
            Self::Chrom(value) => Some(variant.chrom() == value),
            Self::Region { chrom, start, end } => {
                Some(variant.chrom() == chrom && variant.pos() >= *start && variant.pos() <= *end)
            }
            Self::IdIn(values) => Some(values.contains(variant.id())),
            Self::Snp => Some(is_snp(variant)),
            Self::Biallelic => Some(is_biallelic(variant)),
            Self::Qual { min, max } => Some(variant.qual().is_some_and(|qual| {
                min.is_none_or(|threshold| qual >= threshold)
                    && max.is_none_or(|threshold| qual <= threshold)
            })),
            Self::Maf { .. } | Self::Mac { .. } | Self::MissingRate { .. } | Self::Polymorphic => {
                None
            }
        }
    }

    fn evaluate<V: VariantMetadataView + ?Sized>(
        &self,
        variant: &V,
        stats: Option<&VariantStats>,
    ) -> bool {
        match self {
            Self::Chrom(_)
            | Self::Region { .. }
            | Self::IdIn(_)
            | Self::Snp
            | Self::Biallelic
            | Self::Qual { .. } => self.metadata_decision(variant) == Some(true),
            Self::Maf { .. } | Self::Mac { .. } | Self::MissingRate { .. } | Self::Polymorphic => {
                stats.and_then(|stats| self.evaluate_from_stats(stats)) == Some(true)
            }
        }
    }

    fn requires_genotype_stats(&self) -> bool {
        self.genotype_filter_plan().is_some()
    }

    fn is_genotype_stats_only(&self) -> bool {
        self.requires_genotype_stats()
    }

    fn evaluate_genotype_stats(&self, stats: &VariantStats) -> bool {
        self.evaluate_from_stats(stats).unwrap_or(false)
    }

    fn evaluate_from_stats(&self, stats: &VariantStats) -> Option<bool> {
        self.genotype_filter_plan()
            .and_then(|plan| plan.evaluate_stats(stats))
    }

    fn genotype_filter_plan(&self) -> Option<GenotypeFilterPlan> {
        match self {
            Self::Maf { min, max } => Some(GenotypeFilterPlan::MafRange {
                min: *min,
                max: *max,
            }),
            Self::Mac { min, max } => Some(GenotypeFilterPlan::MacRange {
                min: *min,
                max: *max,
            }),
            Self::MissingRate { max } => Some(GenotypeFilterPlan::MissingRateMax { max: *max }),
            Self::Polymorphic => Some(GenotypeFilterPlan::Polymorphic),
            Self::Chrom(_)
            | Self::Region { .. }
            | Self::IdIn(_)
            | Self::Snp
            | Self::Biallelic
            | Self::Qual { .. } => None,
        }
    }

    fn collect_genotype_conjunction(&self, plan: &mut GenotypeFilterConjunction) {
        match self {
            Self::Polymorphic => plan.polymorphic = true,
            Self::Mac { min, max } => {
                plan.mac_min = max_option(plan.mac_min, *min);
                plan.mac_max = min_option(plan.mac_max, *max);
            }
            Self::Maf { min, max } => {
                plan.maf_min = max_f32_option(plan.maf_min, *min);
                plan.maf_max = min_f32_option(plan.maf_max, *max);
            }
            Self::MissingRate { max } => {
                plan.missing_rate_max = min_f32_option(plan.missing_rate_max, Some(*max));
            }
            Self::Chrom(_)
            | Self::Region { .. }
            | Self::IdIn(_)
            | Self::Snp
            | Self::Biallelic
            | Self::Qual { .. } => {}
        }
    }

    fn and_combine(&self, other: &Self) -> PredicateCombine {
        use PredicateCombine::{AlwaysFalse, Combined, Unchanged};
        match (self, other) {
            (left, right) if left == right => Combined(left.clone()),
            (Self::Chrom(left), Self::Chrom(right)) => {
                if left == right {
                    Combined(self.clone())
                } else {
                    AlwaysFalse
                }
            }
            (
                Self::Region {
                    chrom: left_chrom,
                    start: left_start,
                    end: left_end,
                },
                Self::Region {
                    chrom: right_chrom,
                    start: right_start,
                    end: right_end,
                },
            ) => {
                if left_chrom != right_chrom {
                    return AlwaysFalse;
                }
                let start = (*left_start).max(*right_start);
                let end = (*left_end).min(*right_end);
                if start > end {
                    AlwaysFalse
                } else {
                    Combined(Self::Region {
                        chrom: left_chrom.clone(),
                        start,
                        end,
                    })
                }
            }
            (
                Self::Chrom(chrom),
                Self::Region {
                    chrom: region_chrom,
                    ..
                },
            )
            | (
                Self::Region {
                    chrom: region_chrom,
                    ..
                },
                Self::Chrom(chrom),
            ) => {
                if chrom == region_chrom {
                    Combined(match (self, other) {
                        (Self::Region { .. }, _) => self.clone(),
                        (_, Self::Region { .. }) => other.clone(),
                        _ => unreachable!("matched chrom-region pair"),
                    })
                } else {
                    AlwaysFalse
                }
            }
            (Self::IdIn(left), Self::IdIn(right)) => {
                let values = left.intersection(right).cloned().collect::<BTreeSet<_>>();
                if values.is_empty() {
                    AlwaysFalse
                } else {
                    Combined(Self::IdIn(values))
                }
            }
            (
                Self::Qual {
                    min: left_min,
                    max: left_max,
                },
                Self::Qual {
                    min: right_min,
                    max: right_max,
                },
            ) => combine_f32_range(*left_min, *left_max, *right_min, *right_max)
                .map_or(AlwaysFalse, |(min, max)| Combined(Self::Qual { min, max })),
            (
                Self::Maf {
                    min: left_min,
                    max: left_max,
                },
                Self::Maf {
                    min: right_min,
                    max: right_max,
                },
            ) => combine_f32_range(*left_min, *left_max, *right_min, *right_max)
                .map_or(AlwaysFalse, |(min, max)| Combined(Self::Maf { min, max })),
            (
                Self::Mac {
                    min: left_min,
                    max: left_max,
                },
                Self::Mac {
                    min: right_min,
                    max: right_max,
                },
            ) => combine_u32_range(*left_min, *left_max, *right_min, *right_max)
                .map_or(AlwaysFalse, |(min, max)| Combined(Self::Mac { min, max })),
            (Self::MissingRate { max: left }, Self::MissingRate { max: right }) => {
                Combined(Self::MissingRate {
                    max: (*left).min(*right),
                })
            }
            (Self::Snp, Self::Biallelic) | (Self::Biallelic, Self::Snp) => Combined(Self::Snp),
            _ => Unchanged,
        }
    }

    fn or_combine(&self, other: &Self) -> Option<Predicate> {
        match (self, other) {
            (left, right) if left == right => Some(left.clone()),
            (Self::IdIn(left), Self::IdIn(right)) => {
                let mut values = left.clone();
                values.extend(right.iter().cloned());
                Some(Self::IdIn(values))
            }
            _ => None,
        }
    }
}

fn combine_f32_range(
    left_min: Option<f32>,
    left_max: Option<f32>,
    right_min: Option<f32>,
    right_max: Option<f32>,
) -> Option<(Option<f32>, Option<f32>)> {
    let min = max_f32_option(left_min, right_min);
    let max = min_f32_option(left_max, right_max);
    if min.zip(max).is_some_and(|(min, max)| min > max) {
        None
    } else {
        Some((min, max))
    }
}

fn max_f32_option(left: Option<f32>, right: Option<f32>) -> Option<f32> {
    match (left, right) {
        (Some(left), Some(right)) => Some(left.max(right)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn min_f32_option(left: Option<f32>, right: Option<f32>) -> Option<f32> {
    match (left, right) {
        (Some(left), Some(right)) => Some(left.min(right)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn combine_u32_range(
    left_min: Option<u32>,
    left_max: Option<u32>,
    right_min: Option<u32>,
    right_max: Option<u32>,
) -> Option<(Option<u32>, Option<u32>)> {
    let min = max_option(left_min, right_min);
    let max = min_option(left_max, right_max);
    if min.zip(max).is_some_and(|(min, max)| min > max) {
        None
    } else {
        Some((min, max))
    }
}

fn max_option<T: Ord + Copy>(left: Option<T>, right: Option<T>) -> Option<T> {
    match (left, right) {
        (Some(left), Some(right)) => Some(left.max(right)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

fn min_option<T: Ord + Copy>(left: Option<T>, right: Option<T>) -> Option<T> {
    match (left, right) {
        (Some(left), Some(right)) => Some(left.min(right)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

/// Compute frequency and missingness statistics from sparse missing indices.
///
/// `missing_indices` must be sorted and unique. The corresponding `values`
/// entries are ignored, allowing decoders to use any placeholder for missing
/// calls without materializing a full boolean mask.
pub fn compute_dosage_variant_stats(
    values: &[f32],
    missing_indices: &[usize],
) -> Result<VariantStats, GenoioError> {
    validate_missing_indices(values.len(), missing_indices)?;

    let mut allele_count = 0.0_f64;
    let mut called_count = 0_u64;
    let mut missing_cursor = 0_usize;
    for (index, value) in values.iter().enumerate() {
        if missing_indices
            .get(missing_cursor)
            .is_some_and(|&missing_index| missing_index == index)
        {
            missing_cursor += 1;
            continue;
        }
        if !(0.0..=2.0).contains(value) {
            return Err(GenoioError::invalid_source(
                "<filter>",
                format!("dosage statistics require values in [0, 2]; observed {value}"),
            ));
        }
        allele_count += f64::from(*value);
        called_count += 1;
    }

    let missing_count = u64::try_from(missing_indices.len()).map_err(|_| {
        GenoioError::invalid_source("<filter>", "missing genotype count is out of range")
    })?;
    variant_stats_from_dosage_count(allele_count, called_count, missing_count)
}

/// Return true when called dosage values contain both alleles.
pub fn is_dosage_polymorphic(
    values: &[f32],
    missing_indices: &[usize],
) -> Result<bool, GenoioError> {
    validate_missing_indices(values.len(), missing_indices)?;

    let mut allele_count = 0.0_f64;
    let mut called_count = 0_u64;
    let mut missing_cursor = 0_usize;
    for (index, value) in values.iter().enumerate() {
        if missing_indices
            .get(missing_cursor)
            .is_some_and(|&missing_index| missing_index == index)
        {
            missing_cursor += 1;
            continue;
        }
        if !(0.0..=2.0).contains(value) {
            return Err(GenoioError::invalid_source(
                "<filter>",
                format!("dosage statistics require values in [0, 2]; observed {value}"),
            ));
        }
        allele_count += f64::from(*value);
        called_count += 1;
        if allele_count > 0.0 && allele_count < 2.0 * called_count as f64 {
            return Ok(true);
        }
    }
    Ok(false)
}

fn validate_missing_indices(
    values_len: usize,
    missing_indices: &[usize],
) -> Result<(), GenoioError> {
    let mut previous = None;
    for &index in missing_indices {
        if index >= values_len {
            return Err(GenoioError::invalid_source(
                "<filter>",
                "missing genotype index is outside variant values",
            ));
        }
        if previous.is_some_and(|previous| index <= previous) {
            return Err(GenoioError::invalid_source(
                "<filter>",
                "missing genotype indices must be sorted and unique",
            ));
        }
        previous = Some(index);
    }
    Ok(())
}

/// Compute variant statistics from hard-call category counts.
///
/// Counts are kept as `u64` while accumulating and narrowed only after overflow
/// checks, so large cohorts fail with a metadata error instead of wrapping.
pub fn variant_stats_from_counts(
    hom_ref_count: u64,
    het_count: u64,
    hom_alt_count: u64,
    missing_count: u64,
) -> Result<VariantStats, GenoioError> {
    let called_count = hom_ref_count
        .checked_add(het_count)
        .and_then(|count| count.checked_add(hom_alt_count))
        .ok_or_else(|| {
            GenoioError::invalid_source(
                "<filter>",
                "called genotype count exceeds supported metadata range",
            )
        })?;
    let total = called_count.checked_add(missing_count).ok_or_else(|| {
        GenoioError::invalid_source(
            "<filter>",
            "genotype count exceeds supported metadata range",
        )
    })?;
    let n_called = u32::try_from(called_count).map_err(|_| {
        GenoioError::invalid_source(
            "<filter>",
            "called genotype count exceeds supported metadata range",
        )
    })?;

    let missing_rate = if total == 0 {
        0.0
    } else {
        missing_count as f64 / total as f64
    };
    if n_called == 0 {
        return Ok(VariantStats {
            af: None,
            maf: None,
            mac: None,
            missing_rate,
            n_called,
            polymorphic: false,
        });
    }

    let allele_count = het_count
        .checked_add(hom_alt_count.checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source("<filter>", "allele count exceeds supported metadata range")
        })?)
        .ok_or_else(|| {
            GenoioError::invalid_source("<filter>", "allele count exceeds supported metadata range")
        })?;
    let called_alleles = 2_u64 * u64::from(n_called);
    let af = allele_count as f64 / called_alleles as f64;
    let maf = af.min(1.0 - af);
    let mac = allele_count.min(called_alleles - allele_count);
    let mac = u32::try_from(mac).map_err(|_| {
        GenoioError::invalid_source(
            "<filter>",
            "minor allele count exceeds supported metadata range",
        )
    })?;
    Ok(VariantStats {
        af: Some(af),
        maf: Some(maf),
        mac: Some(f64::from(mac)),
        missing_rate,
        n_called,
        polymorphic: mac > 0,
    })
}

fn variant_stats_from_dosage_count(
    allele_count: f64,
    called_count: u64,
    missing_count: u64,
) -> Result<VariantStats, GenoioError> {
    if !allele_count.is_finite() || allele_count < 0.0 {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "allele dosage count is outside the supported range",
        ));
    }
    let total = called_count.checked_add(missing_count).ok_or_else(|| {
        GenoioError::invalid_source(
            "<filter>",
            "genotype count exceeds supported metadata range",
        )
    })?;
    let n_called = u32::try_from(called_count).map_err(|_| {
        GenoioError::invalid_source(
            "<filter>",
            "called genotype count exceeds supported metadata range",
        )
    })?;

    let missing_rate = if total == 0 {
        0.0
    } else {
        missing_count as f64 / total as f64
    };
    if n_called == 0 {
        return Ok(VariantStats {
            af: None,
            maf: None,
            mac: None,
            missing_rate,
            n_called,
            polymorphic: false,
        });
    }

    let called_alleles = 2.0 * f64::from(n_called);
    if allele_count > called_alleles {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "allele dosage count exceeds called allele count",
        ));
    }
    let af = allele_count / called_alleles;
    let maf = af.min(1.0 - af);
    let mac = allele_count.min(called_alleles - allele_count);
    Ok(VariantStats {
        af: Some(af),
        maf: Some(maf),
        mac: Some(mac),
        missing_rate,
        n_called,
        polymorphic: mac > 0.0,
    })
}

/// Attach computed genotype statistics to variant metadata.
pub fn attach_variant_stats(variant: &mut VariantRecord, stats: VariantStats) {
    variant.af = stats.af.map(|value| value as f32);
    variant.maf = stats.maf.map(|value| value as f32);
    // Public variant metadata keeps MAC integer-valued. Dosage filters still
    // evaluate fractional MAC internally via VariantStats.
    variant.mac = stats.mac.and_then(exact_u32_from_f64);
    variant.missing_rate = Some(stats.missing_rate as f32);
    variant.n_called = Some(stats.n_called);
}

fn exact_u32_from_f64(value: f64) -> Option<u32> {
    if value.is_finite() && value.fract() == 0.0 && value >= 0.0 && value <= f64::from(u32::MAX) {
        Some(value as u32)
    } else {
        None
    }
}

fn is_snp(variant: &(impl VariantMetadataView + ?Sized)) -> bool {
    is_biallelic(variant) && variant.a0().len() == 1 && variant.a1().len() == 1
}

fn is_biallelic(variant: &(impl VariantMetadataView + ?Sized)) -> bool {
    variant
        .alt_allele()
        .is_none_or(|alt_allele| !alt_allele.contains(','))
}

fn params_object(params: &Value) -> Result<&serde_json::Map<String, Value>, GenoioError> {
    match params {
        Value::Object(object) => Ok(object),
        _ => Err(GenoioError::invalid_source(
            "<filter>",
            "predicate params must be a JSON object",
        )),
    }
}

fn expect_no_params(params: &Value) -> Result<(), GenoioError> {
    let object = params_object(params)?;
    if object.is_empty() {
        Ok(())
    } else {
        Err(GenoioError::invalid_source(
            "<filter>",
            "predicate does not accept parameters",
        ))
    }
}

fn required_string(params: &Value, key: &str) -> Result<String, GenoioError> {
    match params_object(params)?.get(key) {
        Some(Value::String(value)) if !value.is_empty() => Ok(value.clone()),
        _ => Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} must be a non-empty string"),
        )),
    }
}

fn required_string_set(params: &Value, key: &str) -> Result<BTreeSet<String>, GenoioError> {
    match params_object(params)?.get(key) {
        Some(Value::Array(values)) => {
            let mut set = BTreeSet::new();
            for value in values {
                let Value::String(text) = value else {
                    return Err(GenoioError::invalid_source(
                        "<filter>",
                        format!("predicate parameter {key:?} must contain only strings"),
                    ));
                };
                if !set.insert(text.clone()) {
                    return Err(GenoioError::invalid_source(
                        "<filter>",
                        format!("predicate parameter {key:?} must not contain duplicates"),
                    ));
                }
            }
            Ok(set)
        }
        _ => Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} must be a string array"),
        )),
    }
}

fn optional_rate(params: &Value, key: &str) -> Result<Option<f32>, GenoioError> {
    match params_object(params)?.get(key) {
        Some(value) => Ok(Some(value_to_rate(key, value)?)),
        None => Ok(None),
    }
}

fn required_rate(params: &Value, key: &str) -> Result<f32, GenoioError> {
    match optional_rate(params, key)? {
        Some(value) => Ok(value),
        None => Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} is required"),
        )),
    }
}

fn value_to_rate(key: &str, value: &Value) -> Result<f32, GenoioError> {
    let Some(number) = value.as_f64() else {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} must be numeric"),
        ));
    };
    if !(0.0..=1.0).contains(&number) {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} must be between 0 and 1"),
        ));
    }
    Ok(number as f32)
}

fn optional_nonnegative_f32(params: &Value, key: &str) -> Result<Option<f32>, GenoioError> {
    let Some(value) = params.get(key) else {
        return Ok(None);
    };
    let number = value.as_f64().ok_or_else(|| {
        GenoioError::invalid_source(
            "<filter>",
            format!("{key} must be a non-negative finite number"),
        )
    })?;
    if !number.is_finite() || number < 0.0 {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("{key} must be a non-negative finite number"),
        ));
    }
    Ok(Some(number as f32))
}

fn optional_u32(params: &Value, key: &str) -> Result<Option<u32>, GenoioError> {
    match params_object(params)?.get(key) {
        Some(value) => {
            let Some(number) = value.as_u64() else {
                return Err(GenoioError::invalid_source(
                    "<filter>",
                    format!("predicate parameter {key:?} must be a non-negative integer"),
                ));
            };
            Ok(Some(u32::try_from(number).map_err(|_| {
                GenoioError::invalid_source(
                    "<filter>",
                    format!("predicate parameter {key:?} is out of range"),
                )
            })?))
        }
        None => Ok(None),
    }
}

fn validate_range<T: PartialOrd>(
    name: &str,
    min: Option<T>,
    max: Option<T>,
) -> Result<(), GenoioError> {
    if min.is_none() && max.is_none() {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("{name} predicate requires at least one threshold"),
        ));
    }
    if min.zip(max).is_some_and(|(min, max)| min > max) {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("{name} predicate min must be <= max"),
        ));
    }
    Ok(())
}

fn parse_region(value: &str) -> Result<(String, u32, u32), GenoioError> {
    let Some((chrom, coordinates)) = value.split_once(':') else {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "invalid region syntax; expected chrom:start-end",
        ));
    };
    let Some((start_text, end_text)) = coordinates.split_once('-') else {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "invalid region syntax; expected chrom:start-end",
        ));
    };
    let start = start_text.parse::<u32>().map_err(|error| {
        GenoioError::invalid_source(
            "<filter>",
            format!("invalid region start coordinate: {error}"),
        )
    })?;
    let end = end_text.parse::<u32>().map_err(|error| {
        GenoioError::invalid_source(
            "<filter>",
            format!("invalid region end coordinate: {error}"),
        )
    })?;
    if chrom.is_empty() || start == 0 || end < start {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "invalid region coordinates; expected 1-based start <= end",
        ));
    }
    Ok((chrom.to_string(), start, end))
}
