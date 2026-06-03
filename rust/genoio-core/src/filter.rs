// pattern: Functional Core

use std::collections::BTreeSet;

use serde::Deserialize;
use serde_json::Value;

use crate::{MetadataError, VariantRecord};

#[derive(Debug, Clone, PartialEq)]
pub struct VariantFilter {
    expr: Expr,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionPredicate {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VariantStats {
    pub af: Option<f64>,
    pub maf: Option<f64>,
    pub mac: Option<u32>,
    pub missing_rate: f64,
    pub n_called: u32,
    pub polymorphic: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantWindow {
    pub start: usize,
    pub len: usize,
}

impl VariantWindow {
    pub fn contains(self, retained_index: usize) -> bool {
        retained_index >= self.start && retained_index < self.start.saturating_add(self.len)
    }
}

#[derive(Debug, Clone, PartialEq)]
enum Expr {
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
    pub fn from_json_value(value: Value) -> Result<Self, MetadataError> {
        let raw: RawExpr = serde_json::from_value(value).map_err(|error| {
            MetadataError::parse("<filter>", format!("invalid filter IR: {error}"))
        })?;
        Ok(Self {
            expr: Expr::from_raw(raw)?,
        })
    }

    pub fn metadata_decision(&self, variant: &VariantRecord) -> Option<bool> {
        self.expr.metadata_decision(variant)
    }

    pub fn evaluate(&self, variant: &VariantRecord, stats: Option<&VariantStats>) -> bool {
        self.expr.evaluate(variant, stats)
    }

    pub fn requires_genotype_stats(&self) -> bool {
        self.expr.requires_genotype_stats()
    }

    pub fn has_region_predicate(&self) -> bool {
        self.expr.has_region_predicate()
    }

    pub fn concrete_region_pushdown(&self) -> Option<RegionPredicate> {
        self.expr.concrete_region_pushdown()
    }
}

impl Expr {
    fn from_raw(raw: RawExpr) -> Result<Self, MetadataError> {
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

    fn metadata_decision(&self, variant: &VariantRecord) -> Option<bool> {
        match self {
            Self::Predicate(predicate) => predicate.metadata_decision(variant),
            Self::And(left, right) => match (
                left.metadata_decision(variant),
                right.metadata_decision(variant),
            ) {
                (Some(false), _) | (_, Some(false)) => Some(false),
                (Some(true), Some(true)) => Some(true),
                _ => None,
            },
            Self::Or(left, right) => match (
                left.metadata_decision(variant),
                right.metadata_decision(variant),
            ) {
                (Some(true), _) | (_, Some(true)) => Some(true),
                (Some(false), Some(false)) => Some(false),
                _ => None,
            },
            Self::Not(expr) => expr.metadata_decision(variant).map(|decision| !decision),
        }
    }

    fn evaluate(&self, variant: &VariantRecord, stats: Option<&VariantStats>) -> bool {
        match self {
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
            Self::Predicate(predicate) => predicate.requires_genotype_stats(),
            Self::And(left, right) | Self::Or(left, right) => {
                left.requires_genotype_stats() || right.requires_genotype_stats()
            }
            Self::Not(expr) => expr.requires_genotype_stats(),
        }
    }

    fn has_region_predicate(&self) -> bool {
        match self {
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
}

impl Predicate {
    fn from_raw(name: &str, params: Value) -> Result<Self, MetadataError> {
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
            other => Err(MetadataError::parse(
                "<filter>",
                format!("unknown predicate name: {other}"),
            )),
        }
    }

    fn metadata_decision(&self, variant: &VariantRecord) -> Option<bool> {
        match self {
            Self::Chrom(value) => Some(variant.chrom == *value),
            Self::Region { chrom, start, end } => {
                Some(variant.chrom == *chrom && variant.pos >= *start && variant.pos <= *end)
            }
            Self::IdIn(values) => Some(values.contains(&variant.id)),
            Self::Snp => Some(is_snp(variant)),
            Self::Biallelic => Some(is_biallelic(variant)),
            Self::Qual { min, max } => Some(variant.qual.is_some_and(|qual| {
                min.is_none_or(|threshold| qual >= threshold)
                    && max.is_none_or(|threshold| qual <= threshold)
            })),
            Self::Maf { .. } | Self::Mac { .. } | Self::MissingRate { .. } | Self::Polymorphic => {
                None
            }
        }
    }

    fn evaluate(&self, variant: &VariantRecord, stats: Option<&VariantStats>) -> bool {
        match self {
            Self::Chrom(_)
            | Self::Region { .. }
            | Self::IdIn(_)
            | Self::Snp
            | Self::Biallelic
            | Self::Qual { .. } => self.metadata_decision(variant) == Some(true),
            Self::Maf { min, max } => stats.and_then(|stats| stats.maf).is_some_and(|maf| {
                min.is_none_or(|threshold| maf >= f64::from(threshold))
                    && max.is_none_or(|threshold| maf <= f64::from(threshold))
            }),
            Self::Mac { min, max } => stats.and_then(|stats| stats.mac).is_some_and(|mac| {
                min.is_none_or(|threshold| mac >= threshold)
                    && max.is_none_or(|threshold| mac <= threshold)
            }),
            Self::MissingRate { max } => {
                stats.is_some_and(|stats| stats.missing_rate <= f64::from(*max))
            }
            Self::Polymorphic => stats.is_some_and(|stats| stats.polymorphic),
        }
    }

    fn requires_genotype_stats(&self) -> bool {
        matches!(
            self,
            Self::Maf { .. } | Self::Mac { .. } | Self::MissingRate { .. } | Self::Polymorphic
        )
    }
}

pub fn compute_variant_stats(
    values: &[f32],
    missing_mask: &[bool],
) -> Result<VariantStats, MetadataError> {
    if values.len() != missing_mask.len() {
        return Err(MetadataError::parse(
            "<filter>",
            "variant values and missing mask lengths differ",
        ));
    }

    let mut allele_count = 0_u64;
    let mut n_called = 0_u64;
    for (value, missing) in values.iter().zip(missing_mask) {
        if *missing {
            continue;
        }
        allele_count += discrete_allele_count(*value)?;
        n_called += 1;
    }
    let n_called = u32::try_from(n_called).map_err(|_| {
        MetadataError::parse(
            "<filter>",
            "called genotype count exceeds supported metadata range",
        )
    })?;

    let total = values.len();
    let missing_count = total - usize::try_from(n_called).unwrap_or(usize::MAX);
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

    let called_alleles = 2_u64 * u64::from(n_called);
    let af = allele_count as f64 / called_alleles as f64;
    let maf = af.min(1.0 - af);
    let mac = allele_count.min(called_alleles - allele_count);
    let mac = u32::try_from(mac).map_err(|_| {
        MetadataError::parse(
            "<filter>",
            "minor allele count exceeds supported metadata range",
        )
    })?;
    Ok(VariantStats {
        af: Some(af),
        maf: Some(maf),
        mac: Some(mac),
        missing_rate,
        n_called,
        polymorphic: mac > 0,
    })
}

pub fn attach_variant_stats(variant: &mut VariantRecord, stats: VariantStats) {
    variant.af = stats.af.map(|value| value as f32);
    variant.maf = stats.maf.map(|value| value as f32);
    variant.mac = stats.mac;
    variant.missing_rate = Some(stats.missing_rate as f32);
    variant.n_called = Some(stats.n_called);
}

fn discrete_allele_count(value: f32) -> Result<u64, MetadataError> {
    match value {
        0.0 => Ok(0),
        1.0 => Ok(1),
        2.0 => Ok(2),
        other => Err(MetadataError::parse(
            "<filter>",
            format!("genotype statistics require discrete 0/1/2 values; observed {other}"),
        )),
    }
}

fn is_snp(variant: &VariantRecord) -> bool {
    is_biallelic(variant) && variant.a0.len() == 1 && variant.a1.len() == 1
}

fn is_biallelic(variant: &VariantRecord) -> bool {
    variant
        .alt_allele
        .as_ref()
        .is_none_or(|alt_allele| !alt_allele.contains(','))
}

fn params_object(params: &Value) -> Result<&serde_json::Map<String, Value>, MetadataError> {
    match params {
        Value::Object(object) => Ok(object),
        _ => Err(MetadataError::parse(
            "<filter>",
            "predicate params must be a JSON object",
        )),
    }
}

fn expect_no_params(params: &Value) -> Result<(), MetadataError> {
    let object = params_object(params)?;
    if object.is_empty() {
        Ok(())
    } else {
        Err(MetadataError::parse(
            "<filter>",
            "predicate does not accept parameters",
        ))
    }
}

fn required_string(params: &Value, key: &str) -> Result<String, MetadataError> {
    match params_object(params)?.get(key) {
        Some(Value::String(value)) if !value.is_empty() => Ok(value.clone()),
        _ => Err(MetadataError::parse(
            "<filter>",
            format!("predicate parameter {key:?} must be a non-empty string"),
        )),
    }
}

fn required_string_set(params: &Value, key: &str) -> Result<BTreeSet<String>, MetadataError> {
    match params_object(params)?.get(key) {
        Some(Value::Array(values)) => {
            let mut set = BTreeSet::new();
            for value in values {
                let Value::String(text) = value else {
                    return Err(MetadataError::parse(
                        "<filter>",
                        format!("predicate parameter {key:?} must contain only strings"),
                    ));
                };
                if !set.insert(text.clone()) {
                    return Err(MetadataError::parse(
                        "<filter>",
                        format!("predicate parameter {key:?} must not contain duplicates"),
                    ));
                }
            }
            Ok(set)
        }
        _ => Err(MetadataError::parse(
            "<filter>",
            format!("predicate parameter {key:?} must be a string array"),
        )),
    }
}

fn optional_rate(params: &Value, key: &str) -> Result<Option<f32>, MetadataError> {
    match params_object(params)?.get(key) {
        Some(value) => Ok(Some(value_to_rate(key, value)?)),
        None => Ok(None),
    }
}

fn required_rate(params: &Value, key: &str) -> Result<f32, MetadataError> {
    match optional_rate(params, key)? {
        Some(value) => Ok(value),
        None => Err(MetadataError::parse(
            "<filter>",
            format!("predicate parameter {key:?} is required"),
        )),
    }
}

fn value_to_rate(key: &str, value: &Value) -> Result<f32, MetadataError> {
    let Some(number) = value.as_f64() else {
        return Err(MetadataError::parse(
            "<filter>",
            format!("predicate parameter {key:?} must be numeric"),
        ));
    };
    if !(0.0..=1.0).contains(&number) {
        return Err(MetadataError::parse(
            "<filter>",
            format!("predicate parameter {key:?} must be between 0 and 1"),
        ));
    }
    Ok(number as f32)
}

fn optional_nonnegative_f32(params: &Value, key: &str) -> Result<Option<f32>, MetadataError> {
    let Some(value) = params.get(key) else {
        return Ok(None);
    };
    let number = value.as_f64().ok_or_else(|| {
        MetadataError::parse(
            "<filter>",
            format!("{key} must be a non-negative finite number"),
        )
    })?;
    if !number.is_finite() || number < 0.0 {
        return Err(MetadataError::parse(
            "<filter>",
            format!("{key} must be a non-negative finite number"),
        ));
    }
    Ok(Some(number as f32))
}

fn optional_u32(params: &Value, key: &str) -> Result<Option<u32>, MetadataError> {
    match params_object(params)?.get(key) {
        Some(value) => {
            let Some(number) = value.as_u64() else {
                return Err(MetadataError::parse(
                    "<filter>",
                    format!("predicate parameter {key:?} must be a non-negative integer"),
                ));
            };
            Ok(Some(u32::try_from(number).map_err(|_| {
                MetadataError::parse(
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
) -> Result<(), MetadataError> {
    if min.is_none() && max.is_none() {
        return Err(MetadataError::parse(
            "<filter>",
            format!("{name} predicate requires at least one threshold"),
        ));
    }
    if min.zip(max).is_some_and(|(min, max)| min > max) {
        return Err(MetadataError::parse(
            "<filter>",
            format!("{name} predicate min must be <= max"),
        ));
    }
    Ok(())
}

fn parse_region(value: &str) -> Result<(String, u32, u32), MetadataError> {
    let Some((chrom, coordinates)) = value.split_once(':') else {
        return Err(MetadataError::parse(
            "<filter>",
            "invalid region syntax; expected chrom:start-end",
        ));
    };
    let Some((start_text, end_text)) = coordinates.split_once('-') else {
        return Err(MetadataError::parse(
            "<filter>",
            "invalid region syntax; expected chrom:start-end",
        ));
    };
    let start = start_text.parse::<u32>().map_err(|error| {
        MetadataError::parse(
            "<filter>",
            format!("invalid region start coordinate: {error}"),
        )
    })?;
    let end = end_text.parse::<u32>().map_err(|error| {
        MetadataError::parse(
            "<filter>",
            format!("invalid region end coordinate: {error}"),
        )
    })?;
    if chrom.is_empty() || start == 0 || end < start {
        return Err(MetadataError::parse(
            "<filter>",
            "invalid region coordinates; expected 1-based start <= end",
        ));
    }
    Ok((chrom.to_string(), start, end))
}
