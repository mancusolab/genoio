// pattern: Functional Core
//! Dosage filter evaluation helpers shared by genotype readers.
//!
//! Matrix-only reads can often decide genotype-stat filters from dosage counts
//! without building full `VariantStats`. Metadata-returning reads still attach
//! complete stats to retained variants.

use genoio_core::{
    compute_dosage_variant_stats, is_dosage_polymorphic, GenoioError, GenotypeFilterConjunction,
    GenotypeFilterPlan, VariantFilter, VariantRecord, VariantStats,
};

use crate::error::Result;

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct DosageFilterCounts {
    pub(crate) allele_count: f64,
    pub(crate) called_count: u64,
    pub(crate) missing_count: u64,
}

impl DosageFilterCounts {
    pub(crate) fn record_called_dosage(&mut self, value: f32) {
        self.allele_count += f64::from(value);
        self.called_count += 1;
    }

    pub(crate) fn evaluate_plan(self, plan: GenotypeFilterPlan) -> Result<Option<bool>> {
        match plan {
            GenotypeFilterPlan::Generic => Ok(None),
            GenotypeFilterPlan::Polymorphic => Ok(Some(self.is_polymorphic()?)),
            GenotypeFilterPlan::MacRange { min, max } => Ok(Some(self.mac_in_range(min, max)?)),
            GenotypeFilterPlan::MafRange { min, max } => Ok(Some(self.maf_in_range(min, max)?)),
            GenotypeFilterPlan::MissingRateMax { max } => {
                Ok(Some(self.missing_rate()? <= f64::from(max)))
            }
            GenotypeFilterPlan::Conjunction(plan) => Ok(Some(self.evaluate_conjunction(plan)?)),
        }
    }

    fn evaluate_conjunction(self, plan: GenotypeFilterConjunction) -> Result<bool> {
        if plan.polymorphic && !self.is_polymorphic()? {
            return Ok(false);
        }
        if (plan.mac_min.is_some() || plan.mac_max.is_some())
            && !self.mac_in_range(plan.mac_min, plan.mac_max)?
        {
            return Ok(false);
        }
        if (plan.maf_min.is_some() || plan.maf_max.is_some())
            && !self.maf_in_range(plan.maf_min, plan.maf_max)?
        {
            return Ok(false);
        }
        if let Some(max) = plan.missing_rate_max {
            if self.missing_rate()? > f64::from(max) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn called_alleles(self) -> Result<Option<f64>> {
        if self.called_count == 0 {
            return Ok(None);
        }
        let called_count = u32::try_from(self.called_count).map_err(|_| {
            GenoioError::invalid_source(
                "<filter>",
                "called genotype count exceeds supported metadata range",
            )
        })?;
        Ok(Some(2.0 * f64::from(called_count)))
    }

    fn total_count(self) -> Result<u64> {
        self.called_count
            .checked_add(self.missing_count)
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    "<filter>",
                    "genotype count exceeds supported metadata range",
                )
            })
    }

    fn minor_allele_count(self) -> Result<Option<f64>> {
        let Some(called_alleles) = self.called_alleles()? else {
            return Ok(None);
        };
        Ok(Some(
            self.allele_count.min(called_alleles - self.allele_count),
        ))
    }

    fn missing_rate(self) -> Result<f64> {
        let total = self.total_count()?;
        if total == 0 {
            Ok(0.0)
        } else {
            Ok(self.missing_count as f64 / total as f64)
        }
    }

    fn is_polymorphic(self) -> Result<bool> {
        Ok(self.minor_allele_count()?.is_some_and(|mac| mac > 0.0))
    }

    fn mac_in_range(self, min: Option<u32>, max: Option<u32>) -> Result<bool> {
        let Some(mac) = self.minor_allele_count()? else {
            return Ok(false);
        };
        Ok(min.is_none_or(|threshold| mac >= f64::from(threshold))
            && max.is_none_or(|threshold| mac <= f64::from(threshold)))
    }

    fn maf_in_range(self, min: Option<f32>, max: Option<f32>) -> Result<bool> {
        let Some(mac) = self.minor_allele_count()? else {
            return Ok(false);
        };
        let Some(called_alleles) = self.called_alleles()? else {
            return Ok(false);
        };
        let maf = mac / called_alleles;
        Ok(min.is_none_or(|threshold| maf >= f64::from(threshold))
            && max.is_none_or(|threshold| maf <= f64::from(threshold)))
    }
}

fn dosage_counts_for_filter(values: &[f32], missing: &[bool]) -> Result<DosageFilterCounts> {
    if values.len() != missing.len() {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "variant values and missing mask lengths differ",
        ));
    }

    let mut counts = DosageFilterCounts::default();
    for (value, is_missing) in values.iter().zip(missing) {
        if *is_missing {
            counts.missing_count += 1;
            continue;
        }
        if !(0.0..=2.0).contains(value) {
            return Err(GenoioError::invalid_source(
                "<filter>",
                format!("dosage statistics require values in [0, 2]; observed {value}"),
            ));
        }
        counts.record_called_dosage(*value);
    }
    Ok(counts)
}

pub(crate) fn evaluate_dosage_filter(
    values: &[f32],
    missing: &[bool],
    filter: &VariantFilter,
    variant: &VariantRecord,
    require_stats: bool,
) -> Result<(bool, Option<VariantStats>)> {
    let plan = filter.genotype_filter_plan();
    if !require_stats {
        // Matrix-only reads only need the retain/drop decision. The caller has
        // already run metadata partial evaluation, so compiled genotype plans
        // can bypass `VariantStats` construction for common dosage predicates.
        if matches!(plan, GenotypeFilterPlan::Polymorphic) {
            return Ok((is_dosage_polymorphic(values, missing)?, None));
        }
        if let Some(retain) = dosage_counts_for_filter(values, missing)?.evaluate_plan(plan)? {
            return Ok((retain, None));
        }
    }

    let stats = compute_dosage_variant_stats(values, missing)?;
    Ok((filter.evaluate(variant, Some(&stats)), Some(stats)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn test_variant() -> VariantRecord {
        VariantRecord {
            chrom: "1".to_string(),
            pos: 10,
            id: "rs1".to_string(),
            a0: "A".to_string(),
            a1: "G".to_string(),
            ref_allele: Some("A".to_string()),
            alt_allele: Some("G".to_string()),
            source_a0: "A".to_string(),
            source_a1: "G".to_string(),
            flipped: false,
            qual: None,
            af: None,
            maf: None,
            mac: None,
            missing_rate: None,
            n_called: None,
        }
    }

    fn genotype_filter(name: &str, params: serde_json::Value) -> VariantFilter {
        VariantFilter::from_json_value(json!({
            "op": "predicate",
            "name": name,
            "params": params,
        }))
        .unwrap()
    }

    fn dosage_fixture() -> ([f32; 4], [bool; 4]) {
        ([0.0, 1.0, 2.0, 2.0], [false, false, false, true])
    }

    #[test]
    fn dosage_filter_plan_evaluates_mac_maf_and_missing_rate() {
        let (values, missing) = dosage_fixture();
        let counts = dosage_counts_for_filter(&values, &missing).unwrap();

        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MacRange {
                    min: Some(3),
                    max: Some(3),
                })
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MafRange {
                    min: Some(0.49),
                    max: Some(0.51),
                })
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MissingRateMax { max: 0.20 })
                .unwrap(),
            Some(false)
        );
    }

    #[test]
    fn dosage_filter_plan_evaluates_conjunctions() {
        let (values, missing) = dosage_fixture();
        let counts = dosage_counts_for_filter(&values, &missing).unwrap();

        let passing = GenotypeFilterConjunction {
            polymorphic: true,
            mac_min: Some(2),
            mac_max: Some(4),
            maf_min: Some(0.4),
            maf_max: Some(0.6),
            missing_rate_max: Some(0.3),
        };
        let failing = GenotypeFilterConjunction {
            missing_rate_max: Some(0.2),
            ..passing
        };

        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::Conjunction(passing))
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::Conjunction(failing))
                .unwrap(),
            Some(false)
        );
    }

    #[test]
    fn matrix_only_dosage_filter_returns_decision_without_stats() {
        let (values, missing) = dosage_fixture();
        let filter = genotype_filter("maf", json!({ "min": 0.49, "max": 0.51 }));

        let (retain, stats) =
            evaluate_dosage_filter(&values, &missing, &filter, &test_variant(), false).unwrap();

        assert!(retain);
        assert_eq!(stats, None);
    }

    #[test]
    fn metadata_output_dosage_filter_keeps_stats() {
        let (values, missing) = dosage_fixture();
        let filter = genotype_filter("maf", json!({ "min": 0.49, "max": 0.51 }));

        let (retain, stats) =
            evaluate_dosage_filter(&values, &missing, &filter, &test_variant(), true).unwrap();

        assert!(retain);
        assert_eq!(stats.unwrap().missing_rate, 0.25);
    }

    #[test]
    fn dosage_filter_counts_match_variant_stats_thresholds() {
        let (values, missing) = dosage_fixture();
        let counts = dosage_counts_for_filter(&values, &missing).unwrap();
        let stats = compute_dosage_variant_stats(&values, &missing).unwrap();

        assert_eq!(counts.is_polymorphic().unwrap(), stats.polymorphic);
        assert_eq!(counts.missing_rate().unwrap(), stats.missing_rate);
        assert_eq!(counts.minor_allele_count().unwrap(), stats.mac);
    }

    #[test]
    fn dosage_filter_counts_handle_monomorphic_and_missing_variants() {
        let values = [0.0, 0.0, 0.0];
        let missing = [false, true, false];
        let counts = dosage_counts_for_filter(&values, &missing).unwrap();

        assert_eq!(counts.called_count, 2);
        assert_eq!(counts.missing_count, 1);
        assert_eq!(counts.allele_count, 0.0);
        assert!(!counts.is_polymorphic().unwrap());
        assert_eq!(counts.minor_allele_count().unwrap(), Some(0.0));
        assert!((counts.missing_rate().unwrap() - (1.0 / 3.0)).abs() < f64::EPSILON);
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::Polymorphic)
                .unwrap(),
            Some(false)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MacRange {
                    min: Some(1),
                    max: None,
                })
                .unwrap(),
            Some(false)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MafRange {
                    min: Some(0.01),
                    max: None,
                })
                .unwrap(),
            Some(false)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MissingRateMax { max: 0.34 })
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MissingRateMax { max: 0.30 })
                .unwrap(),
            Some(false)
        );
    }

    #[test]
    fn dosage_filter_counts_handle_all_missing_variants() {
        let values = [0.0, 1.0, 2.0];
        let missing = [true, true, true];
        let counts = dosage_counts_for_filter(&values, &missing).unwrap();

        assert_eq!(counts.called_count, 0);
        assert_eq!(counts.missing_count, 3);
        assert_eq!(counts.minor_allele_count().unwrap(), None);
        assert_eq!(counts.missing_rate().unwrap(), 1.0);
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::Polymorphic)
                .unwrap(),
            Some(false)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MacRange {
                    min: Some(0),
                    max: None,
                })
                .unwrap(),
            Some(false)
        );
        assert_eq!(
            counts
                .evaluate_plan(GenotypeFilterPlan::MafRange {
                    min: Some(0.0),
                    max: Some(1.0),
                })
                .unwrap(),
            Some(false)
        );
    }

    #[test]
    fn dosage_filter_counts_reject_invalid_inputs() {
        let length_error = dosage_counts_for_filter(&[0.0, 1.0], &[false])
            .expect_err("mismatched value and missing lengths should fail");
        assert!(length_error
            .to_string()
            .contains("variant values and missing mask lengths differ"));

        let high_value =
            dosage_counts_for_filter(&[2.1], &[false]).expect_err("dosage above two should fail");
        assert!(high_value.to_string().contains("values in [0, 2]"));

        let low_value =
            dosage_counts_for_filter(&[-0.1], &[false]).expect_err("negative dosage should fail");
        assert!(low_value.to_string().contains("values in [0, 2]"));
    }
}
