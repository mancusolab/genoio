// pattern: Functional Core

use genoio_core::{
    variant_stats_from_counts, GenoioError, GenotypeFilterConjunction, GenotypeFilterPlan,
    VariantFilter, VariantMetadataView, VariantStats,
};

use crate::error::Result;

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct HardcallCounts {
    pub(super) hom_ref: u64,
    pub(super) het: u64,
    pub(super) hom_alt: u64,
    pub(super) missing: u64,
}

impl HardcallCounts {
    pub(crate) fn record_hom_ref(&mut self) {
        self.hom_ref += 1;
    }

    pub(crate) fn record_het(&mut self) {
        self.het += 1;
    }

    pub(crate) fn record_hom_alt(&mut self) {
        self.hom_alt += 1;
    }

    pub(crate) fn record_missing(&mut self) {
        self.missing += 1;
    }

    pub(crate) fn evaluate_plan(self, plan: GenotypeFilterPlan) -> Result<Option<bool>> {
        match plan {
            GenotypeFilterPlan::Generic => Ok(None),
            GenotypeFilterPlan::Polymorphic => Ok(Some(self.is_polymorphic()?)),
            GenotypeFilterPlan::MacRange { min, max } => Ok(Some(self.mac_in_range(min, max)?)),
            GenotypeFilterPlan::MafRange { min, max } => Ok(Some(self.maf_in_range(min, max)?)),
            GenotypeFilterPlan::MissingRateMax { max } => {
                Ok(Some(self.missing_rate() <= f64::from(max)))
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
        if plan
            .missing_rate_max
            .is_some_and(|max| self.missing_rate() > f64::from(max))
        {
            return Ok(false);
        }
        Ok(true)
    }

    fn called_count(self) -> Result<u64> {
        self.hom_ref
            .checked_add(self.het)
            .and_then(|count| count.checked_add(self.hom_alt))
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    "<filter>",
                    "called genotype count exceeds supported metadata range",
                )
            })
    }

    fn total_count(self) -> Result<u64> {
        self.called_count()?
            .checked_add(self.missing)
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    "<filter>",
                    "genotype count exceeds supported metadata range",
                )
            })
    }

    fn allele_count(self) -> Result<Option<u64>> {
        if self.called_count()? == 0 {
            return Ok(None);
        }
        self.het
            .checked_add(self.hom_alt.checked_mul(2).ok_or_else(|| {
                GenoioError::invalid_source(
                    "<filter>",
                    "allele count exceeds supported metadata range",
                )
            })?)
            .map(Some)
            .ok_or_else(|| {
                GenoioError::invalid_source(
                    "<filter>",
                    "allele count exceeds supported metadata range",
                )
            })
    }

    fn minor_allele_count(self) -> Result<Option<u32>> {
        let Some(allele_count) = self.allele_count()? else {
            return Ok(None);
        };
        let called_alleles = self.called_count()?.checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source("<filter>", "allele count exceeds supported metadata range")
        })?;
        let mac = allele_count.min(called_alleles - allele_count);
        u32::try_from(mac).map(Some).map_err(|_| {
            GenoioError::invalid_source(
                "<filter>",
                "minor allele count exceeds supported metadata range",
            )
        })
    }

    fn missing_rate(self) -> f64 {
        let total = self.total_count().unwrap_or(0);
        if total == 0 {
            0.0
        } else {
            self.missing as f64 / total as f64
        }
    }

    fn is_polymorphic(self) -> Result<bool> {
        Ok(self.minor_allele_count()?.is_some_and(|mac| mac > 0))
    }

    fn mac_in_range(self, min: Option<u32>, max: Option<u32>) -> Result<bool> {
        let Some(mac) = self.minor_allele_count()? else {
            return Ok(false);
        };
        Ok(min.is_none_or(|threshold| mac >= threshold)
            && max.is_none_or(|threshold| mac <= threshold))
    }

    fn maf_in_range(self, min: Option<f32>, max: Option<f32>) -> Result<bool> {
        let Some(mac) = self.minor_allele_count()? else {
            return Ok(false);
        };
        let called_alleles = self.called_count()?.checked_mul(2).ok_or_else(|| {
            GenoioError::invalid_source("<filter>", "allele count exceeds supported metadata range")
        })?;
        let maf = f64::from(mac) / called_alleles as f64;
        Ok(min.is_none_or(|threshold| maf >= f64::from(threshold))
            && max.is_none_or(|threshold| maf <= f64::from(threshold)))
    }

    pub(crate) fn variant_stats(self) -> Result<VariantStats> {
        variant_stats_from_counts(self.hom_ref, self.het, self.hom_alt, self.missing)
    }
}

pub(crate) fn evaluate_hardcall_counts_filter<V: VariantMetadataView + ?Sized>(
    counts: HardcallCounts,
    filter: &VariantFilter,
    filter_plan: GenotypeFilterPlan,
    variant: Option<&V>,
    require_stats: bool,
) -> Result<(bool, Option<VariantStats>)> {
    if !require_stats {
        if let Some(retain) = counts.evaluate_plan(filter_plan)? {
            return Ok((retain, None));
        }
    }

    let stats = counts.variant_stats()?;
    let retain = if let Some(variant) = variant {
        filter.evaluate_view(variant, Some(&stats))
    } else {
        filter.evaluate_genotype_stats(&stats).ok_or_else(|| {
            GenoioError::internal_contract(
                "genotype-stats-only fast path received metadata-dependent filter",
            )
        })?
    };
    Ok((retain, Some(stats)))
}
