// pattern: Functional Core

use crate::{GenoioError, VariantRecord};

use super::model::VariantStats;

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
