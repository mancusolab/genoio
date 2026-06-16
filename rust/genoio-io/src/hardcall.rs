// pattern: Functional Core

use genoio_core::{
    variant_stats_from_counts, GenoioError, GenotypeFilterConjunction, GenotypeFilterPlan,
    VariantStats,
};

use crate::error::Result;

pub(crate) const HARDCALL_BATCH_SIZE: usize = 64;

#[derive(Debug, Clone, Default)]
pub(crate) struct PackedHardcalls {
    words: Vec<u64>,
    sample_ct: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct HardcallCounts {
    hom_ref: u64,
    het: u64,
    hom_alt: u64,
    missing: u64,
}

impl HardcallCounts {
    fn evaluate_plan(self, plan: GenotypeFilterPlan) -> Result<Option<bool>> {
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
}

impl PackedHardcalls {
    const SAMPLES_PER_WORD: usize = 32;
    const BITS_PER_SAMPLE: usize = 2;
    const LOW_BITS: u64 = 0x5555_5555_5555_5555;

    pub(crate) fn resize(&mut self, sample_ct: usize) {
        self.sample_ct = sample_ct;
        self.words
            .resize(sample_ct.div_ceil(Self::SAMPLES_PER_WORD), 0);
        self.mask_unused_slots();
    }

    pub(crate) fn load_pgen_payload(&mut self, payload: &[u8], sample_ct: usize) {
        self.load_payload(payload, sample_ct);
    }

    pub(crate) fn load_plink1_bed_payload(&mut self, payload: &[u8], sample_ct: usize) {
        self.load_payload(payload, sample_ct);
        self.rotate_plink1_to_canonical();
    }

    fn load_payload(&mut self, payload: &[u8], sample_ct: usize) {
        self.resize(sample_ct);
        for (word_index, word) in self.words.iter_mut().enumerate() {
            let start = word_index * 8;
            let end = payload.len().min(start + 8);
            let mut bytes = [0_u8; 8];
            if start < end {
                bytes[..end - start].copy_from_slice(&payload[start..end]);
            }
            *word = u64::from_le_bytes(bytes);
        }
        self.mask_unused_slots();
    }

    pub(crate) fn clear_to(&mut self, category: u8) {
        let category = u64::from(category & 0b11);
        let mut word = 0_u64;
        for slot_index in 0..Self::SAMPLES_PER_WORD {
            word |= category << (slot_index * Self::BITS_PER_SAMPLE);
        }
        self.words.fill(word);
        self.mask_unused_slots();
    }

    pub(crate) fn set(&mut self, sample_index: usize, category: u8) {
        debug_assert!(sample_index < self.sample_ct);
        let word_index = sample_index / Self::SAMPLES_PER_WORD;
        let shift = (sample_index % Self::SAMPLES_PER_WORD) * Self::BITS_PER_SAMPLE;
        let mask = 0b11_u64 << shift;
        self.words[word_index] =
            (self.words[word_index] & !mask) | (u64::from(category & 0b11) << shift);
    }

    pub(crate) fn get(&self, sample_index: usize) -> u8 {
        debug_assert!(sample_index < self.sample_ct);
        let word_index = sample_index / Self::SAMPLES_PER_WORD;
        let shift = (sample_index % Self::SAMPLES_PER_WORD) * Self::BITS_PER_SAMPLE;
        ((self.words[word_index] >> shift) & 0b11) as u8
    }

    pub(crate) fn sample_ct(&self) -> usize {
        self.sample_ct
    }

    pub(crate) fn copy_from(&mut self, other: &Self) {
        self.sample_ct = other.sample_ct;
        self.words.clear();
        self.words.extend_from_slice(&other.words);
    }

    pub(crate) fn invert_0_2(&mut self) {
        for sample_index in 0..self.sample_ct {
            match self.get(sample_index) {
                0 => self.set(sample_index, 2),
                2 => self.set(sample_index, 0),
                _ => {}
            }
        }
    }

    pub(crate) fn expand_selected(
        &self,
        source_indices: &[usize],
        values: &mut Vec<f32>,
        missing: &mut Vec<bool>,
    ) {
        values.clear();
        missing.clear();
        for source_index in source_indices {
            let (value, is_missing) = decode_hardcall_code(self.get(*source_index));
            values.push(value);
            missing.push(is_missing);
        }
    }

    pub(crate) fn stats_for_selected(&self, source_indices: &[usize]) -> Result<VariantStats> {
        let counts = self.counts_for_selected(source_indices)?;
        variant_stats_from_counts(counts.hom_ref, counts.het, counts.hom_alt, counts.missing)
    }

    fn counts_for_selected(&self, source_indices: &[usize]) -> Result<HardcallCounts> {
        let mut counts = HardcallCounts::default();
        for source_index in source_indices {
            if *source_index >= self.sample_ct {
                return Err(GenoioError::invalid_source(
                    "<hardcall>",
                    "selected sample index is outside hard-call sample count",
                ));
            }
            match self.get(*source_index) {
                0 => counts.hom_ref += 1,
                1 => counts.het += 1,
                2 => counts.hom_alt += 1,
                3 => counts.missing += 1,
                _ => unreachable!("two-bit hard-call code should be masked"),
            }
        }
        Ok(counts)
    }

    pub(crate) fn stats_for_selection(
        &self,
        source_indices: &[usize],
        all_samples_selected: bool,
    ) -> Result<VariantStats> {
        if all_samples_selected {
            return self.stats_for_all();
        }
        self.stats_for_selected(source_indices)
    }

    pub(crate) fn evaluate_filter_plan_for_selection(
        &self,
        plan: GenotypeFilterPlan,
        source_indices: &[usize],
        all_samples_selected: bool,
    ) -> Result<Option<bool>> {
        if matches!(plan, GenotypeFilterPlan::Generic) {
            return Ok(None);
        }
        if matches!(plan, GenotypeFilterPlan::Polymorphic) {
            return self
                .is_polymorphic_for_selection(source_indices, all_samples_selected)
                .map(Some);
        }
        let counts = if all_samples_selected {
            self.counts_for_all()
        } else {
            self.counts_for_selected(source_indices)?
        };
        counts.evaluate_plan(plan)
    }

    pub(crate) fn is_polymorphic_for_selection(
        &self,
        source_indices: &[usize],
        all_samples_selected: bool,
    ) -> Result<bool> {
        if all_samples_selected {
            return Ok(self.is_polymorphic_for_all());
        }
        self.is_polymorphic_for_selected(source_indices)
    }

    fn is_polymorphic_for_selected(&self, source_indices: &[usize]) -> Result<bool> {
        let mut has_ref = false;
        let mut has_alt = false;
        for source_index in source_indices {
            if *source_index >= self.sample_ct {
                return Err(GenoioError::invalid_source(
                    "<hardcall>",
                    "selected sample index is outside hard-call sample count",
                ));
            }
            match self.get(*source_index) {
                0 => has_ref = true,
                1 => {
                    has_ref = true;
                    has_alt = true;
                }
                2 => has_alt = true,
                3 => {}
                _ => unreachable!("two-bit hard-call code should be masked"),
            }
            if has_ref && has_alt {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn is_polymorphic_for_all(&self) -> bool {
        let mut has_ref = false;
        let mut has_alt = false;
        for (word_index, word) in self.words.iter().copied().enumerate() {
            let valid_lanes = self.valid_lane_mask(word_index);
            let low_bits = word & valid_lanes;
            let high_bits = (word >> 1) & valid_lanes;

            has_ref |= (!high_bits & valid_lanes) != 0;
            has_alt |= ((low_bits ^ high_bits) & valid_lanes) != 0;
            if has_ref && has_alt {
                return true;
            }
        }
        false
    }

    fn stats_for_all(&self) -> Result<VariantStats> {
        let counts = self.counts_for_all();
        variant_stats_from_counts(counts.hom_ref, counts.het, counts.hom_alt, counts.missing)
    }

    fn counts_for_all(&self) -> HardcallCounts {
        let mut counts = HardcallCounts::default();

        for (word_index, word) in self.words.iter().copied().enumerate() {
            let valid_lanes = self.valid_lane_mask(word_index);
            let low_bits = word & valid_lanes;
            let high_bits = (word >> 1) & valid_lanes;

            counts.hom_ref += ((!low_bits & !high_bits) & valid_lanes).count_ones() as u64;
            counts.het += (low_bits & !high_bits & valid_lanes).count_ones() as u64;
            counts.hom_alt += (!low_bits & high_bits & valid_lanes).count_ones() as u64;
            counts.missing += (low_bits & high_bits & valid_lanes).count_ones() as u64;
        }

        counts
    }

    fn valid_lane_mask(&self, word_index: usize) -> u64 {
        let word_start = word_index * Self::SAMPLES_PER_WORD;
        let remaining = self.sample_ct.saturating_sub(word_start);
        match remaining.min(Self::SAMPLES_PER_WORD) {
            0 => 0,
            Self::SAMPLES_PER_WORD => Self::LOW_BITS,
            used_slots => ((1_u64 << (used_slots * Self::BITS_PER_SAMPLE)) - 1) & Self::LOW_BITS,
        }
    }

    fn rotate_plink1_to_canonical(&mut self) {
        const HIGH_BITS: u64 = 0xaaaa_aaaa_aaaa_aaaa;
        for word in &mut self.words {
            let old_high_bits_in_low_position = (*word >> 1) & Self::LOW_BITS;
            let new_low_bits = (*word & Self::LOW_BITS) ^ old_high_bits_in_low_position;
            let new_high_bits = (!*word) & HIGH_BITS;
            *word = new_low_bits | new_high_bits;
        }
        self.mask_unused_slots();
    }

    fn mask_unused_slots(&mut self) {
        let used_slots = self.sample_ct % Self::SAMPLES_PER_WORD;
        if used_slots == 0 || self.words.is_empty() {
            return;
        }
        let used_bits = used_slots * Self::BITS_PER_SAMPLE;
        let mask = (1_u64 << used_bits) - 1;
        if let Some(last_word) = self.words.last_mut() {
            *last_word &= mask;
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct HardcallBatch {
    variants: Vec<PackedHardcalls>,
    sample_ct: usize,
}

impl HardcallBatch {
    pub(crate) fn new(sample_ct: usize) -> Self {
        Self {
            variants: Vec::with_capacity(HARDCALL_BATCH_SIZE),
            sample_ct,
        }
    }

    pub(crate) fn push(&mut self, packed: &PackedHardcalls) {
        debug_assert_eq!(packed.sample_ct, self.sample_ct);
        let mut copy = PackedHardcalls::default();
        copy.copy_from(packed);
        self.variants.push(copy);
    }

    pub(crate) fn len(&self) -> usize {
        self.variants.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.variants.is_empty()
    }

    pub(crate) fn is_full(&self) -> bool {
        self.variants.len() == HARDCALL_BATCH_SIZE
    }

    pub(crate) fn clear(&mut self) {
        self.variants.clear();
    }

    pub(crate) fn expand_into_sample_major(
        &self,
        source_indices: &[usize],
        variant_start: usize,
        n_variants: usize,
        out_values: &mut [f32],
        out_missing: &mut [bool],
    ) {
        debug_assert!(variant_start + self.variants.len() <= n_variants);
        debug_assert_eq!(out_values.len(), source_indices.len() * n_variants);
        debug_assert_eq!(out_missing.len(), source_indices.len() * n_variants);

        for (sample_index, source_index) in source_indices.iter().copied().enumerate() {
            debug_assert!(source_index < self.sample_ct);
            let row_start = sample_index * n_variants;
            for (batch_variant_index, packed) in self.variants.iter().enumerate() {
                let variant_index = variant_start + batch_variant_index;
                let (value, is_missing) = decode_hardcall_code(packed.get(source_index));
                out_values[row_start + variant_index] = value;
                out_missing[row_start + variant_index] = is_missing;
            }
        }
    }
}

pub(crate) fn decode_hardcall_code(code: u8) -> (f32, bool) {
    match code {
        0b00 => (0.0, false),
        0b01 => (1.0, false),
        0b10 => (2.0, false),
        0b11 => (0.0, true),
        _ => unreachable!("two-bit hard-call code should be masked"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn append_variant_to_sample_major(
        values: &[f32],
        missing: &[bool],
        variant_index: usize,
        n_variants: usize,
        out_values: &mut [f32],
        out_missing: &mut [bool],
    ) {
        debug_assert_eq!(values.len(), missing.len());
        debug_assert!(variant_index < n_variants);
        debug_assert_eq!(out_values.len(), values.len() * n_variants);
        debug_assert_eq!(out_missing.len(), missing.len() * n_variants);

        for (sample_index, (&value, &is_missing)) in values.iter().zip(missing).enumerate() {
            let offset = sample_index * n_variants + variant_index;
            out_values[offset] = value;
            out_missing[offset] = is_missing;
        }
    }

    #[test]
    fn plink1_payload_rotates_to_canonical_hardcall_codes() {
        let mut packed = PackedHardcalls::default();

        packed.load_plink1_bed_payload(&[0b1110_0100], 4);

        assert_eq!(packed.get(0), 2);
        assert_eq!(packed.get(1), 3);
        assert_eq!(packed.get(2), 1);
        assert_eq!(packed.get(3), 0);
    }

    #[test]
    fn packed_hardcalls_round_trip_and_expand_selected() {
        let mut packed = PackedHardcalls::default();
        packed.resize(35);
        packed.clear_to(0);
        packed.set(0, 0);
        packed.set(1, 1);
        packed.set(2, 2);
        packed.set(3, 3);
        packed.set(34, 2);

        assert_eq!(packed.get(0), 0);
        assert_eq!(packed.get(1), 1);
        assert_eq!(packed.get(2), 2);
        assert_eq!(packed.get(3), 3);
        assert_eq!(packed.get(34), 2);

        let mut values = vec![99.0];
        let mut missing = vec![true];
        packed.expand_selected(&[3, 1, 34, 0], &mut values, &mut missing);

        assert_eq!(values, vec![0.0, 1.0, 2.0, 0.0]);
        assert_eq!(missing, vec![true, false, false, false]);
    }

    #[test]
    fn hardcall_batch_expands_like_variant_at_a_time() {
        let sample_ct = 5;
        let source_indices = (0..sample_ct).collect::<Vec<_>>();
        let n_variants = HARDCALL_BATCH_SIZE + 3;
        let mut packed_variants = Vec::with_capacity(n_variants);
        let mut expected_values = vec![0.0; sample_ct * n_variants];
        let mut expected_missing = vec![false; sample_ct * n_variants];
        let mut scratch_values = Vec::new();
        let mut scratch_missing = Vec::new();

        for variant_index in 0..n_variants {
            let mut packed = PackedHardcalls::default();
            packed.resize(sample_ct);
            for sample_index in 0..sample_ct {
                packed.set(sample_index, ((variant_index + sample_index) % 4) as u8);
            }
            packed.expand_selected(&source_indices, &mut scratch_values, &mut scratch_missing);
            append_variant_to_sample_major(
                &scratch_values,
                &scratch_missing,
                variant_index,
                n_variants,
                &mut expected_values,
                &mut expected_missing,
            );
            packed_variants.push(packed);
        }

        let mut batch = HardcallBatch::new(sample_ct);
        let mut actual_values = vec![0.0; sample_ct * n_variants];
        let mut actual_missing = vec![false; sample_ct * n_variants];
        let mut batch_start = 0;
        for packed in &packed_variants {
            batch.push(packed);
            if batch.is_full() {
                batch.expand_into_sample_major(
                    &source_indices,
                    batch_start,
                    n_variants,
                    &mut actual_values,
                    &mut actual_missing,
                );
                batch_start += batch.len();
                batch.clear();
            }
        }
        batch.expand_into_sample_major(
            &source_indices,
            batch_start,
            n_variants,
            &mut actual_values,
            &mut actual_missing,
        );

        assert_eq!(actual_values, expected_values);
        assert_eq!(actual_missing, expected_missing);
    }

    #[test]
    fn packed_hardcalls_copy_and_invert_0_2() {
        let mut source = PackedHardcalls::default();
        source.resize(5);
        source.clear_to(3);
        source.set(0, 0);
        source.set(1, 1);
        source.set(2, 2);

        let mut copy = PackedHardcalls::default();
        copy.copy_from(&source);
        copy.invert_0_2();

        assert_eq!(
            (0..5)
                .map(|sample_index| copy.get(sample_index))
                .collect::<Vec<_>>(),
            vec![2, 1, 0, 3, 3]
        );
        assert_eq!(
            (0..5)
                .map(|sample_index| source.get(sample_index))
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 3, 3]
        );
    }

    #[test]
    fn packed_hardcalls_loads_pgen_payload_and_masks_unused_trailing_slots() {
        let mut packed = PackedHardcalls::default();
        packed.load_pgen_payload(&[0b1110_0100, 0xff], 5);

        assert_eq!(
            (0..5)
                .map(|sample_index| packed.get(sample_index))
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 3, 3]
        );

        packed.resize(8);
        assert_eq!(
            (0..8)
                .map(|sample_index| packed.get(sample_index))
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 3, 3, 0, 0, 0]
        );
    }

    #[test]
    fn packed_hardcalls_stats_for_selected_matches_expanded_stats() {
        let mut packed = PackedHardcalls::default();
        packed.resize(8);
        for (sample_index, category) in [0, 1, 2, 3, 2, 0, 1, 3].into_iter().enumerate() {
            packed.set(sample_index, category);
        }

        for source_indices in [&[0, 1, 2, 3, 4, 5, 6, 7][..], &[7, 3][..], &[][..]] {
            let mut values = Vec::new();
            let mut missing = Vec::new();
            packed.expand_selected(source_indices, &mut values, &mut missing);

            let expected = genoio_core::compute_variant_stats(&values, &missing)
                .expect("expanded stats should compute");
            let actual = packed
                .stats_for_selected(source_indices)
                .expect("packed stats should compute");

            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn packed_hardcalls_detects_polymorphic_for_selection() {
        let mut packed = PackedHardcalls::default();
        packed.resize(8);
        for (sample_index, category) in [0, 0, 3, 0, 2, 3, 1, 0].into_iter().enumerate() {
            packed.set(sample_index, category);
        }

        assert!(packed
            .is_polymorphic_for_selection(&[0, 1, 2, 3, 4, 5, 6, 7], true)
            .expect("valid all-sample selection should evaluate"));
        assert!(packed
            .is_polymorphic_for_selection(&[0, 4], false)
            .expect("valid selected samples should evaluate"));
        assert!(packed
            .is_polymorphic_for_selection(&[1, 6], false)
            .expect("valid selected samples should evaluate"));
        assert!(!packed
            .is_polymorphic_for_selection(&[0, 1, 3, 7], false)
            .expect("valid selected samples should evaluate"));
        assert!(!packed
            .is_polymorphic_for_selection(&[2, 5], false)
            .expect("valid selected samples should evaluate"));
        assert!(packed.is_polymorphic_for_selection(&[8], false).is_err());
    }

    #[test]
    fn packed_hardcalls_evaluate_compiled_filter_plans_for_selection() {
        let mut packed = PackedHardcalls::default();
        packed.resize(8);
        for (sample_index, category) in [0, 1, 2, 3, 2, 0, 1, 3].into_iter().enumerate() {
            packed.set(sample_index, category);
        }

        assert_eq!(
            packed
                .evaluate_filter_plan_for_selection(
                    genoio_core::GenotypeFilterPlan::MacRange {
                        min: Some(2),
                        max: Some(8),
                    },
                    &[0, 1, 2, 3, 4, 5, 6, 7],
                    true,
                )
                .expect("valid plan should evaluate"),
            Some(true)
        );
        assert_eq!(
            packed
                .evaluate_filter_plan_for_selection(
                    genoio_core::GenotypeFilterPlan::MafRange {
                        min: Some(0.4),
                        max: Some(0.5),
                    },
                    &[0, 1, 2, 3, 4, 5, 6, 7],
                    true,
                )
                .expect("valid plan should evaluate"),
            Some(true)
        );
        assert_eq!(
            packed
                .evaluate_filter_plan_for_selection(
                    genoio_core::GenotypeFilterPlan::Conjunction(
                        genoio_core::GenotypeFilterConjunction {
                            polymorphic: true,
                            mac_min: Some(2),
                            mac_max: None,
                            maf_min: None,
                            maf_max: None,
                            missing_rate_max: Some(0.20),
                        },
                    ),
                    &[0, 1, 2, 3, 4, 5, 6, 7],
                    true,
                )
                .expect("valid plan should evaluate"),
            Some(false)
        );
        assert_eq!(
            packed
                .evaluate_filter_plan_for_selection(
                    genoio_core::GenotypeFilterPlan::MissingRateMax { max: 0.0 },
                    &[0, 5],
                    false,
                )
                .expect("valid selected plan should evaluate"),
            Some(true)
        );
        assert_eq!(
            packed
                .evaluate_filter_plan_for_selection(
                    genoio_core::GenotypeFilterPlan::Generic,
                    &[0, 1],
                    false,
                )
                .expect("generic plan should fall back"),
            None
        );
    }

    #[test]
    fn packed_hardcalls_stats_for_all_matches_expanded_stats_with_partial_word() {
        let mut packed = PackedHardcalls::default();
        packed.resize(35);
        for sample_index in 0..packed.sample_ct() {
            packed.set(sample_index, (sample_index % 4) as u8);
        }

        let source_indices = (0..packed.sample_ct()).collect::<Vec<_>>();
        let mut values = Vec::new();
        let mut missing = Vec::new();
        packed.expand_selected(&source_indices, &mut values, &mut missing);

        let expected = genoio_core::compute_variant_stats(&values, &missing)
            .expect("expanded stats should compute");
        let actual = packed
            .stats_for_selection(&source_indices, true)
            .expect("packed stats should compute");

        assert_eq!(actual, expected);
    }
}
