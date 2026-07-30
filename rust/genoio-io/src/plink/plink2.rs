// pattern: Imperative Shell
//! PLINK2 reader facade and shared filter helpers.
//!
//! Public entry points are re-exported from dense, dosage, haplotype, and sparse
//! modules. Shared helpers here keep genotype-stat filter behavior consistent
//! across those read shapes.

use std::path::Path;

use genoio_core::{
    GenoioError, MetadataOutput, SampleMetadataBuffers, SourceCapabilities, VariantFilter,
};

use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
#[cfg(test)]
use crate::hardcall::HardcallBatch as PackedVariantBatch;
#[cfg(test)]
use crate::hardcall::PackedHardcalls as PackedGenotypes;
#[cfg(test)]
use crate::hardcall::HARDCALL_BATCH_SIZE;

mod dense;
mod dosage;
mod haplotype;
mod metadata;
mod pgen;
mod session;
mod source;
mod sparse;

pub(crate) use session::Plink2BlockSession;

#[doc(inline)]
pub use dense::read_plink2_dense_windowed;
#[doc(inline)]
pub use dosage::read_plink2_dosage_dense_windowed;
#[doc(inline)]
pub use haplotype::{
    read_plink2_haplotypes_dense_windowed, read_plink2_haplotypes_dosage_dense_windowed,
    read_plink2_haplotypes_sparse_windowed,
};
#[doc(inline)]
pub use sparse::read_plink2_sparse_windowed;

use metadata::{parse_psam, parse_pvar_metadata};
use pgen::{read_supported_pgen_header, validate_plink2_dimensions};

#[cfg(test)]
const PGEN_PACKED_TRANSPOSE_BATCH: usize = HARDCALL_BATCH_SIZE;

pub(super) fn require_genotype_decision_filter(
    variant_filter: Option<&VariantFilter>,
) -> Result<&VariantFilter> {
    variant_filter.ok_or_else(|| {
        GenoioError::internal_contract("genotype decision requires a variant filter")
    })
}

/// Read PLINK2 metadata with variant metadata staged as columnar buffers.
pub fn read_plink2_metadata(pgen: &Path, pvar: &Path, psam: &Path) -> Result<MetadataOutput> {
    let header = read_supported_pgen_header(pgen)?;
    let samples = parse_psam(psam)?;
    let variants = parse_pvar_metadata(pvar)?;
    validate_plink2_dimensions(pgen, &header, samples.len(), variants.len())?;

    Ok(MetadataOutput {
        samples: SampleMetadataBuffers::from_records(&samples, false)?,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{apply_dense_missing_policy_to_variant, write_sample_major_variant_slot};

    fn stats_from_expanded_hardcalls(
        values: &[f32],
        missing_indices: &[usize],
    ) -> genoio_core::VariantStats {
        let mut hom_ref = 0_u64;
        let mut het = 0_u64;
        let mut hom_alt = 0_u64;
        let mut missing_cursor = 0_usize;
        for (index, value) in values.iter().enumerate() {
            if missing_indices
                .get(missing_cursor)
                .is_some_and(|&missing_index| missing_index == index)
            {
                missing_cursor += 1;
                continue;
            }
            match *value as u8 {
                0 => hom_ref += 1,
                1 => het += 1,
                2 => hom_alt += 1,
                _ => panic!("test hardcall value must be in {{0, 1, 2}}"),
            }
        }
        genoio_core::variant_stats_from_counts(
            hom_ref,
            het,
            hom_alt,
            u64::try_from(missing_indices.len()).expect("test missing count fits in u64"),
        )
        .expect("test hardcall stats should compute")
    }

    fn assert_values_with_nan(actual: &[f32], expected: &[f32]) {
        assert_eq!(actual.len(), expected.len());
        for (&actual, &expected) in actual.iter().zip(expected) {
            if expected.is_nan() {
                assert!(actual.is_nan(), "expected NaN, observed {actual}");
            } else {
                assert_eq!(actual, expected);
            }
        }
    }

    #[test]
    fn packed_genotypes_round_trip_and_expand_selected() {
        let mut packed = PackedGenotypes::default();
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
        let mut missing_indices = vec![99];
        packed.expand_selected(&[3, 1, 34, 0], &mut values, &mut missing_indices);

        assert_eq!(values, vec![0.0, 1.0, 2.0, 0.0]);
        assert_eq!(missing_indices, vec![0]);
    }

    #[test]
    fn packed_variant_batch_expands_like_variant_at_a_time() {
        let sample_ct = 5;
        let source_indices = (0..sample_ct).collect::<Vec<_>>();
        let n_variants = PGEN_PACKED_TRANSPOSE_BATCH + 3;
        let mut packed_variants = Vec::with_capacity(n_variants);
        let mut expected_values = vec![0.0; sample_ct * n_variants];
        let mut scratch_values = Vec::new();
        let mut scratch_missing_indices = Vec::new();

        for variant_index in 0..n_variants {
            let mut packed = PackedGenotypes::default();
            packed.resize(sample_ct);
            for sample_index in 0..sample_ct {
                packed.set(sample_index, ((variant_index + sample_index) % 4) as u8);
            }
            packed.expand_selected(
                &source_indices,
                &mut scratch_values,
                &mut scratch_missing_indices,
            );
            apply_dense_missing_policy_to_variant(
                &mut scratch_values,
                &scratch_missing_indices,
                genoio_core::DenseMissingPolicy::Nan,
            )
            .expect("test missing policy should apply");
            write_sample_major_variant_slot(
                &mut expected_values,
                sample_ct,
                n_variants,
                variant_index,
                &scratch_values,
            )
            .expect("expected dense slot should write");
            packed_variants.push(packed);
        }

        let mut batch = PackedVariantBatch::new(sample_ct);
        let mut actual_values = vec![0.0; sample_ct * n_variants];
        let mut variant_values = Vec::with_capacity(sample_ct);
        let mut missing_indices = Vec::new();
        let mut batch_start = 0;
        for packed in &packed_variants {
            batch.push(packed);
            if batch.is_full() {
                crate::hardcall::flush_hardcall_batch_into_sample_major(
                    &mut batch,
                    &source_indices,
                    &mut batch_start,
                    n_variants,
                    &mut actual_values,
                    genoio_core::DenseMissingPolicy::Nan,
                    &mut variant_values,
                    &mut missing_indices,
                )
                .expect("batch flush should succeed");
            }
        }
        crate::hardcall::flush_hardcall_batch_into_sample_major(
            &mut batch,
            &source_indices,
            &mut batch_start,
            n_variants,
            &mut actual_values,
            genoio_core::DenseMissingPolicy::Nan,
            &mut variant_values,
            &mut missing_indices,
        )
        .expect("batch flush should succeed");

        assert_values_with_nan(&actual_values, &expected_values);
    }

    #[test]
    fn packed_genotypes_copy_and_invert_0_2() {
        let mut source = PackedGenotypes::default();
        source.resize(5);
        source.clear_to(3);
        source.set(0, 0);
        source.set(1, 1);
        source.set(2, 2);

        let mut copy = PackedGenotypes::default();
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
    fn packed_genotypes_loads_pgen_payload_and_masks_unused_trailing_slots() {
        let mut packed = PackedGenotypes::default();
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
    fn packed_genotypes_stats_for_selected_matches_expanded_stats() {
        let mut packed = PackedGenotypes::default();
        packed.resize(8);
        for (sample_index, category) in [0, 1, 2, 3, 2, 0, 1, 3].into_iter().enumerate() {
            packed.set(sample_index, category);
        }

        for source_indices in [&[0, 1, 2, 3, 4, 5, 6, 7][..], &[7, 3][..], &[][..]] {
            let mut values = Vec::new();
            let mut missing_indices = Vec::new();
            packed.expand_selected(source_indices, &mut values, &mut missing_indices);

            let expected = stats_from_expanded_hardcalls(&values, &missing_indices);
            let actual = packed
                .stats_for_selected(source_indices)
                .expect("packed stats should compute");

            assert_eq!(actual, expected);
        }
    }
}
