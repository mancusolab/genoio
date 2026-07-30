// pattern: Imperative Shell
//! Sparse hard-call PLINK2 read orchestration.
//!
//! The sparse path decodes retained hard calls, applies the same filters as the
//! dense path, flips common alternate columns when needed, and emits CSC columns
//! without materializing a dense matrix.

use std::path::Path;

use genoio_core::{
    append_sparse_column, flip_values_to_minor_allele, reject_sparse_missing,
    SampleMetadataBuffers, SparseGenotypeMatrix, VariantFilter, VariantMetadataBuffers,
    VariantRecord, VariantWindow,
};

use crate::error::Result;

use super::metadata::parse_pvar_source_window;
use super::pgen::{open_pgen_payload, read_plink2_variant_values, PgenDecoderState, PgenLayout};
use super::source::{require_pvar, Plink2ReadContext};

/// Source-window row that may omit metadata for matrix-only reads.
struct SourceWindowVariant {
    source_index: usize,
    variant: Option<VariantRecord>,
}

#[inline]
fn append_decoded_sparse_column(
    decoder_state: &mut PgenDecoderState,
    variant: &mut VariantRecord,
    indptr: &mut Vec<i32>,
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
) -> Result<()> {
    reject_sparse_missing(!decoder_state.missing_indices.is_empty())?;
    flip_values_to_minor_allele(&mut decoder_state.values, variant);
    append_sparse_column(indptr, indices, data, &decoder_state.values)
}

#[inline]
fn append_decoded_sparse_column_without_variant_metadata(
    decoder_state: &mut PgenDecoderState,
    indptr: &mut Vec<i32>,
    indices: &mut Vec<i32>,
    data: &mut Vec<f32>,
) -> Result<()> {
    reject_sparse_missing(!decoder_state.missing_indices.is_empty())?;
    // Matrix-only source windows must preserve sparse minor-allele orientation
    // without mutating metadata that the caller did not request.
    flip_values_to_minor_allele_without_metadata(&mut decoder_state.values);
    append_sparse_column(indptr, indices, data, &decoder_state.values)
}

fn flip_values_to_minor_allele_without_metadata(values: &mut [f32]) {
    let a1_count = values.iter().sum::<f32>();
    let a0_count = 2.0 * values.len() as f32 - a1_count;
    if a1_count <= a0_count {
        return;
    }
    for value in values {
        *value = 2.0 - *value;
    }
}

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors sparse read options plus metadata return choices"
)]
pub fn read_plink2_sparse_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    // See the dense fast path: unfiltered windows can be interpreted directly
    // in source coordinates, but filtered windows cannot.
    if let (None, Some(window)) = (variant_filter, variant_window) {
        return read_plink2_sparse_source_window(
            pgen,
            pvar,
            psam,
            requested_samples,
            window,
            return_samples,
            return_variants,
        );
    }

    let output = super::session::read_windowed(
        pgen,
        pvar,
        psam,
        crate::blocks::BlockReadOptions {
            matrix_kind: crate::blocks::MatrixKind::Genotype,
            sparse: true,
            requested_samples: requested_samples.map(<[String]>::to_vec),
            variant_filter: variant_filter.cloned(),
            dosage_source: crate::blocks::DosageSource::Hardcall,
            missing_policy: genoio_core::DenseMissingPolicy::Raise,
            return_samples,
            return_variants,
        },
        variant_window,
    )?;
    let crate::blocks::BlockOutput::Sparse(output) = output else {
        return Err(genoio_core::GenoioError::internal_contract(
            "PLINK2 sparse hardcall session returned dense output",
        ));
    };
    Ok(output)
}

fn read_plink2_sparse_source_window(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    window: VariantWindow,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    let decode_variant_ct = window.start.saturating_add(window.len);
    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected: _,
    } = Plink2ReadContext::new_prefix(pgen, psam, requested_samples, decode_variant_ct)?;
    let mut diagnostics = selection.diagnostics;
    let n_variants = header
        .variant_ct
        .saturating_sub(window.start)
        .min(window.len);
    let window_variants: Vec<SourceWindowVariant> = if return_variants {
        parse_pvar_source_window(pvar, window, header.variant_ct)?
            .into_iter()
            .map(|(source_index, variant)| SourceWindowVariant {
                source_index,
                variant: Some(variant),
            })
            .collect()
    } else {
        require_pvar(pvar)?;
        (window.start..window.start + n_variants)
            .map(|source_index| SourceWindowVariant {
                source_index,
                variant: None,
            })
            .collect()
    };
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let output_variant_capacity = n_variants;
    let mut indptr = Vec::with_capacity(output_variant_capacity + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants =
        return_variants.then(|| VariantMetadataBuffers::with_capacity(output_variant_capacity));

    match header.layout {
        PgenLayout::FixedWidth
        | PgenLayout::FixedWidthDosage
        | PgenLayout::FixedWidthPhasedDosage => {
            // Fixed-width records can be decoded by direct source index.
            for source_variant in window_variants {
                read_plink2_variant_values(
                    pgen,
                    &mut file,
                    &header,
                    source_variant.source_index,
                    &selection.source_indices,
                    &mut decoder_state,
                )?;
                if let Some(mut variant) = source_variant.variant {
                    append_decoded_sparse_column(
                        &mut decoder_state,
                        &mut variant,
                        &mut indptr,
                        &mut indices,
                        &mut data,
                    )?;
                    if let Some(variants) = variants.as_mut() {
                        variants.push_record(&variant)?;
                    }
                } else {
                    append_decoded_sparse_column_without_variant_metadata(
                        &mut decoder_state,
                        &mut indptr,
                        &mut indices,
                        &mut data,
                    )?;
                }
            }
        }
        PgenLayout::VariableWidth => {
            let mut window_iter = window_variants.into_iter().peekable();
            let prefix_end = header.record_types.len();
            // Preserve LD state exactly as dense reads do, then append only
            // requested variants to sparse columns.
            for variant_index in 0..prefix_end {
                read_plink2_variant_values(
                    pgen,
                    &mut file,
                    &header,
                    variant_index,
                    &selection.source_indices,
                    &mut decoder_state,
                )?;
                if window_iter
                    .peek()
                    .is_some_and(|source_variant| source_variant.source_index == variant_index)
                {
                    if let Some(source_variant) = window_iter.next() {
                        if let Some(mut variant) = source_variant.variant {
                            append_decoded_sparse_column(
                                &mut decoder_state,
                                &mut variant,
                                &mut indptr,
                                &mut indices,
                                &mut data,
                            )?;
                            if let Some(variants) = variants.as_mut() {
                                variants.push_record(&variant)?;
                            }
                        } else {
                            append_decoded_sparse_column_without_variant_metadata(
                                &mut decoder_state,
                                &mut indptr,
                                &mut indices,
                                &mut data,
                            )?;
                        }
                    }
                }
            }
        }
    }

    let n_variants = indptr.len().saturating_sub(1);
    diagnostics.candidate_variants = n_variants;
    diagnostics.retained_variants = n_variants;
    let samples =
        SampleMetadataBuffers::optional_from_records(&selection.samples, return_samples, false)?;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        samples,
        variants,
        diagnostics,
    )
}
