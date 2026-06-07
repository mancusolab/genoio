// pattern: Imperative Shell

use std::fs;
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, compute_dosage_variant_stats,
    flip_values_to_minor_allele, reject_sparse_missing_values, select_samples_source_order,
    DenseGenotypeMatrix, DenseSampleSelection, GenoioError, MetadataOutput, PartialFilterDecision,
    SampleRecord, SourceCapabilities, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};

use crate::error::Result;
use crate::hardcall::HardcallBatch as PackedVariantBatch;
#[cfg(test)]
use crate::hardcall::{PackedHardcalls as PackedGenotypes, HARDCALL_BATCH_SIZE};
use crate::matrix::{
    finish_dense_matrix, finish_variant_major_dense_matrix, DenseMatrixParts,
    VariantMajorDenseParts,
};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

mod metadata;
mod pgen;
use metadata::{parse_psam, parse_pvar, parse_pvar_source_window, PvarRecordReader};
use pgen::{
    decode_plink2_haplotype_dosage_aux, decode_plink2_haplotype_hardcall_aux, open_pgen_payload,
    read_fixed_width_variant_packed_sequential, read_plink2_variant_dosage,
    read_plink2_variant_haplotype_dosage_track, read_plink2_variant_haplotype_main_track,
    read_plink2_variant_packed, read_plink2_variant_values, read_supported_pgen_header,
    read_supported_pgen_header_prefix, seek_fixed_width_variant_record, validate_plink2_dimensions,
    validate_plink2_sample_count, PgenDecoderState, PgenHaplotypeDecodeState, PgenHeader,
    PgenLayout,
};

#[cfg(test)]
const PGEN_PACKED_TRANSPOSE_BATCH: usize = HARDCALL_BATCH_SIZE;

#[cfg(test)]
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

fn flush_packed_variant_batch(
    batch: &mut PackedVariantBatch,
    source_indices: &[usize],
    batch_start: &mut usize,
    n_variants: usize,
    values: &mut [f32],
    missing_mask: &mut [bool],
) {
    if batch.is_empty() {
        return;
    }
    batch.expand_into_sample_major(
        source_indices,
        *batch_start,
        n_variants,
        values,
        missing_mask,
    );
    *batch_start += batch.len();
    batch.clear();
}

fn expand_selected_samples_to_haplotypes(selection: &DenseSampleSelection) -> Vec<SampleRecord> {
    let mut haplotype_samples = Vec::with_capacity(selection.samples.len() * 2);
    for (sample, &source_index) in selection.samples.iter().zip(&selection.source_indices) {
        for haplotype_index in 0..2 {
            let mut haplotype_sample = sample.clone();
            haplotype_sample.source_sample_index = Some(source_index);
            haplotype_sample.haplotype_index = Some(haplotype_index);
            haplotype_samples.push(haplotype_sample);
        }
    }
    haplotype_samples
}

/// Read PLINK2 sample and variant metadata without returning genotypes.
pub fn read_plink2_metadata(pgen: &Path, pvar: &Path, psam: &Path) -> Result<MetadataOutput> {
    let header = read_supported_pgen_header(pgen)?;
    let samples = parse_psam(psam)?;
    let variants = parse_pvar(pvar)?;
    validate_plink2_dimensions(pgen, &header, samples.len(), variants.len())?;

    Ok(MetadataOutput {
        samples,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}

/// Read all retained PLINK2 hard-call genotypes as a dense matrix.
pub fn read_plink2_dense(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_plink2_dense_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        None,
        false,
    )
}

/// Read retained PLINK2 hard calls as a dense matrix over an optional block window.
pub fn read_plink2_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    // With no variant filter, retained order is identical to source order.
    // This lets block reads avoid full PVAR parsing and full variable-width
    // header validation. Filtered reads use the slower complete path because
    // retained-window membership depends on evaluating earlier variants.
    if let (None, Some(window)) = (variant_filter, variant_window) {
        if matrix_only {
            return read_plink2_dense_matrix_only_source_window(
                pgen,
                psam,
                requested_samples,
                window,
            );
        }
        return read_plink2_dense_source_window(pgen, pvar, psam, requested_samples, window);
    }

    let header = read_supported_pgen_header(pgen)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        fs::metadata(pvar).map_err(|source| GenoioError::Io {
            path: pvar.to_path_buf(),
            source,
        })?;
        diagnostics.retained_variants = 0;
        return finish_dense_matrix(
            DenseMatrixParts {
                n_samples: selection.samples.len(),
                n_variants: 0,
                values: Vec::new(),
                missing_mask: Vec::new(),
                samples: selection.samples,
                variants: Vec::new(),
                diagnostics,
            },
            matrix_only,
        );
    }
    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let output_variant_capacity = variant_window.map_or(header.variant_ct, |window| {
        window.len.min(header.variant_ct)
    });
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut variant_major_missing =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let mut output_variant_count = 0_usize;
    let requires_sequential_decode = matches!(header.layout, PgenLayout::VariableWidth);
    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let mut decoded_packed = false;
        if requires_sequential_decode {
            read_plink2_variant_packed(
                pgen,
                &mut file,
                &header,
                variant_index,
                &mut decoder_state,
            )?;
            decoded_packed = true;
        }
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => {
                stopped_after_window = true;
                break;
            }
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        if !decoded_packed {
            read_plink2_variant_packed(
                pgen,
                &mut file,
                &header,
                variant_index,
                &mut decoder_state,
            )?;
        }
        let stats = if needs_genotype_decision {
            Some(
                decoder_state
                    .packed
                    .stats_for_selected(&selection.source_indices)?,
            )
        } else {
            None
        };
        if needs_genotype_decision {
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        decoder_state.packed.expand_selected(
            &selection.source_indices,
            &mut decoder_state.values,
            &mut decoder_state.missing,
        );
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&decoder_state.values);
        variant_major_missing.extend_from_slice(&decoder_state.missing);
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = selection.samples.len();
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

/// Read retained PLINK2 unphased biallelic dosages as a dense matrix.
pub fn read_plink2_dosage_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let header = read_supported_pgen_header(pgen)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        fs::metadata(pvar).map_err(|source| GenoioError::Io {
            path: pvar.to_path_buf(),
            source,
        })?;
        diagnostics.retained_variants = 0;
        return finish_dense_matrix(
            DenseMatrixParts {
                n_samples: selection.samples.len(),
                n_variants: 0,
                values: Vec::new(),
                missing_mask: Vec::new(),
                samples: selection.samples,
                variants: Vec::new(),
                diagnostics,
            },
            matrix_only,
        );
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let output_variant_capacity = variant_window.map_or(header.variant_ct, |window| {
        window.len.min(header.variant_ct)
    });
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut variant_major_missing =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let mut output_variant_count = 0_usize;

    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => {
                stopped_after_window = true;
                break;
            }
        }

        read_plink2_variant_dosage(
            pgen,
            &mut file,
            &header,
            variant_index,
            &selection.source_indices,
            &mut decoder_state,
        )?;
        if matches!(partial_decision, PartialFilterDecision::NeedGenotypes) {
            let stats =
                compute_dosage_variant_stats(&decoder_state.values, &decoder_state.missing)?;
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
            attach_variant_stats(&mut variant, stats);
        }
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&decoder_state.values);
        variant_major_missing.extend_from_slice(&decoder_state.missing);
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = selection.samples.len();
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

/// Read retained explicit-phased PLINK2 hard calls as dense haplotype rows.
pub fn read_plink2_haplotypes_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let header = read_supported_pgen_header(pgen)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics.clone();
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        fs::metadata(pvar).map_err(|source| GenoioError::Io {
            path: pvar.to_path_buf(),
            source,
        })?;
        diagnostics.retained_variants = 0;
        return finish_dense_matrix(
            DenseMatrixParts {
                n_samples: haplotype_samples.len(),
                n_variants: 0,
                values: Vec::new(),
                missing_mask: Vec::new(),
                samples: haplotype_samples,
                variants: Vec::new(),
                diagnostics,
            },
            matrix_only,
        );
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let mut haplotype_state = PgenHaplotypeDecodeState::default();
    let output_variant_capacity = variant_window.map_or(header.variant_ct, |window| {
        window.len.min(header.variant_ct)
    });
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut variant_major_missing = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let mut output_variant_count = 0_usize;

    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let main_track_cursor = read_plink2_variant_haplotype_main_track(
            pgen,
            &mut file,
            &header,
            variant_index,
            &mut decoder_state,
        )?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => {
                stopped_after_window = true;
                break;
            }
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        if needs_genotype_decision {
            let stats = decoder_state
                .packed
                .stats_for_selected(&selection.source_indices)?;
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
            attach_variant_stats(&mut variant, stats);
        }
        decode_plink2_haplotype_hardcall_aux(
            pgen,
            &header,
            variant_index,
            main_track_cursor,
            &selection.source_indices,
            &decoder_state,
            &mut haplotype_state,
        )?;
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&haplotype_state.selected_haplotype_values);
        variant_major_missing.extend_from_slice(&haplotype_state.selected_haplotype_missing);
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = n_haplotypes;
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: haplotype_samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

/// Read retained explicit-phased PLINK2 dosages as dense haplotype rows.
pub fn read_plink2_haplotypes_dosage_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let header = read_supported_pgen_header(pgen)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics.clone();
    let haplotype_samples = expand_selected_samples_to_haplotypes(&selection);
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        fs::metadata(pvar).map_err(|source| GenoioError::Io {
            path: pvar.to_path_buf(),
            source,
        })?;
        diagnostics.retained_variants = 0;
        return finish_dense_matrix(
            DenseMatrixParts {
                n_samples: haplotype_samples.len(),
                n_variants: 0,
                values: Vec::new(),
                missing_mask: Vec::new(),
                samples: haplotype_samples,
                variants: Vec::new(),
                diagnostics,
            },
            matrix_only,
        );
    }

    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());
    let mut haplotype_state = PgenHaplotypeDecodeState::default();
    let output_variant_capacity = variant_window.map_or(header.variant_ct, |window| {
        window.len.min(header.variant_ct)
    });
    let n_haplotypes = selection.samples.len() * 2;
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut variant_major_missing = Vec::with_capacity(n_haplotypes * output_variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let mut output_variant_count = 0_usize;

    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let main_track_cursor = read_plink2_variant_haplotype_dosage_track(
            pgen,
            &mut file,
            &header,
            variant_index,
            &mut decoder_state,
        )?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => {
                stopped_after_window = true;
                break;
            }
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        if needs_genotype_decision {
            decode_plink2_haplotype_dosage_aux(
                pgen,
                &header,
                variant_index,
                main_track_cursor,
                &selection.source_indices,
                &decoder_state,
                &mut haplotype_state,
            )?;
            let stats = compute_dosage_variant_stats(
                &haplotype_state.selected_collapsed_values,
                &haplotype_state.selected_collapsed_missing,
            )?;
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
            attach_variant_stats(&mut variant, stats);
        } else {
            decode_plink2_haplotype_dosage_aux(
                pgen,
                &header,
                variant_index,
                main_track_cursor,
                &selection.source_indices,
                &decoder_state,
                &mut haplotype_state,
            )?;
        }
        if !matrix_only {
            variants.push(variant);
        }
        variant_major_values.extend_from_slice(&haplotype_state.selected_haplotype_values);
        variant_major_missing.extend_from_slice(&haplotype_state.selected_haplotype_missing);
        output_variant_count += 1;
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = n_haplotypes;
    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: haplotype_samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

/// Read all retained explicit-phased PLINK2 hard calls as sparse haplotype CSC.
pub fn read_plink2_haplotypes_sparse(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_plink2_haplotypes_sparse_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        None,
    )
}

/// Read retained explicit-phased PLINK2 hard calls as sparse haplotype CSC.
pub fn read_plink2_haplotypes_sparse_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    let dense = read_plink2_haplotypes_dense_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        false,
    )?;
    dense_haplotype_hardcalls_to_sparse(dense)
}

fn dense_haplotype_hardcalls_to_sparse(dense: DenseGenotypeMatrix) -> Result<SparseGenotypeMatrix> {
    reject_sparse_missing_values(&dense.missing_mask)?;
    let mut indptr = Vec::with_capacity(dense.n_variants + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();

    for col in 0..dense.n_variants {
        for row in 0..dense.n_samples {
            let value = dense.values[row * dense.n_variants + col];
            if value != 0.0 {
                indices.push(row);
                data.push(value);
            }
        }
        indptr.push(indices.len());
    }

    SparseGenotypeMatrix::new(
        dense.n_samples,
        dense.n_variants,
        indptr,
        indices,
        data,
        dense.samples,
        dense.variants,
        dense.diagnostics,
    )
}

/// Read all retained PLINK2 hard-call genotypes as sparse CSC.
pub fn read_plink2_sparse(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_plink2_sparse_windowed(pgen, pvar, psam, requested_samples, variant_filter, None)
}

/// Read retained PLINK2 hard calls as sparse CSC over an optional block window.
pub fn read_plink2_sparse_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    // See the dense fast path: unfiltered windows can be interpreted directly
    // in source coordinates, but filtered windows cannot.
    if let (None, Some(window)) = (variant_filter, variant_window) {
        return read_plink2_sparse_source_window(pgen, pvar, psam, requested_samples, window);
    }

    let header = read_supported_pgen_header(pgen)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        fs::metadata(pvar).map_err(|source| GenoioError::Io {
            path: pvar.to_path_buf(),
            source,
        })?;
        diagnostics.retained_variants = 0;
        return SparseGenotypeMatrix::new(
            selection.samples.len(),
            0,
            vec![0],
            Vec::new(),
            Vec::new(),
            selection.samples,
            Vec::new(),
            diagnostics,
        );
    }
    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let output_variant_capacity = variant_window.map_or(header.variant_ct, |window| {
        window.len.min(header.variant_ct)
    });
    let mut indptr = Vec::with_capacity(output_variant_capacity + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let requires_sequential_decode = matches!(header.layout, PgenLayout::VariableWidth);
    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        let mut decoded_packed = false;
        if requires_sequential_decode {
            read_plink2_variant_packed(
                pgen,
                &mut file,
                &header,
                variant_index,
                &mut decoder_state,
            )?;
            decoded_packed = true;
        }
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => {
                stopped_after_window = true;
                break;
            }
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        if !decoded_packed {
            read_plink2_variant_packed(
                pgen,
                &mut file,
                &header,
                variant_index,
                &mut decoder_state,
            )?;
        }
        let stats = if needs_genotype_decision {
            Some(
                decoder_state
                    .packed
                    .stats_for_selected(&selection.source_indices)?,
            )
        } else {
            None
        };
        if needs_genotype_decision {
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => {
                    stopped_after_window = true;
                    break;
                }
            }
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        decoder_state.packed.expand_selected(
            &selection.source_indices,
            &mut decoder_state.values,
            &mut decoder_state.missing,
        );
        reject_sparse_missing_values(&decoder_state.missing)?;
        flip_values_to_minor_allele(&mut decoder_state.values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &decoder_state.values);
        variants.push(variant);
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_variants = variants.len();
    diagnostics.retained_variants = n_variants;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        selection.samples,
        variants,
        diagnostics,
    )
}

fn read_plink2_dense_matrix_only_source_window(
    pgen: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    window: VariantWindow,
) -> Result<DenseGenotypeMatrix> {
    let decode_variant_ct = window
        .start
        .checked_add(window.len)
        .ok_or_else(|| GenoioError::invalid_source(pgen, "variant window end is out of range"))?;
    let header = read_supported_pgen_header_prefix(pgen, decode_variant_ct)?;
    if window.start > header.variant_ct {
        return Err(GenoioError::invalid_source(
            pgen,
            format!(
                "variant window start {} exceeds pgen variant count {}",
                window.start, header.variant_ct
            ),
        ));
    }
    let n_variants = window.len.min(header.variant_ct - window.start);
    if let Some(requested_samples) = requested_samples {
        let all_samples = parse_psam(psam)?;
        validate_plink2_sample_count(pgen, &header, all_samples.len())?;
        let selection = select_samples_source_order(&all_samples, Some(requested_samples), pgen)?;
        return decode_plink2_dense_matrix_only_source_window(
            pgen,
            &header,
            window,
            n_variants,
            &selection.source_indices,
            selection.samples.len(),
            selection.diagnostics,
        );
    }

    let n_samples = header.sample_ct;
    let source_indices = (0..n_samples).collect::<Vec<_>>();
    let diagnostics = genoio_core::DenseDiagnostics {
        requested_samples: n_samples,
        retained_samples: n_samples,
        candidate_variants: n_variants,
        retained_variants: n_variants,
        ..genoio_core::DenseDiagnostics::default()
    };

    decode_plink2_dense_matrix_only_source_window(
        pgen,
        &header,
        window,
        n_variants,
        &source_indices,
        n_samples,
        diagnostics,
    )
}

fn decode_plink2_dense_matrix_only_source_window(
    pgen: &Path,
    header: &PgenHeader,
    window: VariantWindow,
    n_variants: usize,
    source_indices: &[usize],
    n_samples: usize,
    mut diagnostics: genoio_core::DenseDiagnostics,
) -> Result<DenseGenotypeMatrix> {
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, n_samples);
    let mut values = vec![0.0; n_samples * n_variants];
    let mut missing_mask = vec![false; n_samples * n_variants];
    let mut packed_batch = PackedVariantBatch::new(header.sample_ct);
    let mut batch_start = 0_usize;

    // Unfiltered source windows know their final retained width up front, so
    // construct the public sample-major buffers directly.
    match header.layout {
        PgenLayout::FixedWidth
        | PgenLayout::FixedWidthDosage
        | PgenLayout::FixedWidthPhasedDosage => {
            if n_variants > 0 {
                seek_fixed_width_variant_record(pgen, &mut file, header, window.start)?;
            }
            for _ in 0..n_variants {
                read_fixed_width_variant_packed_sequential(
                    pgen,
                    &mut file,
                    header,
                    &mut decoder_state,
                )?;
                packed_batch.push(&decoder_state.packed);
                if packed_batch.is_full() {
                    flush_packed_variant_batch(
                        &mut packed_batch,
                        source_indices,
                        &mut batch_start,
                        n_variants,
                        &mut values,
                        &mut missing_mask,
                    );
                }
            }
            flush_packed_variant_batch(
                &mut packed_batch,
                source_indices,
                &mut batch_start,
                n_variants,
                &mut values,
                &mut missing_mask,
            );
        }
        PgenLayout::VariableWidth => {
            let prefix_end = window.start + n_variants;
            for source_variant_index in 0..prefix_end {
                read_plink2_variant_packed(
                    pgen,
                    &mut file,
                    header,
                    source_variant_index,
                    &mut decoder_state,
                )?;
                if source_variant_index >= window.start {
                    packed_batch.push(&decoder_state.packed);
                    if packed_batch.is_full() {
                        flush_packed_variant_batch(
                            &mut packed_batch,
                            source_indices,
                            &mut batch_start,
                            n_variants,
                            &mut values,
                            &mut missing_mask,
                        );
                    }
                }
            }
            flush_packed_variant_batch(
                &mut packed_batch,
                source_indices,
                &mut batch_start,
                n_variants,
                &mut values,
                &mut missing_mask,
            );
        }
    }

    diagnostics.candidate_variants = n_variants;
    diagnostics.retained_variants = n_variants;

    DenseGenotypeMatrix::new_matrix_only(n_samples, n_variants, values, missing_mask, diagnostics)
}

fn read_plink2_dense_source_window(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    window: VariantWindow,
) -> Result<DenseGenotypeMatrix> {
    let decode_variant_ct = window.start.saturating_add(window.len);
    let header = read_supported_pgen_header_prefix(pgen, decode_variant_ct)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    let window_variants = parse_pvar_source_window(pvar, window, header.variant_ct)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let n_variants = window_variants.len();
    let mut variants = Vec::new();
    let mut values = vec![0.0; n_samples * n_variants];
    let mut missing_mask = vec![false; n_samples * n_variants];
    let mut packed_batch = PackedVariantBatch::new(header.sample_ct);
    let mut batch_start = 0_usize;

    // This metadata-bearing source-window path uses the same packed batch
    // expansion as matrix-only windows while preserving metadata alignment.
    match header.layout {
        PgenLayout::FixedWidth
        | PgenLayout::FixedWidthDosage
        | PgenLayout::FixedWidthPhasedDosage => {
            if let Some((first_variant_index, _)) = window_variants.first() {
                seek_fixed_width_variant_record(pgen, &mut file, &header, *first_variant_index)?;
            }
            for (source_variant_index, variant) in window_variants {
                debug_assert!(source_variant_index < header.variant_ct);
                read_fixed_width_variant_packed_sequential(
                    pgen,
                    &mut file,
                    &header,
                    &mut decoder_state,
                )?;
                variants.push(variant);
                packed_batch.push(&decoder_state.packed);
                if packed_batch.is_full() {
                    flush_packed_variant_batch(
                        &mut packed_batch,
                        &selection.source_indices,
                        &mut batch_start,
                        n_variants,
                        &mut values,
                        &mut missing_mask,
                    );
                }
            }
            flush_packed_variant_batch(
                &mut packed_batch,
                &selection.source_indices,
                &mut batch_start,
                n_variants,
                &mut values,
                &mut missing_mask,
            );
        }
        PgenLayout::VariableWidth => {
            let mut window_iter = window_variants.into_iter().peekable();
            let prefix_end = header.record_types.len();
            // Variable-width PGEN can use LD-compressed records that depend on
            // earlier non-LD records. Decode the prefix through the requested
            // window to maintain state, but batch only requested variants.
            for variant_index in 0..prefix_end {
                read_plink2_variant_packed(
                    pgen,
                    &mut file,
                    &header,
                    variant_index,
                    &mut decoder_state,
                )?;
                if window_iter
                    .peek()
                    .is_some_and(|(source_index, _)| *source_index == variant_index)
                {
                    if let Some((_, variant)) = window_iter.next() {
                        variants.push(variant);
                        packed_batch.push(&decoder_state.packed);
                        if packed_batch.is_full() {
                            flush_packed_variant_batch(
                                &mut packed_batch,
                                &selection.source_indices,
                                &mut batch_start,
                                n_variants,
                                &mut values,
                                &mut missing_mask,
                            );
                        }
                    }
                }
            }
            flush_packed_variant_batch(
                &mut packed_batch,
                &selection.source_indices,
                &mut batch_start,
                n_variants,
                &mut values,
                &mut missing_mask,
            );
        }
    }

    debug_assert_eq!(variants.len(), n_variants);
    diagnostics.candidate_variants = n_variants;
    diagnostics.retained_variants = n_variants;

    DenseGenotypeMatrix::new(
        n_samples,
        n_variants,
        values,
        missing_mask,
        selection.samples,
        variants,
        diagnostics,
    )
}

fn read_plink2_sparse_source_window(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    window: VariantWindow,
) -> Result<SparseGenotypeMatrix> {
    let decode_variant_ct = window.start.saturating_add(window.len);
    let header = read_supported_pgen_header_prefix(pgen, decode_variant_ct)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    let window_variants = parse_pvar_source_window(pvar, window, header.variant_ct)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();

    match header.layout {
        PgenLayout::FixedWidth
        | PgenLayout::FixedWidthDosage
        | PgenLayout::FixedWidthPhasedDosage => {
            // Fixed-width records can be decoded by direct source index.
            for (variant_index, mut variant) in window_variants {
                read_plink2_variant_values(
                    pgen,
                    &mut file,
                    &header,
                    variant_index,
                    &selection.source_indices,
                    &mut decoder_state,
                )?;
                reject_sparse_missing_values(&decoder_state.missing)?;
                flip_values_to_minor_allele(&mut decoder_state.values, &mut variant);
                append_sparse_column(&mut indptr, &mut indices, &mut data, &decoder_state.values);
                variants.push(variant);
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
                    .is_some_and(|(source_index, _)| *source_index == variant_index)
                {
                    if let Some((_, mut variant)) = window_iter.next() {
                        reject_sparse_missing_values(&decoder_state.missing)?;
                        flip_values_to_minor_allele(&mut decoder_state.values, &mut variant);
                        append_sparse_column(
                            &mut indptr,
                            &mut indices,
                            &mut data,
                            &decoder_state.values,
                        );
                        variants.push(variant);
                    }
                }
            }
        }
    }

    let n_variants = variants.len();
    diagnostics.candidate_variants = n_variants;
    diagnostics.retained_variants = n_variants;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        selection.samples,
        variants,
        diagnostics,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let mut missing = vec![true];
        packed.expand_selected(&[3, 1, 34, 0], &mut values, &mut missing);

        assert_eq!(values, vec![0.0, 1.0, 2.0, 0.0]);
        assert_eq!(missing, vec![true, false, false, false]);
    }

    #[test]
    fn packed_variant_batch_expands_like_variant_at_a_time() {
        let sample_ct = 5;
        let source_indices = (0..sample_ct).collect::<Vec<_>>();
        let n_variants = PGEN_PACKED_TRANSPOSE_BATCH + 3;
        let mut packed_variants = Vec::with_capacity(n_variants);
        let mut expected_values = vec![0.0; sample_ct * n_variants];
        let mut expected_missing = vec![false; sample_ct * n_variants];
        let mut scratch_values = Vec::new();
        let mut scratch_missing = Vec::new();

        for variant_index in 0..n_variants {
            let mut packed = PackedGenotypes::default();
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

        let mut batch = PackedVariantBatch::new(sample_ct);
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
}
