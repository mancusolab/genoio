// pattern: Imperative Shell
//! Dense hard-call PLINK2 read orchestration.
//!
//! This path decodes PGEN hard calls, evaluates metadata and genotype-stat
//! filters, and writes retained variants into sample-major dense matrices.
//! Matrix-only windows can skip PVAR parsing when filter semantics allow it.

use std::path::Path;

use genoio_core::{
    attach_variant_stats, DenseGenotypeMatrix, DenseSampleSelection, GenoioError,
    GenotypeFilterPlan, PartialFilterDecision, VariantFilter, VariantWindow,
};

use crate::error::Result;
use crate::hardcall::{
    evaluate_packed_hardcall_filter, flush_hardcall_batch_into_sample_major,
    HardcallBatch as PackedVariantBatch,
};
use crate::matrix::{finish_dense_matrix, shrink_sample_major_width, DenseMatrixParts};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use super::metadata::{parse_pvar_source_window, PvarRecordReader};
use super::pgen::{
    open_pgen_payload, read_fixed_width_variant_packed_sequential, read_plink2_variant_packed,
    read_supported_pgen_header_prefix, seek_fixed_width_variant_record, PgenDecoderState,
    PgenHeader, PgenLayout,
};
use super::require_genotype_decision_filter;
use super::source::{
    empty_dense_for_samples, matrix_only_source_window_diagnostics, require_pvar,
    select_samples_for_header, variant_output_capacity, Plink2ReadContext,
};

fn can_skip_pvar_for_matrix_only_genotype_filter(
    matrix_only: bool,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Option<(&VariantFilter, VariantWindow)> {
    if !matrix_only {
        return None;
    }
    let filter = variant_filter?;
    let window = variant_window?;
    (filter.requires_genotype_stats() && filter.is_genotype_stats_only())
        .then_some((filter, window))
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

    if let Some((filter, window)) =
        can_skip_pvar_for_matrix_only_genotype_filter(matrix_only, variant_filter, variant_window)
    {
        // Matrix-only genotype-stat filters need PGEN/PSAM data only. Use
        // prefix PGEN headers so retained-block reads do not pay the full
        // variable-width header parse before decoding the first block.
        require_pvar(pvar)?;
        return read_plink2_dense_matrix_only_genotype_filter(
            pgen,
            psam,
            requested_samples,
            filter,
            window,
        );
    }

    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected,
    } = Plink2ReadContext::new(pgen, psam, requested_samples)?;
    let mut diagnostics = selection.diagnostics.clone();
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        require_pvar(pvar)?;
        return empty_dense_for_samples(selection.samples, diagnostics, matrix_only);
    }
    let mut pvar_reader = PvarRecordReader::new(pvar)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let output_variant_capacity = variant_output_capacity(&header, variant_window);
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let n_samples = selection.samples.len();
    let mut values = vec![0.0; n_samples * output_variant_capacity];
    let mut missing_mask = vec![false; n_samples * output_variant_capacity];
    let mut packed_batch = PackedVariantBatch::new(header.sample_ct);
    let mut batch_start = 0_usize;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut stopped_after_window = false;
    let mut output_variant_count = 0_usize;
    let genotype_filter_plan = variant_filter.map_or(
        GenotypeFilterPlan::Generic,
        VariantFilter::genotype_filter_plan,
    );
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
        let mut stats = None;
        if needs_genotype_decision {
            let filter = require_genotype_decision_filter(variant_filter)?;
            let (retain_variant, computed_stats) = evaluate_packed_hardcall_filter(
                &decoder_state.packed,
                &selection.source_indices,
                all_samples_selected,
                filter,
                genotype_filter_plan,
                Some(&variant),
                !matrix_only,
            )?;
            stats = computed_stats;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
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
        if !matrix_only {
            variants.push(variant);
        }
        packed_batch.push(&decoder_state.packed);
        output_variant_count += 1;
        if packed_batch.is_full() {
            flush_hardcall_batch_into_sample_major(
                &mut packed_batch,
                &selection.source_indices,
                &mut batch_start,
                output_variant_capacity,
                &mut values,
                &mut missing_mask,
            );
        }
        if retention.window_is_satisfied() {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }
    flush_hardcall_batch_into_sample_major(
        &mut packed_batch,
        &selection.source_indices,
        &mut batch_start,
        output_variant_capacity,
        &mut values,
        &mut missing_mask,
    );

    let n_variants = output_variant_count;
    shrink_sample_major_width(&mut values, n_samples, output_variant_capacity, n_variants);
    shrink_sample_major_width(
        &mut missing_mask,
        n_samples,
        output_variant_capacity,
        n_variants,
    );
    diagnostics.retained_variants = n_variants;
    finish_dense_matrix(
        DenseMatrixParts {
            n_samples,
            n_variants,
            values,
            missing_mask,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
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
        let selection = select_samples_for_header(pgen, psam, Some(requested_samples), &header)?;
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
    let diagnostics = matrix_only_source_window_diagnostics(n_samples, n_variants);

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

fn read_plink2_dense_matrix_only_genotype_filter(
    pgen: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    filter: &VariantFilter,
    window: VariantWindow,
) -> Result<DenseGenotypeMatrix> {
    let mut decode_variant_ct = initial_genotype_filter_prefix_len(pgen, window)?;
    loop {
        let Plink2ReadContext {
            header,
            selection,
            all_samples_selected,
        } = Plink2ReadContext::new_prefix(pgen, psam, requested_samples, decode_variant_ct)?;
        let decode_limit = genotype_filter_decode_limit(&header);
        let is_complete = decode_limit >= header.variant_ct;
        let read = decode_plink2_dense_matrix_only_genotype_filter(
            pgen,
            &header,
            selection,
            all_samples_selected,
            filter,
            window,
            decode_limit,
        )?;
        if read.window_satisfied || is_complete {
            return Ok(read.matrix);
        }
        decode_variant_ct = next_genotype_filter_prefix_len(pgen, decode_limit, header.variant_ct)?;
    }
}

struct MatrixOnlyGenotypeFilterRead {
    matrix: DenseGenotypeMatrix,
    window_satisfied: bool,
}

fn initial_genotype_filter_prefix_len(pgen: &Path, window: VariantWindow) -> Result<usize> {
    window
        .start
        .checked_add(window.len)
        .ok_or_else(|| GenoioError::invalid_source(pgen, "variant window end is out of range"))
}

fn next_genotype_filter_prefix_len(
    pgen: &Path,
    current: usize,
    variant_ct: usize,
) -> Result<usize> {
    if current >= variant_ct {
        return Ok(variant_ct);
    }
    let doubled = current
        .checked_mul(2)
        .ok_or_else(|| GenoioError::invalid_source(pgen, "variant window end is out of range"))?;
    Ok(doubled.max(current + 1).min(variant_ct))
}

fn genotype_filter_decode_limit(header: &PgenHeader) -> usize {
    match header.layout {
        PgenLayout::VariableWidth => header.record_types.len(),
        PgenLayout::FixedWidth
        | PgenLayout::FixedWidthDosage
        | PgenLayout::FixedWidthPhasedDosage => header.variant_ct,
    }
}

fn decode_plink2_dense_matrix_only_genotype_filter(
    pgen: &Path,
    header: &PgenHeader,
    selection: DenseSampleSelection,
    all_samples_selected: bool,
    filter: &VariantFilter,
    window: VariantWindow,
    decode_limit: usize,
) -> Result<MatrixOnlyGenotypeFilterRead> {
    let output_variant_capacity = window.len.min(header.variant_ct);
    let n_samples = selection.samples.len();
    if output_variant_capacity == 0 {
        let mut diagnostics = selection.diagnostics;
        diagnostics.retained_variants = 0;
        let matrix = DenseGenotypeMatrix::new_matrix_only(
            n_samples,
            0,
            Vec::new(),
            Vec::new(),
            diagnostics,
        )?;
        return Ok(MatrixOnlyGenotypeFilterRead {
            matrix,
            window_satisfied: true,
        });
    }

    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, n_samples);
    let mut values = vec![0.0; n_samples * output_variant_capacity];
    let mut missing_mask = vec![false; n_samples * output_variant_capacity];
    let mut packed_batch = PackedVariantBatch::new(header.sample_ct);
    let mut batch_start = 0_usize;
    let mut retention = RetainedVariantState::new(Some(window));
    let mut diagnostics = selection.diagnostics;
    let mut output_variant_count = 0_usize;
    let genotype_filter_plan = filter.genotype_filter_plan();

    for variant_index in 0..decode_limit {
        diagnostics.candidate_variants += 1;

        match header.layout {
            PgenLayout::FixedWidth
            | PgenLayout::FixedWidthDosage
            | PgenLayout::FixedWidthPhasedDosage => read_fixed_width_variant_packed_sequential(
                pgen,
                &mut file,
                header,
                &mut decoder_state,
            )?,
            PgenLayout::VariableWidth => {
                read_plink2_variant_packed(
                    pgen,
                    &mut file,
                    header,
                    variant_index,
                    &mut decoder_state,
                )?;
            }
        }

        let (retain_variant, _) = evaluate_packed_hardcall_filter(
            &decoder_state.packed,
            &selection.source_indices,
            all_samples_selected,
            filter,
            genotype_filter_plan,
            None,
            false,
        )?;
        match retention.genotype_decision(retain_variant, &mut diagnostics) {
            RetentionAction::Include => {
                packed_batch.push(&decoder_state.packed);
                output_variant_count += 1;
                if packed_batch.is_full() {
                    flush_hardcall_batch_into_sample_major(
                        &mut packed_batch,
                        &selection.source_indices,
                        &mut batch_start,
                        output_variant_capacity,
                        &mut values,
                        &mut missing_mask,
                    );
                }
            }
            RetentionAction::Skip => {}
            RetentionAction::Stop => break,
        }
        if retention.window_is_satisfied() {
            break;
        }
    }

    flush_hardcall_batch_into_sample_major(
        &mut packed_batch,
        &selection.source_indices,
        &mut batch_start,
        output_variant_capacity,
        &mut values,
        &mut missing_mask,
    );
    shrink_sample_major_width(
        &mut values,
        n_samples,
        output_variant_capacity,
        output_variant_count,
    );
    shrink_sample_major_width(
        &mut missing_mask,
        n_samples,
        output_variant_capacity,
        output_variant_count,
    );
    diagnostics.retained_variants = output_variant_count;

    let matrix = DenseGenotypeMatrix::new_matrix_only(
        n_samples,
        output_variant_count,
        values,
        missing_mask,
        diagnostics,
    )?;
    Ok(MatrixOnlyGenotypeFilterRead {
        matrix,
        window_satisfied: retention.window_is_satisfied(),
    })
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
                    flush_hardcall_batch_into_sample_major(
                        &mut packed_batch,
                        source_indices,
                        &mut batch_start,
                        n_variants,
                        &mut values,
                        &mut missing_mask,
                    );
                }
            }
            flush_hardcall_batch_into_sample_major(
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
                        flush_hardcall_batch_into_sample_major(
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
            flush_hardcall_batch_into_sample_major(
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
    let Plink2ReadContext {
        header,
        selection,
        all_samples_selected: _,
    } = Plink2ReadContext::new_prefix(pgen, psam, requested_samples, decode_variant_ct)?;
    let mut diagnostics = selection.diagnostics;
    let window_variants = parse_pvar_source_window(pvar, window, header.variant_ct)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let n_variants = window_variants.len();
    let mut variants = Vec::with_capacity(n_variants);
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
                    flush_hardcall_batch_into_sample_major(
                        &mut packed_batch,
                        &selection.source_indices,
                        &mut batch_start,
                        n_variants,
                        &mut values,
                        &mut missing_mask,
                    );
                }
            }
            flush_hardcall_batch_into_sample_major(
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
                            flush_hardcall_batch_into_sample_major(
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
            flush_hardcall_batch_into_sample_major(
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
