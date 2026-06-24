// pattern: Imperative Shell
//! PLINK1 BED/BIM/FAM readers.
//!
//! This module coordinates companion-file validation, sample selection, variant
//! filtering, and dense or sparse output assembly. Binary BED decoding and text
//! metadata parsing live in submodules.

use std::fs::{self, File};
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, flip_values_to_minor_allele, reject_sparse_missing,
    select_samples_source_order, DenseDiagnostics, DenseGenotypeMatrix,
    DenseGenotypeMatrixArrowVariants, DenseLayout, DenseMissingPolicy, DenseSampleSelection,
    GenoioError, GenotypeFilterPlan, MetadataArrowOutput, MetadataOutput, PartialFilterDecision,
    SourceCapabilities, SparseGenotypeMatrix, SparseGenotypeMatrixArrowVariants, VariantFilter,
    VariantMetadataArrowBuffers, VariantRecord, VariantWindow,
};

use crate::error::Result;
use crate::hardcall::{
    evaluate_packed_hardcall_filter, flush_hardcall_batch_into_sample_major, HardcallBatch,
};
use crate::matrix::shrink_sample_major_width;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

mod bed;
mod metadata;

use bed::{
    infer_bed_variant_count, open_bed_file, read_plink1_variant_packed,
    read_plink1_variant_packed_sequential, seek_plink1_variant, validate_bed_payload_len,
    Plink1DecoderState,
};
use metadata::{parse_bim, parse_bim_source_window, parse_fam};

fn can_skip_bim_for_matrix_only_genotype_filter(
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

struct Plink1DenseContext<'a> {
    bed: &'a Path,
    bed_file: File,
    n_source_samples: usize,
    n_source_variants: usize,
    bytes_per_variant: usize,
    selection: DenseSampleSelection,
    all_samples_selected: bool,
}

fn dense_arrow_output_to_rows(
    output: DenseGenotypeMatrixArrowVariants,
    context: &'static str,
) -> Result<DenseGenotypeMatrix> {
    output.into_matrix().map_err(|error| {
        GenoioError::internal_contract(format!(
            "PLINK1 {context} Arrow-to-row compatibility conversion failed: {error}"
        ))
    })
}

fn sparse_arrow_output_to_rows(
    output: SparseGenotypeMatrixArrowVariants,
    context: &'static str,
) -> Result<SparseGenotypeMatrix> {
    output.into_matrix().map_err(|error| {
        GenoioError::internal_contract(format!(
            "PLINK1 {context} Arrow-to-row compatibility conversion failed: {error}"
        ))
    })
}

/// Read PLINK1 sample and variant metadata without decoding BED genotypes.
pub fn read_plink1_metadata(bed: &Path, bim: &Path, fam: &Path) -> Result<MetadataOutput> {
    fs::metadata(bed).map_err(|source| GenoioError::Io {
        path: bed.to_path_buf(),
        source,
    })?;
    let samples = parse_fam(fam)?;
    let variants = parse_bim(bim)?;

    Ok(MetadataOutput {
        samples,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}

/// Read PLINK1 metadata with variant metadata staged as Arrow-compatible buffers.
pub fn read_plink1_metadata_arrow(
    bed: &Path,
    bim: &Path,
    fam: &Path,
) -> Result<MetadataArrowOutput> {
    read_plink1_metadata(bed, bim, fam).and_then(MetadataArrowOutput::from_metadata)
}

/// Read all retained PLINK1 genotypes as a dense sample-by-variant matrix.
pub fn read_plink1_dense(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_plink1_dense_windowed_with_missing_policy(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        None,
        DenseMissingPolicy::Nan,
        false,
    )
}

/// Read retained PLINK1 genotypes as a dense matrix over an optional block window.
pub fn read_plink1_dense_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_plink1_dense_windowed_with_missing_policy(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        matrix_only,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "reader boundary keeps PLINK companions and dense policy controls explicit"
)]
pub fn read_plink1_dense_windowed_with_missing_policy(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_plink1_dense_windowed_with_arrow_variants(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        variant_window,
        missing_policy,
        !matrix_only,
        !matrix_only,
    )
    .and_then(|output| dense_arrow_output_to_rows(output, "dense"))
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors dense read options plus metadata return choices"
)]
pub fn read_plink1_dense_windowed_with_arrow_variants(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let matrix_only = !return_samples && !return_variants;
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        return empty_plink1_dense_arrow(
            bed,
            bim,
            fam,
            requested_samples,
            return_samples,
            return_variants,
        );
    }

    if let (None, Some(window)) = (variant_filter, variant_window) {
        return read_plink1_dense_source_window_arrow(
            bed,
            bim,
            fam,
            requested_samples,
            window,
            missing_policy,
            return_samples,
            return_variants,
        );
    }

    let bed_file = open_bed_file(bed)?;

    let all_samples = parse_fam(fam)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
    let all_samples_selected = requested_samples.is_none();
    let n_source_samples = all_samples.len();
    let bytes_per_variant = n_source_samples.div_ceil(4);
    if let Some((filter, window)) =
        can_skip_bim_for_matrix_only_genotype_filter(matrix_only, variant_filter, variant_window)
    {
        // Matrix-only genotype-stat filters need BED hard calls and FAM sample order,
        // but not parsed BIM rows. Still require the companion file to exist.
        fs::metadata(bim).map_err(|source| GenoioError::Io {
            path: bim.to_path_buf(),
            source,
        })?;
        let n_source_variants =
            infer_bed_variant_count(bed, &bed_file, n_source_samples, bytes_per_variant)?;
        let context = Plink1DenseContext {
            bed,
            bed_file,
            n_source_samples,
            n_source_variants,
            bytes_per_variant,
            selection,
            all_samples_selected,
        };
        return read_plink1_dense_matrix_only_genotype_filter_arrow(
            context,
            filter,
            window,
            missing_policy,
        );
    }
    let source_variants = parse_bim(bim)?;
    let n_source_variants = source_variants.len();
    let diagnostics = selection.diagnostics.clone();
    let context = Plink1DenseContext {
        bed,
        bed_file,
        n_source_samples,
        n_source_variants,
        bytes_per_variant,
        selection,
        all_samples_selected,
    };
    read_plink1_dense_with_variants_arrow(
        context,
        source_variants,
        variant_filter,
        variant_window,
        missing_policy,
        return_samples,
        return_variants,
        diagnostics,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "dense PLINK1 loop keeps decoded source state and metadata controls explicit"
)]
fn read_plink1_dense_with_variants_arrow(
    mut context: Plink1DenseContext<'_>,
    source_variants: Vec<VariantRecord>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    mut diagnostics: DenseDiagnostics,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    validate_bed_payload_len(
        context.bed,
        &context.bed_file,
        context.n_source_samples,
        context.n_source_variants,
        context.bytes_per_variant,
    )?;

    let output_variant_capacity = variant_window.map_or(context.n_source_variants, |window| {
        window.len.min(context.n_source_variants)
    });
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let n_samples = context.selection.samples.len();
    let mut values = vec![0.0; n_samples * output_variant_capacity];
    let mut batch = HardcallBatch::new(context.n_source_samples);
    let mut batch_start = 0_usize;
    let mut variant_values = Vec::with_capacity(n_samples);
    let mut missing_indices = Vec::new();
    let mut decoder_state = Plink1DecoderState::new(
        context.n_source_samples,
        context.bytes_per_variant,
        context.selection.samples.len(),
    );
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;
    let genotype_filter_plan = variant_filter.map_or(
        GenotypeFilterPlan::Generic,
        VariantFilter::genotype_filter_plan,
    );
    for (variant_index, mut variant) in source_variants.into_iter().enumerate() {
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        read_plink1_variant_packed(
            context.bed,
            &mut context.bed_file,
            variant_index,
            context.bytes_per_variant,
            context.n_source_samples,
            &mut decoder_state,
        )?;

        let mut stats = None;
        if needs_genotype_decision {
            let filter = variant_filter.ok_or_else(|| {
                GenoioError::internal_contract("genotype decision requires a variant filter")
            })?;
            let (retain_variant, computed_stats) = evaluate_packed_hardcall_filter(
                &decoder_state.packed,
                &context.selection.source_indices,
                context.all_samples_selected,
                filter,
                genotype_filter_plan,
                Some(&variant),
                return_variants,
            )?;
            stats = computed_stats;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        batch.push(&decoder_state.packed);
        output_variant_count += 1;
        if batch.is_full() {
            flush_hardcall_batch_into_sample_major(
                &mut batch,
                &context.selection.source_indices,
                &mut batch_start,
                output_variant_capacity,
                &mut values,
                missing_policy,
                &mut variant_values,
                &mut missing_indices,
            )?;
        }
        if retention.window_is_satisfied() {
            break;
        }
    }
    flush_hardcall_batch_into_sample_major(
        &mut batch,
        &context.selection.source_indices,
        &mut batch_start,
        output_variant_capacity,
        &mut values,
        missing_policy,
        &mut variant_values,
        &mut missing_indices,
    )?;

    let n_variants = output_variant_count;
    shrink_sample_major_width(&mut values, n_samples, output_variant_capacity, n_variants);
    diagnostics.retained_variants = n_variants;
    let samples = if return_samples {
        context.selection.samples
    } else {
        Vec::new()
    };
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        n_variants,
        values,
        DenseLayout::SampleMajor,
        samples,
        variants,
        diagnostics,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "source-window reader keeps companion files and metadata controls explicit"
)]
fn read_plink1_dense_source_window_arrow(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    window: VariantWindow,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let mut bed_file = open_bed_file(bed)?;
    let all_samples = parse_fam(fam)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
    let mut diagnostics = selection.diagnostics;
    let n_source_samples = all_samples.len();
    let bytes_per_variant = n_source_samples.div_ceil(4);
    let n_source_variants =
        infer_bed_variant_count(bed, &bed_file, n_source_samples, bytes_per_variant)?;
    let n_variants = n_source_variants
        .saturating_sub(window.start)
        .min(window.len);
    let variants = if return_variants {
        let source_variants = parse_bim_source_window(bim, window, n_variants)?;
        let mut buffers = VariantMetadataArrowBuffers::with_capacity(source_variants.len());
        for variant in source_variants {
            buffers.push_record(&variant)?;
        }
        Some(buffers)
    } else {
        fs::metadata(bim).map_err(|source| GenoioError::Io {
            path: bim.to_path_buf(),
            source,
        })?;
        None
    };

    let n_samples = selection.samples.len();
    let mut values = vec![0.0; n_samples * n_variants];
    let mut decoder_state =
        Plink1DecoderState::new(n_source_samples, bytes_per_variant, selection.samples.len());
    let mut batch = HardcallBatch::new(n_source_samples);
    let mut batch_start = 0_usize;
    let mut variant_values = Vec::with_capacity(n_samples);
    let mut missing_indices = Vec::new();
    seek_plink1_variant(bed, &mut bed_file, window.start, bytes_per_variant)?;
    for variant_offset in 0..n_variants {
        read_plink1_variant_packed_sequential(
            bed,
            &mut bed_file,
            bytes_per_variant,
            n_source_samples,
            &mut decoder_state,
        )?;
        batch.push(&decoder_state.packed);
        if batch.is_full() {
            flush_hardcall_batch_into_sample_major(
                &mut batch,
                &selection.source_indices,
                &mut batch_start,
                n_variants,
                &mut values,
                missing_policy,
                &mut variant_values,
                &mut missing_indices,
            )?;
        }
        diagnostics.candidate_variants = variant_offset + 1;
    }
    flush_hardcall_batch_into_sample_major(
        &mut batch,
        &selection.source_indices,
        &mut batch_start,
        n_variants,
        &mut values,
        missing_policy,
        &mut variant_values,
        &mut missing_indices,
    )?;

    diagnostics.retained_variants = n_variants;
    let samples = if return_samples {
        selection.samples
    } else {
        Vec::new()
    };
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        n_variants,
        values,
        DenseLayout::SampleMajor,
        samples,
        variants,
        diagnostics,
    )
}

fn read_plink1_dense_matrix_only_genotype_filter_arrow(
    mut context: Plink1DenseContext<'_>,
    filter: &VariantFilter,
    window: VariantWindow,
    missing_policy: DenseMissingPolicy,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let output_variant_capacity = window.len.min(context.n_source_variants);
    let n_samples = context.selection.samples.len();
    if output_variant_capacity == 0 {
        let mut diagnostics = context.selection.diagnostics;
        diagnostics.retained_variants = 0;
        return DenseGenotypeMatrixArrowVariants::new_with_layout(
            n_samples,
            0,
            Vec::new(),
            DenseLayout::SampleMajor,
            Vec::new(),
            None,
            diagnostics,
        );
    }

    let mut values = vec![0.0; n_samples * output_variant_capacity];
    let mut decoder_state = Plink1DecoderState::new(
        context.n_source_samples,
        context.bytes_per_variant,
        context.selection.samples.len(),
    );
    let mut batch = HardcallBatch::new(context.n_source_samples);
    let mut batch_start = 0_usize;
    let mut variant_values = Vec::with_capacity(n_samples);
    let mut missing_indices = Vec::new();
    let mut retention = RetainedVariantState::new(Some(window));
    let mut diagnostics = context.selection.diagnostics;
    let mut output_variant_count = 0_usize;
    let genotype_filter_plan = filter.genotype_filter_plan();

    for _ in 0..context.n_source_variants {
        diagnostics.candidate_variants += 1;
        read_plink1_variant_packed_sequential(
            context.bed,
            &mut context.bed_file,
            context.bytes_per_variant,
            context.n_source_samples,
            &mut decoder_state,
        )?;

        let (retain_variant, _) = evaluate_packed_hardcall_filter(
            &decoder_state.packed,
            &context.selection.source_indices,
            context.all_samples_selected,
            filter,
            genotype_filter_plan,
            None,
            false,
        )?;
        match retention.genotype_decision(retain_variant, &mut diagnostics) {
            RetentionAction::Include => {
                batch.push(&decoder_state.packed);
                output_variant_count += 1;
                if batch.is_full() {
                    flush_hardcall_batch_into_sample_major(
                        &mut batch,
                        &context.selection.source_indices,
                        &mut batch_start,
                        output_variant_capacity,
                        &mut values,
                        missing_policy,
                        &mut variant_values,
                        &mut missing_indices,
                    )?;
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
        &mut batch,
        &context.selection.source_indices,
        &mut batch_start,
        output_variant_capacity,
        &mut values,
        missing_policy,
        &mut variant_values,
        &mut missing_indices,
    )?;
    shrink_sample_major_width(
        &mut values,
        n_samples,
        output_variant_capacity,
        output_variant_count,
    );
    diagnostics.retained_variants = output_variant_count;

    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        output_variant_count,
        values,
        DenseLayout::SampleMajor,
        Vec::new(),
        None,
        diagnostics,
    )
}

/// Read all retained PLINK1 genotypes as a sparse CSC matrix.
pub fn read_plink1_sparse(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_plink1_sparse_windowed(bed, bim, fam, requested_samples, variant_filter, None)
}

/// Read retained PLINK1 genotypes as sparse CSC over an optional block window.
pub fn read_plink1_sparse_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    read_plink1_sparse_windowed_with_arrow_variants(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )
    .and_then(|output| sparse_arrow_output_to_rows(output, "sparse"))
}

#[expect(
    clippy::too_many_arguments,
    reason = "Arrow facade mirrors sparse read options plus metadata return choices"
)]
pub fn read_plink1_sparse_windowed_with_arrow_variants(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        return empty_plink1_sparse_arrow(
            bed,
            bim,
            fam,
            requested_samples,
            return_samples,
            return_variants,
        );
    }

    let mut bed_file = open_bed_file(bed)?;

    let all_samples = parse_fam(fam)?;
    let source_variants = parse_bim(bim)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
    let all_samples_selected = requested_samples.is_none();
    let mut diagnostics = selection.diagnostics;
    let n_source_samples = all_samples.len();
    let n_source_variants = source_variants.len();
    let bytes_per_variant = n_source_samples.div_ceil(4);
    validate_bed_payload_len(
        bed,
        &bed_file,
        n_source_samples,
        n_source_variants,
        bytes_per_variant,
    )?;

    let n_samples = selection.samples.len();
    let output_variant_capacity = variant_window.map_or(n_source_variants, |window| {
        window.len.min(n_source_variants)
    });
    let mut indptr = Vec::with_capacity(output_variant_capacity + 1);
    indptr.push(0);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = return_variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(output_variant_capacity));
    let mut decoder_state =
        Plink1DecoderState::new(n_source_samples, bytes_per_variant, selection.samples.len());
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;
    let genotype_filter_plan = variant_filter.map_or(
        GenotypeFilterPlan::Generic,
        VariantFilter::genotype_filter_plan,
    );
    for (variant_index, mut variant) in source_variants.into_iter().enumerate() {
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        read_plink1_variant_packed(
            bed,
            &mut bed_file,
            variant_index,
            bytes_per_variant,
            n_source_samples,
            &mut decoder_state,
        )?;
        let mut stats = None;
        if needs_genotype_decision {
            let filter = variant_filter.ok_or_else(|| {
                GenoioError::internal_contract("genotype decision requires a variant filter")
            })?;
            let (retain_variant, computed_stats) = evaluate_packed_hardcall_filter(
                &decoder_state.packed,
                &selection.source_indices,
                all_samples_selected,
                filter,
                genotype_filter_plan,
                Some(&variant),
                return_variants,
            )?;
            stats = computed_stats;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        decoder_state.packed.expand_selected(
            &selection.source_indices,
            &mut decoder_state.values,
            &mut decoder_state.missing_indices,
        );
        reject_sparse_missing(!decoder_state.missing_indices.is_empty())?;
        flip_values_to_minor_allele(&mut decoder_state.values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &decoder_state.values);
        output_variant_count += 1;
        if let Some(variants) = variants.as_mut() {
            variants.push_record(&variant)?;
        }
        if retention.window_is_satisfied() {
            break;
        }
    }

    let n_variants = output_variant_count;
    diagnostics.retained_variants = n_variants;
    let samples = if return_samples {
        selection.samples
    } else {
        Vec::new()
    };
    SparseGenotypeMatrixArrowVariants::new(
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

fn empty_plink1_dense_arrow(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    fs::metadata(bed).map_err(|source| GenoioError::Io {
        path: bed.to_path_buf(),
        source,
    })?;
    fs::metadata(bim).map_err(|source| GenoioError::Io {
        path: bim.to_path_buf(),
        source,
    })?;
    let all_samples = parse_fam(fam)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
    let n_samples = selection.samples.len();
    let samples = if return_samples {
        selection.samples
    } else {
        Vec::new()
    };
    let variants = return_variants.then(|| VariantMetadataArrowBuffers::with_capacity(0));
    DenseGenotypeMatrixArrowVariants::new_with_layout(
        n_samples,
        0,
        Vec::new(),
        DenseLayout::SampleMajor,
        samples,
        variants,
        selection.diagnostics,
    )
}

fn empty_plink1_sparse_arrow(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    fs::metadata(bed).map_err(|source| GenoioError::Io {
        path: bed.to_path_buf(),
        source,
    })?;
    fs::metadata(bim).map_err(|source| GenoioError::Io {
        path: bim.to_path_buf(),
        source,
    })?;
    let all_samples = parse_fam(fam)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
    let n_samples = selection.samples.len();
    let samples = if return_samples {
        selection.samples
    } else {
        Vec::new()
    };
    let variants = return_variants.then(|| VariantMetadataArrowBuffers::with_capacity(0));
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    SparseGenotypeMatrixArrowVariants::new(
        n_samples,
        0,
        vec![0],
        Vec::new(),
        Vec::new(),
        samples,
        variants,
        diagnostics,
    )
}
