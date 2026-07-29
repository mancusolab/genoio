// pattern: Imperative Shell
//! PLINK1 BED/BIM/FAM readers.
//!
//! This module coordinates companion-file validation, sample selection, variant
//! filtering, and dense or sparse output assembly. Binary BED decoding and text
//! metadata parsing live in submodules.

use std::fs;
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, flip_values_to_minor_allele, reject_sparse_missing,
    select_samples_source_order, DenseGenotypeMatrix, DenseMissingPolicy, GenoioError,
    GenotypeFilterPlan, MetadataOutput, PartialFilterDecision, SampleMetadataBuffers,
    SourceCapabilities, SparseGenotypeMatrix, VariantFilter, VariantMetadataBuffers, VariantWindow,
};

use crate::error::Result;
use crate::hardcall::evaluate_packed_hardcall_filter;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

mod bed;
mod metadata;
mod session;

pub(crate) use session::Plink1BlockSession;

use bed::{
    open_bed_file, read_plink1_variant_packed, validate_bed_payload_len, Plink1DecoderState,
};
use metadata::{count_bim_records, parse_bim_metadata, parse_fam, BimRecordReader};

/// Read PLINK1 metadata with variant metadata staged as columnar buffers.
pub fn read_plink1_metadata(bed: &Path, bim: &Path, fam: &Path) -> Result<MetadataOutput> {
    fs::metadata(bed).map_err(|source| GenoioError::Io {
        path: bed.to_path_buf(),
        source,
    })?;
    let samples = parse_fam(fam)?;
    let variants = parse_bim_metadata(bim)?;

    Ok(MetadataOutput {
        samples: SampleMetadataBuffers::from_records(&samples, false)?,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors dense read options plus metadata return choices"
)]
pub fn read_plink1_dense_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    let matrix_only = !return_samples && !return_variants;
    let retained_skip = variant_window.map_or(0, |window| window.start);
    let read_bim_records = !(matrix_only
        && variant_window.is_some()
        && variant_filter.is_none_or(|filter| {
            filter.requires_genotype_stats() && filter.is_genotype_stats_only()
        }));
    let options = crate::blocks::BlockReadOptions {
        matrix_kind: crate::blocks::MatrixKind::Genotype,
        sparse: false,
        requested_samples: requested_samples.map(<[String]>::to_vec),
        variant_filter: variant_filter.cloned(),
        dosage_source: crate::blocks::DosageSource::Hardcall,
        missing_policy,
        return_samples,
        return_variants,
    };
    let mut session = Plink1BlockSession::open_windowed(
        bed.to_path_buf(),
        bim.to_path_buf(),
        fam.to_path_buf(),
        options,
        retained_skip,
        read_bim_records,
    )?;
    let block_size = match variant_window {
        Some(window) => window.len,
        None => session.source_record_capacity(),
    };
    match session.next_dense_block(block_size)? {
        Some(matrix) => Ok(matrix),
        None => session.empty_dense_output(),
    }
}

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors sparse read options plus metadata return choices"
)]
pub fn read_plink1_sparse_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        return empty_plink1_sparse(
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
    let n_source_variants = count_bim_records(bim)?;
    let mut source_variants = BimRecordReader::new(bim)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
    let all_samples_selected = requested_samples.is_none();
    let mut diagnostics = selection.diagnostics;
    let n_source_samples = all_samples.len();
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
    let mut variants =
        return_variants.then(|| VariantMetadataBuffers::with_capacity(output_variant_capacity));
    let mut decoder_state =
        Plink1DecoderState::new(n_source_samples, bytes_per_variant, selection.samples.len());
    let mut retention = RetainedVariantState::new(variant_window);
    let mut output_variant_count = 0_usize;
    let genotype_filter_plan = variant_filter.map_or(
        GenotypeFilterPlan::Generic,
        VariantFilter::genotype_filter_plan,
    );
    while let Some((variant_index, mut variant)) = source_variants.next_record()? {
        if variant_index >= n_source_variants {
            return Err(GenoioError::invalid_source(
                bed,
                "bim variant count exceeds bed variant count",
            ));
        }
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
        append_sparse_column(&mut indptr, &mut indices, &mut data, &decoder_state.values)?;
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

fn empty_plink1_sparse(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
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
    let samples =
        SampleMetadataBuffers::optional_from_records(&selection.samples, return_samples, false)?;
    let variants = return_variants.then(|| VariantMetadataBuffers::with_capacity(0));
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    SparseGenotypeMatrix::new(
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
