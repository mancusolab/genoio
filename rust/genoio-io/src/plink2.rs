// pattern: Imperative Shell

use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, flip_values_to_minor_allele,
    reject_sparse_missing_values, select_samples_source_order,
    transpose_variant_major_to_sample_major, DenseGenotypeMatrix, MetadataError, MetadataOutput,
    PartialFilterDecision, SampleRecord, SourceCapabilities, SparseGenotypeMatrix, VariantFilter,
    VariantRecord, VariantWindow,
};

use crate::error::Result;
#[cfg(test)]
use crate::hardcall::HARDCALL_BATCH_SIZE;
use crate::hardcall::{HardcallBatch as PackedVariantBatch, PackedHardcalls as PackedGenotypes};

const PGEN_MAGIC: [u8; 2] = [0x6c, 0x1b];
const PGEN_MODE_FIXED_WIDTH_HARDCALLS: u8 = 0x02;
const PGEN_MODE_FIXED_WIDTH_DOSAGE: u8 = 0x03;
const PGEN_MODE_VARIABLE_WIDTH: u8 = 0x10;
const PGEN_HEADER_LEN: u64 = 12;
const PGEN_VARIANT_BLOCK_SIZE: usize = 65_536;
#[cfg(test)]
const PGEN_PACKED_TRANSPOSE_BATCH: usize = HARDCALL_BATCH_SIZE;

#[derive(Debug, Clone)]
struct PgenHeader {
    layout: PgenLayout,
    variant_ct: usize,
    sample_ct: usize,
    bytes_per_variant: usize,
    record_types: Vec<u8>,
    record_offsets: Vec<u64>,
}

#[derive(Debug, Clone)]
struct PgenDecoderState {
    previous_non_ld_packed: PackedGenotypes,
    has_previous_non_ld: bool,
    record: Vec<u8>,
    packed: PackedGenotypes,
    values: Vec<f32>,
    missing: Vec<bool>,
}

struct DosageOverlayTarget<'a> {
    source_indices: &'a [usize],
    values: &'a mut [f32],
    missing: &'a mut [bool],
}

struct SelectedSampleCursor<'a> {
    source_indices: &'a [usize],
    selected_index: usize,
}

impl<'a> SelectedSampleCursor<'a> {
    fn new(source_indices: &'a [usize]) -> Self {
        Self {
            source_indices,
            selected_index: 0,
        }
    }

    fn selected_index_for(&mut self, source_index: usize) -> Option<usize> {
        // source_indices are stored in PGEN source order, so a forward-only
        // cursor avoids a search for every stored dosage.
        while self
            .source_indices
            .get(self.selected_index)
            .is_some_and(|selected_source_index| *selected_source_index < source_index)
        {
            self.selected_index += 1;
        }
        self.source_indices
            .get(self.selected_index)
            .copied()
            .filter(|selected_source_index| *selected_source_index == source_index)
            .map(|_| self.selected_index)
    }
}

impl PgenDecoderState {
    fn new(sample_ct: usize, selected_sample_ct: usize) -> Self {
        Self {
            previous_non_ld_packed: PackedGenotypes::default(),
            has_previous_non_ld: false,
            record: Vec::with_capacity(sample_ct.div_ceil(4)),
            packed: PackedGenotypes::default(),
            values: Vec::with_capacity(selected_sample_ct),
            missing: Vec::with_capacity(selected_sample_ct),
        }
    }
}

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

#[derive(Debug, Clone)]
enum PgenLayout {
    FixedWidth,
    FixedWidthDosage,
    VariableWidth,
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
        if requested_samples.is_none() && matrix_only {
            return read_plink2_dense_matrix_only_source_window(pgen, window);
        }
        return read_plink2_dense_source_window(pgen, pvar, psam, requested_samples, window);
    }

    let header = read_supported_pgen_header(pgen)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        fs::metadata(pvar).map_err(|source| MetadataError::Io {
            path: pvar.to_path_buf(),
            source,
        })?;
        diagnostics.retained_variants = 0;
        return DenseGenotypeMatrix::new(
            selection.samples.len(),
            0,
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

    let output_variant_capacity = variant_window.map_or(header.variant_ct, |window| {
        window.len.min(header.variant_ct)
    });
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut variant_major_missing =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut retained_index = 0_usize;
    let mut stopped_after_window = false;
    let requires_sequential_decode = matches!(header.layout, PgenLayout::VariableWidth);
    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        diagnostics.candidate_variants += 1;
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
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                        stopped_after_window = true;
                        break;
                    }
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {}
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
        if needs_genotype_decision
            && variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref()))
        {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if needs_genotype_decision {
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                    stopped_after_window = true;
                    break;
                }
                continue;
            }
        }
        decoder_state.packed.expand_selected(
            &selection.source_indices,
            &mut decoder_state.values,
            &mut decoder_state.missing,
        );
        variants.push(variant);
        variant_major_values.extend_from_slice(&decoder_state.values);
        variant_major_missing.extend_from_slice(&decoder_state.missing);
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = selection.samples.len();
    let n_variants = variants.len();
    diagnostics.retained_variants = n_variants;
    let values =
        transpose_variant_major_to_sample_major(&variant_major_values, n_samples, n_variants);
    let missing_mask =
        transpose_variant_major_to_sample_major(&variant_major_missing, n_samples, n_variants);

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

/// Read retained PLINK2 unphased biallelic dosages as a dense matrix.
pub fn read_plink2_dosage_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<DenseGenotypeMatrix> {
    if variant_filter.is_some_and(VariantFilter::requires_genotype_stats) {
        return Err(MetadataError::parse(
            pgen,
            "dosage-backed genotype reads do not support genotype-stat filters yet",
        ));
    }

    let header = read_supported_pgen_header(pgen)?;
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, &header, all_samples.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        fs::metadata(pvar).map_err(|source| MetadataError::Io {
            path: pvar.to_path_buf(),
            source,
        })?;
        diagnostics.retained_variants = 0;
        return DenseGenotypeMatrix::new(
            selection.samples.len(),
            0,
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
    let output_variant_capacity = variant_window.map_or(header.variant_ct, |window| {
        window.len.min(header.variant_ct)
    });
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut variant_major_missing =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut retained_index = 0_usize;
    let mut stopped_after_window = false;

    while let Some((variant_index, variant)) = pvar_reader.next_record()? {
        diagnostics.candidate_variants += 1;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                        stopped_after_window = true;
                        break;
                    }
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {
                return Err(MetadataError::parse(
                    pgen,
                    "dosage-backed genotype reads do not support genotype-stat filters yet",
                ));
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
        variants.push(variant);
        variant_major_values.extend_from_slice(&decoder_state.values);
        variant_major_missing.extend_from_slice(&decoder_state.missing);
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
            stopped_after_window = true;
            break;
        }
    }
    if !stopped_after_window {
        pvar_reader.validate_count(header.variant_ct)?;
    }

    let n_samples = selection.samples.len();
    let n_variants = variants.len();
    diagnostics.retained_variants = n_variants;
    let values =
        transpose_variant_major_to_sample_major(&variant_major_values, n_samples, n_variants);
    let missing_mask =
        transpose_variant_major_to_sample_major(&variant_major_missing, n_samples, n_variants);

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
        fs::metadata(pvar).map_err(|source| MetadataError::Io {
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
    let mut retained_index = 0_usize;
    let mut stopped_after_window = false;
    let requires_sequential_decode = matches!(header.layout, PgenLayout::VariableWidth);
    while let Some((variant_index, mut variant)) = pvar_reader.next_record()? {
        diagnostics.candidate_variants += 1;
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
        match partial_decision {
            PartialFilterDecision::Reject => {
                diagnostics.dropped_metadata_variants += 1;
                continue;
            }
            PartialFilterDecision::Accept => {
                let include_in_window =
                    variant_window.is_none_or(|window| window.contains(retained_index));
                retained_index += 1;
                if !include_in_window {
                    if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                        stopped_after_window = true;
                        break;
                    }
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {}
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
        if needs_genotype_decision
            && variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref()))
        {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if needs_genotype_decision {
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                if variant_window.is_some_and(|window| window.is_past(retained_index)) {
                    stopped_after_window = true;
                    break;
                }
                continue;
            }
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
        if variant_window.is_some_and(|window| window.is_past(retained_index)) {
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
    window: VariantWindow,
) -> Result<DenseGenotypeMatrix> {
    let decode_variant_ct = window
        .start
        .checked_add(window.len)
        .ok_or_else(|| MetadataError::parse(pgen, "variant window end is out of range"))?;
    let header = read_supported_pgen_header_prefix(pgen, decode_variant_ct)?;
    if window.start > header.variant_ct {
        return Err(MetadataError::parse(
            pgen,
            format!(
                "variant window start {} exceeds pgen variant count {}",
                window.start, header.variant_ct
            ),
        ));
    }
    let n_variants = window.len.min(header.variant_ct - window.start);
    let n_samples = header.sample_ct;
    let source_indices = (0..n_samples).collect::<Vec<_>>();
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, n_samples);
    let mut values = vec![0.0; n_samples * n_variants];
    let mut missing_mask = vec![false; n_samples * n_variants];
    let mut packed_batch = PackedVariantBatch::new(header.sample_ct);
    let mut batch_start = 0_usize;

    // Unfiltered source windows know their final retained width up front, so
    // construct the public sample-major buffers directly.
    match header.layout {
        PgenLayout::FixedWidth | PgenLayout::FixedWidthDosage => {
            if n_variants > 0 {
                seek_fixed_width_variant_record(pgen, &mut file, &header, window.start)?;
            }
            for _ in 0..n_variants {
                read_fixed_width_variant_packed_sequential(
                    pgen,
                    &mut file,
                    &header,
                    &mut decoder_state,
                )?;
                packed_batch.push(&decoder_state.packed);
                if packed_batch.is_full() {
                    flush_packed_variant_batch(
                        &mut packed_batch,
                        &source_indices,
                        &mut batch_start,
                        n_variants,
                        &mut values,
                        &mut missing_mask,
                    );
                }
            }
            flush_packed_variant_batch(
                &mut packed_batch,
                &source_indices,
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
                    &header,
                    source_variant_index,
                    &mut decoder_state,
                )?;
                if source_variant_index >= window.start {
                    packed_batch.push(&decoder_state.packed);
                    if packed_batch.is_full() {
                        flush_packed_variant_batch(
                            &mut packed_batch,
                            &source_indices,
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
                &source_indices,
                &mut batch_start,
                n_variants,
                &mut values,
                &mut missing_mask,
            );
        }
    }

    let diagnostics = genoio_core::DenseDiagnostics {
        requested_samples: n_samples,
        retained_samples: n_samples,
        candidate_variants: n_variants,
        retained_variants: n_variants,
        ..genoio_core::DenseDiagnostics::default()
    };

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
        PgenLayout::FixedWidth | PgenLayout::FixedWidthDosage => {
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
                    let (_, variant) = window_iter.next().expect("peeked variant should exist");
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
        PgenLayout::FixedWidth | PgenLayout::FixedWidthDosage => {
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
                    let (_, mut variant) = window_iter.next().expect("peeked variant should exist");
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

fn read_supported_pgen_header(path: &Path) -> Result<PgenHeader> {
    let mut file = File::open(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = [0_u8; PGEN_HEADER_LEN as usize];
    file.read_exact(&mut header)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    if header[0..2] != PGEN_MAGIC {
        return Err(MetadataError::parse(path, "invalid pgen magic bytes"));
    }
    let variant_ct = usize::try_from(u32::from_le_bytes(header[3..7].try_into().unwrap()))
        .map_err(|_| MetadataError::parse(path, "pgen variant count is out of range"))?;
    let sample_ct = usize::try_from(u32::from_le_bytes(header[7..11].try_into().unwrap()))
        .map_err(|_| MetadataError::parse(path, "pgen sample count is out of range"))?;
    let bytes_per_variant = sample_ct.div_ceil(4);
    match header[2] {
        PGEN_MODE_FIXED_WIDTH_HARDCALLS => {
            if header[11] != 0 {
                return Err(MetadataError::parse(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic hardcalls without header extensions are supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(path, &file, variant_ct, bytes_per_variant)?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidth,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_FIXED_WIDTH_DOSAGE => {
            if header[11] != 0 {
                return Err(MetadataError::parse(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic dosage without header extensions is supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(
                path,
                &file,
                variant_ct,
                fixed_width_dosage_record_len(sample_ct),
            )?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidthDosage,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_VARIABLE_WIDTH => {
            let (record_types, record_offsets) =
                read_variable_width_header_body(path, &mut file, variant_ct, header[11])?;
            Ok(PgenHeader {
                layout: PgenLayout::VariableWidth,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types,
                record_offsets,
            })
        }
        mode => Err(MetadataError::parse(
            path,
            format!(
                "unsupported pgen mode 0x{mode:02x}; only fixed-width and variable-width biallelic hardcalls are supported"
            ),
        )),
    }
}

fn read_supported_pgen_header_prefix(
    path: &Path,
    requested_variant_ct: usize,
) -> Result<PgenHeader> {
    let mut file = File::open(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = [0_u8; PGEN_HEADER_LEN as usize];
    file.read_exact(&mut header)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    if header[0..2] != PGEN_MAGIC {
        return Err(MetadataError::parse(path, "invalid pgen magic bytes"));
    }
    let variant_ct = usize::try_from(u32::from_le_bytes(header[3..7].try_into().unwrap()))
        .map_err(|_| MetadataError::parse(path, "pgen variant count is out of range"))?;
    let sample_ct = usize::try_from(u32::from_le_bytes(header[7..11].try_into().unwrap()))
        .map_err(|_| MetadataError::parse(path, "pgen sample count is out of range"))?;
    let bytes_per_variant = sample_ct.div_ceil(4);
    let prefix_variant_ct = requested_variant_ct.min(variant_ct);
    match header[2] {
        PGEN_MODE_FIXED_WIDTH_HARDCALLS => {
            if header[11] != 0 {
                return Err(MetadataError::parse(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic hardcalls without header extensions are supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(path, &file, variant_ct, bytes_per_variant)?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidth,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_FIXED_WIDTH_DOSAGE => {
            if header[11] != 0 {
                return Err(MetadataError::parse(
                    path,
                    "unsupported pgen header flags; only fixed-width biallelic dosage without header extensions is supported",
                ));
            }
            validate_fixed_width_pgen_payload_len(
                path,
                &file,
                variant_ct,
                fixed_width_dosage_record_len(sample_ct),
            )?;
            Ok(PgenHeader {
                layout: PgenLayout::FixedWidthDosage,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types: Vec::new(),
                record_offsets: Vec::new(),
            })
        }
        PGEN_MODE_VARIABLE_WIDTH => {
            let (record_types, record_offsets) = read_variable_width_header_body_prefix(
                path,
                &mut file,
                variant_ct,
                header[11],
                prefix_variant_ct,
            )?;
            Ok(PgenHeader {
                layout: PgenLayout::VariableWidth,
                variant_ct,
                sample_ct,
                bytes_per_variant,
                record_types,
                record_offsets,
            })
        }
        mode => Err(MetadataError::parse(
            path,
            format!(
                "unsupported pgen mode 0x{mode:02x}; only fixed-width and variable-width biallelic hardcalls or unphased dosages are supported"
            ),
        )),
    }
}

fn read_variable_width_header_body(
    path: &Path,
    file: &mut File,
    variant_ct: usize,
    header_format: u8,
) -> Result<(Vec<u8>, Vec<u64>)> {
    let type_length_format = header_format & 0x0f;
    let type_width_bits = match type_length_format {
        0..=3 => 4,
        4..=7 => 8,
        other => {
            return Err(MetadataError::parse(
                path,
                format!("unsupported pgen variant-record type/length format {other}"),
            ));
        }
    };
    let length_width = usize::from((type_length_format & 0x03) + 1);
    let allele_count_format = (header_format >> 4) & 0x03;
    if allele_count_format != 0 {
        return Err(MetadataError::parse(
            path,
            "unsupported pgen allele-count table; multiallelic PGEN decode is not implemented",
        ));
    }

    let block_ct = variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let mut block_offsets = Vec::with_capacity(block_ct);
    for _ in 0..block_ct {
        let mut bytes = [0_u8; 8];
        file.read_exact(&mut bytes)
            .map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        block_offsets.push(u64::from_le_bytes(bytes));
    }

    let mut record_types = Vec::with_capacity(variant_ct);
    let mut record_lengths = Vec::with_capacity(variant_ct);
    for block_index in 0..block_ct {
        let block_variant_ct = block_variant_count(variant_ct, block_index);
        if type_width_bits == 8 {
            let mut types = vec![0_u8; block_variant_ct];
            file.read_exact(&mut types)
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            record_types.extend(types);
        } else {
            let mut packed_types = vec![0_u8; block_variant_ct.div_ceil(2)];
            file.read_exact(&mut packed_types)
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            for variant_in_block in 0..block_variant_ct {
                let byte = packed_types[variant_in_block / 2];
                let record_type = if variant_in_block % 2 == 0 {
                    byte & 0x0f
                } else {
                    byte >> 4
                };
                record_types.push(record_type);
            }
        }

        for _ in 0..block_variant_ct {
            let mut bytes = [0_u8; 4];
            file.read_exact(&mut bytes[..length_width])
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            record_lengths.push(u32::from_le_bytes(bytes));
        }
    }
    for record_type in &record_types {
        validate_supported_variable_record_type(path, *record_type)?;
    }
    let header_end = file.stream_position().map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    let mut record_offsets = Vec::with_capacity(variant_ct + 1);
    for (block_index, block_offset) in block_offsets.iter().enumerate() {
        let block_start = block_index * PGEN_VARIANT_BLOCK_SIZE;
        let block_end =
            (block_start + block_variant_count(variant_ct, block_index)).min(variant_ct);
        let mut offset = *block_offset;
        if record_offsets.len() == block_start {
            if block_index == 0 && offset != header_end {
                return Err(MetadataError::parse(
                    path,
                    "pgen first variant-block offset does not match header length",
                ));
            }
            record_offsets.push(offset);
        } else if record_offsets
            .get(block_start)
            .is_none_or(|expected_offset| *expected_offset != offset)
        {
            return Err(MetadataError::parse(
                path,
                "pgen variant-block offset does not match preceding record lengths",
            ));
        }
        for length in &record_lengths[block_start..block_end] {
            offset = offset
                .checked_add(u64::from(*length))
                .ok_or_else(|| MetadataError::parse(path, "pgen record offset is out of range"))?;
            record_offsets.push(offset);
        }
    }
    if record_offsets.len() != variant_ct + 1 {
        return Err(MetadataError::parse(
            path,
            "pgen variable-width header did not yield one offset per variant",
        ));
    }
    let actual_len = file
        .metadata()
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if record_offsets[variant_ct] > actual_len {
        return Err(MetadataError::parse(
            path,
            "pgen variable-width records extend past end of file",
        ));
    }

    Ok((record_types, record_offsets))
}

fn read_variable_width_header_body_prefix(
    path: &Path,
    file: &mut File,
    variant_ct: usize,
    header_format: u8,
    prefix_variant_ct: usize,
) -> Result<(Vec<u8>, Vec<u64>)> {
    let type_length_format = header_format & 0x0f;
    let type_width_bits = match type_length_format {
        0..=3 => 4,
        4..=7 => 8,
        other => {
            return Err(MetadataError::parse(
                path,
                format!("unsupported pgen variant-record type/length format {other}"),
            ));
        }
    };
    let length_width = usize::from((type_length_format & 0x03) + 1);
    let allele_count_format = (header_format >> 4) & 0x03;
    if allele_count_format != 0 {
        return Err(MetadataError::parse(
            path,
            "unsupported pgen allele-count table; multiallelic PGEN decode is not implemented",
        ));
    }

    let block_ct = variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let mut block_offsets = Vec::with_capacity(block_ct);
    for _ in 0..block_ct {
        let mut bytes = [0_u8; 8];
        file.read_exact(&mut bytes)
            .map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        block_offsets.push(u64::from_le_bytes(bytes));
    }

    let header_end = variable_width_header_end(path, variant_ct, type_width_bits, length_width)?;
    let prefix_block_ct = prefix_variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let mut record_types = Vec::with_capacity(prefix_variant_ct);
    let mut record_offsets = Vec::with_capacity(prefix_variant_ct.saturating_add(1));
    for (block_index, block_offset) in block_offsets
        .iter()
        .take(prefix_block_ct)
        .copied()
        .enumerate()
    {
        let block_variant_ct = block_variant_count(variant_ct, block_index);
        let block_start = block_index * PGEN_VARIANT_BLOCK_SIZE;
        let needed_in_block = prefix_variant_ct
            .saturating_sub(block_start)
            .min(block_variant_ct);
        // Type and length tables are block-grouped in the PGEN header. We
        // still have to skip through unneeded entries in the last touched
        // block so the file cursor reaches the matching length table.
        read_variable_record_type_prefix(
            path,
            file,
            type_width_bits,
            block_variant_ct,
            needed_in_block,
            &mut record_types,
        )?;
        if record_offsets.is_empty() {
            if block_index == 0 && block_offset != header_end {
                return Err(MetadataError::parse(
                    path,
                    "pgen first variant-block offset does not match header length",
                ));
            }
            record_offsets.push(block_offset);
        } else if record_offsets
            .get(block_start)
            .is_none_or(|expected_offset| *expected_offset != block_offset)
        {
            return Err(MetadataError::parse(
                path,
                "pgen variant-block offset does not match preceding record lengths",
            ));
        }
        let mut offset = block_offset;
        for _ in 0..needed_in_block {
            let mut bytes = [0_u8; 4];
            file.read_exact(&mut bytes[..length_width])
                .map_err(|source| MetadataError::Io {
                    path: path.to_path_buf(),
                    source,
                })?;
            offset = offset
                .checked_add(u64::from(u32::from_le_bytes(bytes)))
                .ok_or_else(|| MetadataError::parse(path, "pgen record offset is out of range"))?;
            record_offsets.push(offset);
        }
        let remaining_lengths = block_variant_ct - needed_in_block;
        skip_bytes(path, file, remaining_lengths * length_width)?;
    }
    // Only validate the prefix that may be decoded for this block. Unsupported
    // later records should not prevent first-block reads from succeeding.
    for record_type in &record_types {
        validate_supported_variable_record_type(path, *record_type)?;
    }
    if let Some(prefix_end) = record_offsets.last() {
        let actual_len = file
            .metadata()
            .map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?
            .len();
        if *prefix_end > actual_len {
            return Err(MetadataError::parse(
                path,
                "pgen variable-width records extend past end of file",
            ));
        }
    }
    Ok((record_types, record_offsets))
}

fn variable_width_header_end(
    path: &Path,
    variant_ct: usize,
    type_width_bits: usize,
    length_width: usize,
) -> Result<u64> {
    let block_ct = variant_ct.div_ceil(PGEN_VARIANT_BLOCK_SIZE);
    let block_offsets_len = block_ct
        .checked_mul(8)
        .ok_or_else(|| MetadataError::parse(path, "pgen variable-width header is out of range"))?;
    let mut header_end = PGEN_HEADER_LEN
        .checked_add(u64::try_from(block_offsets_len).map_err(|_| {
            MetadataError::parse(path, "pgen variable-width header is out of range")
        })?)
        .ok_or_else(|| MetadataError::parse(path, "pgen variable-width header is out of range"))?;
    for block_index in 0..block_ct {
        let block_variant_ct = block_variant_count(variant_ct, block_index);
        let type_table_len = if type_width_bits == 8 {
            block_variant_ct
        } else {
            block_variant_ct.div_ceil(2)
        };
        let length_table_len = block_variant_ct.checked_mul(length_width).ok_or_else(|| {
            MetadataError::parse(path, "pgen variable-width header is out of range")
        })?;
        let table_len = type_table_len
            .checked_add(length_table_len)
            .ok_or_else(|| {
                MetadataError::parse(path, "pgen variable-width header is out of range")
            })?;
        header_end = header_end
            .checked_add(u64::try_from(table_len).map_err(|_| {
                MetadataError::parse(path, "pgen variable-width header is out of range")
            })?)
            .ok_or_else(|| {
                MetadataError::parse(path, "pgen variable-width header is out of range")
            })?;
    }
    Ok(header_end)
}

fn read_variable_record_type_prefix(
    path: &Path,
    file: &mut File,
    type_width_bits: usize,
    block_variant_ct: usize,
    needed_in_block: usize,
    record_types: &mut Vec<u8>,
) -> Result<()> {
    if type_width_bits == 8 {
        let mut types = vec![0_u8; needed_in_block];
        file.read_exact(&mut types)
            .map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        record_types.extend(types);
        skip_bytes(path, file, block_variant_ct - needed_in_block)?;
        return Ok(());
    }

    let packed_needed = needed_in_block.div_ceil(2);
    let mut packed_types = vec![0_u8; packed_needed];
    file.read_exact(&mut packed_types)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    for variant_in_block in 0..needed_in_block {
        let byte = packed_types[variant_in_block / 2];
        let record_type = if variant_in_block % 2 == 0 {
            byte & 0x0f
        } else {
            byte >> 4
        };
        record_types.push(record_type);
    }
    // Four-bit type tables pack two variants per byte, so skipping must use
    // packed byte counts rather than raw variant counts.
    skip_bytes(path, file, block_variant_ct.div_ceil(2) - packed_needed)?;
    Ok(())
}

fn skip_bytes(path: &Path, file: &mut File, len: usize) -> Result<()> {
    let offset =
        i64::try_from(len).map_err(|_| MetadataError::parse(path, "pgen skip is out of range"))?;
    file.seek(SeekFrom::Current(offset))
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(())
}

fn block_variant_count(variant_ct: usize, block_index: usize) -> usize {
    let block_start = block_index * PGEN_VARIANT_BLOCK_SIZE;
    (variant_ct - block_start).min(PGEN_VARIANT_BLOCK_SIZE)
}

fn validate_supported_variable_record_type(path: &Path, record_type: u8) -> Result<()> {
    if record_type & 0x08 != 0 {
        return Err(MetadataError::parse(
            path,
            "unsupported pgen multiallelic hard-call patch set",
        ));
    }
    let dosage_bits = (record_type >> 5) & 0x03;
    if dosage_bits != 0 && record_type & 0x10 != 0 {
        return Err(MetadataError::parse(
            path,
            "unsupported pgen hardcall-phase track with dosage",
        ));
    }
    if record_type & 0x80 != 0 {
        return Err(MetadataError::parse(
            path,
            "unsupported pgen phased-dosage track",
        ));
    }
    match record_type & 0x07 {
        0 | 1 | 2 | 3 | 4 | 6 | 7 => Ok(()),
        compression => Err(MetadataError::parse(
            path,
            format!("unsupported pgen main-track compression type {compression}"),
        )),
    }
}

fn fixed_width_dosage_record_len(sample_ct: usize) -> usize {
    sample_ct.div_ceil(4) + sample_ct * 2
}

fn fixed_width_record_len(header: &PgenHeader) -> usize {
    match header.layout {
        PgenLayout::FixedWidth => header.bytes_per_variant,
        PgenLayout::FixedWidthDosage => fixed_width_dosage_record_len(header.sample_ct),
        PgenLayout::VariableWidth => header.bytes_per_variant,
    }
}

fn validate_fixed_width_pgen_payload_len(
    path: &Path,
    file: &File,
    variant_ct: usize,
    bytes_per_record: usize,
) -> Result<()> {
    let payload_len = variant_ct
        .checked_mul(bytes_per_record)
        .ok_or_else(|| MetadataError::parse(path, "pgen payload length is out of range"))?;
    let expected_len = PGEN_HEADER_LEN
        .checked_add(
            u64::try_from(payload_len)
                .map_err(|_| MetadataError::parse(path, "pgen payload length is out of range"))?,
        )
        .ok_or_else(|| MetadataError::parse(path, "pgen payload length is out of range"))?;
    let actual_len = file
        .metadata()
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if actual_len != expected_len {
        return Err(MetadataError::parse(
            path,
            format!("pgen payload length {actual_len} does not match fixed-width header"),
        ));
    }
    Ok(())
}

fn validate_plink2_dimensions(
    path: &Path,
    header: &PgenHeader,
    sample_ct: usize,
    variant_ct: usize,
) -> Result<()> {
    validate_plink2_sample_count(path, header, sample_ct)?;
    if header.variant_ct != variant_ct {
        return Err(MetadataError::parse(
            path,
            format!(
                "pgen variant count {} does not match pvar variant count {variant_ct}",
                header.variant_ct
            ),
        ));
    }
    Ok(())
}

fn validate_plink2_sample_count(path: &Path, header: &PgenHeader, sample_ct: usize) -> Result<()> {
    if header.sample_ct != sample_ct {
        return Err(MetadataError::parse(
            path,
            format!(
                "pgen sample count {} does not match psam sample count {sample_ct}",
                header.sample_ct
            ),
        ));
    }
    Ok(())
}

fn open_pgen_payload(path: &Path) -> Result<File> {
    let mut file = File::open(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    file.seek(SeekFrom::Start(PGEN_HEADER_LEN))
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(file)
}

fn read_plink2_variant_values(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    read_plink2_variant_packed(path, file, header, variant_index, decoder_state)?;
    decoder_state.packed.expand_selected(
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing,
    );
    Ok(())
}

fn read_plink2_variant_dosage(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    match header.layout {
        PgenLayout::FixedWidthDosage => read_fixed_width_dosage_variant_values(
            path,
            file,
            header,
            variant_index,
            source_indices,
            decoder_state,
        ),
        PgenLayout::VariableWidth => read_variable_width_dosage_variant_values(
            path,
            file,
            header,
            variant_index,
            source_indices,
            decoder_state,
        ),
        PgenLayout::FixedWidth => Err(MetadataError::parse(
            path,
            "pgen does not contain dosage values",
        )),
    }
}

fn read_plink2_variant_packed(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    match header.layout {
        PgenLayout::FixedWidth | PgenLayout::FixedWidthDosage => {
            read_fixed_width_variant_packed(path, file, header, variant_index, decoder_state)
        }
        PgenLayout::VariableWidth => {
            read_variable_width_variant_packed(path, file, header, variant_index, decoder_state)
        }
    }
}

fn read_fixed_width_variant_packed(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    seek_fixed_width_variant_record(path, file, header, variant_index)?;
    read_fixed_width_variant_packed_sequential(path, file, header, decoder_state)
}

fn seek_fixed_width_variant_record(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
) -> Result<()> {
    let payload_offset = variant_index
        .checked_mul(fixed_width_record_len(header))
        .ok_or_else(|| MetadataError::parse(path, "pgen variant offset is out of range"))?;
    let offset = PGEN_HEADER_LEN
        .checked_add(
            u64::try_from(payload_offset)
                .map_err(|_| MetadataError::parse(path, "pgen variant offset is out of range"))?,
        )
        .ok_or_else(|| MetadataError::parse(path, "pgen variant offset is out of range"))?;
    file.seek(SeekFrom::Start(offset))
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(())
}

fn read_fixed_width_variant_packed_sequential(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    decoder_state
        .record
        .resize(fixed_width_record_len(header), 0);
    file.read_exact(&mut decoder_state.record)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    decoder_state.packed.load_pgen_payload(
        &decoder_state.record[..header.bytes_per_variant],
        header.sample_ct,
    );
    Ok(())
}

fn read_variable_width_variant_packed(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    let start = header.record_offsets[variant_index];
    let end = header.record_offsets[variant_index + 1];
    let record_len = usize::try_from(
        end.checked_sub(start)
            .ok_or_else(|| MetadataError::parse(path, "pgen record length is out of range"))?,
    )
    .map_err(|_| MetadataError::parse(path, "pgen record length is out of range"))?;
    file.seek(SeekFrom::Start(start))
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    decoder_state.record.resize(record_len, 0);
    file.read_exact(&mut decoder_state.record)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    let record = decoder_state.record.as_slice();
    let record_type = header.record_types[variant_index];
    let compression = record_type & 0x07;
    match compression {
        0 => {
            if record.len() < header.bytes_per_variant {
                return Err(MetadataError::parse(
                    path,
                    "pgen uncompressed record is shorter than expected",
                ));
            }
            decoder_state
                .packed
                .load_pgen_payload(&record[..header.bytes_per_variant], header.sample_ct);
        }
        1 => decode_one_bit_record(path, record, header.sample_ct, &mut decoder_state.packed)?,
        2 | 3 => {
            if !decoder_state.has_previous_non_ld {
                return Err(MetadataError::parse(
                    path,
                    "pgen LD-compressed record appears before any non-LD record",
                ));
            }
            decode_ld_compressed_record(
                path,
                record,
                header.sample_ct,
                &decoder_state.previous_non_ld_packed,
                compression == 3,
                &mut decoder_state.packed,
            )?;
        }
        4 => decode_difflist_record(path, record, header.sample_ct, 0, &mut decoder_state.packed)?,
        6 => decode_difflist_record(path, record, header.sample_ct, 2, &mut decoder_state.packed)?,
        7 => decode_difflist_record(path, record, header.sample_ct, 3, &mut decoder_state.packed)?,
        other => {
            return Err(MetadataError::parse(
                path,
                format!("unsupported pgen main-track compression type {other}"),
            ));
        }
    }
    if decoder_state.packed.sample_ct() != header.sample_ct {
        return Err(MetadataError::parse(
            path,
            "pgen decoded category count does not match sample count",
        ));
    }
    if !matches!(compression, 2 | 3) {
        decoder_state
            .previous_non_ld_packed
            .copy_from(&decoder_state.packed);
        decoder_state.has_previous_non_ld = true;
    }
    Ok(())
}

fn read_fixed_width_dosage_variant_values(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    let record_len = fixed_width_dosage_record_len(header.sample_ct);
    let payload_offset = variant_index
        .checked_mul(record_len)
        .ok_or_else(|| MetadataError::parse(path, "pgen variant offset is out of range"))?;
    let offset = PGEN_HEADER_LEN
        .checked_add(
            u64::try_from(payload_offset)
                .map_err(|_| MetadataError::parse(path, "pgen variant offset is out of range"))?,
        )
        .ok_or_else(|| MetadataError::parse(path, "pgen variant offset is out of range"))?;
    file.seek(SeekFrom::Start(offset))
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    decoder_state.record.resize(record_len, 0);
    file.read_exact(&mut decoder_state.record)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    decoder_state.packed.load_pgen_payload(
        &decoder_state.record[..header.bytes_per_variant],
        header.sample_ct,
    );
    decoder_state.packed.expand_selected(
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing,
    );
    overlay_fixed_width_dosages(
        path,
        &decoder_state.record[header.bytes_per_variant..],
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing,
    )
}

fn read_variable_width_dosage_variant_values(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    source_indices: &[usize],
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    let start = header.record_offsets[variant_index];
    let end = header.record_offsets[variant_index + 1];
    let record_len = usize::try_from(
        end.checked_sub(start)
            .ok_or_else(|| MetadataError::parse(path, "pgen record length is out of range"))?,
    )
    .map_err(|_| MetadataError::parse(path, "pgen record length is out of range"))?;
    file.seek(SeekFrom::Start(start))
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    decoder_state.record.resize(record_len, 0);
    file.read_exact(&mut decoder_state.record)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    let record = decoder_state.record.as_slice();
    let record_type = header.record_types[variant_index];
    let dosage_bits = (record_type >> 5) & 0x03;
    if dosage_bits == 0 {
        return Err(MetadataError::parse(
            path,
            "pgen record does not contain dosage values",
        ));
    }
    let cursor = decode_variable_width_main_track(
        path,
        record,
        record_type,
        header.sample_ct,
        &decoder_state.previous_non_ld_packed,
        decoder_state.has_previous_non_ld,
        &mut decoder_state.packed,
    )?;
    if !matches!(record_type & 0x07, 2 | 3) {
        decoder_state
            .previous_non_ld_packed
            .copy_from(&decoder_state.packed);
        decoder_state.has_previous_non_ld = true;
    }

    decoder_state.packed.expand_selected(
        source_indices,
        &mut decoder_state.values,
        &mut decoder_state.missing,
    );
    overlay_variable_width_dosages(
        path,
        record,
        cursor,
        dosage_bits,
        header.sample_ct,
        DosageOverlayTarget {
            source_indices,
            values: &mut decoder_state.values,
            missing: &mut decoder_state.missing,
        },
    )
}

fn decode_variable_width_main_track(
    path: &Path,
    record: &[u8],
    record_type: u8,
    sample_ct: usize,
    previous_non_ld_packed: &PackedGenotypes,
    has_previous_non_ld: bool,
    packed: &mut PackedGenotypes,
) -> Result<usize> {
    let compression = record_type & 0x07;
    match compression {
        0 => {
            let bytes_per_variant = sample_ct.div_ceil(4);
            if record.len() < bytes_per_variant {
                return Err(MetadataError::parse(
                    path,
                    "pgen uncompressed record is shorter than expected",
                ));
            }
            packed.load_pgen_payload(&record[..bytes_per_variant], sample_ct);
            Ok(bytes_per_variant)
        }
        1 => decode_one_bit_record_with_cursor(path, record, sample_ct, packed),
        2 | 3 => {
            if !has_previous_non_ld {
                return Err(MetadataError::parse(
                    path,
                    "pgen LD-compressed record appears before any non-LD record",
                ));
            }
            let mut cursor = 0;
            let entries = decode_difflist(path, record, &mut cursor, sample_ct, true)?;
            if previous_non_ld_packed.sample_ct() != sample_ct {
                return Err(MetadataError::parse(
                    path,
                    "pgen LD state length does not match sample count",
                ));
            }
            packed.copy_from(previous_non_ld_packed);
            for (sample_index, category) in entries {
                packed.set(sample_index, category);
            }
            if compression == 3 {
                packed.invert_0_2();
            }
            Ok(cursor)
        }
        4 => decode_difflist_record_with_cursor(path, record, sample_ct, 0, packed),
        6 => decode_difflist_record_with_cursor(path, record, sample_ct, 2, packed),
        7 => decode_difflist_record_with_cursor(path, record, sample_ct, 3, packed),
        other => Err(MetadataError::parse(
            path,
            format!("unsupported pgen main-track compression type {other}"),
        )),
    }
}

fn decode_one_bit_record_with_cursor(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    packed: &mut PackedGenotypes,
) -> Result<usize> {
    let common_categories = *record.first().ok_or_else(|| {
        MetadataError::parse(path, "pgen 1-bit record is missing common-category byte")
    })?;
    let (low_category, high_category) = match common_categories {
        1 => (0, 1),
        2 => (0, 2),
        3 => (0, 3),
        5 => (1, 2),
        6 => (1, 3),
        9 => (2, 3),
        other => {
            return Err(MetadataError::parse(
                path,
                format!("invalid pgen 1-bit common-category byte {other}"),
            ));
        }
    };
    let bitarray_len = sample_ct.div_ceil(8);
    if record.len() < 1 + bitarray_len {
        return Err(MetadataError::parse(
            path,
            "pgen 1-bit record is shorter than expected",
        ));
    }
    let bitarray = &record[1..1 + bitarray_len];
    packed.resize(sample_ct);
    packed.clear_to(low_category);
    for sample_index in 0..sample_ct {
        if bit_is_set(bitarray, sample_index) {
            packed.set(sample_index, high_category);
        }
    }
    let mut cursor = 1 + bitarray_len;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    Ok(cursor)
}

fn decode_difflist_record_with_cursor(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    base_category: u8,
    packed: &mut PackedGenotypes,
) -> Result<usize> {
    packed.resize(sample_ct);
    packed.clear_to(base_category);
    let mut cursor = 0;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    Ok(cursor)
}

fn overlay_variable_width_dosages(
    path: &Path,
    record: &[u8],
    mut cursor: usize,
    dosage_bits: u8,
    sample_ct: usize,
    mut target: DosageOverlayTarget<'_>,
) -> Result<()> {
    let mut selected_samples = SelectedSampleCursor::new(target.source_indices);
    match dosage_bits {
        1 => overlay_difflist_dosages(
            path,
            record,
            &mut cursor,
            sample_ct,
            &mut selected_samples,
            &mut target,
        )?,
        2 => {
            let dosage_bytes_len = sample_ct.checked_mul(2).ok_or_else(|| {
                MetadataError::parse(path, "pgen dosage byte count is out of range")
            })?;
            ensure_record_bytes(path, record, cursor, dosage_bytes_len)?;
            for sample_index in 0..sample_ct {
                let byte_index = cursor + sample_index * 2;
                let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
                overlay_selected_pgen_dosage(sample_index, raw, &mut selected_samples, &mut target);
            }
        }
        3 => {
            let bitarray_len = sample_ct.div_ceil(8);
            ensure_record_bytes(path, record, cursor, bitarray_len)?;
            let bitarray = &record[cursor..cursor + bitarray_len];
            cursor += bitarray_len;
            let dosage_ct = (0..sample_ct)
                .filter(|sample_index| bit_is_set(bitarray, *sample_index))
                .count();
            let dosage_bytes_len = dosage_ct.checked_mul(2).ok_or_else(|| {
                MetadataError::parse(path, "pgen dosage byte count is out of range")
            })?;
            ensure_record_bytes(path, record, cursor, dosage_bytes_len)?;
            let mut dosage_index = 0;
            for sample_index in 0..sample_ct {
                if bit_is_set(bitarray, sample_index) {
                    let byte_index = cursor + dosage_index * 2;
                    let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
                    overlay_selected_pgen_dosage(
                        sample_index,
                        raw,
                        &mut selected_samples,
                        &mut target,
                    );
                    dosage_index += 1;
                }
            }
        }
        other => {
            return Err(MetadataError::parse(
                path,
                format!("unsupported pgen dosage track type {other}"),
            ));
        }
    }
    Ok(())
}

fn overlay_difflist_dosages(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    selected_samples: &mut SelectedSampleCursor<'_>,
    target: &mut DosageOverlayTarget<'_>,
) -> Result<()> {
    let list_len = read_base128_varint(path, record, cursor)?;
    if list_len == 0 {
        return Ok(());
    }
    if list_len > sample_ct {
        return Err(MetadataError::parse(
            path,
            "pgen difflist length exceeds sample count",
        ));
    }

    let group_ct = list_len.div_ceil(64);
    let sample_id_width = sample_id_width(sample_ct);
    let mut first_ids = Vec::with_capacity(group_ct);
    for _ in 0..group_ct {
        first_ids.push(read_fixed_width_sample_id(
            path,
            record,
            cursor,
            sample_id_width,
        )?);
    }
    ensure_record_bytes(path, record, *cursor, group_ct.saturating_sub(1))?;
    *cursor += group_ct.saturating_sub(1);

    let deltas_start = *cursor;
    let mut values_start = deltas_start;
    walk_difflist_ids(
        path,
        record,
        &mut values_start,
        sample_ct,
        list_len,
        &first_ids,
        |_, _| {},
    )?;

    let dosage_bytes_len = list_len
        .checked_mul(2)
        .ok_or_else(|| MetadataError::parse(path, "pgen dosage byte count is out of range"))?;
    ensure_record_bytes(path, record, values_start, dosage_bytes_len)?;

    let mut ids_cursor = deltas_start;
    // The values follow the encoded sample IDs, so we first walk IDs to find
    // values_start, then walk them again while overlaying selected samples.
    walk_difflist_ids(
        path,
        record,
        &mut ids_cursor,
        sample_ct,
        list_len,
        &first_ids,
        |sample_index, dosage_index| {
            let byte_index = values_start + dosage_index * 2;
            let raw = u16::from_le_bytes([record[byte_index], record[byte_index + 1]]);
            overlay_selected_pgen_dosage(sample_index, raw, selected_samples, target);
        },
    )?;
    *cursor = values_start + dosage_bytes_len;
    Ok(())
}

fn walk_difflist_ids(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    list_len: usize,
    first_ids: &[usize],
    mut visit: impl FnMut(usize, usize),
) -> Result<()> {
    let mut previous_sample_id = None;
    let mut entry_index = 0;
    for (group_index, first_id) in first_ids.iter().copied().enumerate() {
        let group_len = (list_len - group_index * 64).min(64);
        let mut sample_id = first_id;
        validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
        visit(sample_id, entry_index);
        entry_index += 1;
        for _ in 1..group_len {
            let delta = read_base128_varint(path, record, cursor)?;
            sample_id = sample_id.checked_add(delta).ok_or_else(|| {
                MetadataError::parse(path, "pgen difflist sample id is out of range")
            })?;
            validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
            visit(sample_id, entry_index);
            entry_index += 1;
        }
    }
    Ok(())
}

fn overlay_selected_pgen_dosage(
    source_index: usize,
    raw: u16,
    selected_samples: &mut SelectedSampleCursor<'_>,
    target: &mut DosageOverlayTarget<'_>,
) {
    if let Some(selected_index) = selected_samples.selected_index_for(source_index) {
        apply_pgen_dosage(
            raw,
            &mut target.values[selected_index],
            &mut target.missing[selected_index],
        );
    }
}

fn overlay_fixed_width_dosages(
    path: &Path,
    dosage_bytes: &[u8],
    source_indices: &[usize],
    values: &mut [f32],
    missing: &mut [bool],
) -> Result<()> {
    for (selected_index, source_index) in source_indices.iter().copied().enumerate() {
        let byte_index = source_index
            .checked_mul(2)
            .ok_or_else(|| MetadataError::parse(path, "pgen dosage offset is out of range"))?;
        ensure_record_bytes(path, dosage_bytes, byte_index, 2)?;
        let raw = u16::from_le_bytes([dosage_bytes[byte_index], dosage_bytes[byte_index + 1]]);
        apply_pgen_dosage(
            raw,
            &mut values[selected_index],
            &mut missing[selected_index],
        );
    }
    Ok(())
}

fn apply_pgen_dosage(raw: u16, value: &mut f32, is_missing: &mut bool) {
    if raw == u16::MAX {
        *value = 0.0;
        *is_missing = true;
        return;
    }
    *value = f32::from(raw) * (2.0 / 32768.0);
    *is_missing = false;
}

fn decode_one_bit_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    let common_categories = *record.first().ok_or_else(|| {
        MetadataError::parse(path, "pgen 1-bit record is missing common-category byte")
    })?;
    let (low_category, high_category) = match common_categories {
        1 => (0, 1),
        2 => (0, 2),
        3 => (0, 3),
        5 => (1, 2),
        6 => (1, 3),
        9 => (2, 3),
        other => {
            return Err(MetadataError::parse(
                path,
                format!("invalid pgen 1-bit common-category byte {other}"),
            ));
        }
    };
    let bitarray_len = sample_ct.div_ceil(8);
    if record.len() < 1 + bitarray_len {
        return Err(MetadataError::parse(
            path,
            "pgen 1-bit record is shorter than expected",
        ));
    }
    let bitarray = &record[1..1 + bitarray_len];
    packed.resize(sample_ct);
    packed.clear_to(low_category);
    for sample_index in 0..sample_ct {
        if bit_is_set(bitarray, sample_index) {
            packed.set(sample_index, high_category);
        }
    }
    let mut cursor = 1 + bitarray_len;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    Ok(())
}

fn decode_ld_compressed_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    previous_non_ld_packed: &PackedGenotypes,
    inverted: bool,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    if previous_non_ld_packed.sample_ct() != sample_ct {
        return Err(MetadataError::parse(
            path,
            "pgen LD state length does not match sample count",
        ));
    }
    packed.copy_from(previous_non_ld_packed);
    let mut cursor = 0;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    if inverted {
        packed.invert_0_2();
    }
    Ok(())
}

fn decode_difflist_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    base_category: u8,
    packed: &mut PackedGenotypes,
) -> Result<()> {
    packed.resize(sample_ct);
    packed.clear_to(base_category);
    let mut cursor = 0;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        packed.set(sample_index, category);
    }
    Ok(())
}

fn decode_difflist(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    sample_ct: usize,
    with_values: bool,
) -> Result<Vec<(usize, u8)>> {
    let list_len = read_base128_varint(path, record, cursor)?;
    if list_len == 0 {
        return Ok(Vec::new());
    }
    if list_len > sample_ct {
        return Err(MetadataError::parse(
            path,
            "pgen difflist length exceeds sample count",
        ));
    }
    let group_ct = list_len.div_ceil(64);
    let sample_id_width = sample_id_width(sample_ct);
    let mut first_ids = Vec::with_capacity(group_ct);
    for _ in 0..group_ct {
        first_ids.push(read_fixed_width_sample_id(
            path,
            record,
            cursor,
            sample_id_width,
        )?);
    }
    ensure_record_bytes(path, record, *cursor, group_ct.saturating_sub(1))?;
    *cursor += group_ct.saturating_sub(1);

    let packed_values_start = *cursor;
    if with_values {
        ensure_record_bytes(path, record, *cursor, list_len.div_ceil(4))?;
        *cursor += list_len.div_ceil(4);
    }

    let mut entries = Vec::with_capacity(list_len);
    let mut previous_sample_id = None;
    for (group_index, first_id) in first_ids.into_iter().enumerate() {
        let group_len = (list_len - group_index * 64).min(64);
        let mut sample_id = first_id;
        validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
        entries.push((
            sample_id,
            packed_difflist_value(record, packed_values_start, entries.len(), with_values),
        ));
        for _ in 1..group_len {
            let delta = read_base128_varint(path, record, cursor)?;
            sample_id = sample_id.checked_add(delta).ok_or_else(|| {
                MetadataError::parse(path, "pgen difflist sample id is out of range")
            })?;
            validate_difflist_sample_id(path, sample_id, sample_ct, &mut previous_sample_id)?;
            entries.push((
                sample_id,
                packed_difflist_value(record, packed_values_start, entries.len(), with_values),
            ));
        }
    }
    Ok(entries)
}

fn read_base128_varint(path: &Path, record: &[u8], cursor: &mut usize) -> Result<usize> {
    let mut value = 0_usize;
    let mut shift = 0_u32;
    loop {
        ensure_record_bytes(path, record, *cursor, 1)?;
        let byte = record[*cursor];
        *cursor += 1;
        value |= usize::from(byte & 0x7f)
            .checked_shl(shift)
            .ok_or_else(|| MetadataError::parse(path, "pgen varint is out of range"))?;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
        if shift >= usize::BITS {
            return Err(MetadataError::parse(path, "pgen varint is out of range"));
        }
    }
}

fn read_fixed_width_sample_id(
    path: &Path,
    record: &[u8],
    cursor: &mut usize,
    width: usize,
) -> Result<usize> {
    ensure_record_bytes(path, record, *cursor, width)?;
    let mut value = 0_usize;
    for byte_index in 0..width {
        value |= usize::from(record[*cursor + byte_index]) << (8 * byte_index);
    }
    *cursor += width;
    Ok(value)
}

fn sample_id_width(sample_ct: usize) -> usize {
    if sample_ct <= 1 << 8 {
        1
    } else if sample_ct <= 1 << 16 {
        2
    } else if sample_ct <= 1 << 24 {
        3
    } else {
        4
    }
}

fn packed_difflist_value(record: &[u8], start: usize, index: usize, with_values: bool) -> u8 {
    if !with_values {
        return 0;
    }
    (record[start + index / 4] >> ((index % 4) * 2)) & 0b11
}

fn validate_difflist_sample_id(
    path: &Path,
    sample_id: usize,
    sample_ct: usize,
    previous_sample_id: &mut Option<usize>,
) -> Result<()> {
    if sample_id >= sample_ct {
        return Err(MetadataError::parse(
            path,
            "pgen difflist sample id is outside sample count",
        ));
    }
    if previous_sample_id.is_some_and(|previous| sample_id <= previous) {
        return Err(MetadataError::parse(
            path,
            "pgen difflist sample ids must be strictly increasing",
        ));
    }
    *previous_sample_id = Some(sample_id);
    Ok(())
}

fn ensure_record_bytes(path: &Path, record: &[u8], cursor: usize, len: usize) -> Result<()> {
    if cursor.checked_add(len).is_none_or(|end| end > record.len()) {
        return Err(MetadataError::parse(
            path,
            "pgen record ended before expected data",
        ));
    }
    Ok(())
}

fn bit_is_set(bytes: &[u8], bit_index: usize) -> bool {
    bytes[bit_index / 8] & (1 << (bit_index % 8)) != 0
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

    #[test]
    fn variable_record_helpers_write_packed_genotypes() {
        let path = Path::new("test.pgen");
        let mut packed = PackedGenotypes::default();

        decode_one_bit_record(path, &[2, 0b0000_1010, 0], 4, &mut packed)
            .expect("one-bit record should decode");
        assert_eq!(
            (0..4)
                .map(|sample_index| packed.get(sample_index))
                .collect::<Vec<_>>(),
            vec![0, 2, 0, 2]
        );

        decode_difflist_record(path, &[2, 1, 9, 2], 4, 0, &mut packed)
            .expect("difflist record should decode");
        assert_eq!(
            (0..4)
                .map(|sample_index| packed.get(sample_index))
                .collect::<Vec<_>>(),
            vec![0, 1, 0, 2]
        );

        let mut previous = PackedGenotypes::default();
        previous.resize(4);
        previous.clear_to(0);
        previous.set(1, 1);
        previous.set(2, 2);
        previous.set(3, 3);

        decode_ld_compressed_record(path, &[1, 2, 0], 4, &previous, true, &mut packed)
            .expect("LD-compressed record should decode");
        assert_eq!(
            (0..4)
                .map(|sample_index| packed.get(sample_index))
                .collect::<Vec<_>>(),
            vec![2, 1, 2, 3]
        );
    }

    #[test]
    fn variable_width_dosage_overlay_preserves_hardcall_inferred_values() {
        let path = Path::new("test.pgen");
        let mut record = vec![0b0000_0110];
        record.extend(100_u16.to_le_bytes());
        record.extend(200_u16.to_le_bytes());
        let mut values = vec![0.0, 2.0, 0.0];
        let mut missing = vec![false, false, true];

        overlay_variable_width_dosages(
            path,
            &record,
            0,
            3,
            4,
            DosageOverlayTarget {
                source_indices: &[0, 2, 3],
                values: &mut values,
                missing: &mut missing,
            },
        )
        .expect("dosage overlay should decode");

        assert_eq!(values, vec![0.0, f32::from(200_u16) * (2.0 / 32768.0), 0.0]);
        assert_eq!(missing, vec![false, false, true]);
    }

    #[test]
    fn variable_width_dosage_list_overlay_uses_source_order_without_dense_index() {
        let path = Path::new("test.pgen");
        let mut record = vec![3, 1, 3, 5];
        record.extend(100_u16.to_le_bytes());
        record.extend(200_u16.to_le_bytes());
        record.extend(300_u16.to_le_bytes());
        let mut values = vec![0.0, 1.0];
        let mut missing = vec![false, false];

        overlay_variable_width_dosages(
            path,
            &record,
            0,
            1,
            10,
            DosageOverlayTarget {
                source_indices: &[4, 9],
                values: &mut values,
                missing: &mut missing,
            },
        )
        .expect("dosage-list overlay should decode");

        assert_eq!(
            values,
            vec![
                f32::from(200_u16) * (2.0 / 32768.0),
                f32::from(300_u16) * (2.0 / 32768.0),
            ]
        );
        assert_eq!(missing, vec![false, false]);
    }
}

fn parse_psam(path: &Path) -> Result<Vec<SampleRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = None;
    let mut records = Vec::new();
    for (line_index, line) in contents.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            header = Some(parse_psam_header(trimmed));
            continue;
        }
        let columns = header
            .as_ref()
            .ok_or_else(|| MetadataError::parse(path, "psam header line is required"))?;
        records.push(parse_psam_line(path, line_index + 1, columns, trimmed)?);
    }
    Ok(records)
}

#[derive(Debug, Clone, Copy)]
struct PsamColumns {
    fid: Option<usize>,
    iid: usize,
    father: Option<usize>,
    mother: Option<usize>,
    sex: Option<usize>,
    phenotype: Option<usize>,
}

fn parse_psam_header(line: &str) -> PsamColumns {
    let fields = line
        .trim_start_matches('#')
        .split_whitespace()
        .collect::<Vec<_>>();
    let find = |names: &[&str]| {
        fields
            .iter()
            .position(|field| names.iter().any(|name| field.eq_ignore_ascii_case(name)))
    };
    PsamColumns {
        fid: find(&["FID"]),
        iid: find(&["IID"]).unwrap_or(0),
        father: find(&["PAT", "FATHER"]),
        mother: find(&["MAT", "MOTHER"]),
        sex: find(&["SEX"]),
        phenotype: find(&["PHENO1", "PHENO", "PHENOTYPE"]),
    }
}

fn parse_psam_line(
    path: &Path,
    line_number: usize,
    columns: &PsamColumns,
    line: &str,
) -> Result<SampleRecord> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let required = columns
        .iid
        .max(columns.fid.unwrap_or(0))
        .max(columns.father.unwrap_or(0))
        .max(columns.mother.unwrap_or(0))
        .max(columns.sex.unwrap_or(0))
        .max(columns.phenotype.unwrap_or(0));
    if fields.len() <= required {
        return Err(MetadataError::parse(
            path,
            format!("psam line {line_number} has too few fields"),
        ));
    }
    Ok(SampleRecord {
        fid: columns
            .fid
            .and_then(|index| optional_plink_value(fields[index])),
        iid: fields[columns.iid].to_string(),
        father: columns
            .father
            .and_then(|index| optional_plink_value(fields[index])),
        mother: columns
            .mother
            .and_then(|index| optional_plink_value(fields[index])),
        sex: columns.sex.map(|index| fields[index].to_string()),
        phenotype: columns.phenotype.map(|index| fields[index].to_string()),
        source_sample_index: None,
        haplotype_index: None,
    })
}

fn parse_pvar(path: &Path) -> Result<Vec<VariantRecord>> {
    let mut reader = open_pvar_reader(path)?;
    let mut contents = String::new();
    reader
        .read_to_string(&mut contents)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    let data_lines = contents
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty() && !line.trim_start().starts_with("##"))
        .collect::<Vec<_>>();
    let header_index = data_lines
        .iter()
        .rposition(|(_, line)| line.trim_start().starts_with("#CHROM"));
    let (columns, body_start) = if let Some(header_index) = header_index {
        (
            parse_pvar_header(data_lines[header_index].1)?,
            header_index + 1,
        )
    } else {
        infer_pvar_header(path, data_lines.first().map(|(_, line)| *line))?
    };
    data_lines
        .into_iter()
        .skip(body_start)
        .map(|(index, line)| parse_pvar_line(path, index + 1, &columns, line))
        .collect()
}

fn parse_pvar_source_window(
    path: &Path,
    window: VariantWindow,
    expected_variant_ct: usize,
) -> Result<Vec<(usize, VariantRecord)>> {
    let reader = open_pvar_reader(path)?;
    let mut columns = None;
    let mut body_started = false;
    let mut source_index = 0_usize;
    let window_end = window.start.saturating_add(window.len);
    let mut records = Vec::with_capacity(window.len);

    for (line_index, line_result) in reader.lines().enumerate() {
        let line = line_result.map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with("##") {
            continue;
        }
        if trimmed.starts_with("#CHROM") {
            columns = Some(parse_pvar_header(trimmed)?);
            body_started = true;
            continue;
        }
        if !body_started {
            // Headerless PVAR has PLINK-specific default columns. Infer once
            // from the first data row, then use the same parser as headered
            // PVAR rows.
            let (inferred, _) = infer_pvar_header(path, Some(trimmed))?;
            columns = Some(inferred);
            body_started = true;
        }

        let columns = columns
            .as_ref()
            .expect("pvar columns should be initialized before parsing body rows");
        let variant = parse_pvar_line(path, line_index + 1, columns, trimmed)?;
        if source_index >= window.start && source_index < window_end {
            records.push((source_index, variant));
        }
        source_index += 1;
    }

    if source_index != expected_variant_ct {
        return Err(MetadataError::parse(
            path,
            format!(
                "pvar variant count {source_index} does not match pgen variant count {expected_variant_ct}",
            ),
        ));
    }

    Ok(records)
}

struct PvarRecordReader {
    path: std::path::PathBuf,
    lines: std::iter::Enumerate<std::io::Lines<Box<dyn BufRead>>>,
    columns: Option<PvarColumns>,
    body_started: bool,
    source_index: usize,
}

impl PvarRecordReader {
    fn new(path: &Path) -> Result<Self> {
        Ok(Self {
            path: path.to_path_buf(),
            lines: open_pvar_reader(path)?.lines().enumerate(),
            columns: None,
            body_started: false,
            source_index: 0,
        })
    }

    fn next_record(&mut self) -> Result<Option<(usize, VariantRecord)>> {
        for (line_index, line_result) in self.lines.by_ref() {
            let line = line_result.map_err(|source| MetadataError::Io {
                path: self.path.clone(),
                source,
            })?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with("##") {
                continue;
            }
            if trimmed.starts_with("#CHROM") {
                self.columns = Some(parse_pvar_header(trimmed)?);
                self.body_started = true;
                continue;
            }
            if !self.body_started {
                let (inferred, _) = infer_pvar_header(&self.path, Some(trimmed))?;
                self.columns = Some(inferred);
                self.body_started = true;
            }

            let columns = self
                .columns
                .as_ref()
                .expect("pvar columns should be initialized before parsing body rows");
            let variant = parse_pvar_line(&self.path, line_index + 1, columns, trimmed)?;
            let source_index = self.source_index;
            self.source_index += 1;
            return Ok(Some((source_index, variant)));
        }
        Ok(None)
    }

    fn validate_count(&self, expected_variant_ct: usize) -> Result<()> {
        if self.source_index != expected_variant_ct {
            return Err(MetadataError::parse(
                &self.path,
                format!(
                    "pvar variant count {} does not match pgen variant count {expected_variant_ct}",
                    self.source_index
                ),
            ));
        }
        Ok(())
    }
}

fn open_pvar_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    if path
        .file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.ends_with(".pvar.zst"))
    {
        let decoder =
            zstd::stream::read::Decoder::new(file).map_err(|source| MetadataError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        return Ok(Box::new(BufReader::new(decoder)));
    }
    Ok(Box::new(BufReader::new(file)))
}

#[derive(Debug, Clone, Copy)]
struct PvarColumns {
    chrom: usize,
    pos: usize,
    id: usize,
    ref_allele: usize,
    alt_allele: usize,
    qual: Option<usize>,
}

fn parse_pvar_header(line: &str) -> Result<PvarColumns> {
    let fields = line
        .trim_start_matches('#')
        .split_whitespace()
        .take_while(|field| !field.eq_ignore_ascii_case("FORMAT"))
        .collect::<Vec<_>>();
    let find = |name: &str| {
        fields
            .iter()
            .position(|field| field.eq_ignore_ascii_case(name))
    };
    Ok(PvarColumns {
        chrom: find("CHROM")
            .ok_or_else(|| MetadataError::parse("<pvar>", "pvar header missing #CHROM"))?,
        pos: find("POS")
            .ok_or_else(|| MetadataError::parse("<pvar>", "pvar header missing POS"))?,
        id: find("ID").ok_or_else(|| MetadataError::parse("<pvar>", "pvar header missing ID"))?,
        ref_allele: find("REF")
            .ok_or_else(|| MetadataError::parse("<pvar>", "pvar header missing REF"))?,
        alt_allele: find("ALT")
            .ok_or_else(|| MetadataError::parse("<pvar>", "pvar header missing ALT"))?,
        qual: find("QUAL"),
    })
}

fn infer_pvar_header(path: &Path, first_data_line: Option<&str>) -> Result<(PvarColumns, usize)> {
    let Some(line) = first_data_line else {
        return Ok((
            PvarColumns {
                chrom: 0,
                id: 1,
                pos: 2,
                alt_allele: 3,
                ref_allele: 4,
                qual: None,
            },
            0,
        ));
    };
    let field_count = line.split_whitespace().count();
    match field_count {
        5 => Ok((
            PvarColumns {
                chrom: 0,
                id: 1,
                pos: 2,
                alt_allele: 3,
                ref_allele: 4,
                qual: None,
            },
            0,
        )),
        count if count >= 6 => Ok((
            PvarColumns {
                chrom: 0,
                id: 1,
                pos: 3,
                alt_allele: 4,
                ref_allele: 5,
                qual: None,
            },
            0,
        )),
        _ => Err(MetadataError::parse(
            path,
            "pvar without header must have at least five columns",
        )),
    }
}

fn parse_pvar_line(
    path: &Path,
    line_number: usize,
    columns: &PvarColumns,
    line: &str,
) -> Result<VariantRecord> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let required = columns
        .chrom
        .max(columns.pos)
        .max(columns.id)
        .max(columns.ref_allele)
        .max(columns.alt_allele)
        .max(columns.qual.unwrap_or(0));
    if fields.len() <= required {
        return Err(MetadataError::parse(
            path,
            format!("pvar line {line_number} has too few fields"),
        ));
    }
    let pos = fields[columns.pos].parse::<u32>().map_err(|error| {
        MetadataError::parse(
            path,
            format!("pvar line {line_number} has invalid position: {error}"),
        )
    })?;
    let ref_allele = fields[columns.ref_allele].to_string();
    let alt_allele = fields[columns.alt_allele].to_string();
    let first_alt = alt_allele.split(',').next().unwrap_or("").to_string();
    if first_alt.is_empty() {
        return Err(MetadataError::parse(
            path,
            format!("pvar line {line_number} has empty ALT allele"),
        ));
    }
    let qual = columns
        .qual
        .map(|index| parse_optional_qual(path, line_number, fields[index]))
        .transpose()?
        .flatten();

    Ok(VariantRecord {
        chrom: fields[columns.chrom].to_string(),
        pos,
        id: fields[columns.id].to_string(),
        a0: ref_allele.clone(),
        a1: first_alt.clone(),
        ref_allele: Some(ref_allele.clone()),
        alt_allele: Some(alt_allele),
        source_a0: ref_allele,
        source_a1: first_alt,
        flipped: false,
        qual,
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
}

fn parse_optional_qual(path: &Path, line_number: usize, value: &str) -> Result<Option<f32>> {
    if value == "." {
        return Ok(None);
    }
    let qual = value.parse::<f32>().map_err(|error| {
        MetadataError::parse(
            path,
            format!("pvar line {line_number} has invalid QUAL: {error}"),
        )
    })?;
    if qual.is_finite() {
        Ok(Some(qual))
    } else {
        Ok(None)
    }
}

fn optional_plink_value(value: &str) -> Option<String> {
    if value == "0" || value == "." || value == "NA" {
        None
    } else {
        Some(value.to_string())
    }
}
