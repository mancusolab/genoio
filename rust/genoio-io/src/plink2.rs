// pattern: Imperative Shell

use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, compute_variant_stats, flip_values_to_minor_allele,
    reject_sparse_missing_values, select_samples_source_order,
    transpose_variant_major_to_sample_major, DenseGenotypeMatrix, MetadataError, MetadataOutput,
    PartialFilterDecision, SampleRecord, SourceCapabilities, SparseGenotypeMatrix, VariantFilter,
    VariantRecord, VariantWindow,
};

use crate::error::Result;

const PGEN_MAGIC: [u8; 2] = [0x6c, 0x1b];
const PGEN_MODE_FIXED_WIDTH_HARDCALLS: u8 = 0x02;
const PGEN_MODE_VARIABLE_WIDTH: u8 = 0x10;
const PGEN_HEADER_LEN: u64 = 12;
const PGEN_VARIANT_BLOCK_SIZE: usize = 65_536;

#[derive(Debug, Clone)]
struct PgenHeader {
    layout: PgenLayout,
    variant_ct: usize,
    sample_ct: usize,
    bytes_per_variant: usize,
    record_types: Vec<u8>,
    record_offsets: Vec<u64>,
}

#[derive(Debug, Clone, Default)]
struct PackedGenotypes {
    words: Vec<u64>,
    sample_ct: usize,
}

impl PackedGenotypes {
    const SAMPLES_PER_WORD: usize = 32;
    const BITS_PER_SAMPLE: usize = 2;

    fn resize(&mut self, sample_ct: usize) {
        self.sample_ct = sample_ct;
        self.words.resize(sample_ct.div_ceil(Self::SAMPLES_PER_WORD), 0);
        self.mask_unused_slots();
    }

    fn clear_to(&mut self, category: u8) {
        let category = u64::from(category & 0b11);
        let mut word = 0_u64;
        for slot_index in 0..Self::SAMPLES_PER_WORD {
            word |= category << (slot_index * Self::BITS_PER_SAMPLE);
        }
        self.words.fill(word);
        self.mask_unused_slots();
    }

    fn set(&mut self, sample_index: usize, category: u8) {
        debug_assert!(sample_index < self.sample_ct);
        let word_index = sample_index / Self::SAMPLES_PER_WORD;
        let shift = (sample_index % Self::SAMPLES_PER_WORD) * Self::BITS_PER_SAMPLE;
        let mask = 0b11_u64 << shift;
        self.words[word_index] =
            (self.words[word_index] & !mask) | (u64::from(category & 0b11) << shift);
    }

    fn get(&self, sample_index: usize) -> u8 {
        debug_assert!(sample_index < self.sample_ct);
        let word_index = sample_index / Self::SAMPLES_PER_WORD;
        let shift = (sample_index % Self::SAMPLES_PER_WORD) * Self::BITS_PER_SAMPLE;
        ((self.words[word_index] >> shift) & 0b11) as u8
    }

    fn copy_from(&mut self, other: &Self) {
        self.sample_ct = other.sample_ct;
        self.words.clear();
        self.words.extend_from_slice(&other.words);
    }

    fn invert_0_2(&mut self) {
        for sample_index in 0..self.sample_ct {
            match self.get(sample_index) {
                0 => self.set(sample_index, 2),
                2 => self.set(sample_index, 0),
                _ => {}
            }
        }
    }

    fn expand_selected(
        &self,
        source_indices: &[usize],
        values: &mut Vec<f32>,
        missing: &mut Vec<bool>,
    ) {
        values.clear();
        missing.clear();
        for source_index in source_indices {
            let (value, is_missing) = decode_pgen_code(self.get(*source_index));
            values.push(value);
            missing.push(is_missing);
        }
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
struct PgenDecoderState {
    previous_non_ld_categories: Vec<u8>,
    has_previous_non_ld: bool,
    record: Vec<u8>,
    categories: Vec<u8>,
    values: Vec<f32>,
    missing: Vec<bool>,
}

impl PgenDecoderState {
    fn new(sample_ct: usize, selected_sample_ct: usize) -> Self {
        Self {
            previous_non_ld_categories: Vec::with_capacity(sample_ct),
            has_previous_non_ld: false,
            record: Vec::with_capacity(sample_ct.div_ceil(4)),
            categories: Vec::with_capacity(sample_ct),
            values: Vec::with_capacity(selected_sample_ct),
            missing: Vec::with_capacity(selected_sample_ct),
        }
    }
}

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

#[derive(Debug, Clone)]
enum PgenLayout {
    FixedWidth,
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
    let source_variants = parse_pvar(pvar)?;
    validate_plink2_dimensions(pgen, &header, all_samples.len(), source_variants.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::with_capacity(selection.samples.len() * header.variant_ct);
    let mut variant_major_missing = Vec::with_capacity(selection.samples.len() * header.variant_ct);
    let mut retained_index = 0_usize;
    let requires_sequential_decode = matches!(header.layout, PgenLayout::VariableWidth);
    for (variant_index, mut variant) in source_variants.into_iter().enumerate() {
        diagnostics.candidate_variants += 1;
        let mut decoded_values = false;
        if requires_sequential_decode {
            read_plink2_variant_values(
                pgen,
                &mut file,
                &header,
                variant_index,
                &selection.source_indices,
                &mut decoder_state,
            )?;
            decoded_values = true;
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
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {}
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        if !decoded_values {
            read_plink2_variant_values(
                pgen,
                &mut file,
                &header,
                variant_index,
                &selection.source_indices,
                &mut decoder_state,
            )?;
        }
        let stats = if needs_genotype_decision {
            Some(compute_variant_stats(
                &decoder_state.values,
                &decoder_state.missing,
            )?)
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
                continue;
            }
        }
        variants.push(variant);
        variant_major_values.extend_from_slice(&decoder_state.values);
        variant_major_missing.extend_from_slice(&decoder_state.missing);
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
    let source_variants = parse_pvar(pvar)?;
    validate_plink2_dimensions(pgen, &header, all_samples.len(), source_variants.len())?;
    let selection = select_samples_source_order(&all_samples, requested_samples, pgen)?;
    let mut diagnostics = selection.diagnostics;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();
    let mut retained_index = 0_usize;
    let requires_sequential_decode = matches!(header.layout, PgenLayout::VariableWidth);
    for (variant_index, mut variant) in source_variants.into_iter().enumerate() {
        diagnostics.candidate_variants += 1;
        let mut decoded_values = false;
        if requires_sequential_decode {
            read_plink2_variant_values(
                pgen,
                &mut file,
                &header,
                variant_index,
                &selection.source_indices,
                &mut decoder_state,
            )?;
            decoded_values = true;
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
                    continue;
                }
            }
            PartialFilterDecision::NeedGenotypes => {}
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        if !decoded_values {
            read_plink2_variant_values(
                pgen,
                &mut file,
                &header,
                variant_index,
                &selection.source_indices,
                &mut decoder_state,
            )?;
        }
        let stats = if needs_genotype_decision {
            Some(compute_variant_stats(
                &decoder_state.values,
                &decoder_state.missing,
            )?)
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
                continue;
            }
        }
        reject_sparse_missing_values(&decoder_state.missing)?;
        flip_values_to_minor_allele(&mut decoder_state.values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &decoder_state.values);
        variants.push(variant);
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

    // Unfiltered source windows know their final retained width up front, so
    // construct the public sample-major buffers directly.
    match header.layout {
        PgenLayout::FixedWidth => {
            if n_variants > 0 {
                seek_fixed_width_variant_record(pgen, &mut file, &header, window.start)?;
            }
            for variant_index in 0..n_variants {
                read_fixed_width_variant_categories_sequential(
                    pgen,
                    &mut file,
                    &header,
                    &mut decoder_state,
                )?;
                fill_plink2_variant_values(&source_indices, &mut decoder_state);
                append_variant_to_sample_major(
                    &decoder_state.values,
                    &decoder_state.missing,
                    variant_index,
                    n_variants,
                    &mut values,
                    &mut missing_mask,
                );
            }
        }
        PgenLayout::VariableWidth => {
            let prefix_end = window.start + n_variants;
            for source_variant_index in 0..prefix_end {
                read_plink2_variant_values(
                    pgen,
                    &mut file,
                    &header,
                    source_variant_index,
                    &source_indices,
                    &mut decoder_state,
                )?;
                if source_variant_index >= window.start {
                    append_variant_to_sample_major(
                        &decoder_state.values,
                        &decoder_state.missing,
                        source_variant_index - window.start,
                        n_variants,
                        &mut values,
                        &mut missing_mask,
                    );
                }
            }
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
    let window_variants = parse_pvar_source_window(pvar, window)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let n_variants = window_variants.len();
    let mut variants = Vec::new();
    let mut values = vec![0.0; n_samples * n_variants];
    let mut missing_mask = vec![false; n_samples * n_variants];

    // This metadata-bearing source-window path uses the same direct
    // sample-major construction as matrix-only windows.
    match header.layout {
        PgenLayout::FixedWidth => {
            if let Some((first_variant_index, _)) = window_variants.first() {
                seek_fixed_width_variant_record(pgen, &mut file, &header, *first_variant_index)?;
            }
            for (window_variant_index, (source_variant_index, variant)) in
                window_variants.into_iter().enumerate()
            {
                debug_assert!(source_variant_index < header.variant_ct);
                read_fixed_width_variant_categories_sequential(
                    pgen,
                    &mut file,
                    &header,
                    &mut decoder_state,
                )?;
                fill_plink2_variant_values(&selection.source_indices, &mut decoder_state);
                variants.push(variant);
                append_variant_to_sample_major(
                    &decoder_state.values,
                    &decoder_state.missing,
                    window_variant_index,
                    n_variants,
                    &mut values,
                    &mut missing_mask,
                );
            }
        }
        PgenLayout::VariableWidth => {
            let mut window_iter = window_variants.into_iter().peekable();
            let prefix_end = header.record_types.len();
            // Variable-width PGEN can use LD-compressed records that depend on
            // earlier non-LD records. Decode the prefix through the requested
            // window to maintain state, but append only requested variants.
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
                    let (_, variant) = window_iter.next().expect("peeked variant should exist");
                    let window_variant_index = variants.len();
                    variants.push(variant);
                    append_variant_to_sample_major(
                        &decoder_state.values,
                        &decoder_state.missing,
                        window_variant_index,
                        n_variants,
                        &mut values,
                        &mut missing_mask,
                    );
                }
            }
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
    let window_variants = parse_pvar_source_window(pvar, window)?;
    let mut file = open_pgen_payload(pgen)?;
    let mut decoder_state = PgenDecoderState::new(header.sample_ct, selection.samples.len());

    let n_samples = selection.samples.len();
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();

    match header.layout {
        PgenLayout::FixedWidth => {
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
                "unsupported pgen mode 0x{mode:02x}; only fixed-width and variable-width biallelic hardcalls are supported"
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
            record_offsets.push(block_offset);
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
    Ok((record_types, record_offsets))
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
    if (record_type >> 5) & 0x03 != 0 {
        return Err(MetadataError::parse(path, "unsupported pgen dosage track"));
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

fn validate_fixed_width_pgen_payload_len(
    path: &Path,
    file: &File,
    variant_ct: usize,
    bytes_per_variant: usize,
) -> Result<()> {
    let payload_len = variant_ct
        .checked_mul(bytes_per_variant)
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
    read_plink2_variant_categories(path, file, header, variant_index, decoder_state)?;

    fill_plink2_variant_values(source_indices, decoder_state);
    Ok(())
}

fn fill_plink2_variant_values(source_indices: &[usize], decoder_state: &mut PgenDecoderState) {
    decoder_state.values.clear();
    decoder_state.missing.clear();
    for source_index in source_indices {
        let (value, missing) = decode_pgen_code(decoder_state.categories[*source_index]);
        decoder_state.values.push(value);
        decoder_state.missing.push(missing);
    }
}

fn read_plink2_variant_categories(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    match header.layout {
        PgenLayout::FixedWidth => {
            read_fixed_width_variant_categories(path, file, header, variant_index, decoder_state)
        }
        PgenLayout::VariableWidth => {
            read_variable_width_variant_categories(path, file, header, variant_index, decoder_state)
        }
    }
}

fn read_fixed_width_variant_categories(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    seek_fixed_width_variant_record(path, file, header, variant_index)?;
    read_fixed_width_variant_categories_sequential(path, file, header, decoder_state)
}

fn seek_fixed_width_variant_record(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
) -> Result<()> {
    let payload_offset = variant_index
        .checked_mul(header.bytes_per_variant)
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

fn read_fixed_width_variant_categories_sequential(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    decoder_state.record.resize(header.bytes_per_variant, 0);
    file.read_exact(&mut decoder_state.record)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    decode_packed_categories(
        &decoder_state.record,
        header.sample_ct,
        &mut decoder_state.categories,
    );
    Ok(())
}

fn read_variable_width_variant_categories(
    path: &Path,
    file: &mut File,
    header: &PgenHeader,
    variant_index: usize,
    decoder_state: &mut PgenDecoderState,
) -> Result<()> {
    let start = header.record_offsets[variant_index];
    let end = header.record_offsets[variant_index + 1];
    let record_len = usize::try_from(end - start)
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
            decode_packed_categories(
                &record[..header.bytes_per_variant],
                header.sample_ct,
                &mut decoder_state.categories,
            );
        }
        1 => decode_one_bit_record(
            path,
            record,
            header.sample_ct,
            &mut decoder_state.categories,
        )?,
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
                &decoder_state.previous_non_ld_categories,
                compression == 3,
                &mut decoder_state.categories,
            )?;
        }
        4 => decode_difflist_record(
            path,
            record,
            header.sample_ct,
            0,
            &mut decoder_state.categories,
        )?,
        6 => decode_difflist_record(
            path,
            record,
            header.sample_ct,
            2,
            &mut decoder_state.categories,
        )?,
        7 => decode_difflist_record(
            path,
            record,
            header.sample_ct,
            3,
            &mut decoder_state.categories,
        )?,
        other => {
            return Err(MetadataError::parse(
                path,
                format!("unsupported pgen main-track compression type {other}"),
            ));
        }
    }
    if decoder_state.categories.len() != header.sample_ct {
        return Err(MetadataError::parse(
            path,
            "pgen decoded category count does not match sample count",
        ));
    }
    if !matches!(compression, 2 | 3) {
        decoder_state.previous_non_ld_categories.clear();
        decoder_state
            .previous_non_ld_categories
            .extend_from_slice(&decoder_state.categories);
        decoder_state.has_previous_non_ld = true;
    }
    Ok(())
}

fn decode_packed_categories(payload: &[u8], sample_ct: usize, categories: &mut Vec<u8>) {
    categories.clear();
    for sample_index in 0..sample_ct {
        let byte = payload[sample_index / 4];
        categories.push((byte >> ((sample_index % 4) * 2)) & 0b11);
    }
}

fn decode_one_bit_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    categories: &mut Vec<u8>,
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
    categories.clear();
    categories.resize(sample_ct, low_category);
    for (sample_index, category) in categories.iter_mut().enumerate() {
        if bit_is_set(bitarray, sample_index) {
            *category = high_category;
        }
    }
    let mut cursor = 1 + bitarray_len;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        categories[sample_index] = category;
    }
    Ok(())
}

fn decode_ld_compressed_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    previous_non_ld_categories: &[u8],
    inverted: bool,
    categories: &mut Vec<u8>,
) -> Result<()> {
    if previous_non_ld_categories.len() != sample_ct {
        return Err(MetadataError::parse(
            path,
            "pgen LD state length does not match sample count",
        ));
    }
    categories.clear();
    categories.extend_from_slice(previous_non_ld_categories);
    let mut cursor = 0;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        categories[sample_index] = category;
    }
    if inverted {
        invert_categories(categories);
    }
    Ok(())
}

fn decode_difflist_record(
    path: &Path,
    record: &[u8],
    sample_ct: usize,
    base_category: u8,
    categories: &mut Vec<u8>,
) -> Result<()> {
    categories.clear();
    categories.resize(sample_ct, base_category);
    let mut cursor = 0;
    for (sample_index, category) in decode_difflist(path, record, &mut cursor, sample_ct, true)? {
        categories[sample_index] = category;
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

fn invert_categories(categories: &mut [u8]) {
    for category in categories {
        *category = match *category {
            0 => 2,
            2 => 0,
            other => other,
        };
    }
}

fn decode_pgen_code(code: u8) -> (f32, bool) {
    match code {
        0b00 => (0.0, false),
        0b01 => (1.0, false),
        0b10 => (2.0, false),
        0b11 => (0.0, true),
        _ => unreachable!("two-bit PGEN code should be masked"),
    }
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
            (0..5).map(|sample_index| copy.get(sample_index)).collect::<Vec<_>>(),
            vec![2, 1, 0, 3, 3]
        );
        assert_eq!(
            (0..5).map(|sample_index| source.get(sample_index)).collect::<Vec<_>>(),
            vec![0, 1, 2, 3, 3]
        );
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
    let contents = fs::read_to_string(path).map_err(|source| MetadataError::Io {
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
) -> Result<Vec<(usize, VariantRecord)>> {
    let file = File::open(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let reader = BufReader::new(file);
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
        if source_index >= window_end {
            break;
        }
        if source_index >= window.start {
            let columns = columns
                .as_ref()
                .expect("pvar columns should be initialized before parsing body rows");
            records.push((
                source_index,
                parse_pvar_line(path, line_index + 1, columns, trimmed)?,
            ));
        }
        source_index += 1;
    }

    Ok(records)
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
