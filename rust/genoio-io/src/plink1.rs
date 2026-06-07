// pattern: Imperative Shell

use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, flip_values_to_minor_allele,
    reject_sparse_missing_values, select_samples_source_order, DenseGenotypeMatrix, GenoioError,
    MetadataOutput, PartialFilterDecision, SampleRecord, SourceCapabilities, SparseGenotypeMatrix,
    VariantFilter, VariantRecord, VariantWindow,
};

use crate::error::Result;
use crate::hardcall::{HardcallBatch, PackedHardcalls};
use crate::matrix::{
    empty_dense_matrix, empty_sparse_matrix, finish_variant_major_dense_matrix,
    VariantMajorDenseParts,
};
use crate::plink_common::{optional_plink_value, PLINK1_MISSING_VALUES};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

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

/// Read all retained PLINK1 genotypes as a dense sample-by-variant matrix.
pub fn read_plink1_dense(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_plink1_dense_windowed(bed, bim, fam, requested_samples, variant_filter, None)
}

/// Read retained PLINK1 genotypes as a dense matrix over an optional block window.
pub fn read_plink1_dense_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<DenseGenotypeMatrix> {
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        return empty_plink1_dense(bed, bim, fam, requested_samples);
    }

    if let (None, Some(window)) = (variant_filter, variant_window) {
        return read_plink1_dense_source_window(bed, bim, fam, requested_samples, window);
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

    let output_variant_capacity = variant_window.map_or(n_source_variants, |window| {
        window.len.min(n_source_variants)
    });
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut variant_major_values =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut variant_major_missing =
        Vec::with_capacity(selection.samples.len() * output_variant_capacity);
    let mut decoder_state =
        Plink1DecoderState::new(n_source_samples, bytes_per_variant, selection.samples.len());
    let mut retention = RetainedVariantState::new(variant_window);
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

        let stats = if needs_genotype_decision {
            Some(
                decoder_state
                    .packed
                    .stats_for_selection(&selection.source_indices, all_samples_selected)?,
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
                RetentionAction::Stop => break,
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
        variants.push(variant);
        variant_major_values.extend_from_slice(&decoder_state.values);
        variant_major_missing.extend_from_slice(&decoder_state.missing);
        if retention.window_is_satisfied() {
            break;
        }
    }

    let n_samples = selection.samples.len();
    let n_variants = variants.len();
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
        false,
    )
}

fn read_plink1_dense_source_window(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    window: VariantWindow,
) -> Result<DenseGenotypeMatrix> {
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
    let variants = parse_bim_source_window(bim, window, n_variants)?;

    let n_samples = selection.samples.len();
    let mut values = vec![0.0; n_samples * n_variants];
    let mut missing_mask = vec![false; n_samples * n_variants];
    let mut decoder_state =
        Plink1DecoderState::new(n_source_samples, bytes_per_variant, selection.samples.len());
    let mut batch = HardcallBatch::new(n_source_samples);
    let mut batch_start = 0_usize;
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
            batch.expand_into_sample_major(
                &selection.source_indices,
                batch_start,
                n_variants,
                &mut values,
                &mut missing_mask,
            );
            batch_start += batch.len();
            batch.clear();
        }
        diagnostics.candidate_variants = variant_offset + 1;
    }
    if !batch.is_empty() {
        batch.expand_into_sample_major(
            &selection.source_indices,
            batch_start,
            n_variants,
            &mut values,
            &mut missing_mask,
        );
    }

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
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        return empty_plink1_sparse(bed, bim, fam, requested_samples);
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
    let mut variants = Vec::with_capacity(output_variant_capacity);
    let mut decoder_state =
        Plink1DecoderState::new(n_source_samples, bytes_per_variant, selection.samples.len());
    let mut retention = RetainedVariantState::new(variant_window);
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
        let stats = if needs_genotype_decision {
            Some(
                decoder_state
                    .packed
                    .stats_for_selection(&selection.source_indices, all_samples_selected)?,
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
                RetentionAction::Stop => break,
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
            break;
        }
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

fn empty_plink1_dense(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
) -> Result<DenseGenotypeMatrix> {
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
    empty_dense_matrix(selection.samples, selection.diagnostics)
}

fn empty_plink1_sparse(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
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
    empty_sparse_matrix(selection.samples, selection.diagnostics)
}

fn open_bed_file(path: &Path) -> Result<File> {
    let mut file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = [0_u8; 3];
    file.read_exact(&mut header)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    validate_bed_header(path, &header)?;
    Ok(file)
}

fn validate_bed_header(path: &Path, header: &[u8; 3]) -> Result<()> {
    if header[0] != 0x6c || header[1] != 0x1b {
        return Err(GenoioError::invalid_source(path, "invalid bed magic bytes"));
    }
    if header[2] == 0x00 {
        return Err(GenoioError::invalid_source(
            path,
            "sample-major bed mode is not supported",
        ));
    }
    if header[2] != 0x01 {
        return Err(GenoioError::invalid_source(path, "invalid bed mode byte"));
    }
    Ok(())
}

fn validate_bed_payload_len(
    path: &Path,
    file: &File,
    n_source_samples: usize,
    n_source_variants: usize,
    bytes_per_variant: usize,
) -> Result<()> {
    let expected_len = 3 + n_source_variants * bytes_per_variant;
    let actual_len = file
        .metadata()
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    let expected_len_u64 = u64::try_from(expected_len)
        .map_err(|_| GenoioError::invalid_source(path, "bed payload length is out of range"))?;
    if actual_len != expected_len_u64 {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "bed payload length {actual_len} does not match {n_source_samples} samples and {n_source_variants} variants"
            ),
        ));
    }
    Ok(())
}

fn infer_bed_variant_count(
    path: &Path,
    file: &File,
    n_source_samples: usize,
    bytes_per_variant: usize,
) -> Result<usize> {
    let actual_len = file
        .metadata()
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if actual_len < 3 {
        return Err(GenoioError::invalid_source(
            path,
            "bed payload length is out of range",
        ));
    }
    let payload_len = actual_len - 3;
    let bytes_per_variant_u64 = u64::try_from(bytes_per_variant)
        .map_err(|_| GenoioError::invalid_source(path, "bed payload length is out of range"))?;
    if bytes_per_variant_u64 == 0 || payload_len % bytes_per_variant_u64 != 0 {
        return Err(GenoioError::invalid_source(
            path,
            format!("bed payload length {actual_len} does not match {n_source_samples} samples"),
        ));
    }
    usize::try_from(payload_len / bytes_per_variant_u64)
        .map_err(|_| GenoioError::invalid_source(path, "bed payload length is out of range"))
}

#[derive(Debug, Clone)]
struct Plink1DecoderState {
    payload: Vec<u8>,
    packed: PackedHardcalls,
    values: Vec<f32>,
    missing: Vec<bool>,
}

impl Plink1DecoderState {
    fn new(sample_ct: usize, bytes_per_variant: usize, selected_sample_ct: usize) -> Self {
        Self {
            payload: Vec::with_capacity(bytes_per_variant),
            packed: {
                let mut packed = PackedHardcalls::default();
                packed.resize(sample_ct);
                packed
            },
            values: Vec::with_capacity(selected_sample_ct),
            missing: Vec::with_capacity(selected_sample_ct),
        }
    }
}

fn read_plink1_variant_packed(
    path: &Path,
    file: &mut File,
    variant_index: usize,
    bytes_per_variant: usize,
    sample_ct: usize,
    decoder_state: &mut Plink1DecoderState,
) -> Result<()> {
    seek_plink1_variant(path, file, variant_index, bytes_per_variant)?;
    read_plink1_variant_packed_sequential(path, file, bytes_per_variant, sample_ct, decoder_state)
}

fn seek_plink1_variant(
    path: &Path,
    file: &mut File,
    variant_index: usize,
    bytes_per_variant: usize,
) -> Result<()> {
    let offset = 3 + variant_index * bytes_per_variant;
    file.seek(SeekFrom::Start(u64::try_from(offset).map_err(|_| {
        GenoioError::invalid_source(path, "bed variant offset is out of range")
    })?))
    .map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

fn read_plink1_variant_packed_sequential(
    path: &Path,
    file: &mut File,
    bytes_per_variant: usize,
    sample_ct: usize,
    decoder_state: &mut Plink1DecoderState,
) -> Result<()> {
    decoder_state.payload.resize(bytes_per_variant, 0);
    file.read_exact(&mut decoder_state.payload)
        .map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    decoder_state
        .packed
        .load_plink1_bed_payload(&decoder_state.payload, sample_ct);
    Ok(())
}

fn parse_fam(path: &Path) -> Result<Vec<SampleRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    contents
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty())
        .map(|(index, line)| parse_fam_line(path, index + 1, line))
        .collect()
}

fn parse_fam_line(path: &Path, line_number: usize, line: &str) -> Result<SampleRecord> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if fields.len() < 6 {
        return Err(GenoioError::invalid_source(
            path,
            format!("fam line {line_number} has fewer than six fields"),
        ));
    }

    Ok(SampleRecord {
        fid: Some(fields[0].to_string()),
        iid: fields[1].to_string(),
        father: optional_plink_value(fields[2], PLINK1_MISSING_VALUES),
        mother: optional_plink_value(fields[3], PLINK1_MISSING_VALUES),
        sex: Some(fields[4].to_string()),
        phenotype: Some(fields[5].to_string()),
        source_sample_index: None,
        haplotype_index: None,
    })
}

fn parse_bim(path: &Path) -> Result<Vec<VariantRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    contents
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty())
        .map(|(index, line)| parse_bim_line(path, index + 1, line))
        .collect()
}

fn parse_bim_source_window(
    path: &Path,
    window: VariantWindow,
    expected_records: usize,
) -> Result<Vec<VariantRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut records = Vec::with_capacity(expected_records);
    let end = window.start.saturating_add(expected_records);
    let mut source_index = 0_usize;
    for (line_index, line) in contents.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if source_index >= end {
            break;
        }
        if source_index >= window.start {
            records.push(parse_bim_line(path, line_index + 1, line)?);
        }
        source_index += 1;
    }
    if records.len() != expected_records {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "bim source window contains {} variants but expected {expected_records}",
                records.len()
            ),
        ));
    }
    Ok(records)
}

fn parse_bim_line(path: &Path, line_number: usize, line: &str) -> Result<VariantRecord> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if fields.len() < 6 {
        return Err(GenoioError::invalid_source(
            path,
            format!("bim line {line_number} has fewer than six fields"),
        ));
    }
    let pos = fields[3].parse::<u32>().map_err(|error| {
        GenoioError::invalid_source(
            path,
            format!("bim line {line_number} has invalid position: {error}"),
        )
    })?;
    let a1 = fields[4].to_string();
    let a0 = fields[5].to_string();

    Ok(VariantRecord {
        chrom: fields[0].to_string(),
        pos,
        id: fields[1].to_string(),
        a0: a0.clone(),
        a1: a1.clone(),
        ref_allele: None,
        alt_allele: None,
        source_a0: a0,
        source_a1: a1,
        flipped: false,
        qual: None,
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
}
