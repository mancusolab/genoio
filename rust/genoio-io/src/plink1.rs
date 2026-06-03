// pattern: Imperative Shell

use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use genoio_core::{
    append_sparse_column, attach_variant_stats, compute_variant_stats, flip_values_to_minor_allele,
    reject_sparse_missing_values, select_samples_source_order,
    transpose_variant_major_to_sample_major, DenseGenotypeMatrix, MetadataError, MetadataOutput,
    SampleRecord, SourceCapabilities, SparseGenotypeMatrix, VariantFilter, VariantRecord,
    VariantWindow,
};

use crate::error::Result;

pub fn read_plink1_metadata(bed: &Path, bim: &Path, fam: &Path) -> Result<MetadataOutput> {
    fs::metadata(bed).map_err(|source| MetadataError::Io {
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

pub fn read_plink1_dense(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_plink1_dense_windowed(bed, bim, fam, requested_samples, variant_filter, None)
}

pub fn read_plink1_dense_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<DenseGenotypeMatrix> {
    let mut bed_file = open_bed_file(bed)?;

    let all_samples = parse_fam(fam)?;
    let source_variants = parse_bim(bim)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
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

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::with_capacity(selection.samples.len() * n_source_variants);
    let mut variant_major_missing = Vec::with_capacity(selection.samples.len() * n_source_variants);
    let mut retained_index = 0_usize;
    for (variant_index, mut variant) in source_variants.into_iter().enumerate() {
        diagnostics.candidate_variants += 1;
        if variant_filter.and_then(|filter| filter.metadata_decision(&variant)) == Some(false) {
            diagnostics.dropped_metadata_variants += 1;
            continue;
        }

        let requires_stats = variant_filter.is_some_and(VariantFilter::requires_genotype_stats);
        if !requires_stats {
            if variant_filter.is_some_and(|filter| !filter.evaluate(&variant, None)) {
                diagnostics.dropped_genotype_variants += 1;
                continue;
            }
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                continue;
            }
        }

        let (current_values, current_missing) = read_plink1_variant_values(
            bed,
            &mut bed_file,
            variant_index,
            bytes_per_variant,
            &selection.source_indices,
        )?;

        let stats = if requires_stats {
            Some(compute_variant_stats(&current_values, &current_missing)?)
        } else {
            None
        };
        if variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref())) {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if requires_stats {
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                continue;
            }
        }
        variants.push(variant);
        variant_major_values.extend(current_values);
        variant_major_missing.extend(current_missing);
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

pub fn read_plink1_sparse(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_plink1_sparse_windowed(bed, bim, fam, requested_samples, variant_filter, None)
}

pub fn read_plink1_sparse_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    let mut bed_file = open_bed_file(bed)?;

    let all_samples = parse_fam(fam)?;
    let source_variants = parse_bim(bim)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
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
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();
    let mut retained_index = 0_usize;
    for (variant_index, mut variant) in source_variants.into_iter().enumerate() {
        diagnostics.candidate_variants += 1;
        if variant_filter.and_then(|filter| filter.metadata_decision(&variant)) == Some(false) {
            diagnostics.dropped_metadata_variants += 1;
            continue;
        }

        let requires_stats = variant_filter.is_some_and(VariantFilter::requires_genotype_stats);
        if !requires_stats {
            if variant_filter.is_some_and(|filter| !filter.evaluate(&variant, None)) {
                diagnostics.dropped_genotype_variants += 1;
                continue;
            }
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                continue;
            }
        }

        let (mut current_values, current_missing) = read_plink1_variant_values(
            bed,
            &mut bed_file,
            variant_index,
            bytes_per_variant,
            &selection.source_indices,
        )?;
        let stats = if requires_stats {
            Some(compute_variant_stats(&current_values, &current_missing)?)
        } else {
            None
        };
        if variant_filter.is_some_and(|filter| !filter.evaluate(&variant, stats.as_ref())) {
            diagnostics.dropped_genotype_variants += 1;
            continue;
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        if requires_stats {
            let include_in_window =
                variant_window.is_none_or(|window| window.contains(retained_index));
            retained_index += 1;
            if !include_in_window {
                continue;
            }
        }
        reject_sparse_missing_values(&current_missing)?;
        flip_values_to_minor_allele(&mut current_values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &current_values);
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

fn open_bed_file(path: &Path) -> Result<File> {
    let mut file = File::open(path).map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut header = [0_u8; 3];
    file.read_exact(&mut header)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    validate_bed_header(path, &header)?;
    Ok(file)
}

fn validate_bed_header(path: &Path, header: &[u8; 3]) -> Result<()> {
    if header[0] != 0x6c || header[1] != 0x1b {
        return Err(MetadataError::parse(path, "invalid bed magic bytes"));
    }
    if header[2] == 0x00 {
        return Err(MetadataError::parse(
            path,
            "sample-major bed mode is not supported",
        ));
    }
    if header[2] != 0x01 {
        return Err(MetadataError::parse(path, "invalid bed mode byte"));
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
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?
        .len();
    let expected_len_u64 = u64::try_from(expected_len)
        .map_err(|_| MetadataError::parse(path, "bed payload length is out of range"))?;
    if actual_len != expected_len_u64 {
        return Err(MetadataError::parse(
            path,
            format!(
                "bed payload length {actual_len} does not match {n_source_samples} samples and {n_source_variants} variants"
            ),
        ));
    }
    Ok(())
}

fn read_plink1_variant_values(
    path: &Path,
    file: &mut File,
    variant_index: usize,
    bytes_per_variant: usize,
    source_indices: &[usize],
) -> Result<(Vec<f32>, Vec<bool>)> {
    let offset = 3 + variant_index * bytes_per_variant;
    file.seek(SeekFrom::Start(u64::try_from(offset).map_err(|_| {
        MetadataError::parse(path, "bed variant offset is out of range")
    })?))
    .map_err(|source| MetadataError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let mut payload = vec![0_u8; bytes_per_variant];
    file.read_exact(&mut payload)
        .map_err(|source| MetadataError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    let mut current_values = Vec::with_capacity(source_indices.len());
    let mut current_missing = Vec::with_capacity(source_indices.len());
    for source_index in source_indices {
        let byte = payload[source_index / 4];
        let code = (byte >> ((source_index % 4) * 2)) & 0b11;
        let (value, missing) = decode_plink1_code(code);
        current_values.push(value);
        current_missing.push(missing);
    }
    Ok((current_values, current_missing))
}

fn decode_plink1_code(code: u8) -> (f32, bool) {
    match code {
        0b00 => (2.0, false),
        0b01 => (0.0, true),
        0b10 => (1.0, false),
        0b11 => (0.0, false),
        _ => unreachable!("two-bit PLINK1 code should be masked"),
    }
}

fn parse_fam(path: &Path) -> Result<Vec<SampleRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| MetadataError::Io {
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
        return Err(MetadataError::parse(
            path,
            format!("fam line {line_number} has fewer than six fields"),
        ));
    }

    Ok(SampleRecord {
        fid: Some(fields[0].to_string()),
        iid: fields[1].to_string(),
        father: optional_plink_value(fields[2]),
        mother: optional_plink_value(fields[3]),
        sex: Some(fields[4].to_string()),
        phenotype: Some(fields[5].to_string()),
        source_sample_index: None,
        haplotype_index: None,
    })
}

fn parse_bim(path: &Path) -> Result<Vec<VariantRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| MetadataError::Io {
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

fn parse_bim_line(path: &Path, line_number: usize, line: &str) -> Result<VariantRecord> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if fields.len() < 6 {
        return Err(MetadataError::parse(
            path,
            format!("bim line {line_number} has fewer than six fields"),
        ));
    }
    let pos = fields[3].parse::<u32>().map_err(|error| {
        MetadataError::parse(
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
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
}

fn optional_plink_value(value: &str) -> Option<String> {
    if value == "0" {
        None
    } else {
        Some(value.to_string())
    }
}
