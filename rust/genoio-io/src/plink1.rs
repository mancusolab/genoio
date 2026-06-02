// pattern: Imperative Shell

use std::fs;
use std::path::Path;

use genoio_core::{
    attach_variant_stats, compute_variant_stats, select_samples_source_order,
    sparse_from_dense_minor_flipped, transpose_variant_major_to_sample_major, DenseGenotypeMatrix,
    MetadataError, MetadataOutput, SampleRecord, SourceCapabilities, SparseGenotypeMatrix,
    VariantFilter, VariantRecord,
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
    let bed_bytes = fs::read(bed).map_err(|source| MetadataError::Io {
        path: bed.to_path_buf(),
        source,
    })?;
    validate_bed_header(bed, &bed_bytes)?;

    let all_samples = parse_fam(fam)?;
    let source_variants = parse_bim(bim)?;
    let selection = select_samples_source_order(&all_samples, requested_samples, bed)?;
    let mut diagnostics = selection.diagnostics;
    let n_source_samples = all_samples.len();
    let n_source_variants = source_variants.len();
    let bytes_per_variant = n_source_samples.div_ceil(4);
    let expected_len = 3 + n_source_variants * bytes_per_variant;
    if bed_bytes.len() != expected_len {
        return Err(MetadataError::parse(
            bed,
            format!(
                "bed payload length {} does not match {n_source_samples} samples and {n_source_variants} variants",
                bed_bytes.len()
            ),
        ));
    }

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::with_capacity(selection.samples.len() * n_source_variants);
    let mut variant_major_missing = Vec::with_capacity(selection.samples.len() * n_source_variants);
    for (variant_index, mut variant) in source_variants.into_iter().enumerate() {
        diagnostics.candidate_variants += 1;
        if variant_filter.and_then(|filter| filter.metadata_decision(&variant)) == Some(false) {
            diagnostics.dropped_metadata_variants += 1;
            continue;
        }
        let variant_offset = 3 + variant_index * bytes_per_variant;
        let mut current_values = Vec::with_capacity(selection.source_indices.len());
        let mut current_missing = Vec::with_capacity(selection.source_indices.len());
        for source_index in &selection.source_indices {
            let byte = bed_bytes[variant_offset + source_index / 4];
            let code = (byte >> ((source_index % 4) * 2)) & 0b11;
            let (value, missing) = decode_plink1_code(code);
            current_values.push(value);
            current_missing.push(missing);
        }
        let stats = if variant_filter.is_some_and(VariantFilter::requires_genotype_stats) {
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
    let dense = read_plink1_dense(bed, bim, fam, requested_samples, variant_filter)?;
    Ok(sparse_from_dense_minor_flipped(dense)?)
}

fn validate_bed_header(path: &Path, bed_bytes: &[u8]) -> Result<()> {
    if bed_bytes.len() < 3 {
        return Err(MetadataError::parse(
            path,
            "bed file is shorter than the three magic bytes",
        ));
    }
    if bed_bytes[0] != 0x6c || bed_bytes[1] != 0x1b {
        return Err(MetadataError::parse(path, "invalid bed magic bytes"));
    }
    if bed_bytes[2] == 0x00 {
        return Err(MetadataError::parse(
            path,
            "sample-major bed mode is not supported",
        ));
    }
    if bed_bytes[2] != 0x01 {
        return Err(MetadataError::parse(path, "invalid bed mode byte"));
    }
    Ok(())
}

fn decode_plink1_code(code: u8) -> (f32, bool) {
    match code {
        0b00 => (2.0, false),
        0b01 => (1.0, false),
        0b10 => (0.0, true),
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
