// pattern: Imperative Shell

use std::fs;
use std::path::Path;

use genoio_core::{GenoioError, SampleRecord, VariantRecord, VariantWindow};

use crate::error::Result;

use super::super::common::{optional_plink_value, PLINK1_MISSING_VALUES};

pub(super) fn parse_fam(path: &Path) -> Result<Vec<SampleRecord>> {
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

pub(super) fn parse_bim(path: &Path) -> Result<Vec<VariantRecord>> {
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

pub(super) fn parse_bim_source_window(
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
