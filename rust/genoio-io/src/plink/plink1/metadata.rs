// pattern: Imperative Shell
//! PLINK1 FAM and BIM metadata parsing.
//!
//! The parser is intentionally line-oriented because companion files are small
//! text tables. It normalizes PLINK1 missing tokens while preserving source
//! allele orientation for downstream metadata.

use std::fs::File;
use std::io::{BufRead, BufReader, Lines};
use std::path::Path;

use genoio_core::{GenoioError, SampleRecord, VariantMetadataBuffers, VariantRecord};

use crate::error::Result;

use super::super::common::{optional_plink_value, PLINK1_MISSING_VALUES};

pub(super) fn parse_fam(path: &Path) -> Result<Vec<SampleRecord>> {
    let mut records = Vec::new();
    for (index, line) in open_text_lines(path)?.enumerate() {
        let line = line.map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if line.trim().is_empty() {
            continue;
        }
        records.push(parse_fam_line(path, index + 1, &line)?);
    }
    Ok(records)
}

fn parse_fam_line(path: &Path, line_number: usize, line: &str) -> Result<SampleRecord> {
    let mut fields = line.split_whitespace();
    let fid = required_field(path, "fam", line_number, 6, &mut fields)?;
    let iid = required_field(path, "fam", line_number, 6, &mut fields)?;
    let father = required_field(path, "fam", line_number, 6, &mut fields)?;
    let mother = required_field(path, "fam", line_number, 6, &mut fields)?;
    let sex = required_field(path, "fam", line_number, 6, &mut fields)?;
    let phenotype = required_field(path, "fam", line_number, 6, &mut fields)?;

    Ok(SampleRecord {
        fid: Some(fid.to_string()),
        iid: iid.to_string(),
        father: optional_plink_value(father, PLINK1_MISSING_VALUES),
        mother: optional_plink_value(mother, PLINK1_MISSING_VALUES),
        sex: Some(sex.to_string()),
        phenotype: Some(phenotype.to_string()),
        source_sample_index: None,
        haplotype_index: None,
    })
}

pub(super) fn parse_bim_metadata(path: &Path) -> Result<VariantMetadataBuffers> {
    let mut variants = VariantMetadataBuffers::with_capacity(0);
    let mut reader = BimRecordReader::new(path)?;
    while let Some((_, variant)) = reader.next_record()? {
        variants.push_record(&variant)?;
    }
    Ok(variants)
}

pub(super) fn parse_bim_source_window(
    path: &Path,
    start: usize,
    expected_records: usize,
) -> Result<VariantMetadataBuffers> {
    let mut variants = VariantMetadataBuffers::with_capacity(expected_records);
    let end = start.saturating_add(expected_records);
    let mut reader = BimRecordReader::new(path)?;
    while let Some((source_index, variant)) = reader.next_record()? {
        if source_index >= end {
            break;
        }
        if source_index >= start {
            variants.push_record(&variant)?;
        }
    }
    if variants.len() != expected_records {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "bim source window contains {} variants but expected {expected_records}",
                variants.len()
            ),
        ));
    }
    Ok(variants)
}

pub(super) fn count_bim_records(path: &Path) -> Result<usize> {
    let mut count = 0_usize;
    for line in open_text_lines(path)? {
        let line = line.map_err(|source| GenoioError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if !line.trim().is_empty() {
            count += 1;
        }
    }
    Ok(count)
}

/// Streams non-empty BIM rows as source-indexed variant records.
pub(super) struct BimRecordReader {
    path: std::path::PathBuf,
    lines: std::iter::Enumerate<Lines<BufReader<File>>>,
    source_index: usize,
}

impl BimRecordReader {
    pub(super) fn new(path: &Path) -> Result<Self> {
        Ok(Self {
            path: path.to_path_buf(),
            lines: open_text_lines(path)?.enumerate(),
            source_index: 0,
        })
    }

    pub(super) fn next_record(&mut self) -> Result<Option<(usize, VariantRecord)>> {
        for (line_index, line) in self.lines.by_ref() {
            let line = line.map_err(|source| GenoioError::Io {
                path: self.path.clone(),
                source,
            })?;
            if line.trim().is_empty() {
                continue;
            }
            let variant = parse_bim_line(&self.path, line_index + 1, &line)?;
            let source_index = self.source_index;
            self.source_index += 1;
            return Ok(Some((source_index, variant)));
        }
        Ok(None)
    }
}

fn parse_bim_line(path: &Path, line_number: usize, line: &str) -> Result<VariantRecord> {
    let mut fields = line.split_whitespace();
    let chrom = required_field(path, "bim", line_number, 6, &mut fields)?;
    let id = required_field(path, "bim", line_number, 6, &mut fields)?;
    let _cm = required_field(path, "bim", line_number, 6, &mut fields)?;
    let pos = required_field(path, "bim", line_number, 6, &mut fields)?
        .parse::<u32>()
        .map_err(|error| {
            GenoioError::invalid_source(
                path,
                format!("bim line {line_number} has invalid position: {error}"),
            )
        })?;
    let a1 = required_field(path, "bim", line_number, 6, &mut fields)?.to_string();
    let a0 = required_field(path, "bim", line_number, 6, &mut fields)?.to_string();

    Ok(VariantRecord {
        chrom: chrom.to_string(),
        pos,
        id: id.to_string(),
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

fn open_text_lines(path: &Path) -> Result<Lines<BufReader<File>>> {
    let file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(BufReader::new(file).lines())
}

fn required_field<'a>(
    path: &Path,
    file_kind: &str,
    line_number: usize,
    expected_fields: usize,
    fields: &mut impl Iterator<Item = &'a str>,
) -> Result<&'a str> {
    fields.next().ok_or_else(|| {
        GenoioError::invalid_source(
            path,
            format!("{file_kind} line {line_number} has fewer than {expected_fields} fields"),
        )
    })
}
