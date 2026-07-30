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
#[cfg(test)]
use super::session::Plink1WorkProbe;

pub(super) fn parse_fam(path: &Path) -> Result<Vec<SampleRecord>> {
    #[cfg(test)]
    {
        parse_fam_inner(path, None)
    }
    #[cfg(not(test))]
    {
        parse_fam_inner(path)
    }
}

#[cfg(test)]
pub(super) fn parse_fam_with_probe(
    path: &Path,
    probe: &Plink1WorkProbe,
) -> Result<Vec<SampleRecord>> {
    parse_fam_inner(path, Some(probe))
}

fn parse_fam_inner(
    path: &Path,
    #[cfg(test)] probe: Option<&Plink1WorkProbe>,
) -> Result<Vec<SampleRecord>> {
    let mut records = Vec::new();
    #[cfg(test)]
    let lines = open_text_lines_inner(path, probe.map(|probe| (probe, TextFileKind::Fam)))?;
    #[cfg(not(test))]
    let lines = open_text_lines(path)?;
    for (index, line) in lines.enumerate() {
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

/// Streams non-empty BIM rows as source-indexed variant records.
pub(super) struct BimRecordReader {
    path: std::path::PathBuf,
    lines: std::iter::Enumerate<Lines<BufReader<File>>>,
    source_index: usize,
}

impl BimRecordReader {
    pub(super) fn new(path: &Path) -> Result<Self> {
        #[cfg(test)]
        {
            Self::new_inner(path, None)
        }
        #[cfg(not(test))]
        {
            Self::new_inner(path)
        }
    }

    #[cfg(test)]
    pub(super) fn new_with_probe(path: &Path, probe: &Plink1WorkProbe) -> Result<Self> {
        Self::new_inner(path, Some(probe))
    }

    fn new_inner(path: &Path, #[cfg(test)] probe: Option<&Plink1WorkProbe>) -> Result<Self> {
        #[cfg(test)]
        let lines = open_text_lines_inner(path, probe.map(|probe| (probe, TextFileKind::Bim)))?;
        #[cfg(not(test))]
        let lines = open_text_lines(path)?;
        Ok(Self {
            path: path.to_path_buf(),
            lines: lines.enumerate(),
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

#[cfg(not(test))]
fn open_text_lines(path: &Path) -> Result<Lines<BufReader<File>>> {
    open_text_lines_inner(path)
}

#[cfg(test)]
#[derive(Clone, Copy)]
enum TextFileKind {
    Bim,
    Fam,
}

fn open_text_lines_inner(
    path: &Path,
    #[cfg(test)] probe: Option<(&Plink1WorkProbe, TextFileKind)>,
) -> Result<Lines<BufReader<File>>> {
    let file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    #[cfg(test)]
    if let Some((probe, kind)) = probe {
        match kind {
            TextFileKind::Bim => probe.record_bim_open(),
            TextFileKind::Fam => probe.record_fam_open(),
        }
    }
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
