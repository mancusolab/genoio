// pattern: Imperative Shell
//! PLINK2 PSAM and PVAR metadata parsing.
//!
//! The parser supports plain and zstd-compressed PVAR files, extracts the
//! metadata columns used by core records, and keeps source alleles alongside
//! normalized reference/alternate fields.

use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use genoio_core::{
    GenoioError, SampleRecord, VariantMetadataArrowBuffers, VariantRecord, VariantWindow,
};

use crate::error::Result;
use crate::plink::common::{optional_plink_value, PLINK2_MISSING_VALUES};

pub(super) fn parse_psam(path: &Path) -> Result<Vec<SampleRecord>> {
    let contents = fs::read_to_string(path).map_err(|source| GenoioError::Io {
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
            header = Some(parse_psam_header(path, trimmed)?);
            continue;
        }
        let columns = header
            .as_ref()
            .ok_or_else(|| GenoioError::invalid_source(path, "psam header line is required"))?;
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

fn parse_psam_header(path: &Path, line: &str) -> Result<PsamColumns> {
    let fields = line
        .trim_start_matches('#')
        .split_whitespace()
        .collect::<Vec<_>>();
    let find = |names: &[&str]| {
        fields
            .iter()
            .position(|field| names.iter().any(|name| field.eq_ignore_ascii_case(name)))
    };
    Ok(PsamColumns {
        fid: find(&["FID"]),
        iid: find(&["IID"])
            .ok_or_else(|| GenoioError::invalid_source(path, "psam header missing IID"))?,
        father: find(&["PAT", "FATHER"]),
        mother: find(&["MAT", "MOTHER"]),
        sex: find(&["SEX"]),
        phenotype: find(&["PHENO1", "PHENO", "PHENOTYPE"]),
    })
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
        return Err(GenoioError::invalid_source(
            path,
            format!("psam line {line_number} has too few fields"),
        ));
    }
    Ok(SampleRecord {
        fid: columns
            .fid
            .and_then(|index| optional_plink_value(fields[index], PLINK2_MISSING_VALUES)),
        iid: fields[columns.iid].to_string(),
        father: columns
            .father
            .and_then(|index| optional_plink_value(fields[index], PLINK2_MISSING_VALUES)),
        mother: columns
            .mother
            .and_then(|index| optional_plink_value(fields[index], PLINK2_MISSING_VALUES)),
        sex: columns.sex.map(|index| fields[index].to_string()),
        phenotype: columns.phenotype.map(|index| fields[index].to_string()),
        source_sample_index: None,
        haplotype_index: None,
    })
}

pub(super) fn parse_pvar_arrow(path: &Path) -> Result<VariantMetadataArrowBuffers> {
    let mut reader = open_pvar_reader(path)?;
    let mut contents = String::new();
    reader
        .read_to_string(&mut contents)
        .map_err(|source| GenoioError::Io {
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
    let mut variants =
        VariantMetadataArrowBuffers::with_capacity(data_lines.len().saturating_sub(body_start));
    for (index, line) in data_lines.into_iter().skip(body_start) {
        append_pvar_arrow_line(path, index + 1, &columns, line, &mut variants)?;
    }
    Ok(variants)
}

pub(super) fn parse_pvar_source_window(
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
        let line = line_result.map_err(|source| GenoioError::Io {
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

        let Some(columns) = columns.as_ref() else {
            return Err(GenoioError::invalid_source(
                path,
                "pvar columns were not initialized before parsing body rows",
            ));
        };
        let variant = parse_pvar_line(path, line_index + 1, columns, trimmed)?;
        if source_index >= window.start && source_index < window_end {
            records.push((source_index, variant));
        }
        source_index += 1;
        if source_index >= window_end {
            break;
        }
    }

    let required_rows = window_end.min(expected_variant_ct);
    if source_index < required_rows {
        return Err(GenoioError::invalid_source(
            path,
            format!("pvar variant count {source_index} is shorter than requested source window"),
        ));
    }

    Ok(records)
}

fn append_pvar_arrow_line(
    path: &Path,
    line_number: usize,
    columns: &PvarColumns,
    line: &str,
    variants: &mut VariantMetadataArrowBuffers,
) -> Result<()> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    let required = columns
        .chrom
        .max(columns.pos)
        .max(columns.id)
        .max(columns.ref_allele)
        .max(columns.alt_allele)
        .max(columns.qual.unwrap_or(0));
    if fields.len() <= required {
        return Err(GenoioError::invalid_source(
            path,
            format!("pvar line {line_number} has too few fields"),
        ));
    }
    let pos = fields[columns.pos].parse::<i64>().map_err(|error| {
        GenoioError::invalid_source(
            path,
            format!("pvar line {line_number} has invalid position: {error}"),
        )
    })?;
    let ref_allele = fields[columns.ref_allele];
    let alt_allele = fields[columns.alt_allele];
    let first_alt = alt_allele.split(',').next().unwrap_or("");
    if first_alt.is_empty() {
        return Err(GenoioError::invalid_source(
            path,
            format!("pvar line {line_number} has empty ALT allele"),
        ));
    }
    let qual = columns
        .qual
        .map(|index| parse_optional_qual(path, line_number, fields[index]))
        .transpose()?
        .flatten();

    variants.chroms.append_value(fields[columns.chrom])?;
    variants.positions.push(pos);
    variants.ids.append_value(fields[columns.id])?;
    variants.a0s.append_value(ref_allele)?;
    variants.a1s.append_value(first_alt)?;
    variants.ref_alleles.push(Some(ref_allele.to_string()));
    variants.alt_alleles.push(Some(alt_allele.to_string()));
    variants.source_a0s.append_value(ref_allele)?;
    variants.source_a1s.append_value(first_alt)?;
    variants.flipped.push(false);
    variants.quals.push(qual);
    variants.afs.push(None);
    variants.mafs.push(None);
    variants.macs.push(None);
    variants.missing_rates.push(None);
    variants.n_called.push(None);
    Ok(())
}

pub(super) struct PvarRecordReader {
    path: PathBuf,
    lines: std::iter::Enumerate<std::io::Lines<Box<dyn BufRead>>>,
    columns: Option<PvarColumns>,
    body_started: bool,
    source_index: usize,
}

impl PvarRecordReader {
    pub(super) fn new(path: &Path) -> Result<Self> {
        Ok(Self {
            path: path.to_path_buf(),
            lines: open_pvar_reader(path)?.lines().enumerate(),
            columns: None,
            body_started: false,
            source_index: 0,
        })
    }

    pub(super) fn next_record(&mut self) -> Result<Option<(usize, VariantRecord)>> {
        for (line_index, line_result) in self.lines.by_ref() {
            let line = line_result.map_err(|source| GenoioError::Io {
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

            let Some(columns) = self.columns.as_ref() else {
                return Err(GenoioError::invalid_source(
                    &self.path,
                    "pvar columns were not initialized before parsing body rows",
                ));
            };
            let variant = parse_pvar_line(&self.path, line_index + 1, columns, trimmed)?;
            let source_index = self.source_index;
            self.source_index += 1;
            return Ok(Some((source_index, variant)));
        }
        Ok(None)
    }

    pub(super) fn validate_count(&self, expected_variant_ct: usize) -> Result<()> {
        if self.source_index != expected_variant_ct {
            return Err(GenoioError::invalid_source(
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
    let file = File::open(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    if path
        .file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.ends_with(".pvar.zst"))
    {
        let decoder = zstd::stream::read::Decoder::new(file).map_err(|source| GenoioError::Io {
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
            .ok_or_else(|| GenoioError::invalid_source("<pvar>", "pvar header missing #CHROM"))?,
        pos: find("POS")
            .ok_or_else(|| GenoioError::invalid_source("<pvar>", "pvar header missing POS"))?,
        id: find("ID")
            .ok_or_else(|| GenoioError::invalid_source("<pvar>", "pvar header missing ID"))?,
        ref_allele: find("REF")
            .ok_or_else(|| GenoioError::invalid_source("<pvar>", "pvar header missing REF"))?,
        alt_allele: find("ALT")
            .ok_or_else(|| GenoioError::invalid_source("<pvar>", "pvar header missing ALT"))?,
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
        _ => Err(GenoioError::invalid_source(
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
        return Err(GenoioError::invalid_source(
            path,
            format!("pvar line {line_number} has too few fields"),
        ));
    }
    let pos = fields[columns.pos].parse::<u32>().map_err(|error| {
        GenoioError::invalid_source(
            path,
            format!("pvar line {line_number} has invalid position: {error}"),
        )
    })?;
    let ref_allele = fields[columns.ref_allele].to_string();
    let alt_allele = fields[columns.alt_allele].to_string();
    let first_alt = alt_allele.split(',').next().unwrap_or("").to_string();
    if first_alt.is_empty() {
        return Err(GenoioError::invalid_source(
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
        GenoioError::invalid_source(
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
