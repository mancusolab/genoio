// pattern: Imperative Shell

use std::fs;
use std::path::Path;

use genoio_core::{
    MetadataError, MetadataOutput, SampleRecord, SourceCapabilities, VariantRecord,
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
        MetadataError::parse(path, format!("bim line {line_number} has invalid position: {error}"))
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
    })
}

fn optional_plink_value(value: &str) -> Option<String> {
    if value == "0" {
        None
    } else {
        Some(value.to_string())
    }
}
