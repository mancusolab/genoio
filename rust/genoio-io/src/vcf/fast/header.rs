//! Minimal VCF header handling for the lazy fast path.
//!
//! The fast path needs sample IDs but not structured declarations. Parsing only
//! the final #CHROM line keeps metadata and GT reads permissive for real-world
//! headers that strict parsers can reject.

use std::io::BufRead;
use std::path::Path;

use genoio_core::{GenoioError, SampleRecord};

use crate::error::Result;

pub(super) fn read_sample_records_from_header(
    path: &Path,
    reader: &mut dyn BufRead,
) -> Result<Vec<SampleRecord>> {
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line).map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf fast header read error: {error}"))
        })?;
        if n == 0 {
            return Err(GenoioError::invalid_source(
                path,
                "vcf header is missing #CHROM line",
            ));
        }
        // noodles' full header parser is intentionally strict. Parse the final
        // header line directly because metadata and GT reads only need source
        // sample order here.
        if line.starts_with("#CHROM\t") {
            return sample_records_from_chrom_header(line.trim_end_matches(['\r', '\n']));
        }
        if !line.starts_with("##") {
            return Err(GenoioError::invalid_source(
                path,
                "vcf header is missing #CHROM line before records",
            ));
        }
    }
}

fn sample_records_from_chrom_header(line: &str) -> Result<Vec<SampleRecord>> {
    let mut columns = line.split('\t');
    let expected = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    ];

    for expected_column in expected {
        if columns.next() != Some(expected_column) {
            return Err(GenoioError::invalid_source(
                "<vcf header>",
                "vcf #CHROM header does not have the required fixed columns",
            ));
        }
    }

    let Some(format_column) = columns.next() else {
        return Ok(Vec::new());
    };
    if format_column != "FORMAT" {
        return Err(GenoioError::invalid_source(
            "<vcf header>",
            "vcf #CHROM header has sample columns without FORMAT",
        ));
    }

    Ok(columns
        .map(|sample| SampleRecord {
            fid: None,
            iid: sample.to_string(),
            father: None,
            mother: None,
            sex: None,
            phenotype: None,
            source_sample_index: None,
            haplotype_index: None,
        })
        .collect())
}
