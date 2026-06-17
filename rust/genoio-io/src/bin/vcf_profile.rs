use std::env;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use genoio_core::{VariantFilter, VariantWindow};
use serde_json::json;

#[derive(Debug)]
struct Args {
    path: PathBuf,
    max_variants: usize,
    repeats: usize,
    threads: Option<usize>,
    region: Option<String>,
    samples: Option<Vec<String>>,
    scenario: Scenario,
    filter: FilterMode,
}

#[derive(Debug, Clone, Copy)]
enum Scenario {
    Dense,
    Dosage,
    Sparse,
    HaploDense,
    HaploSparse,
}

#[derive(Debug, Clone, Copy)]
enum FilterMode {
    Biallelic,
    Snp,
    QualMin(f64),
    MafMax(f64),
    None,
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), String> {
    let args = parse_args(env::args().skip(1))?;
    let filter = build_filter(args.filter, args.region.as_deref())?;
    let filter_ref = filter.as_ref();
    let window = Some(VariantWindow {
        start: 0,
        len: args.max_variants,
    });

    for repeat in 0..args.repeats {
        let start = Instant::now();
        let summary = match args.scenario {
            Scenario::Dense => {
                let matrix = genoio_io::read_vcf_dense_windowed_with_threads(
                    &args.path,
                    args.samples.as_deref(),
                    filter_ref,
                    window,
                    false,
                    args.threads,
                )
                .map_err(|error| format!("failed to read dense VCF: {error}"))?;
                ProfileSummary {
                    rows: matrix.n_samples,
                    cols: matrix.n_variants,
                    checksum: matrix.values.iter().copied().map(f64::from).sum(),
                    missing: matrix.missing_mask.iter().filter(|value| **value).count(),
                }
            }
            Scenario::Dosage => {
                let matrix = genoio_io::read_vcf_dosage_dense_windowed_with_threads(
                    &args.path,
                    args.samples.as_deref(),
                    filter_ref,
                    window,
                    false,
                    args.threads,
                )
                .map_err(|error| format!("failed to read VCF dosage: {error}"))?;
                ProfileSummary {
                    rows: matrix.n_samples,
                    cols: matrix.n_variants,
                    checksum: matrix.values.iter().copied().map(f64::from).sum(),
                    missing: matrix.missing_mask.iter().filter(|value| **value).count(),
                }
            }
            Scenario::Sparse => {
                let matrix = genoio_io::read_vcf_sparse_windowed_with_threads(
                    &args.path,
                    args.samples.as_deref(),
                    filter_ref,
                    window,
                    args.threads,
                )
                .map_err(|error| format!("failed to read sparse VCF: {error}"))?;
                ProfileSummary {
                    rows: matrix.n_rows,
                    cols: matrix.n_cols,
                    checksum: matrix.data.iter().copied().map(f64::from).sum(),
                    missing: 0,
                }
            }
            Scenario::HaploDense => {
                let matrix = genoio_io::read_vcf_haplotypes_dense_windowed_with_threads(
                    &args.path,
                    args.samples.as_deref(),
                    filter_ref,
                    window,
                    false,
                    args.threads,
                )
                .map_err(|error| format!("failed to read dense haplotype VCF: {error}"))?;
                ProfileSummary {
                    rows: matrix.n_samples,
                    cols: matrix.n_variants,
                    checksum: matrix.values.iter().copied().map(f64::from).sum(),
                    missing: matrix.missing_mask.iter().filter(|value| **value).count(),
                }
            }
            Scenario::HaploSparse => {
                let matrix = genoio_io::read_vcf_haplotypes_sparse_windowed_with_threads(
                    &args.path,
                    args.samples.as_deref(),
                    filter_ref,
                    window,
                    args.threads,
                )
                .map_err(|error| format!("failed to read sparse haplotype VCF: {error}"))?;
                ProfileSummary {
                    rows: matrix.n_rows,
                    cols: matrix.n_cols,
                    checksum: matrix.data.iter().copied().map(f64::from).sum(),
                    missing: 0,
                }
            }
        };
        let elapsed = start.elapsed();
        println!(
            "repeat={} scenario={:?} threads={} shape=({}, {}) checksum={:.6} missing={} elapsed={:.6}s",
            repeat + 1,
            args.scenario,
            args.threads
                .map(|value| value.to_string())
                .unwrap_or_else(|| "none".to_string()),
            summary.rows,
            summary.cols,
            summary.checksum,
            summary.missing,
            elapsed.as_secs_f64(),
        );
    }

    Ok(())
}

#[derive(Debug)]
struct ProfileSummary {
    rows: usize,
    cols: usize,
    checksum: f64,
    missing: usize,
}

fn build_filter(mode: FilterMode, region: Option<&str>) -> Result<Option<VariantFilter>, String> {
    let base = match mode {
        FilterMode::Biallelic => Some(json!({
            "op": "predicate",
            "name": "biallelic",
            "params": {},
        })),
        FilterMode::Snp => Some(json!({
            "op": "predicate",
            "name": "snp",
            "params": {},
        })),
        FilterMode::QualMin(min) => Some(json!({
            "op": "predicate",
            "name": "qual",
            "params": {"min": min},
        })),
        FilterMode::MafMax(max) => Some(json!({
            "op": "predicate",
            "name": "maf",
            "params": {"max": max},
        })),
        FilterMode::None => None,
    };
    let region = region.map(|value| {
        json!({
            "op": "predicate",
            "name": "region",
            "params": {"value": value},
        })
    });
    let filter = match (base, region) {
        (Some(left), Some(right)) => Some(json!({
            "op": "and",
            "left": left,
            "right": right,
        })),
        (Some(filter), None) | (None, Some(filter)) => Some(filter),
        (None, None) => None,
    };
    filter
        .map(VariantFilter::from_json_value)
        .transpose()
        .map_err(|error| format!("failed to construct variant filter: {error}"))
}

fn parse_args(args: impl Iterator<Item = String>) -> Result<Args, String> {
    let mut path = None;
    let mut max_variants = 10_000;
    let mut repeats = 5;
    let mut threads = None;
    let mut region = None;
    let mut samples = None;
    let mut scenario = Scenario::Dense;
    let mut filter = FilterMode::Biallelic;
    let mut iter = args.peekable();

    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--max-variants" => {
                max_variants = parse_positive_usize("--max-variants", iter.next())?;
            }
            "--repeats" => {
                repeats = parse_positive_usize("--repeats", iter.next())?;
            }
            "--threads" => {
                threads = Some(parse_positive_usize("--threads", iter.next())?);
            }
            "--region" => {
                region = Some(parse_region_arg(iter.next())?);
            }
            "--samples" => {
                samples = Some(parse_samples_arg(iter.next())?);
            }
            "--scenario" => {
                scenario = parse_scenario(iter.next())?;
            }
            "--filter" => {
                filter = parse_filter(iter.next())?;
            }
            "--help" | "-h" => return Err(usage()),
            value if value.starts_with('-') => {
                return Err(format!("unknown argument: {value}\n\n{}", usage()))
            }
            value => {
                if path.is_some() {
                    return Err(format!(
                        "unexpected positional argument: {value}\n\n{}",
                        usage()
                    ));
                }
                path = Some(PathBuf::from(value));
            }
        }
    }

    let path = path.ok_or_else(usage)?;
    Ok(Args {
        path,
        max_variants,
        repeats,
        threads,
        region,
        samples,
        scenario,
        filter,
    })
}

fn parse_positive_usize(name: &str, value: Option<String>) -> Result<usize, String> {
    let value = value.ok_or_else(|| format!("{name} requires a value\n\n{}", usage()))?;
    let parsed = value
        .parse::<usize>()
        .map_err(|error| format!("{name} must be a positive integer: {error}"))?;
    if parsed == 0 {
        return Err(format!("{name} must be a positive integer"));
    }
    Ok(parsed)
}

fn parse_scenario(value: Option<String>) -> Result<Scenario, String> {
    match value.as_deref() {
        Some("dense") => Ok(Scenario::Dense),
        Some("dosage") => Ok(Scenario::Dosage),
        Some("sparse") => Ok(Scenario::Sparse),
        Some("haplo-dense") => Ok(Scenario::HaploDense),
        Some("haplo-sparse") => Ok(Scenario::HaploSparse),
        Some(other) => Err(format!(
            "--scenario must be one of dense, dosage, sparse, haplo-dense, haplo-sparse; got {other:?}"
        )),
        None => Err(format!("--scenario requires a value\n\n{}", usage())),
    }
}

fn parse_region_arg(value: Option<String>) -> Result<String, String> {
    let value = value.ok_or_else(|| format!("--region requires a value\n\n{}", usage()))?;
    if value.split_once(':').is_none() || value.split_once('-').is_none() {
        return Err(format!(
            "--region must use chrom:start-end coordinates; got {value:?}"
        ));
    }
    Ok(value)
}

fn parse_samples_arg(value: Option<String>) -> Result<Vec<String>, String> {
    let value = value.ok_or_else(|| format!("--samples requires a value\n\n{}", usage()))?;
    let samples = value
        .split(',')
        .filter(|sample| !sample.is_empty())
        .map(str::to_string)
        .collect::<Vec<_>>();
    if samples.is_empty() {
        return Err("--samples must contain at least one sample ID".to_string());
    }
    Ok(samples)
}

fn parse_filter(value: Option<String>) -> Result<FilterMode, String> {
    let value = value.ok_or_else(|| format!("--filter requires a value\n\n{}", usage()))?;
    match value.as_str() {
        "biallelic" => Ok(FilterMode::Biallelic),
        "snp" => Ok(FilterMode::Snp),
        "none" => Ok(FilterMode::None),
        value if value.starts_with("qual-min:") => {
            let threshold = parse_filter_f64("--filter qual-min", value.strip_prefix("qual-min:"))?;
            Ok(FilterMode::QualMin(threshold))
        }
        value if value.starts_with("maf-max:") => {
            let threshold = parse_filter_f64("--filter maf-max", value.strip_prefix("maf-max:"))?;
            Ok(FilterMode::MafMax(threshold))
        }
        other => Err(format!(
            "--filter must be one of biallelic, snp, none, qual-min:<N>, maf-max:<N>; got {other:?}"
        )),
    }
}

fn parse_filter_f64(name: &str, value: Option<&str>) -> Result<f64, String> {
    let value = value.ok_or_else(|| format!("{name} requires a value"))?;
    value
        .parse::<f64>()
        .map_err(|error| format!("{name} must be a number: {error}"))
}

fn usage() -> String {
    "usage: vcf_profile <vcf-path> [--max-variants N] [--repeats N] [--threads N] [--region chrom:start-end] [--samples S1,S2] [--scenario dense|dosage|sparse|haplo-dense|haplo-sparse] [--filter biallelic|snp|none|qual-min:N|maf-max:N]".to_string()
}
