// pattern: Imperative Shell

use std::fs;
use std::path::{Path, PathBuf};

const PGEN_DOSAGE_TOLERANCE: f32 = 2.0 / 32768.0;

mod common;

use common::unique_dir;

fn write_vcf(path: &Path) {
    fs::write(
        path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/0\t./.\t0/0
1\t30\trs3\tG\tA\t.\tPASS\t.\tGT\t./.\t./.\t./.
",
    )
    .expect("vcf fixture should be written");
}

fn write_ds_vcf(path: &Path) {
    fs::write(
        path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected alternate allele dosage\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DS\t0/0:0.2\t0/1:1.4\t1/1:1.8
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT:DS\t0/0:0\t0/0:.\t0/1:0.7
",
    )
    .expect("dosage vcf fixture should be written");
}

fn write_fixed_width_plink2(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let pgen = dir.join("tiny.pgen");
    let pvar = dir.join("tiny.pvar");
    let psam = dir.join("tiny.psam");
    fs::write(
        &pgen,
        [
            0x6c, 0x1b, 0x02, 0x04, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x2c, 0x11,
            0x06, 0x3f,
        ],
    )
    .expect("pgen fixture should be written");
    fs::write(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 C T
1 30 rs3 G A
1 40 rs4 T C
",
    )
    .expect("pvar fixture should be written");
    fs::write(
        &psam,
        "\
#IID
S1
S2
S3
",
    )
    .expect("psam fixture should be written");
    (pgen, pvar, psam)
}

fn scaled_dosage(value: f32) -> [u8; 2] {
    let raw = (value / 2.0 * 32768.0).round() as u16;
    raw.to_le_bytes()
}

fn assert_pgen_dosage_close(observed: f32, expected: f32) {
    assert!(
        (observed - expected).abs() <= PGEN_DOSAGE_TOLERANCE,
        "observed dosage {observed} differs from expected {expected}"
    );
}

fn write_fixed_width_plink2_dosage(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let pgen = dir.join("dosage.pgen");
    let pvar = dir.join("dosage.pvar");
    let psam = dir.join("dosage.psam");
    let mut pgen_bytes = vec![0x6c, 0x1b, 0x03];
    pgen_bytes.extend(2_u32.to_le_bytes());
    pgen_bytes.extend(3_u32.to_le_bytes());
    pgen_bytes.push(0);
    pgen_bytes.push(0x24);
    pgen_bytes.extend(scaled_dosage(0.2));
    pgen_bytes.extend(scaled_dosage(1.4));
    pgen_bytes.extend(scaled_dosage(1.8));
    pgen_bytes.push(0x0c);
    pgen_bytes.extend(scaled_dosage(0.0));
    pgen_bytes.extend(u16::MAX.to_le_bytes());
    pgen_bytes.extend(scaled_dosage(0.7));
    fs::write(&pgen, pgen_bytes).expect("dosage pgen fixture should be written");
    fs::write(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 C T
",
    )
    .expect("pvar fixture should be written");
    fs::write(
        &psam,
        "\
#IID
S1
S2
S3
",
    )
    .expect("psam fixture should be written");
    (pgen, pvar, psam)
}

fn variable_width_pgen(record_types: &[u8], records: &[&[u8]], n_samples: u32) -> Vec<u8> {
    let n_variants = u32::try_from(records.len()).expect("test variant count fits u32");
    let header_len = 12 + 8 + record_types.len() + records.len();
    let mut bytes = vec![0x6c, 0x1b, 0x10];
    bytes.extend(n_variants.to_le_bytes());
    bytes.extend(n_samples.to_le_bytes());
    bytes.push(0x04);
    bytes.extend(
        u64::try_from(header_len)
            .expect("test header length fits u64")
            .to_le_bytes(),
    );
    bytes.extend(record_types);
    bytes.extend(
        records
            .iter()
            .map(|record| u8::try_from(record.len()).expect("test record length fits one byte")),
    );
    for record in records {
        bytes.extend(*record);
    }
    bytes
}

fn write_variable_width_plink2(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let pgen = dir.join("variable.pgen");
    let pvar = dir.join("variable.pvar");
    let psam = dir.join("variable.psam");
    let pgen_bytes = variable_width_pgen(
        &[0, 4, 1, 2, 3],
        &[&[0xe4], &[2, 1, 9, 2], &[2, 5, 0], &[1, 1, 3], &[0]],
        4,
    );
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    fs::write(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 A G
1 30 rs3 A G
1 40 rs4 A G
1 50 rs5 A G
",
    )
    .expect("pvar fixture should be written");
    fs::write(
        &psam,
        "\
#IID
S1
S2
S3
S4
",
    )
    .expect("psam fixture should be written");
    (pgen, pvar, psam)
}

#[test]
fn filter_genotype_stats_plink2_dosage_uses_fractional_mac() {
    let dir = unique_dir("plink2-dosage-filter-genotype");
    let (pgen, pvar, psam) = write_fixed_width_plink2_dosage(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "mac",
        "params": {"max": 2}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_plink2_dosage_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        None,
        false,
    )
    .expect("plink2 dosage should filter");

    assert_eq!(dense.n_variants, 1);
    assert_eq!(dense.variants[0].id, "rs2");
    assert_eq!(dense.values[..2], [0.0, 0.0]);
    assert_pgen_dosage_close(dense.values[2], 0.7);
    assert_eq!(dense.missing_mask, vec![false, true, false]);
    assert!(dense.variants[0]
        .af
        .is_some_and(|af| (af - 0.175).abs() <= PGEN_DOSAGE_TOLERANCE));
    assert!(dense.variants[0]
        .maf
        .is_some_and(|maf| (maf - 0.175).abs() <= PGEN_DOSAGE_TOLERANCE));
    assert_eq!(dense.variants[0].mac, None);
    assert_eq!(dense.variants[0].missing_rate, Some(1.0 / 3.0));
    assert_eq!(dense.variants[0].n_called, Some(2));
    assert_eq!(dense.diagnostics.candidate_variants, 2);
    assert_eq!(dense.diagnostics.retained_variants, 1);
    assert_eq!(dense.diagnostics.dropped_genotype_variants, 1);
}

#[test]
fn filter_genotype_stats_vcf_dosage_uses_fractional_mac() {
    let dir = unique_dir("vcf-dosage-filter-genotype");
    let path = dir.join("dosage.vcf");
    write_ds_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "mac",
        "params": {"max": 2}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_vcf_dosage_dense_windowed(&path, None, Some(&filter), None, false)
        .expect("vcf dosage should filter");

    assert_eq!(dense.n_variants, 1);
    assert_eq!(dense.variants[0].id, "rs2");
    assert_eq!(dense.values, vec![0.0, 0.0, 0.7]);
    assert_eq!(dense.missing_mask, vec![false, true, false]);
    assert_eq!(dense.variants[0].af, Some(0.175));
    assert_eq!(dense.variants[0].maf, Some(0.175));
    assert_eq!(dense.variants[0].mac, None);
    assert_eq!(dense.variants[0].missing_rate, Some(1.0 / 3.0));
    assert_eq!(dense.variants[0].n_called, Some(2));
    assert_eq!(dense.diagnostics.candidate_variants, 2);
    assert_eq!(dense.diagnostics.retained_variants, 1);
    assert_eq!(dense.diagnostics.dropped_genotype_variants, 1);

    let matrix_only =
        genoio_io::read_vcf_dosage_dense_windowed(&path, None, Some(&filter), None, true)
            .expect("matrix-only vcf dosage should filter");

    assert_eq!(matrix_only.n_samples, 3);
    assert_eq!(matrix_only.n_variants, 1);
    assert!(matrix_only.samples.is_empty());
    assert!(matrix_only.variants.is_empty());
    assert_eq!(matrix_only.values, vec![0.0, 0.0, 0.7]);
    assert_eq!(matrix_only.missing_mask, vec![false, true, false]);
}

#[test]
fn filter_genotype_stats_use_called_genotypes_before_missing_imputation() {
    let dir = unique_dir("vcf-filter-genotype");
    let path = dir.join("tiny.vcf");
    write_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "maf", "params": {"min": 0.2}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.5}}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_vcf_dense(&path, None, Some(&filter)).expect("vcf should filter");

    assert_eq!(dense.n_variants, 1);
    assert_eq!(dense.variants[0].id, "rs1");
    assert_eq!(dense.values, vec![0.0, 1.0, 2.0]);
    assert_eq!(dense.variants[0].af, Some(0.5));
    assert_eq!(dense.variants[0].maf, Some(0.5));
    assert_eq!(dense.variants[0].mac, Some(3));
    assert_eq!(dense.variants[0].missing_rate, Some(0.0));
    assert_eq!(dense.variants[0].n_called, Some(3));
    assert_eq!(dense.diagnostics.candidate_variants, 3);
    assert_eq!(dense.diagnostics.retained_variants, 1);
    assert_eq!(dense.diagnostics.dropped_genotype_variants, 2);

    let matrix_only = genoio_io::read_vcf_dense_windowed(&path, None, Some(&filter), None, true)
        .expect("matrix-only vcf GT should filter");

    assert_eq!(matrix_only.n_samples, 3);
    assert_eq!(matrix_only.n_variants, 1);
    assert!(matrix_only.samples.is_empty());
    assert!(matrix_only.variants.is_empty());
    assert_eq!(matrix_only.values, vec![0.0, 1.0, 2.0]);
    assert_eq!(matrix_only.missing_mask, vec![false, false, false]);
}

#[test]
fn filter_genotype_stats_plink2_match_expanded_stats_and_attach_metadata() {
    let dir = unique_dir("plink2-filter-genotype");
    let (pgen, pvar, psam) = write_fixed_width_plink2(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "maf", "params": {"min": 0.2}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.5}}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, Some(&filter))
        .expect("plink2 should filter");

    assert_eq!(dense.n_variants, 3);
    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs1", "rs2", "rs3"]
    );
    assert_eq!(
        dense.values,
        vec![0.0, 1.0, 2.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0]
    );
    assert_eq!(dense.variants[0].af, Some(0.5));
    assert_eq!(dense.variants[0].maf, Some(0.5));
    assert_eq!(dense.variants[0].mac, Some(2));
    assert_eq!(dense.variants[0].missing_rate, Some(1.0 / 3.0));
    assert_eq!(dense.variants[0].n_called, Some(2));
    assert_eq!(dense.variants[1].af, Some(1.0 / 3.0));
    assert_eq!(dense.variants[1].maf, Some(1.0 / 3.0));
    assert_eq!(dense.variants[1].mac, Some(2));
    assert_eq!(dense.variants[1].missing_rate, Some(0.0));
    assert_eq!(dense.variants[1].n_called, Some(3));
    assert_eq!(dense.variants[2].af, Some(0.5));
    assert_eq!(dense.variants[2].maf, Some(0.5));
    assert_eq!(dense.variants[2].mac, Some(3));
    assert_eq!(dense.variants[2].missing_rate, Some(0.0));
    assert_eq!(dense.variants[2].n_called, Some(3));
    assert_eq!(dense.diagnostics.candidate_variants, 4);
    assert_eq!(dense.diagnostics.retained_variants, 3);
    assert_eq!(dense.diagnostics.dropped_genotype_variants, 1);
}

#[test]
fn filter_genotype_stats_plink2_variable_width_selected_samples_attach_stats() {
    let dir = unique_dir("plink2-variable-filter-genotype");
    let (pgen, pvar, psam) = write_variable_width_plink2(&dir);
    let keep = vec!["S4".to_string(), "S2".to_string()];
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "maf", "params": {"min": 0.2}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.5}}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, Some(&keep), Some(&filter))
        .expect("variable-width plink2 should filter");

    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S2", "S4"]
    );
    assert_eq!(dense.n_variants, 2);
    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs1", "rs2"]
    );
    assert_eq!(dense.values, vec![1.0, 1.0, 0.0, 2.0]);
    assert_eq!(dense.missing_mask, vec![false, false, true, false]);
    assert_eq!(dense.variants[0].af, Some(0.5));
    assert_eq!(dense.variants[0].maf, Some(0.5));
    assert_eq!(dense.variants[0].mac, Some(1));
    assert_eq!(dense.variants[0].missing_rate, Some(0.5));
    assert_eq!(dense.variants[0].n_called, Some(1));
    assert_eq!(dense.variants[1].af, Some(0.75));
    assert_eq!(dense.variants[1].maf, Some(0.25));
    assert_eq!(dense.variants[1].mac, Some(1));
    assert_eq!(dense.variants[1].missing_rate, Some(0.0));
    assert_eq!(dense.variants[1].n_called, Some(2));
    assert_eq!(dense.diagnostics.candidate_variants, 5);
    assert_eq!(dense.diagnostics.retained_variants, 2);
    assert_eq!(dense.diagnostics.dropped_genotype_variants, 3);
}

#[test]
fn filter_genotype_stats_plink2_sparse_keeps_dense_filter_semantics() {
    let dir = unique_dir("plink2-sparse-filter-genotype");
    let (pgen, pvar, psam) = write_fixed_width_plink2(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "polymorphic", "params": {}},
        "right": {"op": "predicate", "name": "missing_rate", "params": {"max": 0.0}}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, Some(&filter))
        .expect("dense plink2 should filter");
    let sparse = genoio_io::read_plink2_sparse(&pgen, &pvar, &psam, None, Some(&filter))
        .expect("sparse plink2 should filter");

    assert_eq!(
        sparse
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>()
    );
    assert_eq!(sparse.diagnostics.candidate_variants, 4);
    assert_eq!(sparse.diagnostics.retained_variants, 2);
    assert_eq!(sparse.diagnostics.dropped_genotype_variants, 2);
}
