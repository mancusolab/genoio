// pattern: Imperative Shell

use std::fs;
use std::path::{Path, PathBuf};

use genoio_core::{DenseLayout, VariantWindow};

mod common;

use common::dense::assert_values_with_nan;
use common::legacy as legacy_io;
use common::unique_dir;
use common::vcf_arrow as genoio_io;
use common::vcf_arrow::{variant_ids, variants};

fn write_text(path: &Path, contents: &str) {
    fs::write(path, contents).expect("test fixture should be written");
}

fn write_vcf(path: &Path) {
    write_text(
        path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1
2\t30\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/0
1\t40\trs4\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/0
",
    );
}

fn write_vcf_with_invalid_record_after_first(path: &Path) {
    write_text(
        path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t20\tbad\tC\tT,G\t.\tPASS\t.\tGT\t0/1\t1/2
",
    );
}

fn write_vcf_with_metadata_accepted_invalid_first_record(path: &Path) {
    write_text(
        path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\tmetadata_passes\tA\tG,T\t30\tPASS\t.\tGT\t0/0\t1/2
1\t20\tstats_pass\tC\tT\t10\tPASS\t.\tGT\t0/1\t1/1
",
    );
}

fn write_plink_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let bed = dir.join("tiny.bed");
    let bim = dir.join("tiny.bim");
    let fam = dir.join("tiny.fam");
    fs::write(&bed, [0x6c, 0x1b, 0x01, 0x04, 0x0d, 0x03, 0x00])
        .expect("bed fixture should be written");
    write_text(
        &bim,
        "\
1 rs1 0 10 G A
1 rs2 0 20 T C
2 rs3 0 30 A G
1 rs4 0 40 C T
",
    );
    write_text(
        &fam,
        "\
F1 S1 0 0 1 -9
F1 S2 0 0 2 -9
",
    );
    (bed, bim, fam)
}

fn write_plink_source_window_stop_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let bed = dir.join("tiny.bed");
    let bim = dir.join("tiny.bim");
    let fam = dir.join("tiny.fam");
    fs::write(&bed, [0x6c, 0x1b, 0x01, 0x04, 0x0d, 0x03]).expect("bed fixture should be written");
    write_text(
        &bim,
        "\
1 rs1 0 10 A G
1 rs2 0 20 C T
1 bad malformed
",
    );
    write_text(
        &fam,
        "\
F1 S1 0 0 1 -9
F1 S2 0 0 2 -9
",
    );
    (bed, bim, fam)
}

fn fixed_width_pgen(records: &[u8], n_samples: u32, n_variants: u32) -> Vec<u8> {
    let mut bytes = vec![0x6c, 0x1b, 0x02];
    bytes.extend(n_variants.to_le_bytes());
    bytes.extend(n_samples.to_le_bytes());
    bytes.push(0);
    bytes.extend(records);
    bytes
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

fn write_plink2_filter_window_stop_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let pgen = dir.join("tiny.pgen");
    let pvar = dir.join("tiny.pvar");
    let psam = dir.join("tiny.psam");
    fs::write(&pgen, fixed_width_pgen(&[0x11, 0x11, 0x00], 3, 3))
        .expect("pgen fixture should be written");
    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 C T
1 bad malformed A G
",
    );
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
",
    );
    (pgen, pvar, psam)
}

#[test]
fn vcf_dense_window_stops_after_requested_retained_variants() {
    let dir = unique_dir("vcf-dense-window-stop");
    let path = dir.join("blocks.vcf");
    write_vcf_with_invalid_record_after_first(&path);

    let block = genoio_io::read_vcf_dense_windowed(
        &path,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect("windowed vcf should stop before later invalid records");

    assert_eq!(variant_ids(variants(&block.variants)), vec!["rs1"]);
    assert_eq!(block.values, vec![0.0, 1.0]);
    assert_eq!(block.diagnostics.candidate_variants, 1);
}

#[test]
fn vcf_sparse_window_stops_after_requested_retained_variants() {
    let dir = unique_dir("vcf-sparse-window-stop");
    let path = dir.join("blocks.vcf");
    write_vcf_with_invalid_record_after_first(&path);

    let block = genoio_io::read_vcf_sparse_windowed(
        &path,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
    )
    .expect("windowed sparse vcf should stop before later invalid records");

    assert_eq!(variant_ids(variants(&block.variants)), vec!["rs1"]);
    assert_eq!(variants(&block.variants).len(), 1);
    assert_eq!(block.diagnostics.candidate_variants, 1);
}

#[test]
fn vcf_dense_window_uses_retained_variant_order_after_filters() {
    let dir = unique_dir("vcf-block-window");
    let path = dir.join("blocks.vcf");
    write_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "1"}
    }))
    .expect("filter should parse");

    let block = genoio_io::read_vcf_dense_windowed(
        &path,
        None,
        Some(&filter),
        Some(VariantWindow { start: 1, len: 2 }),
        false,
    )
    .expect("windowed vcf should decode");

    assert_eq!(variant_ids(variants(&block.variants)), vec!["rs2", "rs4"]);
    assert_eq!(block.values, vec![1.0, 2.0, 0.0, 0.0]);
    assert_eq!(block.layout, DenseLayout::VariantMajor);
}

#[test]
fn vcf_dense_window_skips_pre_window_metadata_accepted_genotypes() {
    let dir = unique_dir("vcf-partial-filter-window");
    let path = dir.join("blocks.vcf");
    write_vcf_with_metadata_accepted_invalid_first_record(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "or",
        "left": {"op": "predicate", "name": "qual", "params": {"min": 20.0}},
        "right": {"op": "predicate", "name": "maf", "params": {"min": 0.1}}
    }))
    .expect("filter should parse");

    let block = genoio_io::read_vcf_dense_windowed(
        &path,
        None,
        Some(&filter),
        Some(VariantWindow { start: 1, len: 1 }),
        false,
    )
    .expect("windowed vcf should skip pre-window metadata-accepted genotypes");

    assert_eq!(variant_ids(variants(&block.variants)), vec!["stats_pass"]);
    assert_eq!(block.values, vec![1.0, 2.0]);
    assert_eq!(block.diagnostics.candidate_variants, 2);
    assert_eq!(block.diagnostics.retained_variants, 1);
}

#[test]
fn plink1_dense_window_uses_retained_variant_order_after_filters() {
    let dir = unique_dir("plink-block-window");
    let (bed, bim, fam) = write_plink_fixture(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "1"}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink1_dense_windowed(
        &bed,
        &bim,
        &fam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 1, len: 2 }),
        false,
    )
    .expect("windowed plink should decode");

    assert_eq!(
        block
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs2", "rs4"]
    );
    assert_values_with_nan(&block.values, &[f32::NAN, 2.0, 0.0, 2.0]);
}

#[test]
fn plink1_dense_unfiltered_window_stops_after_requested_source_variants() {
    let dir = unique_dir("plink1-source-window-stop");
    let (bed, bim, fam) = write_plink_source_window_stop_fixture(&dir);

    let block = legacy_io::read_plink1_dense_windowed(
        &bed,
        &bim,
        &fam,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect("unfiltered plink1 source window should stop before later malformed metadata");

    assert_eq!(
        block
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs1"]
    );
    assert_values_with_nan(&block.values, &[2.0, f32::NAN]);
    assert_eq!(block.diagnostics.candidate_variants, 1);
}

#[test]
fn plink1_matrix_only_source_window_skips_bim_rows() {
    let dir = unique_dir("plink1-matrix-only-source-window");
    let (bed, bim, fam) = write_plink_source_window_stop_fixture(&dir);
    write_text(
        &bim,
        "\
malformed
",
    );

    let block = legacy_io::read_plink1_dense_windowed(
        &bed,
        &bim,
        &fam,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect("matrix-only plink1 source window should not parse bim rows");

    assert_eq!(block.n_variants, 1);
    assert!(block.variants.is_empty());
    assert_values_with_nan(&block.values, &[2.0, f32::NAN]);
    assert_eq!(block.diagnostics.candidate_variants, 1);
}

#[test]
fn plink1_dense_impossible_filter_returns_empty_without_parsing_bim_variants() {
    let dir = unique_dir("plink1-impossible-filter");
    let (bed, bim, fam) = write_plink_source_window_stop_fixture(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "2"}},
        "right": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink1_dense_windowed(
        &bed,
        &bim,
        &fam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect("impossible plink1 filter should not parse malformed bim rows");

    assert_eq!(block.n_variants, 0);
    assert_eq!(block.diagnostics.candidate_variants, 0);
}

#[test]
fn plink1_matrix_only_genotype_filter_window_skips_bim_rows() {
    let dir = unique_dir("plink1-genotype-filter-skips-bim");
    let (bed, bim, fam) = write_plink_source_window_stop_fixture(&dir);
    write_text(
        &bim,
        "\
malformed
",
    );
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink1_dense_windowed(
        &bed,
        &bim,
        &fam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect("matrix-only genotype filter should not parse bim rows");

    assert_eq!(block.n_variants, 1);
    assert!(block.variants.is_empty());
    assert_eq!(block.values, vec![0.0, 2.0]);
    assert_eq!(block.diagnostics.candidate_variants, 3);
}

#[test]
fn plink2_dense_filtered_window_stops_after_requested_retained_variants() {
    let dir = unique_dir("plink2-filtered-window-stop");
    let (pgen, pvar, psam) = write_plink2_filter_window_stop_fixture(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect("windowed filtered plink2 should stop before later malformed metadata");

    assert_eq!(
        block
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs1"]
    );
    assert_eq!(block.values, vec![1.0, 0.0, 1.0]);
    assert_eq!(block.diagnostics.candidate_variants, 1);
}

#[test]
fn plink2_matrix_only_genotype_filter_window_skips_pvar_rows() {
    let dir = unique_dir("plink2-genotype-filter-skips-pvar");
    let (pgen, pvar, psam) = write_plink2_filter_window_stop_fixture(&dir);
    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT
malformed
",
    );
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect("matrix-only genotype filter should not parse pvar rows");

    assert_eq!(block.n_variants, 1);
    assert!(block.variants.is_empty());
    assert_eq!(block.values, vec![1.0, 0.0, 1.0]);
    assert_eq!(block.diagnostics.candidate_variants, 1);
}

#[test]
fn plink2_matrix_only_genotype_filter_extends_variable_width_prefix() {
    let dir = unique_dir("plink2-genotype-filter-prefix-extend");
    let pgen = dir.join("tiny.pgen");
    let pvar = dir.join("tiny.pvar");
    let psam = dir.join("tiny.psam");
    let all_hom_ref = [0x00];
    let all_het = [0x15];
    fs::write(
        &pgen,
        variable_width_pgen(&[0, 0], &[&all_hom_ref, &all_het], 3),
    )
    .expect("pgen fixture should be written");
    write_text(
        &pvar,
        "exists but matrix-only genotype filters do not parse this\n",
    );
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
",
    );
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect("matrix-only genotype filter should extend the prefix when needed");

    assert_eq!(block.n_variants, 1);
    assert_eq!(block.values, vec![1.0, 1.0, 1.0]);
    assert_eq!(block.diagnostics.candidate_variants, 2);
    assert_eq!(block.diagnostics.dropped_genotype_variants, 1);
}

#[test]
fn plink2_matrix_only_genotype_filter_prefix_ignores_later_unsupported_records() {
    let dir = unique_dir("plink2-genotype-filter-prefix-ignore-later");
    let pgen = dir.join("tiny.pgen");
    let pvar = dir.join("tiny.pvar");
    let psam = dir.join("tiny.psam");
    let all_het = [0x15];
    fs::write(&pgen, variable_width_pgen(&[0, 5], &[&all_het, &[0]], 3))
        .expect("pgen fixture should be written");
    write_text(
        &pvar,
        "exists but matrix-only genotype filters do not parse this\n",
    );
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
",
    );
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect("matrix-only genotype filter should not validate unused later records");

    assert_eq!(block.n_variants, 1);
    assert_eq!(block.values, vec![1.0, 1.0, 1.0]);
    assert_eq!(block.diagnostics.candidate_variants, 1);
}

#[test]
fn plink2_sparse_filtered_window_stops_after_requested_retained_variants() {
    let dir = unique_dir("plink2-sparse-filtered-window-stop");
    let (pgen, pvar, psam) = write_plink2_filter_window_stop_fixture(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink2_sparse_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
    )
    .expect("windowed sparse filtered plink2 should stop before later malformed metadata");

    assert_eq!(
        block
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs1"]
    );
    assert_eq!(block.diagnostics.candidate_variants, 1);
}

#[test]
fn plink2_dense_impossible_filter_returns_empty_without_parsing_pvar_variants() {
    let dir = unique_dir("plink2-impossible-filter");
    let (pgen, pvar, psam) = write_plink2_filter_window_stop_fixture(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "2"}},
        "right": {"op": "predicate", "name": "region", "params": {"value": "1:10-20"}}
    }))
    .expect("filter should parse");

    let block = legacy_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect("impossible plink2 filter should not parse malformed pvar rows");

    assert_eq!(block.n_variants, 0);
    assert_eq!(block.diagnostics.candidate_variants, 0);
}
