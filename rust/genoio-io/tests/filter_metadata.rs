// pattern: Imperative Shell

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{self, Format, Header, Writer};

mod common;

use common::unique_dir;

fn write_vcf(path: &Path) {
    fs::write(
        path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1
2\t30\tindel1\tAT\tA\t.\tPASS\t.\tGT\t0/1
",
    )
    .expect("vcf fixture should be written");
}

fn write_indexed_vcf(path: &Path) {
    let mut header = Header::new();
    header.push_record(br#"##fileformat=VCFv4.2"#);
    header.push_record(br#"##contig=<ID=1>"#);
    header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    header.push_sample(b"S1");

    {
        let mut writer =
            Writer::from_path(path, &header, false, Format::Vcf).expect("vcf writer should open");
        for (id, pos) in [("rs10", 10_i64), ("rs20", 20), ("rs30", 30), ("rs40", 40)] {
            let mut record = writer.empty_record();
            let rid = writer
                .header()
                .name2rid(b"1")
                .expect("contig should resolve");
            record.set_rid(Some(rid));
            record.set_pos(pos - 1);
            record.set_id(id.as_bytes()).expect("id should be set");
            record
                .set_alleles(&[b"A", b"G"])
                .expect("alleles should be set");
            record
                .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
                .expect("genotype should be set");
            writer.write(&record).expect("record should be written");
        }
    }

    build_tabix_index(path);
}

fn write_bgzf_file(path: &Path, contents: &str) {
    let file = fs::File::create(path).expect("test fixture should be created");
    let mut writer = noodles_bgzf::io::Writer::new(file);
    writer
        .write_all(contents.as_bytes())
        .expect("test fixture should be compressed");
}

fn build_tabix_index(path: &Path) {
    let index_path = tabix_index_path(path);
    bcf::index::build(path, Some(index_path.as_path()), 1, bcf::index::Type::Tbx)
        .expect("tabix index should build");
}

fn tabix_index_path(path: &Path) -> PathBuf {
    PathBuf::from(format!("{}.tbi", path.to_string_lossy()))
}

fn csc_to_dense(sparse: &genoio_core::SparseGenotypeMatrix) -> Vec<f32> {
    let mut dense = vec![0.0; sparse.n_rows * sparse.n_cols];
    for col in 0..sparse.n_cols {
        for offset in sparse.indptr[col]..sparse.indptr[col + 1] {
            let row = sparse.indices[offset];
            dense[row * sparse.n_cols + col] = sparse.data[offset];
        }
    }
    dense
}

fn write_compressed_unindexed_vcf(path: &Path) {
    let mut header = Header::new();
    header.push_record(br#"##fileformat=VCFv4.2"#);
    header.push_record(br#"##contig=<ID=1>"#);
    header.push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    header.push_sample(b"S1");

    let mut writer =
        Writer::from_path(path, &header, false, Format::Vcf).expect("vcf writer should open");
    let mut record = writer.empty_record();
    let rid = writer
        .header()
        .name2rid(b"1")
        .expect("contig should resolve");
    record.set_rid(Some(rid));
    record.set_pos(9);
    record.set_id(b"rs10").expect("id should be set");
    record
        .set_alleles(&[b"A", b"G"])
        .expect("alleles should be set");
    record
        .push_genotypes(&[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)])
        .expect("genotype should be set");
    writer.write(&record).expect("record should be written");
}

#[test]
fn indexed_vcf_region_dosage_uses_permissive_fast_path() {
    let dir = unique_dir("vcf-filter-indexed-dosage-fast");
    let path = dir.join("indexed.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected dosage\"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs10\tA\tG\t.\tPASS\t.\tDS\t0.1\t1.1
1\t20\trs20\tA\tG\t.\tPASS\t.\tDS\t0.2\t1.2
1\t30\trs30\tA\tG\t.\tPASS\t.\tDS\t0.3\t1.3
",
    );
    build_tabix_index(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:20-20"}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_vcf_dosage_dense_windowed(&path, None, Some(&filter), None, false)
        .expect("indexed dosage VCF should decode through permissive fast path");

    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs20"]
    );
    assert_eq!(dense.values, vec![0.2, 1.2]);
}

#[test]
fn indexed_vcf_region_sparse_uses_permissive_fast_path() {
    let dir = unique_dir("vcf-filter-indexed-sparse-fast");
    let path = dir.join("indexed.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t20\trs20\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1
1\t30\trs30\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
",
    );
    build_tabix_index(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:20-20"}
    }))
    .expect("filter should parse");

    let sparse = genoio_io::read_vcf_sparse(&path, None, Some(&filter))
        .expect("indexed sparse VCF should decode through permissive fast path");

    assert_eq!(
        sparse
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs20"]
    );
    assert_eq!(csc_to_dense(&sparse), vec![1.0, 0.0]);
    assert!(sparse.variants[0].flipped);
}

#[test]
fn indexed_vcf_region_haplotypes_use_permissive_fast_path() {
    let dir = unique_dir("vcf-filter-indexed-haplo-fast");
    let path = dir.join("indexed.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1
1\t20\trs20\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0
1\t30\trs30\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1
",
    );
    build_tabix_index(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:20-20"}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_vcf_haplotypes_dense(&path, None, Some(&filter))
        .expect("indexed haplotype dense VCF should decode through permissive fast path");
    let sparse = genoio_io::read_vcf_haplotypes_sparse(&path, None, Some(&filter))
        .expect("indexed haplotype sparse VCF should decode through permissive fast path");

    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs20"]
    );
    assert_eq!(dense.values, vec![0.0, 1.0, 1.0, 0.0]);
    assert_eq!(csc_to_dense(&sparse), dense.values);
}

#[test]
fn vcf_metadata_filters_retain_expected_variants() {
    let dir = unique_dir("vcf-filter-metadata");
    let path = dir.join("tiny.vcf");
    write_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
        "right": {"op": "predicate", "name": "snp", "params": {}}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_vcf_dense(&path, None, Some(&filter)).expect("vcf should filter");

    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs1", "rs2"]
    );
    assert_eq!(dense.diagnostics.candidate_variants, 3);
    assert_eq!(dense.diagnostics.retained_variants, 2);
    assert_eq!(dense.diagnostics.dropped_metadata_variants, 1);
}

#[test]
fn region_filter_includes_start_and_end_positions_once() {
    let dir = unique_dir("vcf-filter-region");
    let path = dir.join("tiny.vcf");
    write_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:10-20"}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_vcf_dense(&path, None, Some(&filter)).expect("vcf should filter");

    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| (variant.id.as_str(), variant.pos))
            .collect::<Vec<_>>(),
        vec![("rs1", 10), ("rs2", 20)]
    );
}

#[test]
fn indexed_vcf_region_filter_fetches_exact_start_and_end_positions() {
    let dir = unique_dir("vcf-filter-indexed-region");
    let path = dir.join("indexed.vcf.gz");
    write_indexed_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:20-30"}
    }))
    .expect("filter should parse");

    let dense =
        genoio_io::read_vcf_dense(&path, None, Some(&filter)).expect("indexed vcf should filter");

    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| (variant.id.as_str(), variant.pos))
            .collect::<Vec<_>>(),
        vec![("rs20", 20), ("rs30", 30)]
    );
    assert_eq!(dense.diagnostics.candidate_variants, 2);
}

#[test]
fn indexed_vcf_region_uses_tabix_reference_names() {
    let dir = unique_dir("vcf-filter-indexed-fast-region");
    let path = dir.join("indexed.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/0
1\t20\trs20\tA\tG\t.\tPASS\t.\tGT\t0/1
1\t21\tbad_outside_region\tA\tG\t.\tPASS\t.\tGT\t0/3
1\t30\trs30\tA\tG\t.\tPASS\t.\tGT\t1/1
",
    );
    build_tabix_index(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:20-20"}
    }))
    .expect("filter should parse");

    let dense =
        genoio_io::read_vcf_dense(&path, None, Some(&filter)).expect("indexed VCF should decode");

    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| (variant.id.as_str(), variant.pos))
            .collect::<Vec<_>>(),
        vec![("rs20", 20)]
    );
    assert_eq!(dense.values, vec![1.0]);
    assert_eq!(dense.diagnostics.candidate_variants, 1);
}

#[test]
fn indexed_vcf_region_absent_contig_returns_empty_dense_matrix() {
    let dir = unique_dir("vcf-filter-indexed-absent-contig");
    let path = dir.join("indexed.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/0
1\t20\trs20\tA\tG\t.\tPASS\t.\tGT\t0/1
",
    );
    build_tabix_index(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "2:1-10"}
    }))
    .expect("filter should parse");

    let dense =
        genoio_io::read_vcf_dense(&path, None, Some(&filter)).expect("indexed VCF should decode");

    assert_eq!(dense.n_samples, 1);
    assert_eq!(dense.n_variants, 0);
    assert!(dense.values.is_empty());
    assert!(dense.variants.is_empty());
}

#[test]
fn indexed_vcf_region_sample_filter_preserves_source_order() {
    let dir = unique_dir("vcf-filter-indexed-samples");
    let path = dir.join("indexed.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs20\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1
",
    );
    build_tabix_index(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:10-20"}
    }))
    .expect("filter should parse");
    let samples = vec!["S3".to_string(), "S1".to_string()];

    let dense = genoio_io::read_vcf_dense(&path, Some(&samples), Some(&filter))
        .expect("indexed VCF should decode");

    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs10", "rs20"]
    );
    assert_eq!(dense.values, vec![0.0, 1.0, 2.0, 2.0]);
}

#[test]
fn indexed_vcf_region_falls_back_when_noodles_rejects_header() {
    let dir = unique_dir("vcf-filter-indexed-header-fallback");
    let path = dir.join("indexed.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/0
1\t20\trs20\tA\tG\t.\tPASS\t.\tGT\t0/1
",
    );
    build_tabix_index(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:20-20"}
    }))
    .expect("filter should parse");

    let dense =
        genoio_io::read_vcf_dense(&path, None, Some(&filter)).expect("indexed VCF should decode");

    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs20"]
    );
    assert_eq!(dense.values, vec![1.0]);
}

#[test]
fn compressed_vcf_region_filter_requires_index() {
    let dir = unique_dir("vcf-filter-unindexed-region");
    let path = dir.join("unindexed.vcf.gz");
    write_compressed_unindexed_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "region",
        "params": {"value": "1:10-10"}
    }))
    .expect("filter should parse");

    let error = genoio_io::read_vcf_dense(&path, None, Some(&filter))
        .expect_err("unindexed compressed region filter should fail");

    assert!(error.to_string().contains("requires an index"));
}

#[test]
fn compressed_vcf_non_pushdown_region_filter_falls_back_to_full_scan() {
    let dir = unique_dir("vcf-filter-unindexed-non-pushdown-region");
    let path = dir.join("unindexed.vcf.gz");
    write_compressed_unindexed_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "not",
        "expr": {
            "op": "predicate",
            "name": "region",
            "params": {"value": "1:20-30"}
        }
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_vcf_dense(&path, None, Some(&filter))
        .expect("non-pushdown region filter should full-scan unindexed compressed vcf");

    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs10"]
    );
    assert_eq!(dense.diagnostics.candidate_variants, 1);
}
