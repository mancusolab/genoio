use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{self, Format, Header, Writer};

fn unique_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock should be after unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("genoio-{name}-{nanos}"));
    fs::create_dir(&dir).expect("test temp dir should be created");
    dir
}

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

    let index_path = PathBuf::from(format!("{}.tbi", path.to_string_lossy()));
    bcf::index::build(path, Some(index_path.as_path()), 1, bcf::index::Type::Tbx)
        .expect("tabix index should build");
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
