use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

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
