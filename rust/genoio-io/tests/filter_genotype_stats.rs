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
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/0\t./.\t0/0
1\t30\trs3\tG\tA\t.\tPASS\t.\tGT\t./.\t./.\t./.
",
    )
    .expect("vcf fixture should be written");
}

#[test]
fn genotype_filters_use_called_genotypes_before_missing_imputation() {
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
}
