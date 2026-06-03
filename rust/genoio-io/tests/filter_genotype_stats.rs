// pattern: Imperative Shell

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
