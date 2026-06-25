// pattern: Imperative Shell

use std::fs;
use std::io::Write;

mod common;

use common::unique_dir;
use common::vcf_arrow as genoio_io;
use common::vcf_arrow::{
    dense_values_sample_major_arrow as dense_values_sample_major, sparse_values_dense_arrow,
    variant_a0, variant_a1, variant_ids, variants,
};

fn phased_vcf() -> String {
    "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|0
"
    .to_string()
}

fn mixed_phase_vcf() -> String {
    "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t1/1\t0|0
"
    .to_string()
}

fn mixed_phase_stat_filter_vcf() -> String {
    "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs_phased\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0
1\t20\trs_unphased_monomorphic\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0
"
    .to_string()
}

fn phased_alt_major_vcf() -> String {
    "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs_alt_major\tA\tG\t.\tPASS\t.\tGT\t1|1\t1|0
"
    .to_string()
}

fn write_bgzf_file(path: &std::path::Path, contents: &str) {
    let file = fs::File::create(path).expect("test fixture should be created");
    let mut writer = noodles_bgzf::io::Writer::new(file);
    writer
        .write_all(contents.as_bytes())
        .expect("test fixture should be compressed");
}

#[test]
fn phased_vcf_haplotype_dense_counts_a1_by_sample_haplotype_rows() {
    let dir = unique_dir("vcf-haplo-dense");
    let path = dir.join("phased.vcf");
    fs::write(&path, phased_vcf()).expect("fixture should be written");

    let haplotypes =
        genoio_io::read_vcf_haplotypes_dense(&path, None, None).expect("haplotypes should decode");

    assert_eq!(haplotypes.n_samples, 4);
    assert_eq!(haplotypes.n_variants, 2);
    assert_eq!(
        dense_values_sample_major(&haplotypes),
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        haplotypes
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S1", "S2", "S2"]
    );
}

#[test]
fn phased_vcf_haplotype_dense_matrix_only_omits_metadata() {
    let dir = unique_dir("vcf-haplo-dense-matrix-only");
    let path = dir.join("phased.vcf");
    fs::write(&path, phased_vcf()).expect("fixture should be written");

    let haplotypes = genoio_io::read_vcf_haplotypes_dense_windowed(&path, None, None, None, true)
        .expect("matrix-only haplotypes should decode");

    assert_eq!(haplotypes.n_samples, 4);
    assert_eq!(haplotypes.n_variants, 2);
    assert_eq!(
        dense_values_sample_major(&haplotypes),
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    );
    assert!(haplotypes.samples.is_none());
    assert!(haplotypes.variants.is_none());
}

#[test]
fn compressed_vcf_haplotype_dense_uses_text_backend_semantics() {
    let dir = unique_dir("vcf-haplo-dense-text-compressed");
    let path = dir.join("phased.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DP\t0|1:7\t1|0:8\t1|1:9
1\t20\trs2\tC\tT\t.\tPASS\t.\tDP:GT\t5:1|0\t6:0|1\t7:0|0
",
    );
    let samples = vec!["S3".to_string(), "S1".to_string()];

    let haplotypes = genoio_io::read_vcf_haplotypes_dense(&path, Some(&samples), None)
        .expect("compressed haplotypes should decode");

    assert_eq!(haplotypes.n_samples, 4);
    assert_eq!(haplotypes.n_variants, 2);
    assert_eq!(
        dense_values_sample_major(&haplotypes),
        vec![0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
    );
    assert_eq!(
        haplotypes
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| (sample.iid.as_str(), sample.haplotype_index))
            .collect::<Vec<_>>(),
        vec![
            ("S1", Some(0)),
            ("S1", Some(1)),
            ("S3", Some(0)),
            ("S3", Some(1)),
        ]
    );
    assert_eq!(
        variant_ids(variants(&haplotypes.variants)),
        vec!["rs1", "rs2"]
    );
}

#[test]
fn filtered_haplotype_samples_preserve_source_sample_index() {
    let dir = unique_dir("vcf-haplo-filtered-samples");
    let path = dir.join("phased.vcf");
    fs::write(&path, phased_vcf()).expect("fixture should be written");
    let samples = vec!["S2".to_string()];

    let haplotypes = genoio_io::read_vcf_haplotypes_dense(&path, Some(&samples), None)
        .expect("haplotypes should decode");

    assert_eq!(haplotypes.n_samples, 2);
    assert_eq!(
        haplotypes
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S2", "S2"]
    );
    assert_eq!(
        haplotypes
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.source_sample_index)
            .collect::<Vec<_>>(),
        vec![Some(1), Some(1)]
    );
    assert_eq!(
        haplotypes
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.haplotype_index)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(1)]
    );
}

#[test]
fn phased_vcf_haplotype_sparse_reconstructs_dense_values() {
    let dir = unique_dir("vcf-haplo-sparse");
    let path = dir.join("phased.vcf");
    fs::write(&path, phased_vcf()).expect("fixture should be written");

    let sparse =
        genoio_io::read_vcf_haplotypes_sparse(&path, None, None).expect("haplotypes should decode");

    assert_eq!(sparse.n_rows, 4);
    assert_eq!(sparse.n_cols, 2);
    assert_eq!(
        sparse_values_dense_arrow(&sparse),
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    );
}

#[test]
fn phased_vcf_haplotype_sparse_flips_common_alt_allele() {
    let dir = unique_dir("vcf-haplo-sparse-flip");
    let path = dir.join("phased.vcf");
    fs::write(&path, phased_alt_major_vcf()).expect("fixture should be written");

    let sparse =
        genoio_io::read_vcf_haplotypes_sparse(&path, None, None).expect("haplotypes should decode");

    assert_eq!(sparse.indptr, vec![0, 1]);
    assert_eq!(sparse.indices, vec![3]);
    assert_eq!(sparse.data, vec![1.0]);
    let variants = variants(&sparse.variants);
    assert!(variants.flipped[0]);
    assert_eq!(variant_a0(variants, 0), "G");
    assert_eq!(variant_a1(variants, 0), "A");
}

#[test]
fn compressed_vcf_haplotype_sparse_windowed_matches_existing_semantics() {
    let dir = unique_dir("vcf-haplo-sparse-text-compressed");
    let path = dir.join("phased.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DP\t1|1:7\t1|0:8\t0|1:9
1\t20\trs2\tC\tT\t.\tPASS\t.\tDP:GT\t5:1|0\t6:0|1\t7:0|0
1\t30\trs3\tC\tA\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0
",
    );
    let samples = vec!["S3".to_string(), "S1".to_string()];
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("MAF filter should parse");

    let sparse = genoio_io::read_vcf_haplotypes_sparse_windowed(
        &path,
        Some(&samples),
        Some(&filter),
        Some(genoio_core::VariantWindow { start: 0, len: 2 }),
    )
    .expect("compressed haplotype sparse VCF should decode");

    assert_eq!(sparse.n_rows, 4);
    assert_eq!(sparse.n_cols, 2);
    assert_eq!(sparse.indptr, vec![0, 1, 2]);
    assert_eq!(sparse.indices, vec![2, 0]);
    assert_eq!(sparse.data, vec![1.0, 1.0]);
    assert_eq!(
        sparse
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| (sample.iid.as_str(), sample.haplotype_index))
            .collect::<Vec<_>>(),
        vec![
            ("S1", Some(0)),
            ("S1", Some(1)),
            ("S3", Some(0)),
            ("S3", Some(1)),
        ]
    );
    assert_eq!(variant_ids(variants(&sparse.variants)), vec!["rs1", "rs2"]);
    assert_eq!(variants(&sparse.variants).flipped, vec![true, false]);
}

#[test]
fn threaded_haplotype_reads_match_unthreaded_reads() {
    let dir = unique_dir("vcf-haplo-threaded");
    let path = dir.join("phased.vcf.gz");
    let file = fs::File::create(&path).expect("fixture should be created");
    let mut writer = noodles_bgzf::io::Writer::new(file);
    writer
        .write_all(phased_vcf().as_bytes())
        .expect("fixture should be compressed");
    drop(writer);

    let dense = genoio_io::read_vcf_haplotypes_dense_windowed(&path, None, None, None, false)
        .expect("unthreaded dense haplotypes should decode");
    let threaded_dense = genoio_io::read_vcf_haplotypes_dense_windowed_with_threads(
        &path,
        None,
        None,
        None,
        false,
        Some(2),
    )
    .expect("threaded dense haplotypes should decode");
    assert_eq!(threaded_dense.values, dense.values);
    assert_eq!(threaded_dense.n_samples, dense.n_samples);
    assert_eq!(threaded_dense.n_variants, dense.n_variants);

    let sparse = genoio_io::read_vcf_haplotypes_sparse_windowed(&path, None, None, None)
        .expect("unthreaded sparse haplotypes should decode");
    let threaded_sparse = genoio_io::read_vcf_haplotypes_sparse_windowed_with_threads(
        &path,
        None,
        None,
        None,
        Some(2),
    )
    .expect("threaded sparse haplotypes should decode");
    assert_eq!(threaded_sparse.data, sparse.data);
    assert_eq!(threaded_sparse.indices, sparse.indices);
    assert_eq!(threaded_sparse.indptr, sparse.indptr);
    assert_eq!(threaded_sparse.n_rows, sparse.n_rows);
    assert_eq!(threaded_sparse.n_cols, sparse.n_cols);
}

#[test]
fn haplotype_decode_rejects_unphased_retained_genotype() {
    let dir = unique_dir("vcf-haplo-unphased");
    let path = dir.join("mixed.vcf");
    fs::write(&path, mixed_phase_vcf()).expect("fixture should be written");

    let error = genoio_io::read_vcf_haplotypes_dense(&path, None, None)
        .expect_err("unphased retained genotype should fail");

    assert!(error.to_string().contains("unphased"));
}

#[test]
fn haplotype_stat_filter_drops_unphased_genotype_before_separator_check() {
    let dir = unique_dir("vcf-haplo-filter-unphased");
    let path = dir.join("mixed-stat.vcf");
    fs::write(&path, mixed_phase_stat_filter_vcf()).expect("fixture should be written");
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("filter should parse");

    let haplotypes = genoio_io::read_vcf_haplotypes_dense(&path, None, Some(&filter))
        .expect("unphased dropped genotype should not fail haplotype decode");

    assert_eq!(
        variant_ids(variants(&haplotypes.variants)),
        vec!["rs_phased"]
    );
    assert_eq!(
        dense_values_sample_major(&haplotypes),
        vec![0.0, 1.0, 1.0, 0.0]
    );
}

#[test]
fn haplotype_matrix_only_stat_filter_drops_unphased_genotype_before_separator_check() {
    let dir = unique_dir("vcf-haplo-filter-unphased-matrix-only");
    let path = dir.join("mixed-stat.vcf");
    fs::write(&path, mixed_phase_stat_filter_vcf()).expect("fixture should be written");
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("filter should parse");

    let haplotypes =
        genoio_io::read_vcf_haplotypes_dense_windowed(&path, None, Some(&filter), None, true)
            .expect("matrix-only unphased dropped genotype should not fail haplotype decode");

    assert_eq!(haplotypes.n_samples, 4);
    assert_eq!(haplotypes.n_variants, 1);
    assert!(haplotypes.samples.is_none());
    assert!(haplotypes.variants.is_none());
    assert_eq!(
        dense_values_sample_major(&haplotypes),
        vec![0.0, 1.0, 1.0, 0.0]
    );
}
