use std::fs;
use std::io::Write;
use std::path::Path;

mod common;

use common::unique_dir;

fn write_file(path: &Path, contents: &str) {
    fs::write(path, contents).expect("test fixture should be written");
}

fn write_bgzf_file(path: &Path, contents: &str) {
    let file = fs::File::create(path).expect("test fixture should be created");
    let mut writer = noodles_bgzf::io::Writer::new(file);
    writer
        .write_all(contents.as_bytes())
        .expect("test fixture should be compressed");
}

#[test]
fn vcf_dense_values_count_a1_in_sample_by_variant_shape() {
    let dir = unique_dir("vcf-dense-values");
    let path = dir.join("tiny.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1\t0/0\t1|1
",
    );

    let dense = genoio_io::read_vcf_dense(&path, None, None).expect("dense vcf should decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(dense.values, vec![0.0, 1.0, 1.0, 0.0, 2.0, 2.0]);
    assert_eq!(dense.missing_mask, vec![false; 6]);
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S2", "S3"]
    );
    assert_eq!(
        dense
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs1", "rs2"]
    );
}

#[test]
fn vcf_dense_matrix_only_omits_metadata() {
    let dir = unique_dir("vcf-dense-matrix-only");
    let path = dir.join("tiny.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1\t0/0\t1|1
",
    );

    let dense = genoio_io::read_vcf_dense_windowed(&path, None, None, None, true)
        .expect("matrix-only dense vcf should decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(dense.values, vec![0.0, 1.0, 1.0, 0.0, 2.0, 2.0]);
    assert_eq!(dense.missing_mask, vec![false; 6]);
    assert!(dense.samples.is_empty());
    assert!(dense.variants.is_empty());
}

#[test]
fn vcf_dosage_dense_matrix_only_omits_metadata() {
    let dir = unique_dir("vcf-dosage-dense-matrix-only");
    let path = dir.join("dosage.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected alternate allele dosage\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DS\t0/0:0.2\t0/1:1.4\t1/1:1.8
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT:DS\t0/0:0\t0/0:.\t0/1:0.7
",
    );

    let dense = genoio_io::read_vcf_dosage_dense_windowed(&path, None, None, None, true)
        .expect("matrix-only vcf dosage should decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(dense.values, vec![0.2, 0.0, 1.4, 0.0, 1.8, 0.7]);
    assert_eq!(
        dense.missing_mask,
        vec![false, false, false, true, false, false]
    );
    assert!(dense.samples.is_empty());
    assert!(dense.variants.is_empty());
}

#[test]
fn compressed_vcf_matrix_only_uses_fast_path_semantics() {
    let dir = unique_dir("vcf-dense-fast-compressed");
    let path = dir.join("tiny.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DP\t0/0:7\t0/1:8\t1/1:9
1\t20\trs2\tC\tT\t.\tPASS\t.\tDP:GT\t5:0/1\t6:0/0\t7:./.
1\t30\trs3\tC\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0
",
    );
    let samples = vec!["S3".to_string(), "S1".to_string()];
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("MAF filter should parse");

    let dense = genoio_io::read_vcf_dense_windowed(
        &path,
        Some(&samples),
        Some(&filter),
        Some(genoio_core::VariantWindow { start: 0, len: 2 }),
        true,
    )
    .expect("compressed matrix-only VCF should decode");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(dense.values, vec![0.0, 1.0, 2.0, 0.0]);
    assert_eq!(dense.missing_mask, vec![false, false, false, true]);
    assert!(dense.samples.is_empty());
    assert!(dense.variants.is_empty());
}

#[test]
fn vcf_dense_sample_subset_preserves_source_order() {
    let dir = unique_dir("vcf-dense-sample-subset");
    let path = dir.join("subset.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0
",
    );
    let samples = vec!["S3".to_string(), "S1".to_string()];

    let dense =
        genoio_io::read_vcf_dense(&path, Some(&samples), None).expect("dense vcf should decode");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(dense.values, vec![0.0, 1.0, 2.0, 0.0]);
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
            .samples
            .iter()
            .map(|sample| sample.source_sample_index)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(2)]
    );
}

#[test]
fn vcf_dense_marks_missing_gt_calls() {
    let dir = unique_dir("vcf-dense-missing");
    let path = dir.join("missing.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t./.
",
    );

    let dense = genoio_io::read_vcf_dense(&path, None, None).expect("dense vcf should decode");

    assert_eq!(dense.values, vec![1.0, 0.0]);
    assert_eq!(dense.missing_mask, vec![false, true]);
}

#[test]
fn vcf_dense_contract_validates_shape_and_metadata_lengths() {
    let sample = genoio_core::SampleRecord {
        fid: None,
        iid: "S1".to_string(),
        father: None,
        mother: None,
        sex: None,
        phenotype: None,
        source_sample_index: None,
        haplotype_index: None,
    };
    let variant = genoio_core::VariantRecord {
        chrom: "1".to_string(),
        pos: 10,
        id: "rs1".to_string(),
        a0: "A".to_string(),
        a1: "G".to_string(),
        ref_allele: Some("A".to_string()),
        alt_allele: Some("G".to_string()),
        source_a0: "A".to_string(),
        source_a1: "G".to_string(),
        flipped: false,
        qual: None,
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    };

    let result = genoio_core::DenseGenotypeMatrix::new(
        1,
        2,
        vec![0.0],
        vec![false],
        vec![sample],
        vec![variant],
        genoio_core::DenseDiagnostics::default(),
    );

    assert!(result.is_err());
}

#[test]
fn vcf_dense_rejects_multiallelic_genotype_states() {
    let dir = unique_dir("vcf-dense-multiallelic");
    let path = dir.join("multi.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG,T\t.\tPASS\t.\tGT\t1/2
",
    );

    let error =
        genoio_io::read_vcf_dense(&path, None, None).expect_err("multiallelic GT should fail");

    assert!(error.to_string().contains("multiallelic"));
}

#[test]
fn vcf_dense_rejects_non_diploid_gt_calls() {
    let dir = unique_dir("vcf-dense-haploid");
    let path = dir.join("haploid.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t1
",
    );

    let error = genoio_io::read_vcf_dense(&path, None, None).expect_err("haploid GT should fail");

    assert!(error.to_string().contains("non-diploid GT"));
}

#[test]
fn vcf_dense_rejects_multi_alt_records_even_when_gt_uses_first_alt() {
    let dir = unique_dir("vcf-dense-multi-alt-record");
    let path = dir.join("multi-alt.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG,T\t.\tPASS\t.\tGT\t0/1
",
    );

    let error =
        genoio_io::read_vcf_dense(&path, None, None).expect_err("multi-ALT records should fail");

    assert!(error.to_string().contains("multi-ALT"));
}
