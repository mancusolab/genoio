// pattern: Imperative Shell

use std::fs;
use std::io::Write;
use std::path::Path;

mod common;

use common::dense::assert_values_with_nan;
use common::unique_dir;
use common::vcf_output as genoio_io;
use common::vcf_output::{
    dense_missing_sample_major_output as dense_missing_sample_major,
    dense_values_sample_major_output as dense_values_sample_major, variant_ids, variants,
};

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

fn contradictory_chrom_filter() -> genoio_core::VariantFilter {
    genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "and",
        "left": {"op": "predicate", "name": "chrom", "params": {"value": "1"}},
        "right": {"op": "predicate", "name": "chrom", "params": {"value": "2"}}
    }))
    .expect("contradictory chrom filter should parse")
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
    assert_eq!(
        dense_values_sample_major(&dense),
        vec![0.0, 1.0, 1.0, 0.0, 2.0, 2.0]
    );
    assert_eq!(dense_missing_sample_major(&dense), vec![false; 6]);
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S2", "S3"]
    );
    assert_eq!(variant_ids(variants(&dense.variants)), vec!["rs1", "rs2"]);
}

#[test]
fn plain_vcf_dense_uses_permissive_text_header_path() {
    let dir = unique_dir("vcf-dense-text-plain");
    let path = dir.join("tiny.vcf");
    write_file(
        &path,
        "\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DP\t0/0:7\t0/1:8\t1/1:9
1\t20\trs2\tC\tT\t.\tPASS\t.\tDP:GT\t5:0/1\t6:0/0\t7:./.
",
    );
    let samples = vec!["S3".to_string(), "S1".to_string()];

    let dense = genoio_io::read_vcf_dense(&path, Some(&samples), None)
        .expect("plain VCF should decode through the permissive text backend");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_values_with_nan(
        &dense_values_sample_major(&dense),
        &[0.0, 1.0, 2.0, f32::NAN],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, true]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
}

#[test]
fn plain_vcf_always_false_filter_uses_text_empty_path() {
    let dir = unique_dir("vcf-dense-text-plain-empty");
    let path = dir.join("tiny.vcf");
    write_file(
        &path,
        "\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
this record is intentionally not valid VCF and should not be decoded
",
    );
    let filter = contradictory_chrom_filter();

    let dense = genoio_io::read_vcf_dense(&path, None, Some(&filter))
        .expect("always-false VCF filter should only need samples");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 0);
    assert!(dense.values.is_empty());
    assert!(dense_missing_sample_major(&dense).is_empty());
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S2"]
    );
    assert!(variants(&dense.variants).is_empty());
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
    assert_eq!(
        dense_values_sample_major(&dense),
        vec![0.0, 1.0, 1.0, 0.0, 2.0, 2.0]
    );
    assert_eq!(dense_missing_sample_major(&dense), vec![false; 6]);
    assert!(dense.samples.is_none());
    assert!(dense.variants.is_none());
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
    assert_values_with_nan(
        &dense_values_sample_major(&dense),
        &[0.2, 0.0, 1.4, f32::NAN, 1.8, 0.7],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, true, false, false]
    );
    assert!(dense.samples.is_none());
    assert!(dense.variants.is_none());
}

#[test]
fn compressed_vcf_dosage_dense_uses_text_backend_semantics() {
    let dir = unique_dir("vcf-dosage-text-compressed");
    let path = dir.join("dosage.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected alternate allele dosage\"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tDS\t0.2\t1.4\t1.8
1\t20\trs2\tC\tT\t.\tPASS\t.\tDS\t0\t.\t0.7
",
    );
    let samples = vec!["S3".to_string(), "S1".to_string()];

    let dense = genoio_io::read_vcf_dosage_dense_windowed(&path, Some(&samples), None, None, false)
        .expect("compressed dosage VCF should decode");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(dense_values_sample_major(&dense), vec![0.2, 0.0, 1.8, 0.7]);
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, false]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_eq!(variant_ids(variants(&dense.variants)), vec!["rs1", "rs2"]);
}

#[test]
fn threaded_compressed_vcf_dosage_uses_text_backend_semantics() {
    let dir = unique_dir("vcf-dosage-threaded-text-compressed");
    let path = dir.join("dosage.vcf.gz");
    // The malformed FORMAT header is intentional: it proves threaded compressed
    // reads still use the permissive noodles text backend.
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected alternate allele dosage\"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tDS\t0.2\t1.4\t1.8
1\t20\trs2\tC\tT\t.\tPASS\t.\tDS\t0\t.\t0.7
",
    );
    let samples = vec!["S3".to_string(), "S1".to_string()];

    let dense = genoio_io::read_vcf_dosage_dense_windowed_with_threads(
        &path,
        Some(&samples),
        None,
        None,
        false,
        Some(2),
    )
    .expect("threaded compressed dosage VCF should decode through noodles text backend");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(dense_values_sample_major(&dense), vec![0.2, 0.0, 1.8, 0.7]);
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, false]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
}

#[test]
fn threaded_compressed_vcf_rejects_zero_threads() {
    let dir = unique_dir("vcf-threaded-zero");
    let path = dir.join("tiny.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1
",
    );

    let error =
        genoio_io::read_vcf_dense_windowed_with_threads(&path, None, None, None, false, Some(0))
            .expect_err("zero VCF thread count should fail");

    assert!(error
        .to_string()
        .contains("vcf thread count must be greater than zero"));
}

#[test]
fn plain_vcf_rejects_threaded_reads() {
    let dir = unique_dir("vcf-threaded-plain");
    let path = dir.join("tiny.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1
",
    );

    let error =
        genoio_io::read_vcf_dense_windowed_with_threads(&path, None, None, None, false, Some(2))
            .expect_err("plain VCF should reject explicit thread count");

    assert!(error
        .to_string()
        .contains("threaded reads are only supported for compressed VCF"));
}

#[test]
fn compressed_vcf_dosage_rejects_invalid_ds_values() {
    let cases = [
        ("GT\t0/0", "FORMAT/DS"),
        ("DS\t0.1,0.2", "multiple values"),
        ("DS\tNaN", "finite value in [0, 2]"),
        ("DS\tinf", "finite value in [0, 2]"),
        ("DS\t2.1", "finite value in [0, 2]"),
    ];

    for (format_and_sample, expected) in cases {
        let dir = unique_dir("vcf-dosage-text-invalid");
        let path = dir.join("dosage.vcf.gz");
        write_bgzf_file(
            &path,
            &format!(
                "\
##fileformat=VCFv4.2
##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected alternate allele dosage\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t10\trs1\tA\tG\t.\tPASS\t.\t{format_and_sample}
"
            ),
        );

        let error = genoio_io::read_vcf_dosage_dense_windowed(&path, None, None, None, false)
            .expect_err("invalid compressed DS should fail");

        assert!(
            error.to_string().contains(expected),
            "expected error containing {expected:?}, got {error}"
        );
    }
}

#[test]
fn compressed_vcf_matrix_only_uses_text_backend_semantics() {
    let dir = unique_dir("vcf-dense-text-compressed");
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
    assert_values_with_nan(
        &dense_values_sample_major(&dense),
        &[0.0, 1.0, 2.0, f32::NAN],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, true]
    );
    assert!(dense.samples.is_none());
    assert!(dense.variants.is_none());
}

#[test]
fn compressed_vcf_dense_with_metadata_uses_text_backend_semantics() {
    let dir = unique_dir("vcf-dense-text-metadata");
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
",
    );
    let samples = vec!["S3".to_string(), "S1".to_string()];

    let dense =
        genoio_io::read_vcf_dense(&path, Some(&samples), None).expect("dense VCF should decode");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_values_with_nan(
        &dense_values_sample_major(&dense),
        &[0.0, 1.0, 2.0, f32::NAN],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, true]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_eq!(variant_ids(variants(&dense.variants)), vec!["rs1", "rs2"]);
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
    assert_eq!(dense_values_sample_major(&dense), vec![0.0, 1.0, 2.0, 0.0]);
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.source_sample_index)
            .collect::<Vec<_>>(),
        vec![None, None]
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

    assert_values_with_nan(&dense.values, &[1.0, f32::NAN]);
    assert_eq!(dense_missing_sample_major(&dense), vec![false, true]);
}

#[test]
fn vcf_dense_contract_validates_shape_and_metadata_lengths() {
    let result = genoio_core::DenseGenotypeMatrix::new_with_layout(
        1,
        2,
        vec![0.0],
        genoio_core::DenseLayout::SampleMajor,
        None,
        None,
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
