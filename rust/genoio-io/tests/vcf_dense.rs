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

fn write_file(path: &Path, contents: &str) {
    fs::write(path, contents).expect("test fixture should be written");
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
