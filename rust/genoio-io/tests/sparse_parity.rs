use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

mod common;

use common::unique_dir;

fn write_text(path: &Path, contents: &str) {
    fs::write(path, contents).expect("test fixture should be written");
}

fn write_bgzf_file(path: &Path, contents: &str) {
    let file = fs::File::create(path).expect("test fixture should be created");
    let mut writer = noodles_bgzf::io::Writer::new(file);
    writer
        .write_all(contents.as_bytes())
        .expect("test fixture should be compressed");
}

fn write_vcf(path: &Path, body: &str) {
    write_text(
        path,
        &format!(
            "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
{body}"
        ),
    );
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

fn write_plink_fixture(dir: &Path, bed_bytes: &[u8]) -> (PathBuf, PathBuf, PathBuf) {
    let bed = dir.join("tiny.bed");
    let bim = dir.join("tiny.bim");
    let fam = dir.join("tiny.fam");
    fs::write(&bed, bed_bytes).expect("bed fixture should be written");
    write_text(
        &bim,
        "\
1 rs1 0 10 G A
1 rs2 0 20 T C
",
    );
    write_text(
        &fam,
        "\
F1 S1 0 0 1 -9
F1 S2 0 0 2 -9
F1 S3 0 0 0 -9
",
    );
    (bed, bim, fam)
}

#[test]
fn vcf_sparse_reconstructs_dense_when_no_missing_calls() {
    let dir = unique_dir("vcf-sparse-parity");
    let path = dir.join("tiny.vcf");
    write_vcf(
        &path,
        "\
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1
",
    );

    let dense = genoio_io::read_vcf_dense(&path, None, None).expect("dense vcf should decode");
    let sparse = genoio_io::read_vcf_sparse(&path, None, None).expect("sparse vcf should decode");

    assert_eq!(sparse.n_rows, dense.n_samples);
    assert_eq!(sparse.n_cols, dense.n_variants);
    assert_eq!(csc_to_dense(&sparse), dense.values);
}

#[test]
fn compressed_vcf_sparse_windowed_matches_existing_sparse_semantics() {
    let dir = unique_dir("vcf-sparse-fast-compressed");
    let path = dir.join("tiny.vcf.gz");
    write_bgzf_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT:DP\t0/0:7\t0/1:8\t1/1:9
1\t20\trs2\tC\tT\t.\tPASS\t.\tDP:GT\t5:0/1\t6:0/0\t7:1/1
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

    let sparse = genoio_io::read_vcf_sparse_windowed(
        &path,
        Some(&samples),
        Some(&filter),
        Some(genoio_core::VariantWindow { start: 0, len: 2 }),
    )
    .expect("compressed sparse VCF should decode");

    assert_eq!(sparse.n_rows, 2);
    assert_eq!(sparse.n_cols, 2);
    assert_eq!(sparse.indptr, vec![0, 1, 2]);
    assert_eq!(sparse.indices, vec![1, 0]);
    assert_eq!(sparse.data, vec![2.0, 1.0]);
    assert_eq!(
        sparse
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_eq!(
        sparse
            .variants
            .iter()
            .map(|variant| (variant.id.as_str(), variant.flipped))
            .collect::<Vec<_>>(),
        vec![("rs1", false), ("rs2", true)]
    );
}

#[test]
fn plink1_sparse_reconstructs_dense_when_no_missing_calls() {
    let dir = unique_dir("plink1-sparse-parity");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x0b, 0x2c]);

    let dense = genoio_io::read_plink1_dense(&bed, &bim, &fam, None, None)
        .expect("dense plink should decode");
    let sparse = genoio_io::read_plink1_sparse(&bed, &bim, &fam, None, None)
        .expect("sparse plink should decode");

    assert_eq!(csc_to_dense(&sparse), dense.values);
}

#[test]
fn sparse_reads_fail_when_retained_missing_calls_are_present() {
    let dir = unique_dir("vcf-sparse-missing");
    let path = dir.join("missing.vcf");
    write_vcf(
        &path,
        "\
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t./.\t1/1
",
    );

    let error = genoio_io::read_vcf_sparse(&path, None, None)
        .expect_err("missing sparse genotypes should fail");

    assert!(error.to_string().contains("sparse missing values"));
}

#[test]
fn sparse_reads_flip_common_minor_allele_columns_by_default() {
    let dir = unique_dir("vcf-sparse-flip");
    let path = dir.join("flip.vcf");
    write_vcf(
        &path,
        "\
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t1/1\t1/1\t0/1
",
    );

    let dense = genoio_io::read_vcf_dense(&path, None, None).expect("dense vcf should decode");
    let sparse = genoio_io::read_vcf_sparse(&path, None, None).expect("sparse vcf should decode");

    assert_eq!(dense.values, vec![2.0, 2.0, 1.0]);
    assert_eq!(csc_to_dense(&sparse), vec![0.0, 0.0, 1.0]);
    assert!(sparse.variants[0].flipped);
    assert_eq!(sparse.variants[0].a0, "G");
    assert_eq!(sparse.variants[0].a1, "A");
    assert!(!dense.variants[0].flipped);
}
