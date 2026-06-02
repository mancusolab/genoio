// pattern: Imperative Shell

use std::fs;
use std::path::PathBuf;

fn unique_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "genoio-{name}-{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("test temp dir should be created");
    dir
}

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
        haplotypes.values,
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(
        haplotypes
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S1", "S2", "S2"]
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
        csc_to_dense(&sparse),
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    );
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
