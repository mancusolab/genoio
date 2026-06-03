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

fn write_text(path: &Path, contents: &str) {
    fs::write(path, contents).expect("test fixture should be written");
}

fn write_plink2_fixture(dir: &Path, pgen_bytes: &[u8]) -> (PathBuf, PathBuf, PathBuf) {
    let pgen = dir.join("tiny.pgen");
    let pvar = dir.join("tiny.pvar");
    let psam = dir.join("tiny.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
2 30 rs3 G A 50
",
    );
    write_text(
        &psam,
        "\
#FID IID PAT MAT SEX PHENO
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
",
    );
    (pgen, pvar, psam)
}

fn fixed_width_pgen(records: &[u8], n_samples: u32, n_variants: u32) -> Vec<u8> {
    let mut bytes = vec![0x6c, 0x1b, 0x02];
    bytes.extend(n_variants.to_le_bytes());
    bytes.extend(n_samples.to_le_bytes());
    bytes.push(0);
    bytes.extend(records);
    bytes
}

fn variable_width_pgen(record_types: &[u8], records: &[&[u8]], n_samples: u32) -> Vec<u8> {
    let n_variants = u32::try_from(records.len()).expect("test variant count fits u32");
    let header_len = 12 + 8 + record_types.len() + records.len();
    let mut bytes = vec![0x6c, 0x1b, 0x10];
    bytes.extend(n_variants.to_le_bytes());
    bytes.extend(n_samples.to_le_bytes());
    bytes.push(0x04);
    bytes.extend(
        u64::try_from(header_len)
            .expect("test header length fits u64")
            .to_le_bytes(),
    );
    bytes.extend(record_types);
    bytes.extend(
        records
            .iter()
            .map(|record| u8::try_from(record.len()).expect("test record length fits one byte")),
    );
    for record in records {
        bytes.extend(*record);
    }
    bytes
}

#[test]
fn plink2_dense_decodes_fixed_width_unphased_biallelic_hardcalls() {
    let dir = unique_dir("plink2-dense-values");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    let dense =
        genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None).expect("pgen should decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 3);
    assert_eq!(
        dense.values,
        vec![0.0, 1.0, 2.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0]
    );
    assert_eq!(
        dense.missing_mask,
        vec![false, false, false, true, false, false, false, false, false]
    );
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S2", "S3"]
    );
    assert_eq!(dense.variants[0].ref_allele.as_deref(), Some("A"));
    assert_eq!(dense.variants[0].alt_allele.as_deref(), Some("G"));
    assert_eq!(dense.variants[0].qual, Some(30.0));
}

#[test]
fn plink2_dense_filters_samples_in_source_order() {
    let dir = unique_dir("plink2-dense-samples");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);
    let keep = vec!["S3".to_string(), "S1".to_string()];

    let dense = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, Some(&keep), None)
        .expect("pgen should filter samples");

    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_eq!(dense.values, vec![0.0, 1.0, 2.0, 2.0, 1.0, 0.0]);
}

#[test]
fn plink2_dense_rejects_unsupported_pgen_modes() {
    let dir = unique_dir("plink2-dense-unsupported-mode");
    let mut pgen_bytes = fixed_width_pgen(&[0x00], 1, 1);
    pgen_bytes[2] = 0x03;
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    let error = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect_err("unsupported pgen mode should fail");

    assert!(error.to_string().contains("unsupported pgen mode"));
}

#[test]
fn plink2_dense_decodes_variable_width_hardcall_records() {
    let dir = unique_dir("plink2-dense-variable-width");
    let pgen_bytes = variable_width_pgen(
        &[0, 4, 1, 2, 3],
        &[&[0xe4], &[2, 1, 9, 2], &[2, 5, 0], &[1, 1, 3], &[0]],
        4,
    );
    let pgen = dir.join("variable.pgen");
    let pvar = dir.join("variable.pvar");
    let psam = dir.join("variable.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 A G
1 30 rs3 A G
1 40 rs4 A G
1 50 rs5 A G
",
    );
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
S4
",
    );

    let dense = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect("variable-width pgen should decode");

    assert_eq!(
        dense.values,
        vec![
            0.0, 0.0, 2.0, 2.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0,
            0.0, 0.0, 2.0,
        ]
    );
    assert_eq!(
        dense.missing_mask,
        vec![
            false, false, false, false, false, false, false, false, true, false, false, false,
            false, false, false, true, false, false, false, false,
        ]
    );
}

#[test]
fn plink2_dense_rejects_non_increasing_difflist_sample_ids() {
    let dir = unique_dir("plink2-dense-bad-difflist");
    let pgen_bytes = variable_width_pgen(&[4], &[&[2, 1, 5, 0]], 4);
    let pgen = dir.join("bad_difflist.pgen");
    let pvar = dir.join("bad_difflist.pvar");
    let psam = dir.join("bad_difflist.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
S4
",
    );

    let error = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect_err("duplicate difflist sample id should fail");

    assert!(error.to_string().contains("strictly increasing"));
}

#[test]
fn plink2_dense_rejects_variable_width_block_offset_mismatch() {
    let dir = unique_dir("plink2-dense-bad-block-offset");
    let mut pgen_bytes = variable_width_pgen(&[0], &[&[0x00]], 4);
    pgen_bytes[12] = pgen_bytes[12].saturating_add(1);
    let pgen = dir.join("bad_offset.pgen");
    let pvar = dir.join("bad_offset.pvar");
    let psam = dir.join("bad_offset.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
S4
",
    );

    let error = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect_err("bad block offset should fail");

    assert!(error.to_string().contains("block offset"));
}

#[test]
fn plink2_metadata_reads_psam_and_pvar() {
    let dir = unique_dir("plink2-metadata");
    let pgen_bytes = fixed_width_pgen(&[0x00, 0x00, 0x00], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    let metadata =
        genoio_io::read_plink2_metadata(&pgen, &pvar, &psam).expect("metadata should decode");

    assert_eq!(metadata.samples.len(), 3);
    assert_eq!(metadata.variants.len(), 3);
    assert!(metadata.capabilities.supports_geno);
    assert!(!metadata.capabilities.supports_haplo);
}
