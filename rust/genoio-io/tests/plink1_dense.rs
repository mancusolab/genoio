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
2 indel1 0 30 A AT
",
    );
    write_text(
        &fam,
        "\
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
",
    );
    (bed, bim, fam)
}

#[test]
fn plink1_dense_decodes_variant_major_bed_to_sample_by_variant_matrix() {
    let dir = unique_dir("plink1-dense-values");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x07, 0x2d, 0x38]);

    let dense =
        genoio_io::read_plink1_dense(&bed, &bim, &fam, None, None).expect("plink1 should decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 3);
    assert_eq!(
        dense.values,
        vec![0.0, 1.0, 2.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0]
    );
    assert_eq!(
        dense.missing_mask,
        vec![false, false, false, false, false, true, false, true, false]
    );
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S2", "S3"]
    );
}

#[test]
fn plink1_dense_rejects_invalid_magic_bytes() {
    let dir = unique_dir("plink1-dense-bad-magic");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x00, 0x1b, 0x01, 0x00, 0x00, 0x00]);

    let error = genoio_io::read_plink1_dense(&bed, &bim, &fam, None, None)
        .expect_err("bad magic should fail");

    assert!(error.to_string().contains("magic"));
}

#[test]
fn plink1_dense_rejects_sample_major_mode() {
    let dir = unique_dir("plink1-dense-sample-major");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x00, 0x00, 0x00, 0x00]);

    let error = genoio_io::read_plink1_dense(&bed, &bim, &fam, None, None)
        .expect_err("sample-major should fail");

    assert!(error.to_string().contains("sample-major"));
}

#[test]
fn plink1_dense_filters_samples_in_source_order() {
    let dir = unique_dir("plink1-dense-sample-filter");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x07, 0x2d, 0x38]);

    let keep = vec!["S3".to_string(), "S1".to_string()];
    let dense = genoio_io::read_plink1_dense(&bed, &bim, &fam, Some(&keep), None)
        .expect("plink1 should filter");

    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_eq!(dense.values, vec![0.0, 1.0, 2.0, 2.0, 0.0, 0.0]);
    assert_eq!(dense.diagnostics.requested_samples, 2);
    assert_eq!(dense.diagnostics.retained_samples, 2);
    assert_eq!(dense.diagnostics.missing_samples, 0);
}
