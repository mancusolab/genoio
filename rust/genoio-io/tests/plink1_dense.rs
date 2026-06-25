// pattern: Imperative Shell

use std::fs;
use std::path::{Path, PathBuf};

mod common;

use common::dense::assert_values_with_nan;
use common::plink_arrow as genoio_io;
use common::plink_arrow::dense_missing_sample_major_arrow as dense_missing_sample_major;
use common::unique_dir;

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

fn write_plink_code_table_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let bed = dir.join("codes.bed");
    let bim = dir.join("codes.bim");
    let fam = dir.join("codes.fam");
    // Low-order two-bit sample codes: 00, 01, 10, 11.
    fs::write(&bed, [0x6c, 0x1b, 0x01, 0xe4]).expect("bed fixture should be written");
    write_text(&bim, "1 rs_codes 0 10 G A\n");
    write_text(
        &fam,
        "\
F1 S1 0 0 1 -9
F1 S2 0 0 1 -9
F1 S3 0 0 1 -9
F1 S4 0 0 1 -9
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
    assert_values_with_nan(
        &dense.values,
        &[0.0, f32::NAN, 2.0, f32::NAN, 0.0, 1.0, 2.0, 1.0, 0.0],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, true, false, true, false, false, false, false, false]
    );
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
}

#[test]
fn plink1_dense_matrix_only_omits_metadata() {
    let dir = unique_dir("plink1-dense-matrix-only");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x07, 0x2d, 0x38]);

    let dense = genoio_io::read_plink1_dense_windowed(&bed, &bim, &fam, None, None, None, true)
        .expect("plink1 matrix-only read should decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 3);
    assert!(dense.samples.is_none());
    assert!(dense.variants.is_none());
    assert_values_with_nan(
        &dense.values,
        &[0.0, f32::NAN, 2.0, f32::NAN, 0.0, 1.0, 2.0, 1.0, 0.0],
    );
}

#[test]
fn plink1_dense_decodes_official_two_bit_code_table() {
    let dir = unique_dir("plink1-dense-code-table");
    let (bed, bim, fam) = write_plink_code_table_fixture(&dir);

    let dense =
        genoio_io::read_plink1_dense(&bed, &bim, &fam, None, None).expect("plink1 should decode");

    assert_values_with_nan(&dense.values, &[2.0, f32::NAN, 1.0, 0.0]);
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, true, false, false]
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
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_values_with_nan(&dense.values, &[0.0, f32::NAN, 2.0, 2.0, 1.0, 0.0]);
    assert_eq!(dense.diagnostics.requested_samples, 2);
    assert_eq!(dense.diagnostics.retained_samples, 2);
    assert_eq!(dense.diagnostics.missing_samples, 0);
}
