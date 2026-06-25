// pattern: Imperative Shell

use std::fs;
use std::path::{Path, PathBuf};

use genoio_core::{SparseGenotypeMatrixArrowVariants, VariantWindow};

mod common;

use common::dense::assert_values_with_nan;
use common::plink_arrow as genoio_io;
use common::plink_arrow::{
    dense_missing_sample_major_arrow as dense_missing_sample_major,
    dense_values_sample_major_arrow as dense_values_sample_major, variant_alt_allele, variant_id,
    variant_ids, variant_ref_allele, variants,
};
use common::{unique_dir, TestDir};

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

fn write_plink2_fixture_with_variants(
    dir: &Path,
    pgen_bytes: &[u8],
    variants: &str,
) -> (PathBuf, PathBuf, PathBuf) {
    let pgen = dir.join("tiny.pgen");
    let pvar = dir.join("tiny.pvar");
    let psam = dir.join("tiny.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(&pvar, variants);
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

fn explicit_phased_hardcall_pgen(records: &[&[u8]], n_samples: u32) -> Vec<u8> {
    // PGEN variable-width vrtype 0x10: uncompressed biallelic main track
    // followed by auxiliary track #2 hardcall phase bits.
    variable_width_pgen(&vec![0x10; records.len()], records, n_samples)
}

fn explicit_phased_dosage_pgen(records: &[&[u8]], n_samples: u32) -> Vec<u8> {
    // PGEN variable-width vrtype 0xc0: uncompressed biallelic main track,
    // full dosage track #4, and full explicit phased-dosage track #8.
    variable_width_pgen(&vec![0xc0; records.len()], records, n_samples)
}

fn fixed_width_phased_dosage_pgen(records: &[&[u8]], n_samples: u32) -> Vec<u8> {
    // PGEN fixed-width mode 0x04 stores each variant as an effective vrtype
    // 0xc0 record: hardcall main track, full dosage track, and full explicit
    // phased-dosage track, without a per-variant vrtype byte.
    let n_variants = u32::try_from(records.len()).expect("test variant count fits u32");
    let mut bytes = vec![0x6c, 0x1b, 0x04];
    bytes.extend(n_variants.to_le_bytes());
    bytes.extend(n_samples.to_le_bytes());
    bytes.push(0);
    for record in records {
        bytes.extend(*record);
    }
    bytes
}

fn explicit_phased_hardcall_ld_pgen() -> Vec<u8> {
    // rs1 is a non-LD explicit-phased hardcall. rs2 is vrtype 0x12:
    // LD-compressed main track plus hardcall phase bits, so windowed reads
    // must decode rs1 before retaining rs2.
    let record_1 = [0x21, 0x00];
    let record_2 = [0x02, 0x01, 0x0d, 0x01, 0x02];
    variable_width_pgen(&[0x10, 0x12], &[&record_1, &record_2], 3)
}

fn explicit_phased_dosage_ld_pgen() -> Vec<u8> {
    // rs1 is a non-LD explicit phased dosage record. rs2 is vrtype 0xc2:
    // LD-compressed main track followed by full dosage and phased-dosage
    // auxiliary tracks.
    let mut record_1 = vec![0x25];
    record_1.extend(scaled_dosage(1.0));
    record_1.extend(scaled_dosage(0.5));
    record_1.extend(scaled_dosage(2.0));
    record_1.extend(scaled_phase_delta(0.25, 0.75));
    record_1.extend(scaled_phase_delta(0.0, 0.5));
    record_1.extend(scaled_phase_delta(1.0, 1.0));

    let mut record_2 = vec![0x03, 0x00, 0x00, 0x01, 0x01];
    record_2.extend(scaled_dosage(0.0));
    record_2.extend(scaled_dosage(0.2));
    record_2.extend(scaled_dosage(0.4));
    record_2.extend(scaled_phase_delta(0.0, 0.0));
    record_2.extend(scaled_phase_delta(0.1, 0.1));
    record_2.extend(scaled_phase_delta(0.2, 0.2));

    variable_width_pgen(&[0xc0, 0xc2], &[&record_1, &record_2], 3)
}

fn scaled_dosage(value: f32) -> [u8; 2] {
    let raw = (value / 2.0 * 32768.0).round() as u16;
    raw.to_le_bytes()
}

fn scaled_phase_delta(left: f32, right: f32) -> [u8; 2] {
    let raw = ((left - right) * 16384.0).round() as i16;
    raw.to_le_bytes()
}

fn expected_pgen_dosage(raw: [u8; 2]) -> f32 {
    f32::from(u16::from_le_bytes(raw)) * (2.0 / 32768.0)
}

fn expected_pgen_phase_delta(raw: [u8; 2]) -> f32 {
    f32::from(i16::from_le_bytes(raw)) / 16384.0
}

fn expected_pgen_haplotype_dosages(left: f32, right: f32) -> (f32, f32) {
    let total = expected_pgen_dosage(scaled_dosage(left + right));
    let delta = expected_pgen_phase_delta(scaled_phase_delta(left, right));
    ((total + delta) * 0.5, (total - delta) * 0.5)
}

fn csc_to_dense(sparse: &SparseGenotypeMatrixArrowVariants) -> Vec<f32> {
    let mut dense = vec![0.0; sparse.n_rows * sparse.n_cols];
    for col in 0..sparse.n_cols {
        let start = usize::try_from(sparse.indptr[col]).expect("sparse pointer is nonnegative");
        let end = usize::try_from(sparse.indptr[col + 1]).expect("sparse pointer is nonnegative");
        for offset in start..end {
            let row = usize::try_from(sparse.indices[offset]).expect("sparse row is nonnegative");
            dense[row * sparse.n_cols + col] = sparse.data[offset];
        }
    }
    dense
}

fn write_bad_variable_width_block_offset_fixture(
    name: &str,
) -> (TestDir, PathBuf, PathBuf, PathBuf) {
    let dir = unique_dir(name);
    let mut pgen_bytes = variable_width_pgen(&[0], &[&[0x00]], 4);
    pgen_bytes[12] = pgen_bytes[12].saturating_sub(1);
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
    (dir, pgen, pvar, psam)
}

fn variable_width_two_block_pgen_with_bad_second_offset() -> Vec<u8> {
    let first_block_variant_ct = 65_536_usize;
    let n_variants = first_block_variant_ct + 1;
    let header_len = 12 + 16 + first_block_variant_ct + first_block_variant_ct + 1 + 1;
    let second_block_offset = header_len + first_block_variant_ct - 1;
    let mut bytes = vec![0x6c, 0x1b, 0x10];
    bytes.extend(
        u32::try_from(n_variants)
            .expect("test variant count fits u32")
            .to_le_bytes(),
    );
    bytes.extend(4_u32.to_le_bytes());
    bytes.push(0x04);
    bytes.extend(
        u64::try_from(header_len)
            .expect("test header length fits u64")
            .to_le_bytes(),
    );
    bytes.extend(
        u64::try_from(second_block_offset)
            .expect("test block offset fits u64")
            .to_le_bytes(),
    );
    bytes.extend(std::iter::repeat_n(0_u8, first_block_variant_ct));
    bytes.extend(std::iter::repeat_n(1_u8, first_block_variant_ct));
    bytes.push(0);
    bytes.push(1);
    bytes.extend(std::iter::repeat_n(0_u8, n_variants));
    bytes
}

fn sample_major_window<T: Copy>(
    values: &[T],
    n_samples: usize,
    n_variants: usize,
    start: usize,
    len: usize,
) -> Vec<T> {
    let mut window = Vec::with_capacity(n_samples * len);
    for sample_index in 0..n_samples {
        let row_start = sample_index * n_variants + start;
        window.extend_from_slice(&values[row_start..row_start + len]);
    }
    window
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
    assert_values_with_nan(
        &dense.values,
        &[0.0, 1.0, 2.0, f32::NAN, 0.0, 1.0, 2.0, 1.0, 0.0],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
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
    let variant_metadata = variants(&dense.variants);
    assert_eq!(variant_ref_allele(variant_metadata, 0), Some("A"));
    assert_eq!(variant_alt_allele(variant_metadata, 0), Some("G"));
    assert_eq!(variant_metadata.quals[0], Some(30.0));
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
    assert_values_with_nan(&dense.values, &[0.0, 1.0, 2.0, 2.0, 1.0, 0.0]);
}

#[test]
fn plink2_dense_rejects_unsupported_pgen_modes() {
    let dir = unique_dir("plink2-dense-unsupported-mode");
    let mut pgen_bytes = fixed_width_pgen(&[0x00], 1, 1);
    pgen_bytes[2] = 0x05;
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

    assert_values_with_nan(
        &dense.values,
        &[
            0.0,
            0.0,
            2.0,
            2.0,
            0.0,
            1.0,
            1.0,
            0.0,
            f32::NAN,
            2.0,
            2.0,
            0.0,
            2.0,
            2.0,
            0.0,
            f32::NAN,
            2.0,
            0.0,
            0.0,
            2.0,
        ],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![
            false, false, false, false, false, false, false, false, true, false, false, false,
            false, false, false, true, false, false, false, false,
        ]
    );
}

#[test]
fn plink2_dosage_dense_decodes_variable_width_full_dosage_records() {
    let dir = unique_dir("plink2-dosage-variable-width");
    let mut record_1 = vec![0x24];
    record_1.extend(scaled_dosage(0.2));
    record_1.extend(scaled_dosage(1.4));
    record_1.extend(scaled_dosage(1.8));
    let mut record_2 = vec![0x0c];
    record_2.extend(scaled_dosage(0.0));
    record_2.extend(u16::MAX.to_le_bytes());
    record_2.extend(scaled_dosage(0.7));
    let mut record_3 = vec![0x06];
    record_3.extend(scaled_dosage(2.0));
    record_3.extend(scaled_dosage(1.0));
    record_3.extend(scaled_dosage(0.0));
    let pgen_bytes =
        variable_width_pgen(&[0x40, 0x40, 0x40], &[&record_1, &record_2, &record_3], 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    let dense =
        genoio_io::read_plink2_dosage_dense_windowed(&pgen, &pvar, &psam, None, None, None, false)
            .expect("variable-width dosage pgen should decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 3);
    let scale = 2.0 / 32768.0;
    let expected = vec![
        f32::from(u16::from_le_bytes(scaled_dosage(0.2))) * scale,
        0.0,
        2.0,
        f32::from(u16::from_le_bytes(scaled_dosage(1.4))) * scale,
        f32::NAN,
        1.0,
        f32::from(u16::from_le_bytes(scaled_dosage(1.8))) * scale,
        f32::from(u16::from_le_bytes(scaled_dosage(0.7))) * scale,
        0.0,
    ];
    assert_values_with_nan(&dense_values_sample_major(&dense), &expected);
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, false, true, false, false, false, false]
    );
}

#[test]
fn plink2_haplotype_dense_decodes_explicit_phased_hardcalls() {
    let dir = unique_dir("plink2-haplo-hardcall");
    // v1: S1 0|1, S2 0/0, S3 1/1. v2: S1 1|0, S2 0|1, S3 missing.
    let record_1 = [0x21, 0x00];
    let record_2 = [0x35, 0x02];
    let pgen_bytes = explicit_phased_hardcall_pgen(&[&record_1, &record_2], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
",
    );

    let dense = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect("explicit phased hardcall pgen should decode");

    assert_eq!(dense.n_samples, 6);
    assert_eq!(dense.n_variants, 2);
    assert_values_with_nan(
        &dense_values_sample_major(&dense),
        &[
            0.0,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            f32::NAN,
            1.0,
            f32::NAN,
        ],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, false, false, false, false, false, false, true, false, true]
    );
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S1", "S2", "S2", "S3", "S3"]
    );
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.source_sample_index)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(0), Some(1), Some(1), Some(2), Some(2)]
    );
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.haplotype_index)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(1), Some(0), Some(1), Some(0), Some(1)]
    );
}

#[test]
fn plink2_haplotype_dense_sample_filter_uses_source_order_and_haplotype_order() {
    let dir = unique_dir("plink2-haplo-hardcall-sample-filter");
    let record_1 = [0x21, 0x00];
    let record_2 = [0x35, 0x02];
    let pgen_bytes = explicit_phased_hardcall_pgen(&[&record_1, &record_2], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 C T
",
    );
    let keep = vec!["S3".to_string(), "S1".to_string()];

    let dense = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        None,
        false,
    )
    .expect("explicit phased hardcall pgen should filter samples");

    assert_eq!(dense.n_samples, 4);
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S1", "S3", "S3"]
    );
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.source_sample_index)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(0), Some(2), Some(2)]
    );
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.haplotype_index)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(1), Some(0), Some(1)]
    );
    assert_values_with_nan(
        &dense_values_sample_major(&dense),
        &[0.0, 1.0, 1.0, 0.0, 1.0, f32::NAN, 1.0, f32::NAN],
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, false, false, true, false, true]
    );
}

#[test]
fn plink2_haplotype_dense_ld_compressed_window_matches_full_read_slice() {
    let dir = unique_dir("plink2-haplo-hardcall-ld-window");
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &explicit_phased_hardcall_ld_pgen(),
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
",
    );

    let full = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect("full explicit phased hardcall pgen should decode");
    let window = VariantWindow { start: 1, len: 1 };
    let block = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(window),
        false,
    )
    .expect("window beginning at LD-compressed phased hardcall should decode");

    assert_eq!(block.n_variants, 1);
    assert_eq!(variant_id(variants(&block.variants), 0), "rs2");
    let full_values = dense_values_sample_major(&full);
    let full_missing = dense_missing_sample_major(&full);
    assert_values_with_nan(
        &dense_values_sample_major(&block),
        &sample_major_window(&full_values, full.n_samples, full.n_variants, 1, 1),
    );
    assert_eq!(
        dense_missing_sample_major(&block),
        sample_major_window(&full_missing, full.n_samples, full.n_variants, 1, 1)
    );
}

#[test]
fn plink2_haplotype_sparse_reconstructs_dense_hardcalls() {
    let dir = unique_dir("plink2-haplo-sparse-hardcall");
    let record_1 = [0x21, 0x00];
    let record_2 = [0x35, 0x02];
    let pgen_bytes = explicit_phased_hardcall_pgen(&[&record_1, &record_2], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
",
    );
    let keep = vec!["S1".to_string(), "S2".to_string()];

    let dense = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        None,
        false,
    )
    .expect("dense explicit phased hardcall pgen should decode");
    let sparse = genoio_io::read_plink2_haplotypes_sparse_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        None,
    )
    .expect("sparse explicit phased hardcall pgen should decode");

    assert_eq!(sparse.n_rows, dense.n_samples);
    assert_eq!(sparse.n_cols, dense.n_variants);
    assert_eq!(csc_to_dense(&sparse), dense_values_sample_major(&dense));
    assert_eq!(
        sparse
            .samples
            .iter()
            .map(|sample| (sample.iid.as_str(), sample.haplotype_index))
            .collect::<Vec<_>>(),
        dense
            .samples
            .iter()
            .map(|sample| (sample.iid.as_str(), sample.haplotype_index))
            .collect::<Vec<_>>()
    );
    assert_eq!(variant_ids(variants(&sparse.variants)), vec!["rs1", "rs2"]);
}

#[test]
fn plink2_haplotype_sparse_rejects_retained_missing_hardcalls() {
    let dir = unique_dir("plink2-haplo-sparse-hardcall-missing");
    let record_1 = [0x21, 0x00];
    let record_2 = [0x35, 0x02];
    let pgen_bytes = explicit_phased_hardcall_pgen(&[&record_1, &record_2], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
",
    );

    let error =
        genoio_io::read_plink2_haplotypes_sparse_windowed(&pgen, &pvar, &psam, None, None, None)
            .expect_err("sparse retained missing haplotypes should fail");

    assert!(error.to_string().contains("sparse missing values"));
}

#[test]
fn plink2_haplotype_sparse_ld_compressed_window_matches_dense_slice() {
    let dir = unique_dir("plink2-haplo-sparse-hardcall-ld-window");
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &explicit_phased_hardcall_ld_pgen(),
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
",
    );
    let window = VariantWindow { start: 1, len: 1 };
    let keep = vec!["S1".to_string(), "S2".to_string()];

    let dense = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        Some(window),
        false,
    )
    .expect("dense LD-compressed explicit phased hardcall pgen should decode");
    let sparse = genoio_io::read_plink2_haplotypes_sparse_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        Some(window),
    )
    .expect("sparse LD-compressed explicit phased hardcall pgen should decode");

    assert_eq!(variant_id(variants(&sparse.variants), 0), "rs2");
    assert_eq!(csc_to_dense(&sparse), dense_values_sample_major(&dense));
}

#[test]
fn plink2_haplotype_dense_sample_filter_ignores_unselected_unphased_het() {
    let dir = unique_dir("plink2-haplo-hardcall-unselected-unphased");
    let pgen_bytes = variable_width_pgen(&[0x10], &[&[0x15, 0x0d, 0x02]], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
",
    );
    let keep = vec!["S2".to_string(), "S3".to_string()];

    let dense = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        None,
        false,
    )
    .expect("unselected unphased heterozygote should not reject sample-filtered read");

    assert_eq!(dense.values, vec![0.0, 1.0, 1.0, 0.0]);
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S2", "S2", "S3", "S3"]
    );
}

#[test]
fn plink2_haplotype_dense_unphased_retained_record_fails() {
    let dir = unique_dir("plink2-haplo-hardcall-unphased");
    let pgen_bytes = variable_width_pgen(&[0], &[&[0x21]], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );

    let error = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect_err("unphased retained hardcall record should fail");

    assert!(error.to_string().contains("unphased"));
}

#[test]
fn plink2_haplotype_dosage_dense_decodes_explicit_phased_dosage() {
    let dir = unique_dir("plink2-haplo-dosage");
    let (s1_l, s1_r) = expected_pgen_haplotype_dosages(0.25, 0.75);
    let (s2_l, s2_r) = expected_pgen_haplotype_dosages(0.0, 0.5);
    let (s3_l, s3_r) = expected_pgen_haplotype_dosages(1.0, 1.0);
    let mut record = vec![0x25];
    record.extend(scaled_dosage(1.0));
    record.extend(scaled_dosage(0.5));
    record.extend(scaled_dosage(2.0));
    record.extend(scaled_phase_delta(0.25, 0.75));
    record.extend(scaled_phase_delta(0.0, 0.5));
    record.extend(scaled_phase_delta(1.0, 1.0));
    let pgen_bytes = explicit_phased_dosage_pgen(&[&record], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );

    let dense = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect("explicit phased dosage pgen should decode");

    assert_eq!(dense.n_samples, 6);
    assert_eq!(dense.n_variants, 1);
    assert_eq!(dense.values, vec![s1_l, s1_r, s2_l, s2_r, s3_l, s3_r]);
    assert_eq!(dense_missing_sample_major(&dense), vec![false; 6]);
}

#[test]
fn plink2_haplotype_dosage_dense_decodes_fixed_width_phased_dosage() {
    let dir = unique_dir("plink2-haplo-fixed-width-phased-dosage");
    let (s1_l, s1_r) = expected_pgen_haplotype_dosages(0.25, 0.75);
    let (s2_l, s2_r) = expected_pgen_haplotype_dosages(0.0, 0.5);
    let (s3_l, s3_r) = expected_pgen_haplotype_dosages(1.0, 1.0);
    let mut record = vec![0x25];
    record.extend(scaled_dosage(1.0));
    record.extend(scaled_dosage(0.5));
    record.extend(scaled_dosage(2.0));
    record.extend(scaled_phase_delta(0.25, 0.75));
    record.extend(scaled_phase_delta(0.0, 0.5));
    record.extend(scaled_phase_delta(1.0, 1.0));
    let pgen_bytes = fixed_width_phased_dosage_pgen(&[&record], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );

    let dense = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect("fixed-width explicit phased dosage pgen should decode");

    assert_eq!(dense.n_samples, 6);
    assert_eq!(dense.n_variants, 1);
    assert_eq!(dense.values, vec![s1_l, s1_r, s2_l, s2_r, s3_l, s3_r]);
    assert_eq!(dense_missing_sample_major(&dense), vec![false; 6]);
}

#[test]
fn plink2_haplotype_dosage_missing_values_emit_nan() {
    let dir = unique_dir("plink2-haplo-dosage-missing");
    let (s1_l, s1_r) = expected_pgen_haplotype_dosages(0.25, 0.75);
    let (s3_l, s3_r) = expected_pgen_haplotype_dosages(1.0, 1.0);
    let mut record = vec![0x25];
    record.extend(scaled_dosage(1.0));
    record.extend(u16::MAX.to_le_bytes());
    record.extend(scaled_dosage(2.0));
    record.extend(scaled_phase_delta(0.25, 0.75));
    record.extend(i16::MIN.to_le_bytes());
    record.extend(scaled_phase_delta(1.0, 1.0));
    let pgen_bytes = explicit_phased_dosage_pgen(&[&record], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );

    let dense = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect("explicit phased dosage pgen should decode missing values");

    assert_values_with_nan(&dense.values, &[s1_l, s1_r, f32::NAN, f32::NAN, s3_l, s3_r]);
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, true, true, false, false]
    );
}

fn invalid_phased_dosage_record(dosage_raw: u16, phase_raw: i16) -> Vec<u8> {
    let mut record = vec![0x25];
    record.extend(dosage_raw.to_le_bytes());
    record.extend(scaled_dosage(0.0));
    record.extend(scaled_dosage(0.0));
    record.extend(phase_raw.to_le_bytes());
    record.extend(scaled_phase_delta(0.0, 0.0));
    record.extend(scaled_phase_delta(0.0, 0.0));
    record
}

fn assert_invalid_phased_dosage_record(record: &[u8], expected: &str) {
    let dir = unique_dir("plink2-haplo-dosage-invalid");
    let pgen_bytes = explicit_phased_dosage_pgen(&[record], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );

    let error = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect_err("invalid phased dosage pgen should fail");

    assert!(
        error.to_string().contains(expected),
        "expected error containing {expected:?}, got {error}"
    );
}

#[test]
fn plink2_haplotype_dosage_rejects_out_of_range_total_dosage_raw_values() {
    for dosage_raw in [32_769, 65_534] {
        let record = invalid_phased_dosage_record(dosage_raw, 0);

        assert_invalid_phased_dosage_record(&record, "dosage");
    }
}

#[test]
fn plink2_haplotype_dosage_rejects_out_of_range_phase_raw_values() {
    for phase_raw in [16_385, -16_385] {
        let record = invalid_phased_dosage_record(16_384, phase_raw);

        assert_invalid_phased_dosage_record(&record, "phase");
    }
}

#[test]
fn plink2_haplotype_dosage_rejects_out_of_range_haplotype_components() {
    let record = invalid_phased_dosage_record(0, 16_384);

    assert_invalid_phased_dosage_record(&record, "haplotype");
}

#[test]
fn plink2_haplotype_dosage_genotype_stat_filters_use_collapsed_diploid_dosage() {
    let dir = unique_dir("plink2-haplo-dosage-filter");
    let (v2_s1_l, v2_s1_r) = expected_pgen_haplotype_dosages(0.0, 0.0);
    let (v2_s2_l, v2_s2_r) = expected_pgen_haplotype_dosages(0.1, 0.1);
    let (v2_s3_l, v2_s3_r) = expected_pgen_haplotype_dosages(0.2, 0.2);
    let mut record_1 = vec![0x25];
    record_1.extend(scaled_dosage(1.0));
    record_1.extend(scaled_dosage(0.5));
    record_1.extend(scaled_dosage(2.0));
    record_1.extend(scaled_phase_delta(0.25, 0.75));
    record_1.extend(scaled_phase_delta(0.0, 0.5));
    record_1.extend(scaled_phase_delta(1.0, 1.0));
    let mut record_2 = vec![0x00];
    record_2.extend(scaled_dosage(0.0));
    record_2.extend(scaled_dosage(0.2));
    record_2.extend(scaled_dosage(0.4));
    record_2.extend(scaled_phase_delta(0.0, 0.0));
    record_2.extend(scaled_phase_delta(0.1, 0.1));
    record_2.extend(scaled_phase_delta(0.2, 0.2));
    let pgen_bytes = explicit_phased_dosage_pgen(&[&record_1, &record_2], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 C T
",
    );
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "mac",
        "params": {"max": 1}
    }))
    .expect("filter should parse");

    let dense = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        None,
        false,
    )
    .expect("explicit phased dosage pgen should filter on collapsed dosages");

    assert_eq!(dense.n_samples, 6);
    assert_eq!(dense.n_variants, 1);
    let variant_metadata = variants(&dense.variants);
    assert_eq!(variant_id(variant_metadata, 0), "rs2");
    assert_eq!(
        dense.values,
        vec![v2_s1_l, v2_s1_r, v2_s2_l, v2_s2_r, v2_s3_l, v2_s3_r]
    );
    assert_eq!(dense.diagnostics.dropped_genotype_variants, 1);
    assert!(variant_metadata.afs[0].is_some_and(|af| (af - 0.1).abs() <= 2.0 / 32768.0));
}

#[test]
fn plink2_haplotype_dosage_ld_compressed_window_matches_full_read_slice() {
    let dir = unique_dir("plink2-haplo-dosage-ld-window");
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &explicit_phased_dosage_ld_pgen(),
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
",
    );

    let full = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect("full explicit phased dosage pgen should decode");
    let window = VariantWindow { start: 1, len: 1 };
    let block = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(window),
        false,
    )
    .expect("window beginning at LD-compressed phased dosage should decode");

    assert_eq!(block.n_variants, 1);
    assert_eq!(variant_id(variants(&block.variants), 0), "rs2");
    let full_values = dense_values_sample_major(&full);
    let full_missing = dense_missing_sample_major(&full);
    assert_eq!(
        dense_values_sample_major(&block),
        sample_major_window(&full_values, full.n_samples, full.n_variants, 1, 1)
    );
    assert_eq!(
        dense_missing_sample_major(&block),
        sample_major_window(&full_missing, full.n_samples, full.n_variants, 1, 1)
    );
}

#[test]
fn plink2_haplotype_dosage_unphased_retained_record_fails() {
    let dir = unique_dir("plink2-haplo-dosage-unphased");
    let mut record = vec![0x24];
    record.extend(scaled_dosage(0.2));
    record.extend(scaled_dosage(1.4));
    record.extend(scaled_dosage(1.8));
    let pgen_bytes = variable_width_pgen(&[0x40], &[&record], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );

    let error = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect_err("unphased retained dosage record should fail");

    assert!(error.to_string().contains("phased dosage"));
}

#[test]
fn plink2_dosage_dense_still_decodes_unphased_dosage_records() {
    let dir = unique_dir("plink2-dosage-unphased-regression");
    let mut record = vec![0x24];
    record.extend(scaled_dosage(0.2));
    record.extend(scaled_dosage(1.4));
    record.extend(scaled_dosage(1.8));
    let pgen_bytes = variable_width_pgen(&[0x40], &[&record], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
",
    );

    let dense =
        genoio_io::read_plink2_dosage_dense_windowed(&pgen, &pvar, &psam, None, None, None, false)
            .expect("unphased dosage pgen should still decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 1);
    assert_eq!(
        dense.values,
        vec![
            expected_pgen_dosage(scaled_dosage(0.2)),
            expected_pgen_dosage(scaled_dosage(1.4)),
            expected_pgen_dosage(scaled_dosage(1.8)),
        ]
    );
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false]
    );
}

#[test]
fn plink2_dense_metadata_window_skips_malformed_pvar_records_after_window() {
    let dir = unique_dir("plink2-window-pvar-prefix");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11], 3, 2);
    let pgen = dir.join("tiny.pgen");
    let pvar = dir.join("tiny.pvar");
    let psam = dir.join("tiny.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 bad rs2 A G
",
    );
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
",
    );

    let dense = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect("metadata-bearing window should not parse later pvar records");

    assert_eq!(dense.n_variants, 1);
    assert_eq!(variant_id(variants(&dense.variants), 0), "rs1");
    assert_values_with_nan(&dense.values, &[0.0, f32::NAN, 2.0]);

    let matrix_only = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect("matrix-only first window may skip later pvar records");
    assert_eq!(matrix_only.n_variants, 1);
    assert!(matrix_only.variants.is_none());
    assert_values_with_nan(&matrix_only.values, &[0.0, f32::NAN, 2.0]);
}

#[test]
fn plink2_metadata_source_windows_do_not_require_later_pvar_records() {
    let dir = unique_dir("plink2-window-pvar-missing-later");
    let pgen_bytes = fixed_width_pgen(&[0x00, 0x00, 0x00], 3, 3);
    let pgen = dir.join("tiny.pgen");
    let pvar = dir.join("tiny.pvar");
    let psam = dir.join("tiny.psam");
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
",
    );
    let window = VariantWindow { start: 0, len: 1 };

    let dense =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), false)
            .expect("metadata dense window should not require later pvar records");
    assert_eq!(dense.n_variants, 1);
    assert_eq!(variant_id(variants(&dense.variants), 0), "rs1");

    let sparse =
        genoio_io::read_plink2_sparse_windowed(&pgen, &pvar, &psam, None, None, Some(window))
            .expect("metadata sparse window should not require later pvar records");
    assert_eq!(sparse.n_cols, 1);
    assert_eq!(variant_id(variants(&sparse.variants), 0), "rs1");
}

#[test]
fn plink2_dense_window_aligns_variant_metadata_with_source_window() {
    let dir = unique_dir("plink2-window-metadata");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    let dense = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 1, len: 2 }),
        false,
    )
    .expect("window should decode");

    assert_eq!(dense.n_samples, 3);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(variant_ids(variants(&dense.variants)), vec!["rs2", "rs3"]);
    assert_eq!(dense.values, vec![1.0, 2.0, 0.0, 1.0, 1.0, 0.0]);
    assert_eq!(
        dense_missing_sample_major(&dense),
        vec![false, false, false, false, false, false]
    );
}

#[test]
fn plink2_dense_window_does_not_validate_variable_records_after_window() {
    let dir = unique_dir("plink2-window-pgen-prefix");
    let pgen_bytes = variable_width_pgen(&[0, 5], &[&[0xe4], &[0]], 4);
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

    let dense = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect("first window should not validate later pgen records");

    assert_eq!(dense.n_variants, 1);
    assert_eq!(variant_id(variants(&dense.variants), 0), "rs1");
    assert_values_with_nan(&dense.values, &[0.0, 1.0, 2.0, f32::NAN]);
}

#[test]
fn plink2_dense_matrix_only_window_skips_malformed_metadata() {
    let dir = unique_dir("plink2-matrix-only-malformed-metadata");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    write_text(&pvar, "#CHROM POS ID REF ALT\n1 bad rs1 A G\n");
    let genoio_error = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect_err("metadata-bearing window should parse and reject malformed pvar");
    assert!(genoio_error.to_string().contains("invalid position"));

    let matrix_only = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect("matrix-only window should not parse malformed pvar");
    assert_eq!(matrix_only.n_samples, 3);
    assert_eq!(matrix_only.n_variants, 1);
    assert!(matrix_only.samples.is_empty());
    assert!(matrix_only.variants.is_none());
    assert_values_with_nan(&matrix_only.values, &[0.0, f32::NAN, 2.0]);

    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 C T
1 30 rs3 G A
",
    );
    write_text(&psam, "#IID\n");
    let genoio_error = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        false,
    )
    .expect_err("metadata-bearing window should validate malformed psam dimensions");
    assert!(genoio_error.to_string().contains("sample count"));

    let matrix_only = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 1, len: 1 }),
        true,
    )
    .expect("matrix-only window should not parse malformed psam");
    assert_eq!(matrix_only.n_samples, 3);
    assert_eq!(matrix_only.n_variants, 1);
    assert!(matrix_only.samples.is_empty());
    assert!(matrix_only.variants.is_none());
    assert_eq!(matrix_only.values, vec![1.0, 0.0, 1.0]);
}

#[test]
fn plink2_dense_window_sample_filter_rejects_malformed_psam_even_with_matrix_only_flag() {
    let dir = unique_dir("plink2-window-sample-filter-malformed-psam");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);
    write_text(&psam, "#FID IID\nF1\n");
    let keep = vec!["S1".to_string()];

    let error = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect_err("sample-filtered window should parse and reject malformed psam");

    assert!(error.to_string().contains("too few fields"));
}

#[test]
fn plink2_dense_window_variant_filter_rejects_malformed_pvar_even_with_matrix_only_flag() {
    let dir = unique_dir("plink2-window-variant-filter-malformed-pvar");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);
    write_text(&pvar, "#CHROM POS ID REF ALT\n1 bad rs1 A G\n");
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "1"}
    }))
    .expect("filter should parse");

    let error = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 0, len: 1 }),
        true,
    )
    .expect_err("variant-filtered window should parse and reject malformed pvar");

    assert!(error.to_string().contains("invalid position"));
}

#[test]
fn plink2_dense_matrix_only_window_matches_metadata_source_window_values() {
    let dir = unique_dir("plink2-matrix-only-source-window-values");
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

    let window = VariantWindow { start: 2, len: 2 };
    let metadata_bearing =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), false)
            .expect("metadata-bearing source window should decode");
    let matrix_only =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), true)
            .expect("matrix-only source window should decode");

    assert_eq!(matrix_only.n_samples, metadata_bearing.n_samples);
    assert_eq!(matrix_only.n_variants, metadata_bearing.n_variants);
    assert_values_with_nan(&matrix_only.values, &metadata_bearing.values);
    assert_eq!(
        dense_missing_sample_major(&matrix_only),
        dense_missing_sample_major(&metadata_bearing)
    );
    assert!(matrix_only.samples.is_empty());
    assert!(matrix_only.variants.is_none());
}

#[test]
fn plink2_dense_fixed_width_source_window_matches_full_read_slice() {
    let dir = unique_dir("plink2-fixed-window-matches-full");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11, 0x06, 0x3f], 3, 4);
    let pgen = dir.join("fixed.pgen");
    let pvar = dir.join("fixed.pvar");
    let psam = dir.join("fixed.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(
        &pvar,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
1 20 rs2 A G
1 30 rs3 A G
1 40 rs4 A G
",
    );
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
",
    );

    let full = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect("full fixed-width pgen should decode");
    let window = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 1, len: 2 }),
        false,
    )
    .expect("fixed-width source window should decode");
    let matrix_only = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        None,
        None,
        Some(VariantWindow { start: 1, len: 2 }),
        true,
    )
    .expect("fixed-width matrix-only source window should decode");

    assert_eq!(window.n_samples, 3);
    assert_eq!(window.n_variants, 2);
    assert_eq!(
        window.values,
        vec![
            full.values[1],
            full.values[2],
            full.values[5],
            full.values[6],
            full.values[9],
            full.values[10]
        ]
    );
    assert_eq!(
        dense_missing_sample_major(&window),
        vec![
            dense_missing_sample_major(&full)[1],
            dense_missing_sample_major(&full)[2],
            dense_missing_sample_major(&full)[5],
            dense_missing_sample_major(&full)[6],
            dense_missing_sample_major(&full)[9],
            dense_missing_sample_major(&full)[10]
        ]
    );
    assert_values_with_nan(&matrix_only.values, &window.values);
    assert_eq!(
        dense_missing_sample_major(&matrix_only),
        dense_missing_sample_major(&window)
    );
}

#[test]
fn plink2_dense_fixed_width_source_window_crosses_packed_batch_boundary() {
    let dir = unique_dir("plink2-fixed-window-batch-boundary");
    let n_variants = 70;
    let records = (0..n_variants)
        .map(|variant_index| match variant_index % 4 {
            0 => 0x2c,
            1 => 0x11,
            2 => 0x06,
            _ => 0x3f,
        })
        .collect::<Vec<_>>();
    let pgen_bytes = fixed_width_pgen(&records, 3, n_variants);
    let pgen = dir.join("fixed.pgen");
    let pvar = dir.join("fixed.pvar");
    let psam = dir.join("fixed.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    let pvar_body = (0..n_variants)
        .map(|index| format!("1 {} rs{} A G\n", 10 + index, index + 1))
        .collect::<String>();
    write_text(&pvar, &format!("#CHROM POS ID REF ALT\n{pvar_body}"));
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
",
    );

    let full = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect("full fixed-width pgen should decode");
    let window = VariantWindow { start: 1, len: 66 };
    let metadata_bearing =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), false)
            .expect("metadata-bearing batch-spanning source window should decode");
    let matrix_only =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), true)
            .expect("matrix-only batch-spanning source window should decode");

    let expected_values = sample_major_window(
        &full.values,
        full.n_samples,
        full.n_variants,
        window.start,
        window.len,
    );
    let expected_missing = sample_major_window(
        &dense_missing_sample_major(&full),
        full.n_samples,
        full.n_variants,
        window.start,
        window.len,
    );
    assert_values_with_nan(&metadata_bearing.values, &expected_values);
    assert_eq!(
        dense_missing_sample_major(&metadata_bearing),
        expected_missing
    );
    assert_values_with_nan(&matrix_only.values, &expected_values);
    assert_eq!(dense_missing_sample_major(&matrix_only), expected_missing);
    let variant_ids = variant_ids(variants(&metadata_bearing.variants));
    assert_eq!(variant_ids.first().copied(), Some("rs2"));
    assert_eq!(variant_ids.last().copied(), Some("rs67"));
}

#[test]
fn plink2_dense_variable_width_source_window_matches_full_read_slice() {
    let dir = unique_dir("plink2-variable-window-matches-full");
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

    let full = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect("full variable-width pgen should decode");
    let window = VariantWindow { start: 1, len: 3 };
    let metadata_bearing =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), false)
            .expect("variable-width source window should decode");
    let matrix_only =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), true)
            .expect("variable-width matrix-only source window should decode");
    let keep = vec!["S4".to_string(), "S2".to_string()];
    let filtered_full = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, Some(&keep), None)
        .expect("sample-filtered full variable-width pgen should decode");
    let filtered_window = genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        Some(window),
        false,
    )
    .expect("sample-filtered variable-width source window should decode");

    let expected_values = sample_major_window(
        &full.values,
        full.n_samples,
        full.n_variants,
        window.start,
        window.len,
    );
    let expected_missing = sample_major_window(
        &dense_missing_sample_major(&full),
        full.n_samples,
        full.n_variants,
        window.start,
        window.len,
    );
    assert_values_with_nan(&metadata_bearing.values, &expected_values);
    assert_eq!(
        dense_missing_sample_major(&metadata_bearing),
        expected_missing
    );
    assert_values_with_nan(&matrix_only.values, &expected_values);
    assert_eq!(dense_missing_sample_major(&matrix_only), expected_missing);

    let expected_filtered_values = sample_major_window(
        &filtered_full.values,
        filtered_full.n_samples,
        filtered_full.n_variants,
        window.start,
        window.len,
    );
    let expected_filtered_missing = sample_major_window(
        &dense_missing_sample_major(&filtered_full),
        filtered_full.n_samples,
        filtered_full.n_variants,
        window.start,
        window.len,
    );
    assert_values_with_nan(&filtered_window.values, &expected_filtered_values);
    assert_eq!(
        dense_missing_sample_major(&filtered_window),
        expected_filtered_missing
    );
    assert_eq!(
        filtered_window
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S2", "S4"]
    );
}

#[test]
fn plink2_dense_rejects_truncated_one_bit_variable_record() {
    let dir = unique_dir("plink2-dense-truncated-one-bit");
    let pgen_bytes = variable_width_pgen(&[1], &[&[2]], 9);
    let pgen = dir.join("bad_one_bit.pgen");
    let pvar = dir.join("bad_one_bit.pvar");
    let psam = dir.join("bad_one_bit.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(&pvar, "#CHROM POS ID REF ALT\n1 10 rs1 A G\n");
    write_text(
        &psam,
        "\
#IID
S1
S2
S3
S4
S5
S6
S7
S8
S9
",
    );

    let error = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect_err("truncated 1-bit pgen record should fail");

    assert!(error.to_string().contains("1-bit record is shorter"));
}

#[test]
fn plink2_dense_rejects_truncated_fixed_width_records() {
    let dir = unique_dir("plink2-dense-truncated-fixed");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    let error = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect_err("truncated fixed-width pgen should fail");

    assert!(error.to_string().contains("truncated"));
}

#[test]
fn plink2_dense_rejects_unsupported_variable_width_compression() {
    let dir = unique_dir("plink2-dense-unsupported-compression");
    let pgen_bytes = variable_width_pgen(&[5], &[&[0]], 4);
    let pgen = dir.join("bad_compression.pgen");
    let pvar = dir.join("bad_compression.pvar");
    let psam = dir.join("bad_compression.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(&pvar, "#CHROM POS ID REF ALT\n1 10 rs1 A G\n");
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
        .expect_err("unsupported variable-width compression should fail");

    assert!(error
        .to_string()
        .contains("unsupported pgen main-track compression type 5"));
}

#[test]
fn plink2_dense_rejects_ld_compressed_record_before_non_ld_state() {
    let dir = unique_dir("plink2-dense-ld-before-state");
    let pgen_bytes = variable_width_pgen(&[2], &[&[0]], 4);
    let pgen = dir.join("ld_before_state.pgen");
    let pvar = dir.join("ld_before_state.pvar");
    let psam = dir.join("ld_before_state.psam");
    fs::write(&pgen, pgen_bytes).expect("pgen fixture should be written");
    write_text(&pvar, "#CHROM POS ID REF ALT\n1 10 rs1 A G\n");
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
        .expect_err("LD-compressed first record should fail");

    assert!(error
        .to_string()
        .contains("LD-compressed record appears before any non-LD record"));
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
fn plink2_dense_matrix_only_source_window_rejects_variable_width_block_offset_mismatch() {
    let (_dir, pgen, pvar, psam) = write_bad_variable_width_block_offset_fixture(
        "plink2-dense-window-bad-block-offset-matrix-only",
    );
    let window = VariantWindow { start: 0, len: 1 };

    let error =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), true)
            .expect_err("matrix-only source window should reject bad block offset");

    assert!(error.to_string().contains("block offset"));
}

#[test]
fn plink2_dense_metadata_source_window_rejects_variable_width_block_offset_mismatch() {
    let (_dir, pgen, pvar, psam) = write_bad_variable_width_block_offset_fixture(
        "plink2-dense-window-bad-block-offset-metadata",
    );
    let window = VariantWindow { start: 0, len: 1 };

    let error =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), false)
            .expect_err("metadata source window should reject bad block offset");

    assert!(error.to_string().contains("block offset"));
}

#[test]
fn plink2_dense_source_window_rejects_variable_width_block_offset_mismatch_after_block_boundary() {
    let dir = unique_dir("plink2-dense-window-bad-second-block-offset");
    let pgen = dir.join("bad_second_offset.pgen");
    fs::write(
        &pgen,
        variable_width_two_block_pgen_with_bad_second_offset(),
    )
    .expect("pgen fixture should be written");
    let pvar = dir.join("unused.pvar");
    let psam = dir.join("unused.psam");
    let window = VariantWindow {
        start: 65_536,
        len: 1,
    };

    let error =
        genoio_io::read_plink2_dense_windowed(&pgen, &pvar, &psam, None, None, Some(window), true)
            .expect_err("source window crossing block boundary should reject bad block offset");

    assert!(error.to_string().contains("block offset"));
}

#[test]
fn plink2_metadata_reads_psam_and_pvar() {
    let dir = unique_dir("plink2-metadata");
    let pgen_bytes = fixed_width_pgen(&[0x00, 0x00, 0x00], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    let metadata =
        genoio_io::read_plink2_metadata_arrow(&pgen, &pvar, &psam).expect("metadata should decode");

    assert_eq!(metadata.samples.len(), 3);
    assert_eq!(metadata.variants.len(), 3);
    assert!(metadata.capabilities.supports_geno);
    assert!(!metadata.capabilities.supports_haplo);
}

#[test]
fn plink2_metadata_reads_zstd_compressed_pvar() {
    let dir = unique_dir("plink2-compressed-pvar");
    let pgen_bytes = fixed_width_pgen(&[0x00, 0x00, 0x00], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);
    let pvar_zst = pvar.with_extension("pvar.zst");
    let pvar_contents = fs::read(&pvar).expect("pvar fixture should be readable");
    let compressed = zstd::stream::encode_all(&pvar_contents[..], 0)
        .expect("pvar fixture should compress as zstd");
    fs::write(&pvar_zst, compressed).expect("compressed pvar fixture should be written");
    fs::remove_file(&pvar).expect("uncompressed pvar fixture should be removed");

    let metadata = genoio_io::read_plink2_metadata_arrow(&pgen, &pvar_zst, &psam)
        .expect("metadata should decode");

    assert_eq!(metadata.variants.len(), 3);
    assert_eq!(variant_id(&metadata.variants, 0), "rs1");
}
