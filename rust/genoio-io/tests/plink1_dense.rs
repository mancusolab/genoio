// pattern: Imperative Shell

use std::fs;
use std::path::{Path, PathBuf};

use ::genoio_io::{
    BlockOutput, BlockReadOptions, BlockReader, BlockSource, DosageSource, MatrixKind,
};
use genoio_core::{DenseGenotypeMatrix, DenseMissingPolicy, GenoioError, VariantFilter};

mod common;

use common::dense::assert_values_with_nan;
use common::plink_output as genoio_io;
use common::plink_output::{
    dense_missing_sample_major_output as dense_missing_sample_major,
    dense_values_sample_major_output, variant_ids, variants,
};
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

fn plink1_block_options(
    requested_samples: Option<Vec<String>>,
    variant_filter: Option<VariantFilter>,
    missing_policy: DenseMissingPolicy,
    return_metadata: bool,
) -> BlockReadOptions {
    BlockReadOptions {
        matrix_kind: MatrixKind::Genotype,
        sparse: false,
        requested_samples,
        variant_filter,
        dosage_source: DosageSource::Hardcall,
        missing_policy,
        return_samples: return_metadata,
        return_variants: return_metadata,
    }
}

fn collect_plink1_dense_blocks(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    options: BlockReadOptions,
    block_size: usize,
) -> Vec<DenseGenotypeMatrix> {
    let mut reader = BlockReader::open(
        BlockSource::Plink1 {
            bed: bed.to_path_buf(),
            bim: bim.to_path_buf(),
            fam: fam.to_path_buf(),
        },
        options,
        block_size,
    )
    .expect("persistent plink1 reader should open");
    let mut blocks = Vec::new();
    while let Some(output) = reader
        .next_block()
        .expect("persistent plink1 block should decode")
    {
        let BlockOutput::Dense(matrix) = output else {
            panic!("plink1 dense reader should return dense blocks");
        };
        blocks.push(matrix);
    }
    assert!(reader
        .next_block()
        .expect("persistent plink1 EOF should be sticky")
        .is_none());
    blocks
}

fn concatenate_dense_blocks_sample_major(blocks: &[DenseGenotypeMatrix]) -> Vec<f32> {
    let Some(first) = blocks.first() else {
        return Vec::new();
    };
    let block_values = blocks
        .iter()
        .map(dense_values_sample_major_output)
        .collect::<Vec<_>>();
    let mut values = Vec::with_capacity(
        first.n_samples * blocks.iter().map(|block| block.n_variants).sum::<usize>(),
    );
    for sample_index in 0..first.n_samples {
        for (block, block_values) in blocks.iter().zip(&block_values) {
            let start = sample_index * block.n_variants;
            values.extend_from_slice(&block_values[start..start + block.n_variants]);
        }
    }
    values
}

fn concatenate_block_variant_ids(blocks: &[DenseGenotypeMatrix]) -> Vec<String> {
    blocks
        .iter()
        .flat_map(|block| {
            variant_ids(variants(&block.variants))
                .into_iter()
                .map(str::to_owned)
        })
        .collect()
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

#[test]
fn pbr_rust_plink1_001_dense_blocks_match_whole_read_at_exact_and_partial_boundaries() {
    let dir = unique_dir("pbr-plink1-dense-boundaries");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x07, 0x2d, 0x38]);
    let expected = ::genoio_io::read_plink1_dense_windowed(
        &bed,
        &bim,
        &fam,
        None,
        None,
        None,
        DenseMissingPolicy::Nan,
        true,
        true,
    )
    .expect("whole plink1 read should decode");

    for (block_size, expected_widths) in [(3, vec![3]), (2, vec![2, 1])] {
        let blocks = collect_plink1_dense_blocks(
            &bed,
            &bim,
            &fam,
            plink1_block_options(None, None, DenseMissingPolicy::Nan, true),
            block_size,
        );

        assert_eq!(
            blocks
                .iter()
                .map(|block| block.n_variants)
                .collect::<Vec<_>>(),
            expected_widths
        );
        assert_values_with_nan(
            &concatenate_dense_blocks_sample_major(&blocks),
            &dense_values_sample_major_output(&expected),
        );
        assert_eq!(
            concatenate_block_variant_ids(&blocks),
            variant_ids(variants(&expected.variants))
        );
        assert!(blocks.iter().all(|block| block.samples == expected.samples));
    }
}

#[test]
fn pbr_rust_plink1_001_dense_blocks_preserve_samples_filters_metadata_and_imputation() {
    let dir = unique_dir("pbr-plink1-dense-filters");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x07, 0x2d, 0x38]);
    let requested_samples = vec!["S3".to_owned(), "S1".to_owned()];
    let filter = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("genotype-stat filter should parse");
    let expected = ::genoio_io::read_plink1_dense_windowed(
        &bed,
        &bim,
        &fam,
        Some(&requested_samples),
        Some(&filter),
        None,
        DenseMissingPolicy::Impute,
        true,
        true,
    )
    .expect("filtered whole plink1 read should decode");

    let blocks = collect_plink1_dense_blocks(
        &bed,
        &bim,
        &fam,
        plink1_block_options(
            Some(requested_samples),
            Some(filter),
            DenseMissingPolicy::Impute,
            true,
        ),
        1,
    );

    assert_eq!(
        concatenate_dense_blocks_sample_major(&blocks),
        dense_values_sample_major_output(&expected)
    );
    assert_eq!(
        concatenate_block_variant_ids(&blocks),
        variant_ids(variants(&expected.variants))
    );
    assert!(blocks.iter().all(|block| block.samples == expected.samples));
}

#[test]
fn pbr_rust_plink1_001_dense_blocks_handle_matrix_only_all_filtered_and_missing_rejection() {
    let dir = unique_dir("pbr-plink1-dense-empty-missing");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x07, 0x2d, 0x38]);
    let all_filtered = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "X"}
    }))
    .expect("metadata filter should parse");

    let matrix_only = collect_plink1_dense_blocks(
        &bed,
        &bim,
        &fam,
        plink1_block_options(None, None, DenseMissingPolicy::Nan, false),
        2,
    );
    assert!(matrix_only
        .iter()
        .all(|block| block.samples.is_none() && block.variants.is_none()));

    let empty = collect_plink1_dense_blocks(
        &bed,
        &bim,
        &fam,
        plink1_block_options(None, Some(all_filtered), DenseMissingPolicy::Nan, true),
        2,
    );
    assert!(empty.is_empty());

    let mut reader = BlockReader::open(
        BlockSource::Plink1 { bed, bim, fam },
        plink1_block_options(None, None, DenseMissingPolicy::Raise, true),
        2,
    )
    .expect("persistent plink1 reader should open");
    let error = reader
        .next_block()
        .expect_err("retained missing hard calls should be rejected");
    assert!(matches!(error, GenoioError::MissingData { .. }));
}

#[test]
fn pbr_rust_plink1_001_rejects_unsupported_plink1_block_representations() {
    let source = BlockSource::Plink1 {
        bed: PathBuf::from("unused.bed"),
        bim: PathBuf::from("unused.bim"),
        fam: PathBuf::from("unused.fam"),
    };
    let options = plink1_block_options(None, None, DenseMissingPolicy::Nan, true);

    let dosage_error = BlockReader::open(
        source.clone(),
        BlockReadOptions {
            dosage_source: DosageSource::Dosage,
            ..options.clone()
        },
        2,
    )
    .err()
    .expect("plink1 dosage blocks should be rejected before source opening");
    let haplotype_error = BlockReader::open(
        source,
        BlockReadOptions {
            matrix_kind: MatrixKind::Haplotype,
            ..options
        },
        2,
    )
    .err()
    .expect("plink1 haplotype blocks should be rejected before source opening");

    assert!(matches!(
        dosage_error,
        GenoioError::UnsupportedRepresentation { .. }
    ));
    assert!(matches!(
        haplotype_error,
        GenoioError::UnsupportedRepresentation { .. }
    ));
}
