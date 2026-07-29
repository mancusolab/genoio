// pattern: Imperative Shell

use std::fs;
use std::path::{Path, PathBuf};

use ::genoio_io::{
    BlockOutput, BlockReadOptions, BlockReader, BlockSource, DosageSource, MatrixKind,
};
use genoio_core::{
    DenseGenotypeMatrix, DenseMissingPolicy, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};

mod common;

use common::dense::assert_values_with_nan;
use common::plink_output as genoio_io;
use common::plink_output::{
    dense_missing_sample_major_output as dense_missing_sample_major,
    dense_values_sample_major_output as dense_values_sample_major, variant_a0, variant_a1,
    variant_id, variant_ids, variants,
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
    const VARIANT_BLOCK_SIZE: usize = 65_536;

    assert_eq!(
        record_types.len(),
        records.len(),
        "test record types and payloads must align"
    );
    let n_variants = u32::try_from(records.len()).expect("test variant count fits u32");
    let block_ct = records.len().div_ceil(VARIANT_BLOCK_SIZE);
    let header_len = 12 + block_ct * 8 + record_types.len() + records.len();
    let mut bytes = vec![0x6c, 0x1b, 0x10];
    bytes.extend(n_variants.to_le_bytes());
    bytes.extend(n_samples.to_le_bytes());
    bytes.push(0x04);
    let mut block_offset = header_len;
    for block_records in records.chunks(VARIANT_BLOCK_SIZE) {
        bytes.extend(
            u64::try_from(block_offset)
                .expect("test block offset fits u64")
                .to_le_bytes(),
        );
        block_offset += block_records
            .iter()
            .map(|record| record.len())
            .sum::<usize>();
    }
    for (block_types, block_records) in record_types
        .chunks(VARIANT_BLOCK_SIZE)
        .zip(records.chunks(VARIANT_BLOCK_SIZE))
    {
        bytes.extend(block_types);
        bytes.extend(
            block_records.iter().map(|record| {
                u8::try_from(record.len()).expect("test record length fits one byte")
            }),
        );
    }
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

fn fixed_width_dosage_pgen(records: &[&[u8]], n_samples: u32) -> Vec<u8> {
    // PGEN fixed-width mode 0x03 stores each variant's hardcall main track
    // followed by one full unphased dosage value per sample.
    let n_variants = u32::try_from(records.len()).expect("test variant count fits u32");
    let mut bytes = vec![0x6c, 0x1b, 0x03];
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

fn variable_phased_dosage_existence_pgen() -> Vec<u8> {
    // Type 1: one sparse dosage at S2, implicitly phased by the 0|1
    // hardcall. A total dosage of 1.25 therefore maps to (0.25, 1.0).
    let mut implicit_difflist = vec![0x24, 0x00, 1, 1];
    implicit_difflist.extend(scaled_dosage(1.25));

    // Type 3: one S2 dosage selected by the existence bitarray, with the
    // explicit phase-existence bit mapping the following delta to that dosage.
    let mut mapped_bitarray = vec![0x24, 0x00, 0b0000_0010];
    mapped_bitarray.extend(scaled_dosage(1.0));
    mapped_bitarray.push(0b0000_0001);
    mapped_bitarray.extend(scaled_phase_delta(0.25, 0.75));

    // Type 2: full dosage and full phase arrays. S2 uses the representation's
    // paired missing sentinels; S1/S3 have valid endpoint dosages.
    let mut full = vec![0x24];
    full.extend(scaled_dosage(0.0));
    full.extend(u16::MAX.to_le_bytes());
    full.extend(scaled_dosage(2.0));
    full.extend(scaled_phase_delta(0.0, 0.0));
    full.extend(i16::MIN.to_le_bytes());
    full.extend(scaled_phase_delta(1.0, 1.0));

    variable_width_pgen(
        &[0x30, 0xf0, 0xc0],
        &[&implicit_difflist, &mapped_bitarray, &full],
        3,
    )
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

fn csc_to_dense(sparse: &SparseGenotypeMatrix) -> Vec<f32> {
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

fn plink2_block_options(
    sparse: bool,
    requested_samples: Option<Vec<String>>,
    variant_filter: Option<VariantFilter>,
    missing_policy: DenseMissingPolicy,
    return_metadata: bool,
) -> BlockReadOptions {
    plink2_mode_block_options(
        MatrixKind::Genotype,
        sparse,
        requested_samples,
        variant_filter,
        DosageSource::Hardcall,
        missing_policy,
        return_metadata,
    )
}

fn plink2_mode_block_options(
    matrix_kind: MatrixKind,
    sparse: bool,
    requested_samples: Option<Vec<String>>,
    variant_filter: Option<VariantFilter>,
    dosage_source: DosageSource,
    missing_policy: DenseMissingPolicy,
    return_metadata: bool,
) -> BlockReadOptions {
    BlockReadOptions {
        matrix_kind,
        sparse,
        requested_samples,
        variant_filter,
        dosage_source,
        missing_policy,
        return_samples: return_metadata,
        return_variants: return_metadata,
    }
}

fn collect_plink2_blocks(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    options: BlockReadOptions,
    block_size: usize,
) -> Vec<BlockOutput> {
    let mut reader = BlockReader::open(
        BlockSource::Plink2 {
            pgen: pgen.to_path_buf(),
            pvar: pvar.to_path_buf(),
            psam: psam.to_path_buf(),
        },
        options,
        block_size,
    )
    .expect("persistent plink2 reader should open");
    let mut blocks = Vec::new();
    while let Some(block) = reader
        .next_block()
        .expect("persistent plink2 block should decode")
    {
        blocks.push(block);
    }
    assert!(reader
        .next_block()
        .expect("persistent plink2 EOF should be sticky")
        .is_none());
    blocks
}

fn concatenate_plink2_dense_blocks(blocks: &[DenseGenotypeMatrix]) -> Vec<f32> {
    let Some(first) = blocks.first() else {
        return Vec::new();
    };
    let block_values = blocks
        .iter()
        .map(dense_values_sample_major)
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

fn concatenate_plink2_sparse_blocks(blocks: &[SparseGenotypeMatrix]) -> Vec<f32> {
    let Some(first) = blocks.first() else {
        return Vec::new();
    };
    let block_values = blocks.iter().map(csc_to_dense).collect::<Vec<_>>();
    let mut values =
        Vec::with_capacity(first.n_rows * blocks.iter().map(|block| block.n_cols).sum::<usize>());
    for sample_index in 0..first.n_rows {
        for (block, block_values) in blocks.iter().zip(&block_values) {
            let start = sample_index * block.n_cols;
            values.extend_from_slice(&block_values[start..start + block.n_cols]);
        }
    }
    values
}

fn assert_dense_matrix_matches_oracle(
    actual: &DenseGenotypeMatrix,
    expected: &DenseGenotypeMatrix,
) {
    assert_eq!(actual.n_samples, expected.n_samples);
    assert_eq!(actual.n_variants, expected.n_variants);
    assert_values_with_nan(
        &dense_values_sample_major(actual),
        &dense_values_sample_major(expected),
    );
    assert_eq!(
        dense_missing_sample_major(actual),
        dense_missing_sample_major(expected)
    );
    assert_eq!(actual.samples, expected.samples);
    assert_eq!(actual.variants, expected.variants);
    assert_eq!(actual.diagnostics, expected.diagnostics);
}

#[test]
fn pbr_rust_plink2_001_fixed_width_dense_and_sparse_blocks_match_whole_reads() {
    let dir = unique_dir("pbr-plink2-fixed-block-parity");
    let pgen_bytes = fixed_width_pgen(&[0x24, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);

    let expected_dense = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, None)
        .expect("whole dense plink2 read should decode");
    let dense_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_block_options(false, None, None, DenseMissingPolicy::Nan, true),
        2,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("dense PLINK2 reader returned a sparse block"),
    })
    .collect::<Vec<_>>();
    assert_eq!(
        dense_blocks
            .iter()
            .map(|block| block.n_variants)
            .collect::<Vec<_>>(),
        vec![2, 1]
    );
    assert_values_with_nan(
        &concatenate_plink2_dense_blocks(&dense_blocks),
        &dense_values_sample_major(&expected_dense),
    );
    assert_eq!(
        dense_blocks
            .iter()
            .flat_map(|block| variant_ids(variants(&block.variants)))
            .collect::<Vec<_>>(),
        variant_ids(variants(&expected_dense.variants))
    );

    let expected_sparse = genoio_io::read_plink2_sparse(&pgen, &pvar, &psam, None, None)
        .expect("whole sparse plink2 read should decode");
    let sparse_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_block_options(true, None, None, DenseMissingPolicy::Raise, true),
        2,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Sparse(matrix) => matrix,
        BlockOutput::Dense(_) => panic!("sparse PLINK2 reader returned a dense block"),
    })
    .collect::<Vec<_>>();
    assert_eq!(
        concatenate_plink2_sparse_blocks(&sparse_blocks),
        csc_to_dense(&expected_sparse)
    );
}

#[test]
fn pbr_rust_plink2_001_fixed_width_blocks_preserve_filters_metadata_and_missingness() {
    let dir = unique_dir("pbr-plink2-fixed-block-filters");
    let pgen_bytes = fixed_width_pgen(&[0x2c, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);
    let requested_samples = vec!["S3".to_owned(), "S1".to_owned()];
    let genotype_filter = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("genotype-stat filter should parse");
    let expected = ::genoio_io::read_plink2_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&requested_samples),
        Some(&genotype_filter),
        None,
        DenseMissingPolicy::Impute,
        true,
        true,
    )
    .expect("filtered whole plink2 read should decode");
    let blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_block_options(
            false,
            Some(requested_samples),
            Some(genotype_filter),
            DenseMissingPolicy::Impute,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("dense PLINK2 reader returned a sparse block"),
    })
    .collect::<Vec<_>>();
    assert_eq!(
        concatenate_plink2_dense_blocks(&blocks),
        dense_values_sample_major(&expected)
    );
    assert!(blocks.iter().all(|block| block.samples == expected.samples));

    let matrix_only = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_block_options(false, None, None, DenseMissingPolicy::Nan, false),
        2,
    );
    assert!(matrix_only.iter().all(|block| match block {
        BlockOutput::Dense(matrix) => matrix.samples.is_none() && matrix.variants.is_none(),
        BlockOutput::Sparse(_) => false,
    }));

    let all_filtered = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "X"}
    }))
    .expect("metadata filter should parse");
    assert!(collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_block_options(
            false,
            None,
            Some(all_filtered),
            DenseMissingPolicy::Nan,
            true,
        ),
        2,
    )
    .is_empty());

    let mut missing_reader = BlockReader::open(
        BlockSource::Plink2 { pgen, pvar, psam },
        plink2_block_options(false, None, None, DenseMissingPolicy::Raise, true),
        1,
    )
    .expect("persistent plink2 missingness reader should open");
    assert!(matches!(
        missing_reader
            .next_block()
            .expect_err("retained missing hard calls should be rejected"),
        genoio_core::GenoioError::MissingData { .. }
    ));
}

#[test]
fn pbr_rust_plink2_001_header_and_later_pvar_errors_keep_their_boundaries() {
    let dir = unique_dir("pbr-plink2-error-boundaries");
    let pgen_bytes = fixed_width_pgen(&[0x24, 0x11, 0x06], 3, 3);
    let (pgen, pvar, psam) = write_plink2_fixture(&dir, &pgen_bytes);
    let invalid_pgen = dir.join("invalid.pgen");
    fs::write(&invalid_pgen, [0_u8; 12]).expect("invalid pgen fixture should be written");
    assert!(BlockReader::open(
        BlockSource::Plink2 {
            pgen: invalid_pgen,
            pvar: pvar.clone(),
            psam: psam.clone(),
        },
        plink2_block_options(false, None, None, DenseMissingPolicy::Nan, true),
        1,
    )
    .is_err());

    write_text(
        &pvar,
        "#CHROM POS ID REF ALT\n1 10 rs1 A G\n1 bad rs2 C T\n2 30 rs3 G A\n",
    );
    let mut reader = BlockReader::open(
        BlockSource::Plink2 { pgen, pvar, psam },
        plink2_block_options(false, None, None, DenseMissingPolicy::Nan, true),
        1,
    )
    .expect("later-pvar session should construct");
    assert!(reader
        .next_block()
        .expect("first valid PLINK2 block should decode")
        .is_some());
    assert!(reader
        .next_block()
        .expect_err("malformed later PVAR row should fail when reached")
        .to_string()
        .contains("invalid position"));
}

#[test]
fn pbr_rust_plink2_002_ld_base_survives_a_rejected_record_across_blocks() {
    let dir = unique_dir("pbr-plink2-ld-rejected-base");
    let pgen_bytes = variable_width_pgen(&[0, 0, 2], &[&[0x24], &[0x06], &[0]], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "\
#CHROM POS ID REF ALT
1 10 rs1 A G
2 20 rejected_base C T
1 30 rs3 G A
",
    );
    let filter = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "1"}
    }))
    .expect("metadata filter should parse");

    let expected_dense = genoio_io::read_plink2_dense(&pgen, &pvar, &psam, None, Some(&filter))
        .expect("whole LD-compressed dense read should decode");
    let dense_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_block_options(
            false,
            None,
            Some(filter.clone()),
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("dense PLINK2 reader returned a sparse block"),
    })
    .collect::<Vec<_>>();
    assert_eq!(
        dense_blocks
            .iter()
            .map(|block| block.n_variants)
            .collect::<Vec<_>>(),
        vec![1, 1]
    );
    assert_values_with_nan(
        &concatenate_plink2_dense_blocks(&dense_blocks),
        &dense_values_sample_major(&expected_dense),
    );
    assert_eq!(
        dense_values_sample_major(&dense_blocks[1]),
        vec![2.0, 1.0, 0.0]
    );

    let expected_sparse = genoio_io::read_plink2_sparse(&pgen, &pvar, &psam, None, Some(&filter))
        .expect("whole LD-compressed sparse read should decode");
    let sparse_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_block_options(true, None, Some(filter), DenseMissingPolicy::Raise, true),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Sparse(matrix) => matrix,
        BlockOutput::Dense(_) => panic!("sparse PLINK2 reader returned a dense block"),
    })
    .collect::<Vec<_>>();
    assert_eq!(
        concatenate_plink2_sparse_blocks(&sparse_blocks),
        csc_to_dense(&expected_sparse)
    );
}

#[test]
fn pbr_rust_plink2_002_variable_width_pvar_errors_are_delayed_and_counts_are_validated() {
    let dir = unique_dir("pbr-plink2-variable-pvar-boundaries");
    let pgen_bytes = variable_width_pgen(&[0, 0, 0], &[&[0x24], &[0x11], &[0x06]], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "#CHROM POS ID REF ALT\n1 10 rs1 A G\n1 bad rs2 C T\n1 30 rs3 G A\n",
    );
    let options = plink2_block_options(false, None, None, DenseMissingPolicy::Nan, true);
    let mut malformed = BlockReader::open(
        BlockSource::Plink2 {
            pgen: pgen.clone(),
            pvar: pvar.clone(),
            psam: psam.clone(),
        },
        options.clone(),
        1,
    )
    .expect("malformed-later PVAR session should open");
    assert!(malformed
        .next_block()
        .expect("first valid variable-width block should decode")
        .is_some());
    assert!(malformed
        .next_block()
        .expect_err("malformed later variable-width PVAR row should fail when reached")
        .to_string()
        .contains("invalid position"));

    write_text(&pvar, "#CHROM POS ID REF ALT\n1 10 rs1 A G\n1 20 rs2 C T\n");
    let mut short = BlockReader::open(
        BlockSource::Plink2 {
            pgen: pgen.clone(),
            pvar: pvar.clone(),
            psam: psam.clone(),
        },
        options.clone(),
        2,
    )
    .expect("short-PVAR session should open");
    assert!(short
        .next_block()
        .expect("available short-PVAR prefix should decode")
        .is_some());
    assert!(short
        .next_block()
        .expect_err("short PVAR should fail when the missing row is reached")
        .to_string()
        .contains("fewer"));

    write_text(
        &pvar,
        "#CHROM POS ID REF ALT\n1 10 rs1 A G\n1 20 rs2 C T\n1 30 rs3 G A\n1 40 extra A G\n",
    );
    let mut long = BlockReader::open(BlockSource::Plink2 { pgen, pvar, psam }, options, 3)
        .expect("long-PVAR session should open");
    assert!(long
        .next_block()
        .expect_err("long PVAR should fail at terminal count validation")
        .to_string()
        .contains("exceeds"));
}

#[test]
fn pbr_rust_plink2_002_invalid_variable_width_header_table_fails_during_open() {
    let (_dir, pgen, pvar, psam) =
        write_bad_variable_width_block_offset_fixture("pbr-plink2-session-bad-header-table");
    let error = BlockReader::open(
        BlockSource::Plink2 { pgen, pvar, psam },
        plink2_block_options(false, None, None, DenseMissingPolicy::Nan, true),
        1,
    )
    .expect_err("invalid variable-width header table should fail during session open");

    assert!(error.to_string().contains("block offset"));
}

#[test]
fn pbr_rust_plink2_002_ld_record_at_variant_block_start_fails_during_header_parse() {
    const VARIANT_BLOCK_SIZE: usize = 65_536;

    for compression in [2_u8, 3] {
        let dir = unique_dir(&format!(
            "pbr-plink2-invalid-ld-{compression}-at-variant-block-start"
        ));
        let mut record_types = vec![0_u8; VARIANT_BLOCK_SIZE + 1];
        record_types[VARIANT_BLOCK_SIZE] = compression;
        let records = vec![&[0_u8][..]; VARIANT_BLOCK_SIZE + 1];
        let pgen_bytes = variable_width_pgen(&record_types, &records, 3);
        let (pgen, pvar, psam) =
            write_plink2_fixture_with_variants(&dir, &pgen_bytes, "#CHROM POS ID REF ALT\n");

        let session_error = BlockReader::open(
            BlockSource::Plink2 {
                pgen: pgen.clone(),
                pvar: pvar.clone(),
                psam: psam.clone(),
            },
            plink2_block_options(false, None, None, DenseMissingPolicy::Nan, true),
            1,
        )
        .expect_err("LD at the first record of a later PGEN variant block must fail during open");
        assert!(matches!(
            session_error,
            genoio_core::GenoioError::InvalidSource { .. }
        ));
        assert!(session_error
            .to_string()
            .contains("first record of a variant block"));

        let prefix_error = ::genoio_io::read_plink2_dense_windowed(
            &pgen,
            &pvar,
            &psam,
            None,
            None,
            Some(VariantWindow {
                start: VARIANT_BLOCK_SIZE,
                len: 1,
            }),
            DenseMissingPolicy::Nan,
            false,
            false,
        )
        .expect_err("prefix header parsing must enforce the same variant-block LD boundary");
        assert!(matches!(
            prefix_error,
            genoio_core::GenoioError::InvalidSource { .. }
        ));
        assert!(prefix_error
            .to_string()
            .contains("first record of a variant block"));
    }
}

#[test]
fn pbr_rust_plink2_003_genotype_and_haplotype_dosage_blocks_match_whole_reads() {
    let dir = unique_dir("pbr-plink2-dosage-block-parity");
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &explicit_phased_dosage_ld_pgen(),
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
",
    );

    let expected_genotype =
        genoio_io::read_plink2_dosage_dense_windowed(&pgen, &pvar, &psam, None, None, None, false)
            .expect("whole genotype dosage read should decode");
    let genotype_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_mode_block_options(
            MatrixKind::Genotype,
            false,
            None,
            None,
            DosageSource::Dosage,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("genotype dosage session returned a sparse block"),
    })
    .collect::<Vec<_>>();
    assert_eq!(
        genotype_blocks
            .iter()
            .map(|block| block.n_variants)
            .collect::<Vec<_>>(),
        vec![1, 1]
    );
    assert_values_with_nan(
        &concatenate_plink2_dense_blocks(&genotype_blocks),
        &dense_values_sample_major(&expected_genotype),
    );

    let expected_haplotype = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen, &pvar, &psam, None, None, None, false,
    )
    .expect("whole haplotype dosage read should decode");
    let haplotype_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_mode_block_options(
            MatrixKind::Haplotype,
            false,
            None,
            None,
            DosageSource::Dosage,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("haplotype dosage session returned a sparse block"),
    })
    .collect::<Vec<_>>();
    assert_values_with_nan(
        &concatenate_plink2_dense_blocks(&haplotype_blocks),
        &dense_values_sample_major(&expected_haplotype),
    );
    assert_eq!(
        haplotype_blocks
            .iter()
            .flat_map(|block| variant_ids(variants(&block.variants)))
            .collect::<Vec<_>>(),
        vec!["rs1", "rs2"]
    );
}

fn assert_unselected_dosage_corruption_is_rejected(
    name: &str,
    invalid_pgen_bytes: &[u8],
    valid_pgen_bytes: &[u8],
    matrix_kind: MatrixKind,
    expected_valid_values: &[f32],
) {
    let variants = "#CHROM POS ID REF ALT\n1 10 valid A G\n1 20 invalid C T\n";
    let selected = vec!["S1".to_owned()];
    let invalid_dir = unique_dir(&format!("pbr-plink2-unselected-invalid-{name}"));
    let (invalid_pgen, invalid_pvar, invalid_psam) =
        write_plink2_fixture_with_variants(&invalid_dir, invalid_pgen_bytes, variants);

    let stateless_error = match matrix_kind {
        MatrixKind::Genotype => genoio_io::read_plink2_dosage_dense_windowed(
            &invalid_pgen,
            &invalid_pvar,
            &invalid_psam,
            Some(&selected),
            None,
            None,
            false,
        ),
        MatrixKind::Haplotype => genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
            &invalid_pgen,
            &invalid_pvar,
            &invalid_psam,
            Some(&selected),
            None,
            None,
            false,
        ),
    }
    .expect_err("an invalid unselected dosage entry must fail a stateless read");
    assert!(
        matches!(
            stateless_error,
            genoio_core::GenoioError::InvalidSource { .. }
        ),
        "{name}: {stateless_error}"
    );

    let mut invalid_reader = BlockReader::open(
        BlockSource::Plink2 {
            pgen: invalid_pgen,
            pvar: invalid_pvar,
            psam: invalid_psam,
        },
        plink2_mode_block_options(
            matrix_kind,
            false,
            Some(selected.clone()),
            None,
            DosageSource::Dosage,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .expect("an invalid-later dosage session should open");
    assert!(invalid_reader
        .next_block()
        .expect("the valid first dosage block should decode")
        .is_some());
    let session_error = invalid_reader
        .next_block()
        .expect_err("an invalid unselected dosage entry must fail when its block is requested");
    assert!(
        matches!(
            session_error,
            genoio_core::GenoioError::InvalidSource { .. }
        ),
        "{name}: {session_error}"
    );
    assert!(invalid_reader
        .next_block()
        .expect("a failed dosage session must remain terminal")
        .is_none());

    let valid_dir = unique_dir(&format!("pbr-plink2-unselected-valid-{name}"));
    let (valid_pgen, valid_pvar, valid_psam) =
        write_plink2_fixture_with_variants(&valid_dir, valid_pgen_bytes, variants);
    let stateless_valid = match matrix_kind {
        MatrixKind::Genotype => genoio_io::read_plink2_dosage_dense_windowed(
            &valid_pgen,
            &valid_pvar,
            &valid_psam,
            Some(&selected),
            None,
            None,
            false,
        ),
        MatrixKind::Haplotype => genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
            &valid_pgen,
            &valid_pvar,
            &valid_psam,
            Some(&selected),
            None,
            None,
            false,
        ),
    }
    .expect("the equivalent valid selected-sample read should decode");
    assert_values_with_nan(
        &dense_values_sample_major(&stateless_valid),
        expected_valid_values,
    );

    let valid_blocks = collect_plink2_blocks(
        &valid_pgen,
        &valid_pvar,
        &valid_psam,
        plink2_mode_block_options(
            matrix_kind,
            false,
            Some(selected),
            None,
            DosageSource::Dosage,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("a dosage session returned sparse output"),
    })
    .collect::<Vec<_>>();
    assert_values_with_nan(
        &concatenate_plink2_dense_blocks(&valid_blocks),
        expected_valid_values,
    );
}

#[test]
fn pbr_rust_plink2_003_fixed_genotype_dosage_validates_unselected_samples() {
    let valid_record = {
        let mut record = vec![0x00];
        record.extend(0_u16.to_le_bytes());
        record.extend(0_u16.to_le_bytes());
        record.extend(0_u16.to_le_bytes());
        record
    };
    let invalid_record = {
        let mut record = valid_record.clone();
        record[3..5].copy_from_slice(&32_769_u16.to_le_bytes());
        record
    };
    assert_unselected_dosage_corruption_is_rejected(
        "fixed-genotype",
        &fixed_width_dosage_pgen(&[&valid_record, &invalid_record], 3),
        &fixed_width_dosage_pgen(&[&valid_record, &valid_record], 3),
        MatrixKind::Genotype,
        &[0.0, 0.0],
    );
}

#[test]
fn pbr_rust_plink2_003_fixed_phased_dosage_validates_unselected_samples() {
    let valid_record = {
        let mut record = vec![0x00];
        for _ in 0..3 {
            record.extend(0_u16.to_le_bytes());
        }
        for _ in 0..3 {
            record.extend(0_i16.to_le_bytes());
        }
        record
    };
    let invalid_record = {
        let mut record = valid_record.clone();
        record[3..5].copy_from_slice(&u16::MAX.to_le_bytes());
        record
    };
    assert_unselected_dosage_corruption_is_rejected(
        "fixed-phased",
        &fixed_width_phased_dosage_pgen(&[&valid_record, &invalid_record], 3),
        &fixed_width_phased_dosage_pgen(&[&valid_record, &valid_record], 3),
        MatrixKind::Haplotype,
        &[0.0, 0.0, 0.0, 0.0],
    );
}

#[test]
fn pbr_rust_plink2_003_variable_explicit_phase_validates_unselected_samples() {
    let valid_record = {
        let mut record = vec![0x00];
        for _ in 0..3 {
            record.extend(0_u16.to_le_bytes());
        }
        for _ in 0..3 {
            record.extend(0_i16.to_le_bytes());
        }
        record
    };
    let invalid_record = {
        let mut record = valid_record.clone();
        record[9..11].copy_from_slice(&16_384_i16.to_le_bytes());
        record
    };
    assert_unselected_dosage_corruption_is_rejected(
        "variable-explicit-phase",
        &explicit_phased_dosage_pgen(&[&valid_record, &invalid_record], 3),
        &explicit_phased_dosage_pgen(&[&valid_record, &valid_record], 3),
        MatrixKind::Haplotype,
        &[0.0, 0.0, 0.0, 0.0],
    );
}

#[test]
fn pbr_rust_plink2_003_sparse_explicit_phase_validates_unselected_samples() {
    let valid_record = {
        let mut record = vec![0x00, 0x00, 0b0000_0010];
        record.extend(0_u16.to_le_bytes());
        record.push(0b0000_0001);
        record.extend(0_i16.to_le_bytes());
        record
    };
    let invalid_record = {
        let mut record = valid_record.clone();
        record[6..8].copy_from_slice(&16_384_i16.to_le_bytes());
        record
    };
    assert_unselected_dosage_corruption_is_rejected(
        "sparse-explicit-phase",
        &variable_width_pgen(&[0xf0, 0xf0], &[&valid_record, &invalid_record], 3),
        &variable_width_pgen(&[0xf0, 0xf0], &[&valid_record, &valid_record], 3),
        MatrixKind::Haplotype,
        &[0.0, 0.0, 0.0, 0.0],
    );
}

#[test]
fn pbr_rust_plink2_003_variable_dosage_existence_and_phase_mapping_follow_the_spec() {
    let dir = unique_dir("pbr-plink2-variable-dosage-existence-phase");
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &variable_phased_dosage_existence_pgen(),
        "\
#CHROM POS ID REF ALT
1 10 implicit_type1 A G
1 20 explicit_type3 C T
1 30 full_type2 G A
",
    );
    let selected = vec!["S3".to_owned(), "S2".to_owned()];

    let expected_genotype = [1.25, 1.0, f32::NAN, 2.0, 2.0, 2.0];
    let stateless_genotype = genoio_io::read_plink2_dosage_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&selected),
        None,
        None,
        false,
    )
    .expect("stateless dosage decode should support existence types 1, 2, and 3");
    assert_values_with_nan(
        &dense_values_sample_major(&stateless_genotype),
        &expected_genotype,
    );
    let genotype_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_mode_block_options(
            MatrixKind::Genotype,
            false,
            Some(selected.clone()),
            None,
            DosageSource::Dosage,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("genotype dosage session returned sparse output"),
    })
    .collect::<Vec<_>>();
    assert_values_with_nan(
        &concatenate_plink2_dense_blocks(&genotype_blocks),
        &expected_genotype,
    );

    let expected_haplotype = [
        0.25,
        0.25,
        f32::NAN,
        1.0,
        0.75,
        f32::NAN,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    ];
    let stateless_haplotype = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&selected),
        None,
        None,
        false,
    )
    .expect("stateless haplotype dosage decode should support implicit and mapped phase");
    assert_values_with_nan(
        &dense_values_sample_major(&stateless_haplotype),
        &expected_haplotype,
    );
    let haplotype_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_mode_block_options(
            MatrixKind::Haplotype,
            false,
            Some(selected),
            None,
            DosageSource::Dosage,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("haplotype dosage session returned sparse output"),
    })
    .collect::<Vec<_>>();
    assert_values_with_nan(
        &concatenate_plink2_dense_blocks(&haplotype_blocks),
        &expected_haplotype,
    );
    assert_eq!(
        haplotype_blocks
            .iter()
            .flat_map(|block| variant_ids(variants(&block.variants)))
            .collect::<Vec<_>>(),
        vec!["implicit_type1", "explicit_type3", "full_type2"]
    );

    let complete_filter = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "missing_rate",
        "params": {"max": 0.0}
    }))
    .expect("collapsed-dosage filter should parse");
    let filtered = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&["S2".to_owned(), "S3".to_owned()]),
        Some(&complete_filter),
        None,
        false,
    )
    .expect("stateless implicit/mapped phased dosages should support genotype-stat filters");
    assert_eq!(
        variant_ids(variants(&filtered.variants)),
        vec!["implicit_type1", "explicit_type3"]
    );
    assert_eq!(filtered.diagnostics.dropped_genotype_variants, 1);
}

#[test]
fn pbr_rust_plink2_003_mapped_phase_errors_are_delayed_and_terminal() {
    let dir = unique_dir("pbr-plink2-delayed-mapped-phase-error");
    let mut valid = vec![0x24, 0x00, 1, 1];
    valid.extend(scaled_dosage(1.25));
    let mut malformed = vec![0x24, 0x00, 0b0000_0010];
    malformed.extend(scaled_dosage(1.0));
    malformed.push(0b0000_0001);
    let pgen_bytes = variable_width_pgen(&[0x30, 0xf0], &[&valid, &malformed], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "#CHROM POS ID REF ALT\n1 10 valid A G\n1 20 malformed C T\n",
    );
    let mut reader = BlockReader::open(
        BlockSource::Plink2 { pgen, pvar, psam },
        plink2_mode_block_options(
            MatrixKind::Haplotype,
            false,
            None,
            None,
            DosageSource::Dosage,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .expect("mapped-phase delayed-error session should open");

    assert!(reader
        .next_block()
        .expect("first implicit-phase block should decode")
        .is_some());
    assert!(matches!(
        reader
            .next_block()
            .expect_err("malformed later mapped phase must fail only when reached"),
        genoio_core::GenoioError::InvalidSource { .. }
    ));
    assert!(reader
        .next_block()
        .expect("failed mapped-phase session should remain terminal")
        .is_none());
}

#[test]
fn pbr_rust_plink2_003_invalid_raw_dosage_is_rejected_before_filtering_and_tail_work() {
    let valid_record = {
        let mut record = vec![0x24];
        record.extend(0_u16.to_le_bytes());
        record.extend(16_384_u16.to_le_bytes());
        record.extend(32_768_u16.to_le_bytes());
        record
    };
    let invalid_record = {
        let mut record = vec![0x24];
        record.extend(0_u16.to_le_bytes());
        record.extend(32_769_u16.to_le_bytes());
        record.extend(32_768_u16.to_le_bytes());
        record
    };
    let tail_record = {
        let mut record = vec![0x24];
        record.extend(0_u16.to_le_bytes());
        record.extend(8_192_u16.to_le_bytes());
        record.extend(32_768_u16.to_le_bytes());
        record
    };
    let variant_text = "#CHROM POS ID REF ALT\n1 10 valid A G\n1 20 invalid C T\n1 30 tail G A\n";
    let genotype_filter = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "missing_rate",
        "params": {"max": 1.0}
    }))
    .expect("genotype-stat filter should parse");

    for (name, pgen_bytes) in [
        (
            "variable",
            variable_width_pgen(
                &[0x40, 0x40, 0x40],
                &[&valid_record, &invalid_record, &tail_record],
                3,
            ),
        ),
        (
            "fixed",
            fixed_width_dosage_pgen(&[&valid_record, &invalid_record, &tail_record], 3),
        ),
    ] {
        let dir = unique_dir(&format!("pbr-plink2-invalid-raw-{name}"));
        let (pgen, pvar, psam) =
            write_plink2_fixture_with_variants(&dir, &pgen_bytes, variant_text);

        let stateless_error = genoio_io::read_plink2_dosage_dense_windowed(
            &pgen,
            &pvar,
            &psam,
            None,
            Some(&genotype_filter),
            None,
            false,
        )
        .expect_err("reserved raw dosage must fail before genotype-stat filtering");
        assert!(
            matches!(
                stateless_error,
                genoio_core::GenoioError::InvalidSource { .. }
            ),
            "{name}: {stateless_error}"
        );

        let mut reader = BlockReader::open(
            BlockSource::Plink2 {
                pgen: pgen.clone(),
                pvar: pvar.clone(),
                psam: psam.clone(),
            },
            plink2_mode_block_options(
                MatrixKind::Genotype,
                false,
                None,
                Some(genotype_filter.clone()),
                DosageSource::Dosage,
                DenseMissingPolicy::Nan,
                true,
            ),
            1,
        )
        .expect("invalid-raw delayed-error session should open");
        assert!(reader
            .next_block()
            .expect("valid first dosage block should decode")
            .is_some());
        let error = reader
            .next_block()
            .expect_err("reserved later raw dosage must fail only when reached");
        assert!(
            matches!(error, genoio_core::GenoioError::InvalidSource { .. }),
            "{name}: {error}"
        );
        assert!(reader
            .next_block()
            .expect("failed invalid-raw session should remain terminal")
            .is_none());

        let reject_invalid_metadata = VariantFilter::from_json_value(serde_json::json!({
            "op": "not",
            "expr": {
                "op": "predicate",
                "name": "id_in",
                "params": {"values": ["invalid"]}
            }
        }))
        .expect("metadata filter should parse");
        let blocks = collect_plink2_blocks(
            &pgen,
            &pvar,
            &psam,
            plink2_mode_block_options(
                MatrixKind::Genotype,
                false,
                None,
                Some(reject_invalid_metadata),
                DosageSource::Dosage,
                DenseMissingPolicy::Nan,
                true,
            ),
            1,
        );
        assert_eq!(
            blocks
                .iter()
                .flat_map(|block| match block {
                    BlockOutput::Dense(matrix) => variant_ids(variants(&matrix.variants)),
                    BlockOutput::Sparse(_) => Vec::new(),
                })
                .collect::<Vec<_>>(),
            vec!["valid", "tail"],
            "{name}"
        );
    }
}

#[test]
fn pbr_rust_plink2_003_fixed_width_and_rejected_ld_base_dosage_semantics_are_persistent() {
    let dir = unique_dir("pbr-plink2-fixed-dosage-block-parity");
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
    let mut record_3 = vec![0x30];
    record_3.extend(scaled_dosage(0.4));
    record_3.extend(scaled_dosage(1.0));
    record_3.extend(u16::MAX.to_le_bytes());
    record_3.extend(scaled_phase_delta(0.1, 0.3));
    record_3.extend(scaled_phase_delta(0.5, 0.5));
    record_3.extend(i16::MIN.to_le_bytes());
    let fixed = fixed_width_phased_dosage_pgen(&[&record_1, &record_2, &record_3], 3);
    let (fixed_pgen, fixed_pvar, fixed_psam) = write_plink2_fixture_with_variants(
        &dir,
        &fixed,
        "\
#CHROM POS ID REF ALT QUAL
1 10 rs1 A G 30
1 20 rs2 C T 40
2 30 rs3 G A 50
",
    );
    let requested_samples = vec!["S3".to_owned(), "S1".to_owned()];
    for matrix_kind in [MatrixKind::Genotype, MatrixKind::Haplotype] {
        let expected = match matrix_kind {
            MatrixKind::Genotype => genoio_io::read_plink2_dosage_dense_windowed(
                &fixed_pgen,
                &fixed_pvar,
                &fixed_psam,
                Some(&requested_samples),
                None,
                None,
                false,
            ),
            MatrixKind::Haplotype => genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
                &fixed_pgen,
                &fixed_pvar,
                &fixed_psam,
                Some(&requested_samples),
                None,
                None,
                false,
            ),
        }
        .expect("fixed-width phased-dosage whole read should decode");
        let blocks = collect_plink2_blocks(
            &fixed_pgen,
            &fixed_pvar,
            &fixed_psam,
            plink2_mode_block_options(
                matrix_kind,
                false,
                Some(requested_samples.clone()),
                None,
                DosageSource::Dosage,
                DenseMissingPolicy::Nan,
                true,
            ),
            2,
        )
        .into_iter()
        .map(|block| match block {
            BlockOutput::Dense(matrix) => matrix,
            BlockOutput::Sparse(_) => panic!("fixed-width dosage session returned sparse output"),
        })
        .collect::<Vec<_>>();

        assert_eq!(
            blocks
                .iter()
                .map(|block| block.n_variants)
                .collect::<Vec<_>>(),
            vec![2, 1]
        );
        let actual_values = concatenate_plink2_dense_blocks(&blocks);
        assert_values_with_nan(&actual_values, &dense_values_sample_major(&expected));
        assert_eq!(
            actual_values
                .iter()
                .map(|value| value.is_nan())
                .collect::<Vec<_>>(),
            dense_missing_sample_major(&expected)
        );
        assert!(blocks.iter().all(|block| block.samples == expected.samples));
        assert!(blocks.iter().all(|block| {
            block
                .variants
                .as_ref()
                .is_some_and(|variants| variants.len() == block.n_variants)
        }));
        assert_eq!(
            blocks
                .iter()
                .flat_map(|block| variant_ids(variants(&block.variants)))
                .collect::<Vec<_>>(),
            variant_ids(variants(&expected.variants))
        );
        assert_eq!(
            blocks
                .iter()
                .flat_map(|block| {
                    let metadata = variants(&block.variants);
                    (0..metadata.len())
                        .map(|index| variant_a0(metadata, index))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            (0..variants(&expected.variants).len())
                .map(|index| variant_a0(variants(&expected.variants), index))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            blocks
                .iter()
                .flat_map(|block| {
                    let metadata = variants(&block.variants);
                    (0..metadata.len())
                        .map(|index| variant_a1(metadata, index))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            (0..variants(&expected.variants).len())
                .map(|index| variant_a1(variants(&expected.variants), index))
                .collect::<Vec<_>>()
        );
        assert_eq!(blocks[0].diagnostics.candidate_variants, 2);
        assert_eq!(blocks[0].diagnostics.retained_variants, 2);
        let mut expected_final_diagnostics = expected.diagnostics.clone();
        expected_final_diagnostics.retained_variants = 1;
        assert_eq!(blocks[1].diagnostics, expected_final_diagnostics);
    }

    let rejected_dir = unique_dir("pbr-plink2-rejected-dosage-ld-base");
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &rejected_dir,
        &explicit_phased_dosage_ld_pgen(),
        "#CHROM POS ID REF ALT\n2 10 rejected_base A G\n1 20 retained_ld C T\n",
    );
    let filter = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "1"}
    }))
    .expect("metadata filter should parse");

    for matrix_kind in [MatrixKind::Genotype, MatrixKind::Haplotype] {
        let blocks = collect_plink2_blocks(
            &pgen,
            &pvar,
            &psam,
            plink2_mode_block_options(
                matrix_kind,
                false,
                None,
                Some(filter.clone()),
                DosageSource::Dosage,
                DenseMissingPolicy::Nan,
                true,
            ),
            1,
        );
        assert_eq!(blocks.len(), 1);
        let BlockOutput::Dense(block) = &blocks[0] else {
            panic!("dosage session returned a sparse block");
        };
        assert_eq!(variant_ids(variants(&block.variants)), vec!["retained_ld"]);
    }
}

#[test]
fn pbr_rust_plink2_003_invalid_later_dosage_track_is_not_prefetched() {
    let dir = unique_dir("pbr-plink2-delayed-dosage-error");
    let mut valid = vec![0x24];
    valid.extend(scaled_dosage(0.2));
    valid.extend(scaled_dosage(1.4));
    valid.extend(scaled_dosage(1.8));
    let malformed = [0x24];
    let pgen_bytes = variable_width_pgen(&[0x40, 0x40], &[&valid, &malformed], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "#CHROM POS ID REF ALT\n1 10 rs1 A G\n1 20 rs2 C T\n",
    );
    let mut reader = BlockReader::open(
        BlockSource::Plink2 { pgen, pvar, psam },
        plink2_mode_block_options(
            MatrixKind::Genotype,
            false,
            None,
            None,
            DosageSource::Dosage,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .expect("delayed dosage-error session should open");

    assert!(reader
        .next_block()
        .expect("first valid dosage block should decode")
        .is_some());
    assert!(reader
        .next_block()
        .expect_err("malformed later dosage track should fail only when reached")
        .to_string()
        .contains("dosage"));
}

#[test]
fn pbr_rust_plink2_004_hardcall_haplotype_dense_and_sparse_blocks_match_whole_reads() {
    let dir = unique_dir("pbr-plink2-hardcall-haplotype-block-parity");
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &explicit_phased_hardcall_ld_pgen(),
        "#CHROM POS ID REF ALT\n1 10 rs1 A G\n1 20 rs2 C T\n",
    );
    let keep = vec!["S1".to_owned(), "S2".to_owned()];

    let expected_dense = genoio_io::read_plink2_haplotypes_dense_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        None,
        false,
    )
    .expect("whole hardcall haplotype read should decode");
    let dense_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_mode_block_options(
            MatrixKind::Haplotype,
            false,
            Some(keep.clone()),
            None,
            DosageSource::Hardcall,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Dense(matrix) => matrix,
        BlockOutput::Sparse(_) => panic!("hardcall haplotype session returned a sparse block"),
    })
    .collect::<Vec<_>>();
    assert_values_with_nan(
        &concatenate_plink2_dense_blocks(&dense_blocks),
        &dense_values_sample_major(&expected_dense),
    );
    assert!(dense_blocks
        .iter()
        .all(|block| block.samples == expected_dense.samples));

    let expected_sparse = genoio_io::read_plink2_haplotypes_sparse_windowed(
        &pgen,
        &pvar,
        &psam,
        Some(&keep),
        None,
        None,
    )
    .expect("whole sparse hardcall haplotype read should decode");
    let sparse_blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_mode_block_options(
            MatrixKind::Haplotype,
            true,
            Some(keep),
            None,
            DosageSource::Hardcall,
            DenseMissingPolicy::Raise,
            true,
        ),
        1,
    )
    .into_iter()
    .map(|block| match block {
        BlockOutput::Sparse(matrix) => matrix,
        BlockOutput::Dense(_) => panic!("sparse haplotype session returned a dense block"),
    })
    .collect::<Vec<_>>();
    assert_eq!(
        concatenate_plink2_sparse_blocks(&sparse_blocks),
        csc_to_dense(&expected_sparse)
    );
}

#[test]
fn pbr_rust_plink2_004_hardcall_haplotype_retains_rejected_ld_base_without_prefetch() {
    let dir = unique_dir("pbr-plink2-hardcall-haplotype-rejected-ld-base");
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &explicit_phased_hardcall_ld_pgen(),
        "#CHROM POS ID REF ALT\n2 10 rejected_base A G\n1 20 retained_ld C T\n",
    );
    let filter = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "1"}
    }))
    .expect("metadata filter should parse");
    let blocks = collect_plink2_blocks(
        &pgen,
        &pvar,
        &psam,
        plink2_mode_block_options(
            MatrixKind::Haplotype,
            false,
            None,
            Some(filter),
            DosageSource::Hardcall,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    );

    assert_eq!(blocks.len(), 1);
    let BlockOutput::Dense(block) = &blocks[0] else {
        panic!("hardcall haplotype session returned a sparse block");
    };
    assert_eq!(variant_ids(variants(&block.variants)), vec!["retained_ld"]);
}

#[test]
fn pbr_rust_plink2_004_hardcall_haplotype_errors_remain_record_lazy() {
    let dir = unique_dir("pbr-plink2-hardcall-haplotype-delayed-error");
    let pgen_bytes = variable_width_pgen(&[0x10, 0x00], &[&[0x21, 0x00], &[0x00]], 3);
    let (pgen, pvar, psam) = write_plink2_fixture_with_variants(
        &dir,
        &pgen_bytes,
        "#CHROM POS ID REF ALT\n1 10 valid A G\n1 20 unphased C T\n",
    );
    let mut reader = BlockReader::open(
        BlockSource::Plink2 { pgen, pvar, psam },
        plink2_mode_block_options(
            MatrixKind::Haplotype,
            false,
            None,
            None,
            DosageSource::Hardcall,
            DenseMissingPolicy::Nan,
            true,
        ),
        1,
    )
    .expect("delayed hardcall haplotype-error session should open");

    assert!(reader
        .next_block()
        .expect("first valid hardcall haplotype block should decode")
        .is_some());
    assert!(reader
        .next_block()
        .expect_err("unphased later record should fail only when reached")
        .to_string()
        .contains("unphased"));
}

#[test]
fn pbr_rust_plink2_004_complete_mode_matrix_uses_real_session_branches() {
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Storage {
        FixedHardcall,
        FixedDosage,
        FixedPhasedDosage,
        VariableWidth,
    }

    struct Fixture {
        storage: Storage,
        hardcall: (PathBuf, PathBuf, PathBuf),
        dosage: (PathBuf, PathBuf, PathBuf),
    }

    let variants = "#CHROM POS ID REF ALT QUAL\n1 10 rs1 A G 30\n";
    let fixed_hardcall_dir = unique_dir("pbr-plink2-mode-matrix-fixed-hardcall");
    let fixed_hardcall = write_plink2_fixture_with_variants(
        &fixed_hardcall_dir,
        &fixed_width_pgen(&[0x21], 3, 1),
        variants,
    );

    let mut dosage_record = vec![0x21];
    dosage_record.extend(scaled_dosage(1.0));
    dosage_record.extend(scaled_dosage(0.5));
    dosage_record.extend(scaled_dosage(2.0));
    let fixed_dosage_dir = unique_dir("pbr-plink2-mode-matrix-fixed-dosage");
    let fixed_dosage = write_plink2_fixture_with_variants(
        &fixed_dosage_dir,
        &fixed_width_dosage_pgen(&[&dosage_record], 3),
        variants,
    );

    let mut phased_dosage_record = dosage_record.clone();
    phased_dosage_record.extend(scaled_phase_delta(0.25, 0.75));
    phased_dosage_record.extend(scaled_phase_delta(0.0, 0.5));
    phased_dosage_record.extend(scaled_phase_delta(1.0, 1.0));
    let fixed_phased_dir = unique_dir("pbr-plink2-mode-matrix-fixed-phased-dosage");
    let fixed_phased = write_plink2_fixture_with_variants(
        &fixed_phased_dir,
        &fixed_width_phased_dosage_pgen(&[&phased_dosage_record], 3),
        variants,
    );

    let variable_hardcall_dir = unique_dir("pbr-plink2-mode-matrix-variable-hardcall");
    let variable_hardcall = write_plink2_fixture_with_variants(
        &variable_hardcall_dir,
        &explicit_phased_hardcall_pgen(&[&[0x21, 0x00]], 3),
        variants,
    );
    let variable_dosage_dir = unique_dir("pbr-plink2-mode-matrix-variable-dosage");
    let variable_dosage = write_plink2_fixture_with_variants(
        &variable_dosage_dir,
        &explicit_phased_dosage_pgen(&[&phased_dosage_record], 3),
        variants,
    );

    let fixtures = [
        Fixture {
            storage: Storage::FixedHardcall,
            hardcall: fixed_hardcall.clone(),
            dosage: fixed_hardcall,
        },
        Fixture {
            storage: Storage::FixedDosage,
            hardcall: fixed_dosage.clone(),
            dosage: fixed_dosage,
        },
        Fixture {
            storage: Storage::FixedPhasedDosage,
            hardcall: fixed_phased.clone(),
            dosage: fixed_phased,
        },
        Fixture {
            storage: Storage::VariableWidth,
            hardcall: variable_hardcall,
            dosage: variable_dosage,
        },
    ];

    let mut case_count = 0_usize;
    for fixture in fixtures {
        for matrix_kind in [MatrixKind::Genotype, MatrixKind::Haplotype] {
            for dosage_source in [DosageSource::Hardcall, DosageSource::Dosage] {
                for sparse in [false, true] {
                    case_count += 1;
                    let (pgen, pvar, psam) = match dosage_source {
                        DosageSource::Hardcall => &fixture.hardcall,
                        DosageSource::Dosage => &fixture.dosage,
                    };
                    let options = plink2_mode_block_options(
                        matrix_kind,
                        sparse,
                        None,
                        None,
                        dosage_source,
                        if sparse {
                            DenseMissingPolicy::Raise
                        } else {
                            DenseMissingPolicy::Nan
                        },
                        true,
                    );
                    let case = format!(
                        "{:?}/{matrix_kind:?}/{dosage_source:?}/{}",
                        fixture.storage,
                        if sparse { "sparse" } else { "dense" }
                    );

                    if sparse && dosage_source == DosageSource::Dosage {
                        let error = BlockReader::open(
                            BlockSource::Plink2 {
                                pgen: pgen.clone(),
                                pvar: pvar.clone(),
                                psam: psam.clone(),
                            },
                            options,
                            1,
                        )
                        .expect_err(&format!("{case} must fail at construction"));
                        assert!(
                            matches!(
                                error,
                                genoio_core::GenoioError::UnsupportedRepresentation { .. }
                            ),
                            "{case}: unexpected construction error {error}"
                        );
                        continue;
                    }

                    let supported = match (matrix_kind, dosage_source) {
                        (MatrixKind::Genotype, DosageSource::Hardcall) => true,
                        (MatrixKind::Haplotype, DosageSource::Hardcall) => {
                            fixture.storage == Storage::VariableWidth
                        }
                        (MatrixKind::Genotype, DosageSource::Dosage) => {
                            fixture.storage != Storage::FixedHardcall
                        }
                        (MatrixKind::Haplotype, DosageSource::Dosage) => matches!(
                            fixture.storage,
                            Storage::FixedPhasedDosage | Storage::VariableWidth
                        ),
                    };
                    let mut reader = BlockReader::open(
                        BlockSource::Plink2 {
                            pgen: pgen.clone(),
                            pvar: pvar.clone(),
                            psam: psam.clone(),
                        },
                        options,
                        1,
                    )
                    .unwrap_or_else(|error| {
                        panic!("{case} should pass session construction: {error}")
                    });

                    if !supported {
                        let error = reader
                            .next_block()
                            .expect_err(&format!("{case} must fail at the lazy record boundary"));
                        assert!(
                            matches!(
                                error,
                                genoio_core::GenoioError::UnsupportedRepresentation { .. }
                            ),
                            "{case}: unexpected lazy error {error}"
                        );
                        continue;
                    }

                    let actual = reader
                        .next_block()
                        .unwrap_or_else(|error| panic!("{case} should decode: {error}"))
                        .unwrap_or_else(|| panic!("{case} should yield one block"));
                    match (matrix_kind, dosage_source, sparse, actual) {
                        (
                            MatrixKind::Genotype,
                            DosageSource::Hardcall,
                            false,
                            BlockOutput::Dense(actual),
                        ) => {
                            let expected =
                                genoio_io::read_plink2_dense(pgen, pvar, psam, None, None)
                                    .expect("stateless genotype hardcall oracle should decode");
                            assert_dense_matrix_matches_oracle(&actual, &expected);
                        }
                        (
                            MatrixKind::Genotype,
                            DosageSource::Hardcall,
                            true,
                            BlockOutput::Sparse(actual),
                        ) => {
                            let expected =
                                genoio_io::read_plink2_sparse(pgen, pvar, psam, None, None)
                                    .expect("stateless sparse genotype oracle should decode");
                            assert_eq!(actual, expected, "{case}");
                        }
                        (
                            MatrixKind::Haplotype,
                            DosageSource::Hardcall,
                            false,
                            BlockOutput::Dense(actual),
                        ) => {
                            let expected = genoio_io::read_plink2_haplotypes_dense_windowed(
                                pgen, pvar, psam, None, None, None, false,
                            )
                            .expect("stateless hardcall haplotype oracle should decode");
                            assert_dense_matrix_matches_oracle(&actual, &expected);
                        }
                        (
                            MatrixKind::Haplotype,
                            DosageSource::Hardcall,
                            true,
                            BlockOutput::Sparse(actual),
                        ) => {
                            let expected = genoio_io::read_plink2_haplotypes_sparse_windowed(
                                pgen, pvar, psam, None, None, None,
                            )
                            .expect("stateless sparse haplotype oracle should decode");
                            assert_eq!(actual, expected, "{case}");
                        }
                        (
                            MatrixKind::Genotype,
                            DosageSource::Dosage,
                            false,
                            BlockOutput::Dense(actual),
                        ) => {
                            let expected = genoio_io::read_plink2_dosage_dense_windowed(
                                pgen, pvar, psam, None, None, None, false,
                            )
                            .expect("stateless genotype dosage oracle should decode");
                            assert_dense_matrix_matches_oracle(&actual, &expected);
                        }
                        (
                            MatrixKind::Haplotype,
                            DosageSource::Dosage,
                            false,
                            BlockOutput::Dense(actual),
                        ) => {
                            let expected = genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
                                pgen, pvar, psam, None, None, None, false,
                            )
                            .expect("stateless haplotype dosage oracle should decode");
                            assert_dense_matrix_matches_oracle(&actual, &expected);
                        }
                        other => panic!("{case}: unexpected output branch {other:?}"),
                    }
                    assert!(
                        reader
                            .next_block()
                            .unwrap_or_else(|error| panic!("{case} EOF should be clean: {error}"))
                            .is_none(),
                        "{case}: expected one-record EOF"
                    );
                }
            }
        }
    }
    assert_eq!(case_count, 32);
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
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S2", "S3"]
    );
    let variant_metadata = variants(&dense.variants);
    assert_eq!(variant_a0(variant_metadata, 0), "A");
    assert_eq!(variant_a1(variant_metadata, 0), "G");
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
            .as_ref()
            .expect("sample metadata")
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
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S1", "S2", "S2", "S3", "S3"]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.source_sample_index)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(0), Some(1), Some(1), Some(2), Some(2)]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
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
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S1", "S3", "S3"]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.source_sample_index)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(0), Some(2), Some(2)]
    );
    assert_eq!(
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
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
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| (sample.iid.as_str(), sample.haplotype_index))
            .collect::<Vec<_>>(),
        dense
            .samples
            .as_ref()
            .expect("sample metadata")
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
            .as_ref()
            .expect("sample metadata")
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
    assert!(matrix_only.samples.is_none());
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
    assert!(matrix_only.samples.is_none());
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
    assert!(matrix_only.samples.is_none());
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
            .as_ref()
            .expect("sample metadata")
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
        genoio_io::read_plink2_metadata(&pgen, &pvar, &psam).expect("metadata should decode");

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

    let metadata =
        genoio_io::read_plink2_metadata(&pgen, &pvar_zst, &psam).expect("metadata should decode");

    assert_eq!(metadata.variants.len(), 3);
    assert_eq!(variant_id(&metadata.variants, 0), "rs1");
}
