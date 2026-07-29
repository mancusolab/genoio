// pattern: Imperative Shell

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use ::genoio_io::{
    BlockOutput, BlockReadOptions, BlockReader, BlockSource, DosageSource, MatrixKind,
};
use genoio_core::{DenseMissingPolicy, GenoioError, SparseGenotypeMatrix, VariantFilter};

mod common;

use common::plink_output as plink_io;
use common::unique_dir;
use common::vcf_output as genoio_io;
use common::vcf_output::{
    dense_values_sample_major_output, sparse_values_dense_output, variant_a0, variant_a1,
    variant_ids, variants,
};

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

fn write_three_variant_plink_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let bed = dir.join("three.bed");
    let bim = dir.join("three.bim");
    let fam = dir.join("three.fam");
    fs::write(&bed, [0x6c, 0x1b, 0x01, 0x0b, 0x2c, 0x08]).expect("bed fixture should be written");
    write_text(
        &bim,
        "\
1 rs1 0 10 G A
1 rs2 0 20 T C
2 rs3 0 30 A G
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

fn plink1_sparse_block_options(
    requested_samples: Option<Vec<String>>,
    variant_filter: Option<VariantFilter>,
) -> BlockReadOptions {
    BlockReadOptions {
        matrix_kind: MatrixKind::Genotype,
        sparse: true,
        requested_samples,
        variant_filter,
        dosage_source: DosageSource::Hardcall,
        missing_policy: DenseMissingPolicy::Raise,
        return_samples: true,
        return_variants: true,
    }
}

fn collect_plink1_sparse_blocks(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    options: BlockReadOptions,
    block_size: usize,
) -> Vec<SparseGenotypeMatrix> {
    let mut reader = BlockReader::open(
        BlockSource::Plink1 {
            bed: bed.to_path_buf(),
            bim: bim.to_path_buf(),
            fam: fam.to_path_buf(),
        },
        options,
        block_size,
    )
    .expect("persistent sparse plink1 reader should open");
    let mut blocks = Vec::new();
    while let Some(output) = reader
        .next_block()
        .expect("persistent sparse plink1 block should decode")
    {
        let BlockOutput::Sparse(matrix) = output else {
            panic!("plink1 sparse reader should return sparse blocks");
        };
        blocks.push(matrix);
    }
    blocks
}

fn concatenate_sparse_blocks(blocks: &[SparseGenotypeMatrix]) -> (Vec<i32>, Vec<i32>, Vec<f32>) {
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    for block in blocks {
        let nnz_offset =
            i32::try_from(indices.len()).expect("test sparse nonzero count should fit i32");
        indptr.extend(block.indptr.iter().skip(1).map(|pointer| {
            nnz_offset
                .checked_add(*pointer)
                .expect("test sparse pointer should fit i32")
        }));
        indices.extend_from_slice(&block.indices);
        data.extend_from_slice(&block.data);
    }
    (indptr, indices, data)
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
    assert_eq!(
        sparse_values_dense_output(&sparse),
        dense_values_sample_major_output(&dense)
    );
}

#[test]
fn pbr_rust_textvcf_002_sequential_sparse_genotype_blocks_preserve_csc_parity() {
    let dir = unique_dir("pbr-text-vcf-sparse-blocks");
    let path = dir.join("sparse.vcf");
    write_vcf(
        &path,
        "\
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1
1\t30\trs3\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1
",
    );
    let expected =
        genoio_io::read_vcf_sparse(&path, None, None).expect("whole sparse VCF should decode");
    let mut reader = BlockReader::open(
        BlockSource::Vcf { vcf: path },
        BlockReadOptions {
            matrix_kind: MatrixKind::Genotype,
            sparse: true,
            requested_samples: None,
            variant_filter: None,
            dosage_source: DosageSource::Hardcall,
            missing_policy: DenseMissingPolicy::Raise,
            return_samples: true,
            return_variants: true,
        },
        2,
    )
    .expect("persistent sparse text VCF reader should open");
    let mut blocks = Vec::new();
    while let Some(output) = reader
        .next_block()
        .expect("sparse text VCF block should decode")
    {
        let BlockOutput::Sparse(block) = output else {
            panic!("sparse text VCF session should return sparse output");
        };
        blocks.push(block);
    }
    let (indptr, indices, data) = concatenate_sparse_blocks(&blocks);

    assert_eq!(indptr, expected.indptr);
    assert_eq!(indices, expected.indices);
    assert_eq!(data, expected.data);
    assert_eq!(
        blocks.iter().map(|block| block.n_cols).collect::<Vec<_>>(),
        vec![2, 1]
    );
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
            .as_ref()
            .expect("sample metadata")
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["S1", "S3"]
    );
    assert_eq!(variant_ids(variants(&sparse.variants)), vec!["rs1", "rs2"]);
    let sparse_variants = variants(&sparse.variants);
    assert_eq!(variant_a0(sparse_variants, 0), "A");
    assert_eq!(variant_a1(sparse_variants, 0), "G");
    assert_eq!(variant_a0(sparse_variants, 1), "T");
    assert_eq!(variant_a1(sparse_variants, 1), "C");
}

#[test]
fn plink1_sparse_reconstructs_dense_when_no_missing_calls() {
    let dir = unique_dir("plink1-sparse-parity");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x0b, 0x2c]);

    let dense = plink_io::read_plink1_dense(&bed, &bim, &fam, None, None)
        .expect("dense plink should decode");
    let sparse = plink_io::read_plink1_sparse(&bed, &bim, &fam, None, None)
        .expect("sparse plink should decode");

    assert_eq!(
        sparse_values_dense_output(&sparse),
        dense_values_sample_major_output(&dense)
    );
}

#[test]
fn pbr_rust_plink1_002_sparse_blocks_match_stateless_csc_at_unit_and_partial_widths() {
    let dir = unique_dir("pbr-plink1-sparse-block-parity");
    let (bed, bim, fam) = write_three_variant_plink_fixture(&dir);
    let expected =
        ::genoio_io::read_plink1_sparse_windowed(&bed, &bim, &fam, None, None, None, true, true)
            .expect("whole sparse plink1 read should decode");

    for (block_size, expected_widths) in [(1, vec![1, 1, 1]), (2, vec![2, 1])] {
        let blocks = collect_plink1_sparse_blocks(
            &bed,
            &bim,
            &fam,
            plink1_sparse_block_options(None, None),
            block_size,
        );
        let (indptr, indices, data) = concatenate_sparse_blocks(&blocks);
        let ids = blocks
            .iter()
            .flat_map(|block| {
                variant_ids(variants(&block.variants))
                    .into_iter()
                    .map(str::to_owned)
            })
            .collect::<Vec<_>>();

        assert_eq!(
            blocks.iter().map(|block| block.n_cols).collect::<Vec<_>>(),
            expected_widths
        );
        assert_eq!(indptr, expected.indptr);
        assert_eq!(indices, expected.indices);
        assert_eq!(data, expected.data);
        assert_eq!(ids, variant_ids(variants(&expected.variants)));
        assert!(blocks.iter().all(|block| block.samples == expected.samples));
    }
}

#[test]
fn pbr_rust_plink1_002_sparse_blocks_preserve_sample_metadata_and_genotype_filters() {
    let dir = unique_dir("pbr-plink1-sparse-filter-parity");
    let (bed, bim, fam) = write_three_variant_plink_fixture(&dir);
    let requested_samples = vec!["S3".to_owned(), "S1".to_owned()];
    let filter = VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "maf",
        "params": {"min": 0.1}
    }))
    .expect("genotype-stat filter should parse");
    let expected = ::genoio_io::read_plink1_sparse_windowed(
        &bed,
        &bim,
        &fam,
        Some(&requested_samples),
        Some(&filter),
        None,
        true,
        true,
    )
    .expect("filtered whole sparse plink1 read should decode");
    let blocks = collect_plink1_sparse_blocks(
        &bed,
        &bim,
        &fam,
        plink1_sparse_block_options(Some(requested_samples), Some(filter)),
        1,
    );
    let (indptr, indices, data) = concatenate_sparse_blocks(&blocks);

    assert_eq!(indptr, expected.indptr);
    assert_eq!(indices, expected.indices);
    assert_eq!(data, expected.data);
    assert!(blocks.iter().all(|block| block.samples == expected.samples));
}

#[test]
fn pbr_rust_plink1_002_sparse_blocks_reject_retained_missing_calls() {
    let dir = unique_dir("pbr-plink1-sparse-missing");
    let (bed, bim, fam) = write_plink_fixture(&dir, &[0x6c, 0x1b, 0x01, 0x07, 0x2c]);
    let mut reader = BlockReader::open(
        BlockSource::Plink1 { bed, bim, fam },
        plink1_sparse_block_options(None, None),
        1,
    )
    .expect("persistent sparse plink1 reader should open");

    let error = reader
        .next_block()
        .expect_err("retained sparse missing hard calls should fail");

    assert!(matches!(error, GenoioError::MissingData { .. }));
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
    assert_eq!(sparse_values_dense_output(&sparse), vec![0.0, 0.0, 1.0]);
    let sparse_variants = variants(&sparse.variants);
    assert_eq!(variant_a0(sparse_variants, 0), "G");
    assert_eq!(variant_a1(sparse_variants, 0), "A");
    let dense_variants = variants(&dense.variants);
    assert_eq!(variant_a0(dense_variants, 0), "A");
    assert_eq!(variant_a1(dense_variants, 0), "G");
}
