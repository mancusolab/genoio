// pattern: Imperative Shell

use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::write::ZlibEncoder;
use flate2::Compression;
use genoio_core::{VariantFilter, VariantWindow};
use serde_json::json;

const FLAG_LAYOUT2: u32 = 2 << 2;
const FLAG_SAMPLE_IDENTIFIERS: u32 = 1 << 31;
const FLAG_ZLIB_COMPRESSION: u32 = 1;
const FLAG_ZSTD_COMPRESSION: u32 = 2;
const FLAG_RESERVED_COMPRESSION: u32 = 3;

fn unique_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock should be after unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("genoio-{name}-{nanos}"));
    fs::create_dir(&dir).expect("test temp dir should be created");
    dir
}

fn write_bgen_header(
    writer: &mut impl Write,
    n_samples: u32,
    n_variants: u32,
    flags: u32,
    has_sample_ids: bool,
) -> io::Result<()> {
    let flags = if has_sample_ids {
        flags | FLAG_SAMPLE_IDENTIFIERS
    } else {
        flags
    };

    writer.write_all(&20_u32.to_le_bytes())?;
    writer.write_all(&20_u32.to_le_bytes())?;
    writer.write_all(&n_variants.to_le_bytes())?;
    writer.write_all(&n_samples.to_le_bytes())?;
    writer.write_all(b"bgen")?;
    writer.write_all(&flags.to_le_bytes())
}

fn write_bgen_header_with_free_data(
    writer: &mut impl Write,
    n_samples: u32,
    n_variants: u32,
    flags: u32,
    has_sample_ids: bool,
    free_data: &[u8],
) -> io::Result<()> {
    let flags = if has_sample_ids {
        flags | FLAG_SAMPLE_IDENTIFIERS
    } else {
        flags
    };
    let header_length = 20_u32
        .checked_add(u32::try_from(free_data.len()).expect("free data length should fit u32"))
        .expect("header length should fit u32");

    writer.write_all(&header_length.to_le_bytes())?;
    writer.write_all(&header_length.to_le_bytes())?;
    writer.write_all(&n_variants.to_le_bytes())?;
    writer.write_all(&n_samples.to_le_bytes())?;
    writer.write_all(b"bgen")?;
    writer.write_all(free_data)?;
    writer.write_all(&flags.to_le_bytes())
}

fn write_sample_identifier_block(writer: &mut impl Write, sample_ids: &[&str]) -> io::Result<()> {
    let block_len = 8 + sample_ids
        .iter()
        .map(|sample_id| 2 + sample_id.len())
        .sum::<usize>();
    let block_len = u32::try_from(block_len).expect("sample block length should fit u32");

    writer.write_all(&block_len.to_le_bytes())?;
    writer.write_all(
        &u32::try_from(sample_ids.len())
            .expect("sample count should fit u32")
            .to_le_bytes(),
    )?;
    for sample_id in sample_ids {
        writer.write_all(
            &u16::try_from(sample_id.len())
                .expect("sample id length should fit u16")
                .to_le_bytes(),
        )?;
        writer.write_all(sample_id.as_bytes())?;
    }
    Ok(())
}

fn write_layout2_variant_identifying_data(
    writer: &mut impl Write,
    id: &str,
    rsid: &str,
    chrom: &str,
    pos: u32,
    alleles: &[&str],
) -> io::Result<()> {
    writer.write_all(
        &u16::try_from(id.len())
            .expect("variant id length should fit u16")
            .to_le_bytes(),
    )?;
    writer.write_all(id.as_bytes())?;
    writer.write_all(
        &u16::try_from(rsid.len())
            .expect("variant rsid length should fit u16")
            .to_le_bytes(),
    )?;
    writer.write_all(rsid.as_bytes())?;
    writer.write_all(
        &u16::try_from(chrom.len())
            .expect("variant chromosome length should fit u16")
            .to_le_bytes(),
    )?;
    writer.write_all(chrom.as_bytes())?;
    writer.write_all(&pos.to_le_bytes())?;
    writer.write_all(
        &u16::try_from(alleles.len())
            .expect("allele count should fit u16")
            .to_le_bytes(),
    )?;
    for allele in alleles {
        writer.write_all(
            &u32::try_from(allele.len())
                .expect("allele length should fit u32")
                .to_le_bytes(),
        )?;
        writer.write_all(allele.as_bytes())?;
    }
    Ok(())
}

fn write_empty_layout2_probability_block(
    writer: &mut impl Write,
    n_samples: u32,
    allele_count: u16,
) -> io::Result<()> {
    let sample_ploidies =
        vec![2; usize::try_from(n_samples).expect("sample count should fit usize")];
    let mut block = layout2_probability_block_header_bytes(ProbabilityBlockHeader {
        n_samples,
        allele_count,
        min_ploidy: 2,
        max_ploidy: 2,
        sample_ploidies: &sample_ploidies,
        phased: 0,
        bit_depth: 8,
    });
    let packed_len = usize::try_from(n_samples)
        .expect("sample count should fit usize")
        .checked_mul(2)
        .expect("packed probability byte count should fit usize");
    block.resize(block.len() + packed_len, 0);
    let c = u32::try_from(block.len()).expect("probability block length should fit u32");
    writer.write_all(&c.to_le_bytes())?;
    writer.write_all(&block)
}

fn write_layout2_dosage_probability_block(
    writer: &mut impl Write,
    bit_depth: u8,
    calls: &[Option<(u32, u32)>],
) -> io::Result<()> {
    let sample_ploidies = calls
        .iter()
        .map(|call| if call.is_some() { 2 } else { 0b1000_0010 })
        .collect::<Vec<_>>();
    let mut block = layout2_probability_block_header_bytes(ProbabilityBlockHeader {
        n_samples: u32::try_from(calls.len()).expect("call count should fit u32"),
        allele_count: 2,
        min_ploidy: 2,
        max_ploidy: 2,
        sample_ploidies: &sample_ploidies,
        phased: 0,
        bit_depth,
    });
    append_packed_probabilities(&mut block, bit_depth, calls);

    let c = u32::try_from(block.len()).expect("probability block length should fit u32");
    writer.write_all(&c.to_le_bytes())?;
    writer.write_all(&block)?;
    Ok(())
}

fn write_layout2_probability_block(
    writer: &mut impl Write,
    header: ProbabilityBlockHeader<'_>,
    packed_probabilities: &[u8],
) -> io::Result<()> {
    let mut block = layout2_probability_block_header_bytes(header);
    block.extend_from_slice(packed_probabilities);
    let c = u32::try_from(block.len()).expect("probability block length should fit u32");
    writer.write_all(&c.to_le_bytes())?;
    writer.write_all(&block)?;
    Ok(())
}

fn append_packed_probabilities(output: &mut Vec<u8>, bit_depth: u8, calls: &[Option<(u32, u32)>]) {
    let mut current_byte = 0_u8;
    let mut bits_in_current_byte = 0_u8;
    for &(p_aa, p_ab) in calls.iter().flatten() {
        append_packed_probability_value(
            output,
            &mut current_byte,
            &mut bits_in_current_byte,
            bit_depth,
            p_aa,
        );
        append_packed_probability_value(
            output,
            &mut current_byte,
            &mut bits_in_current_byte,
            bit_depth,
            p_ab,
        );
    }
    if bits_in_current_byte > 0 {
        output.push(current_byte);
    }
}

fn append_packed_probability_value(
    output: &mut Vec<u8>,
    current_byte: &mut u8,
    bits_in_current_byte: &mut u8,
    bit_depth: u8,
    value: u32,
) {
    assert!(bit_depth <= 32);
    assert!(bit_depth == 32 || value < (1_u32 << bit_depth));

    for bit_index in 0..bit_depth {
        let bit = ((value >> bit_index) & 1) as u8;
        *current_byte |= bit << *bits_in_current_byte;
        *bits_in_current_byte += 1;
        if *bits_in_current_byte == 8 {
            output.push(*current_byte);
            *current_byte = 0;
            *bits_in_current_byte = 0;
        }
    }
}

struct ProbabilityBlockHeader<'a> {
    n_samples: u32,
    allele_count: u16,
    min_ploidy: u8,
    max_ploidy: u8,
    sample_ploidies: &'a [u8],
    phased: u8,
    bit_depth: u8,
}

fn write_layout2_probability_block_header(
    writer: &mut impl Write,
    header: ProbabilityBlockHeader<'_>,
) -> io::Result<()> {
    let block = layout2_probability_block_header_bytes(header);
    let c = u32::try_from(block.len()).expect("probability block length should fit u32");
    writer.write_all(&c.to_le_bytes())?;
    writer.write_all(&block)?;
    Ok(())
}

fn layout2_probability_block_header_bytes(header: ProbabilityBlockHeader<'_>) -> Vec<u8> {
    let mut block = Vec::new();
    block.extend_from_slice(&header.n_samples.to_le_bytes());
    block.extend_from_slice(&header.allele_count.to_le_bytes());
    block.extend_from_slice(&header.min_ploidy.to_le_bytes());
    block.extend_from_slice(&header.max_ploidy.to_le_bytes());
    block.extend_from_slice(header.sample_ploidies);
    block.extend_from_slice(&header.phased.to_le_bytes());
    block.extend_from_slice(&header.bit_depth.to_le_bytes());
    block
}

#[derive(Clone, Copy)]
enum TestCompression {
    Zlib,
    Zstd,
}

fn write_compressed_layout2_probability_block(
    writer: &mut impl Write,
    compression: TestCompression,
    header: ProbabilityBlockHeader<'_>,
    packed_probabilities: &[u8],
) -> io::Result<()> {
    let mut decompressed_payload = layout2_probability_block_header_bytes(header);
    decompressed_payload.extend_from_slice(packed_probabilities);
    let compressed_payload = match compression {
        TestCompression::Zlib => {
            let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(&decompressed_payload)?;
            encoder.finish()?
        }
        TestCompression::Zstd => zstd::stream::encode_all(&decompressed_payload[..], 0)?,
    };
    let c = 4_u32
        .checked_add(
            u32::try_from(compressed_payload.len())
                .expect("compressed payload length should fit u32"),
        )
        .expect("compressed block length should fit u32");
    let d = u32::try_from(decompressed_payload.len()).expect("decompressed length should fit u32");
    writer.write_all(&c.to_le_bytes())?;
    writer.write_all(&d.to_le_bytes())?;
    writer.write_all(&compressed_payload)?;
    Ok(())
}

fn write_valid_compressed_probability_block(
    writer: &mut impl Write,
    compression: TestCompression,
) -> io::Result<()> {
    write_compressed_layout2_probability_block(
        writer,
        compression,
        ProbabilityBlockHeader {
            n_samples: 2,
            allele_count: 2,
            min_ploidy: 2,
            max_ploidy: 2,
            sample_ploidies: &[2, 2],
            phased: 0,
            bit_depth: 8,
        },
        &[0, 0, 0, 0],
    )
}

fn write_bgen_fixture(
    path: &Path,
    flags: u32,
    n_variants: u32,
    write_variants: impl FnOnce(&mut Vec<u8>),
) {
    let mut bgen = Vec::new();
    write_bgen_header(&mut bgen, 2, n_variants, flags, true).expect("header should write");
    write_sample_identifier_block(&mut bgen, &["sample_1", "sample_2"])
        .expect("sample block should write");
    let variant_offset = u32::try_from(bgen.len() - 4).expect("variant offset should fit u32");
    bgen[0..4].copy_from_slice(&variant_offset.to_le_bytes());
    write_variants(&mut bgen);
    fs::write(path, bgen).expect("bgen fixture should be written");
}

fn write_two_sample_two_variant_bgen(path: &Path) {
    write_two_sample_two_variant_bgen_with_sample_ids(path, true);
}

fn write_two_sample_two_variant_bgen_without_sample_ids(path: &Path) {
    write_two_sample_two_variant_bgen_with_sample_ids(path, false);
}

fn write_two_sample_two_variant_bgen_with_sample_ids(path: &Path, has_sample_ids: bool) {
    let mut bgen = Vec::new();
    write_bgen_header(&mut bgen, 2, 2, FLAG_LAYOUT2, has_sample_ids).expect("header should write");
    if has_sample_ids {
        write_sample_identifier_block(&mut bgen, &["sample_1", "sample_2"])
            .expect("sample block should write");
    }
    let variant_offset = u32::try_from(bgen.len() - 4).expect("variant offset should fit u32");
    bgen[0..4].copy_from_slice(&variant_offset.to_le_bytes());

    write_layout2_variant_identifying_data(&mut bgen, "var1", "rs1", "1", 10, &["A", "G"])
        .expect("first variant identifying data should write");
    write_empty_layout2_probability_block(&mut bgen, 2, 2)
        .expect("first probability block should write");
    write_layout2_variant_identifying_data(&mut bgen, "var2", "rs2", "2", 20, &["C", "T"])
        .expect("second variant identifying data should write");
    write_empty_layout2_probability_block(&mut bgen, 2, 2)
        .expect("second probability block should write");

    fs::write(path, bgen).expect("bgen fixture should be written");
}

fn write_two_sample_two_variant_dosage_bgen(
    path: &Path,
    bit_depth: u8,
    variant_calls: &[[Option<(u32, u32)>; 2]; 2],
) {
    write_bgen_fixture(path, FLAG_LAYOUT2, 2, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("first variant identifying data should write");
        write_layout2_dosage_probability_block(writer, bit_depth, &variant_calls[0])
            .expect("first dosage probability block should write");
        write_layout2_variant_identifying_data(writer, "var2", "rs2", "2", 20, &["C", "T"])
            .expect("second variant identifying data should write");
        write_layout2_dosage_probability_block(writer, bit_depth, &variant_calls[1])
            .expect("second dosage probability block should write");
    });
}

fn write_three_sample_two_variant_dosage_bgen(
    path: &Path,
    bit_depth: u8,
    variant_calls: &[[Option<(u32, u32)>; 3]; 2],
) {
    let mut bgen = Vec::new();
    write_bgen_header(&mut bgen, 3, 2, FLAG_LAYOUT2, true).expect("header should write");
    write_sample_identifier_block(&mut bgen, &["sample_1", "sample_2", "sample_3"])
        .expect("sample block should write");
    let variant_offset = u32::try_from(bgen.len() - 4).expect("variant offset should fit u32");
    bgen[0..4].copy_from_slice(&variant_offset.to_le_bytes());

    write_layout2_variant_identifying_data(&mut bgen, "var1", "rs1", "1", 10, &["A", "G"])
        .expect("first variant identifying data should write");
    write_layout2_dosage_probability_block(&mut bgen, bit_depth, &variant_calls[0])
        .expect("first dosage probability block should write");
    write_layout2_variant_identifying_data(&mut bgen, "var2", "rs2", "2", 20, &["C", "T"])
        .expect("second variant identifying data should write");
    write_layout2_dosage_probability_block(&mut bgen, bit_depth, &variant_calls[1])
        .expect("second dosage probability block should write");

    fs::write(path, bgen).expect("bgen fixture should be written");
}

fn expected_dosage(bit_depth: u8, p_aa: u32, p_ab: u32) -> f32 {
    let denominator = ((1_u64 << bit_depth) - 1) as f32;
    let p_aa = p_aa as f32 / denominator;
    let p_ab = p_ab as f32 / denominator;
    p_ab + 2.0 * (1.0 - p_aa - p_ab)
}

fn write_sample_file(path: &Path, rows: &[&str]) {
    let mut contents = String::from("ID_1 ID_2 missing\n0 0 0\n");
    for row in rows {
        contents.push_str(row);
        contents.push('\n');
    }
    fs::write(path, contents).expect("sample fixture should be written");
}

fn assert_metadata_error_contains(error: genoio_core::MetadataError, expected: &str) {
    let message = error.to_string();
    assert!(
        message.contains(expected),
        "expected error containing {expected:?}, got {message:?}"
    );
}

fn chrom_filter(chrom: &str) -> VariantFilter {
    VariantFilter::from_json_value(json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": chrom},
    }))
    .expect("chrom filter should parse")
}

#[test]
fn bgen_dosage_dense_decodes_uncompressed_bit_depth_8() {
    let dir = unique_dir("bgen-dosage-bit-depth-8");
    let bgen = dir.join("tiny.bgen");
    let calls = [
        [Some((204, 26)), Some((51, 128))],
        [Some((0, 255)), Some((102, 102))],
    ];
    write_two_sample_two_variant_dosage_bgen(&bgen, 8, &calls);

    let dense = genoio_io::read_bgen_dosage_dense_windowed(&bgen, None, None, None, None)
        .expect("bgen dosage should decode");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(dense.values.len(), 4);
    let expected = vec![
        expected_dosage(8, 204, 26),
        expected_dosage(8, 0, 255),
        expected_dosage(8, 51, 128),
        expected_dosage(8, 102, 102),
    ];
    for (observed, expected) in dense.values.iter().zip(expected) {
        assert!((observed - expected).abs() <= 2.0 / 255.0);
    }
    assert_eq!(dense.missing_mask, vec![false, false, false, false]);
}

#[test]
fn bgen_dosage_dense_decodes_uncompressed_bit_depth_16() {
    let dir = unique_dir("bgen-dosage-bit-depth-16");
    let bgen = dir.join("tiny.bgen");
    let calls = [
        [Some((52_428, 6_554)), Some((13_107, 32_768))],
        [Some((0, 65_535)), Some((26_214, 26_214))],
    ];
    write_two_sample_two_variant_dosage_bgen(&bgen, 16, &calls);

    let dense = genoio_io::read_bgen_dosage_dense_windowed(&bgen, None, None, None, None)
        .expect("bgen dosage should decode");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    let expected = vec![
        expected_dosage(16, 52_428, 6_554),
        expected_dosage(16, 0, 65_535),
        expected_dosage(16, 13_107, 32_768),
        expected_dosage(16, 26_214, 26_214),
    ];
    for (observed, expected) in dense.values.iter().zip(expected) {
        assert!((observed - expected).abs() <= 2.0 / 65_535.0);
    }
    assert_eq!(dense.missing_mask, vec![false, false, false, false]);
}

#[test]
fn bgen_dosage_dense_preserves_missing_sample_calls() {
    let dir = unique_dir("bgen-dosage-missing");
    let bgen = dir.join("tiny.bgen");
    let calls = [[Some((204, 26)), None], [Some((0, 255)), Some((102, 102))]];
    write_two_sample_two_variant_dosage_bgen(&bgen, 8, &calls);

    let dense = genoio_io::read_bgen_dosage_dense_windowed(&bgen, None, None, None, None)
        .expect("bgen dosage should decode");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(
        dense.values,
        vec![
            expected_dosage(8, 204, 26),
            expected_dosage(8, 0, 255),
            0.0,
            expected_dosage(8, 102, 102),
        ]
    );
    assert_eq!(dense.missing_mask, vec![false, false, true, false]);
}

#[test]
fn bgen_dosage_dense_sample_filter_uses_source_order() {
    let dir = unique_dir("bgen-dosage-sample-filter");
    let bgen = dir.join("tiny.bgen");
    let calls = [
        [Some((204, 26)), Some((51, 128)), Some((0, 0))],
        [Some((0, 255)), Some((102, 102)), Some((255, 0))],
    ];
    write_three_sample_two_variant_dosage_bgen(&bgen, 8, &calls);
    let requested_samples = vec!["sample_3".to_string(), "sample_1".to_string()];

    let dense = genoio_io::read_bgen_dosage_dense_windowed(
        &bgen,
        None,
        Some(&requested_samples),
        None,
        None,
    )
    .expect("bgen dosage sample filter should decode");

    assert_eq!(dense.n_samples, 2);
    assert_eq!(dense.n_variants, 2);
    assert_eq!(
        dense
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["sample_1", "sample_3"]
    );
    assert_eq!(
        dense.values,
        vec![
            expected_dosage(8, 204, 26),
            expected_dosage(8, 0, 255),
            expected_dosage(8, 0, 0),
            expected_dosage(8, 255, 0),
        ]
    );
    assert_eq!(dense.missing_mask, vec![false, false, false, false]);
}

#[test]
fn bgen_dosage_dense_skips_metadata_rejected_unsupported_probability_block() {
    let dir = unique_dir("bgen-dosage-skip-metadata-reject");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 2, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("first variant identifying data should write");
        write_layout2_probability_block(
            writer,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 2],
                phased: 1,
                bit_depth: 8,
            },
            &[],
        )
        .expect("unsupported probability block should write");
        write_layout2_variant_identifying_data(writer, "var2", "rs2", "2", 20, &["C", "T"])
            .expect("second variant identifying data should write");
        write_layout2_dosage_probability_block(writer, 8, &[Some((0, 255)), Some((102, 102))])
            .expect("second dosage probability block should write");
    });
    let filter = chrom_filter("2");

    let dense = genoio_io::read_bgen_dosage_dense_windowed(&bgen, None, None, Some(&filter), None)
        .expect("metadata-rejected unsupported probabilities should be skipped");

    assert_eq!(dense.n_variants, 1);
    assert_eq!(dense.variants[0].id, "rs2");
    assert_eq!(dense.values.len(), 2);
}

#[test]
fn bgen_dosage_dense_skips_out_of_window_unsupported_probability_block() {
    let dir = unique_dir("bgen-dosage-skip-window");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 2, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("first variant identifying data should write");
        write_layout2_probability_block(
            writer,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 2],
                phased: 1,
                bit_depth: 8,
            },
            &[],
        )
        .expect("unsupported probability block should write");
        write_layout2_variant_identifying_data(writer, "var2", "rs2", "2", 20, &["C", "T"])
            .expect("second variant identifying data should write");
        write_layout2_dosage_probability_block(writer, 8, &[Some((0, 255)), Some((102, 102))])
            .expect("second dosage probability block should write");
    });

    let dense = genoio_io::read_bgen_dosage_dense_windowed(
        &bgen,
        None,
        None,
        None,
        Some(VariantWindow { start: 1, len: 1 }),
    )
    .expect("out-of-window unsupported probabilities should be skipped");

    assert_eq!(dense.n_variants, 1);
    assert_eq!(dense.variants[0].id, "rs2");
}

#[test]
fn bgen_metadata_reads_header_with_free_data_before_flags() {
    let dir = unique_dir("bgen-header-free-data");
    let bgen = dir.join("tiny.bgen");
    let mut contents = Vec::new();

    write_bgen_header_with_free_data(
        &mut contents,
        2,
        1,
        FLAG_LAYOUT2,
        true,
        &[0xAA, 0xBB, 0xCC, 0xDD],
    )
    .expect("header should write");
    write_sample_identifier_block(&mut contents, &["sample_1", "sample_2"])
        .expect("sample block should write");
    let variant_offset = u32::try_from(contents.len() - 4).expect("variant offset should fit u32");
    contents[0..4].copy_from_slice(&variant_offset.to_le_bytes());
    write_layout2_variant_identifying_data(&mut contents, "var1", "rs1", "1", 10, &["A", "G"])
        .expect("variant identifying data should write");
    write_empty_layout2_probability_block(&mut contents, 2, 2)
        .expect("probability block should write");
    fs::write(&bgen, contents).expect("bgen fixture should be written");

    let metadata = genoio_io::read_bgen_metadata(&bgen, None)
        .expect("bgen metadata with header free data should parse");

    assert_eq!(
        metadata
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["sample_1", "sample_2"]
    );
    assert_eq!(metadata.variants[0].id, "rs1");
}

#[test]
fn bgen_metadata_reads_embedded_sample_ids_and_variant_rows() {
    let dir = unique_dir("bgen-metadata");
    let bgen = dir.join("tiny.bgen");
    write_two_sample_two_variant_bgen(&bgen);

    let metadata = genoio_io::read_bgen_metadata(&bgen, None).expect("bgen metadata should parse");

    assert_eq!(
        metadata
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["sample_1", "sample_2"]
    );
    assert_eq!(
        metadata
            .variants
            .iter()
            .map(|variant| {
                (
                    variant.chrom.as_str(),
                    variant.pos,
                    variant.id.as_str(),
                    variant.a0.as_str(),
                    variant.a1.as_str(),
                )
            })
            .collect::<Vec<_>>(),
        vec![("1", 10, "rs1", "A", "G"), ("2", 20, "rs2", "C", "T")]
    );
    assert_eq!(metadata.variants[0].ref_allele.as_deref(), Some("A"));
    assert_eq!(metadata.variants[0].alt_allele.as_deref(), Some("G"));
    assert_eq!(metadata.variants[0].source_a0, "A");
    assert_eq!(metadata.variants[0].source_a1, "G");
    assert!(!metadata.variants[0].flipped);
    assert_eq!(metadata.variants[0].qual, None);
    assert_eq!(metadata.variants[0].af, None);
    assert_eq!(metadata.variants[0].maf, None);
    assert_eq!(metadata.variants[0].mac, None);
    assert_eq!(metadata.variants[0].missing_rate, None);
    assert_eq!(metadata.variants[0].n_called, None);
    assert!(metadata.capabilities.supports_geno);
    assert!(!metadata.capabilities.supports_haplo);
    assert!(!metadata.capabilities.phased);
}

#[test]
fn bgen_metadata_reads_companion_sample_ids_when_not_embedded() {
    let dir = unique_dir("bgen-companion-sample");
    let bgen = dir.join("tiny.bgen");
    let sample = dir.join("tiny.sample");
    write_two_sample_two_variant_bgen_without_sample_ids(&bgen);
    write_sample_file(&sample, &["sample_a sample_a 0", "sample_b sample_b 0"]);

    let metadata =
        genoio_io::read_bgen_metadata(&bgen, Some(&sample)).expect("metadata should parse");

    assert_eq!(
        metadata
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["sample_a", "sample_b"]
    );
    assert_eq!(metadata.variants.len(), 2);
}

#[test]
fn bgen_metadata_rejects_companion_sample_count_mismatch() {
    let dir = unique_dir("bgen-companion-count-mismatch");
    let bgen = dir.join("tiny.bgen");
    let sample = dir.join("tiny.sample");
    write_two_sample_two_variant_bgen_without_sample_ids(&bgen);
    write_sample_file(&sample, &["sample_a sample_a 0"]);

    let error = genoio_io::read_bgen_metadata(&bgen, Some(&sample))
        .expect_err("sample count mismatch should fail");

    assert_metadata_error_contains(error, "sample count");
}

#[test]
fn bgen_metadata_rejects_duplicate_companion_sample_ids() {
    let dir = unique_dir("bgen-companion-duplicates");
    let bgen = dir.join("tiny.bgen");
    let sample = dir.join("tiny.sample");
    write_two_sample_two_variant_bgen_without_sample_ids(&bgen);
    write_sample_file(&sample, &["sample_a sample_a 0", "sample_a sample_a 0"]);

    let error = genoio_io::read_bgen_metadata(&bgen, Some(&sample))
        .expect_err("duplicate sample ids should fail");

    assert_metadata_error_contains(error, "duplicate sample identifier");
}

#[test]
fn bgen_metadata_rejects_missing_companion_path_when_sample_ids_not_embedded() {
    let dir = unique_dir("bgen-missing-companion");
    let bgen = dir.join("tiny.bgen");
    write_two_sample_two_variant_bgen_without_sample_ids(&bgen);

    let error =
        genoio_io::read_bgen_metadata(&bgen, None).expect_err("missing sample IDs should fail");

    assert_metadata_error_contains(error, "companion sample path");
}

#[test]
fn bgen_metadata_uses_variant_id_when_rsid_is_empty() {
    let dir = unique_dir("bgen-empty-rsid");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_empty_layout2_probability_block(writer, 2, 2)
            .expect("probability block should write");
    });

    let metadata = genoio_io::read_bgen_metadata(&bgen, None).expect("metadata should parse");

    assert_eq!(metadata.variants[0].id, "var1");
}

#[test]
fn bgen_metadata_skips_compressed_probability_blocks_before_next_variant() {
    let dir = unique_dir("bgen-compressed-skip");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2 | FLAG_ZLIB_COMPRESSION, 2, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("first variant identifying data should write");
        write_valid_compressed_probability_block(writer, TestCompression::Zlib)
            .expect("first compressed probability block should write");
        write_layout2_variant_identifying_data(writer, "var2", "rs2", "2", 20, &["C", "T"])
            .expect("second variant identifying data should write");
        write_valid_compressed_probability_block(writer, TestCompression::Zlib)
            .expect("second compressed probability block should write");
    });

    let metadata =
        genoio_io::read_bgen_metadata(&bgen, None).expect("compressed metadata should parse");

    assert_eq!(
        metadata
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs1", "rs2"]
    );
}

#[test]
fn bgen_metadata_rejects_multiallelic_variants() {
    let dir = unique_dir("bgen-multiallelic");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "C", "G"])
            .expect("variant identifying data should write");
        write_empty_layout2_probability_block(writer, 2, 3)
            .expect("probability block should write");
    });

    let error =
        genoio_io::read_bgen_metadata(&bgen, None).expect_err("multiallelic BGEN should fail");

    assert_metadata_error_contains(error, "multiallelic");
}

#[test]
fn bgen_metadata_rejects_phased_layout2_probability_blocks() {
    let dir = unique_dir("bgen-phased-probability-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_layout2_probability_block_header(
            writer,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 2],
                phased: 1,
                bit_depth: 8,
            },
        )
        .expect("phased probability block should write");
    });

    let error = genoio_io::read_bgen_metadata(&bgen, None).expect_err("phased BGEN should fail");

    assert_metadata_error_contains(error, "phased");
}

#[test]
fn bgen_metadata_rejects_variable_ploidy_layout2_probability_blocks() {
    let dir = unique_dir("bgen-variable-ploidy-probability-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_layout2_probability_block_header(
            writer,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 1,
                max_ploidy: 2,
                sample_ploidies: &[1, 2],
                phased: 0,
                bit_depth: 8,
            },
        )
        .expect("variable-ploidy probability block should write");
    });

    let error =
        genoio_io::read_bgen_metadata(&bgen, None).expect_err("variable-ploidy BGEN should fail");

    assert_metadata_error_contains(error, "variable-ploidy");
}

#[test]
fn bgen_metadata_rejects_non_diploid_sample_ploidy_bytes() {
    let dir = unique_dir("bgen-non-diploid-sample-ploidy");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_layout2_probability_block_header(
            writer,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 1],
                phased: 0,
                bit_depth: 8,
            },
        )
        .expect("non-diploid probability block should write");
    });

    let error =
        genoio_io::read_bgen_metadata(&bgen, None).expect_err("non-diploid BGEN should fail");

    assert_metadata_error_contains(error, "variable-ploidy");
}

#[test]
fn bgen_metadata_allows_missing_diploid_sample_ploidy_bytes() {
    let dir = unique_dir("bgen-missing-diploid-sample-ploidy");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_layout2_dosage_probability_block(writer, 8, &[Some((0, 0)), None])
            .expect("missing probability block should write");
    });

    let metadata = genoio_io::read_bgen_metadata(&bgen, None).expect("missing calls should parse");

    assert_eq!(metadata.variants.len(), 1);
}

#[test]
fn bgen_metadata_rejects_zero_bit_depth_probability_blocks() {
    let dir = unique_dir("bgen-zero-bit-depth-probability-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_layout2_probability_block_header(
            writer,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 2],
                phased: 0,
                bit_depth: 0,
            },
        )
        .expect("zero-bit-depth probability block should write");
    });

    let error = genoio_io::read_bgen_metadata(&bgen, None).expect_err("zero bit depth should fail");

    assert_metadata_error_contains(error, "bit depth");
}

#[test]
fn bgen_metadata_rejects_too_large_bit_depth_probability_blocks() {
    let dir = unique_dir("bgen-too-large-bit-depth-probability-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_layout2_probability_block_header(
            writer,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 2],
                phased: 0,
                bit_depth: 33,
            },
        )
        .expect("too-large-bit-depth probability block should write");
    });

    let error =
        genoio_io::read_bgen_metadata(&bgen, None).expect_err("large bit depth should fail");

    assert_metadata_error_contains(error, "bit depth");
}

#[test]
fn bgen_metadata_rejects_truncated_packed_probability_bytes() {
    let dir = unique_dir("bgen-truncated-packed-probabilities");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_layout2_probability_block_header(
            writer,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 2],
                phased: 0,
                bit_depth: 8,
            },
        )
        .expect("truncated probability block should write");
    });

    let error = genoio_io::read_bgen_metadata(&bgen, None)
        .expect_err("truncated packed probabilities should fail");

    assert_metadata_error_contains(error, "truncated");
}

#[test]
fn bgen_metadata_rejects_zlib_compressed_phased_probability_blocks() {
    let dir = unique_dir("bgen-zlib-phased-probability-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2 | FLAG_ZLIB_COMPRESSION, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_compressed_layout2_probability_block(
            writer,
            TestCompression::Zlib,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 2],
                phased: 1,
                bit_depth: 8,
            },
            &[],
        )
        .expect("compressed phased probability block should write");
    });

    let error = genoio_io::read_bgen_metadata(&bgen, None)
        .expect_err("zlib-compressed phased BGEN should fail");

    assert_metadata_error_contains(error, "phased");
}

#[test]
fn bgen_metadata_rejects_zlib_compressed_variable_ploidy_probability_blocks() {
    let dir = unique_dir("bgen-zlib-variable-ploidy-probability-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2 | FLAG_ZLIB_COMPRESSION, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_compressed_layout2_probability_block(
            writer,
            TestCompression::Zlib,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 1,
                max_ploidy: 2,
                sample_ploidies: &[1, 2],
                phased: 0,
                bit_depth: 8,
            },
            &[],
        )
        .expect("compressed variable-ploidy probability block should write");
    });

    let error = genoio_io::read_bgen_metadata(&bgen, None)
        .expect_err("zlib-compressed variable-ploidy BGEN should fail");

    assert_metadata_error_contains(error, "variable-ploidy");
}

#[test]
fn bgen_metadata_rejects_zlib_compressed_non_diploid_sample_ploidy_bytes() {
    let dir = unique_dir("bgen-zlib-non-diploid-sample-ploidy");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2 | FLAG_ZLIB_COMPRESSION, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_compressed_layout2_probability_block(
            writer,
            TestCompression::Zlib,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 1],
                phased: 0,
                bit_depth: 8,
            },
            &[],
        )
        .expect("compressed non-diploid probability block should write");
    });

    let error = genoio_io::read_bgen_metadata(&bgen, None)
        .expect_err("zlib-compressed non-diploid BGEN should fail");

    assert_metadata_error_contains(error, "variable-ploidy");
}

#[test]
fn bgen_metadata_rejects_zstd_compressed_phased_probability_blocks() {
    let dir = unique_dir("bgen-zstd-phased-probability-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2 | FLAG_ZSTD_COMPRESSION, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_compressed_layout2_probability_block(
            writer,
            TestCompression::Zstd,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 2,
                max_ploidy: 2,
                sample_ploidies: &[2, 2],
                phased: 1,
                bit_depth: 8,
            },
            &[],
        )
        .expect("compressed phased probability block should write");
    });

    let error = genoio_io::read_bgen_metadata(&bgen, None)
        .expect_err("zstd-compressed phased BGEN should fail");

    assert_metadata_error_contains(error, "phased");
}

#[test]
fn bgen_metadata_rejects_zstd_compressed_variable_ploidy_probability_blocks() {
    let dir = unique_dir("bgen-zstd-variable-ploidy-probability-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2 | FLAG_ZSTD_COMPRESSION, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        write_compressed_layout2_probability_block(
            writer,
            TestCompression::Zstd,
            ProbabilityBlockHeader {
                n_samples: 2,
                allele_count: 2,
                min_ploidy: 1,
                max_ploidy: 2,
                sample_ploidies: &[1, 2],
                phased: 0,
                bit_depth: 8,
            },
            &[],
        )
        .expect("compressed variable-ploidy probability block should write");
    });

    let error = genoio_io::read_bgen_metadata(&bgen, None)
        .expect_err("zstd-compressed variable-ploidy BGEN should fail");

    assert_metadata_error_contains(error, "variable-ploidy");
}

#[test]
fn bgen_metadata_rejects_compressed_probability_block_shorter_than_length_prefix() {
    let dir = unique_dir("bgen-short-compressed-block");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2 | FLAG_ZLIB_COMPRESSION, 1, |writer| {
        write_layout2_variant_identifying_data(writer, "var1", "rs1", "1", 10, &["A", "G"])
            .expect("variant identifying data should write");
        writer
            .write_all(&3_u32.to_le_bytes())
            .expect("short C should write");
    });

    let error =
        genoio_io::read_bgen_metadata(&bgen, None).expect_err("short compressed block should fail");

    assert_metadata_error_contains(error, "probability block");
}

#[test]
fn bgen_metadata_rejects_unsupported_compression_flag() {
    let dir = unique_dir("bgen-reserved-compression");
    let bgen = dir.join("tiny.bgen");
    write_bgen_fixture(&bgen, FLAG_LAYOUT2 | FLAG_RESERVED_COMPRESSION, 0, |_| {});

    let error =
        genoio_io::read_bgen_metadata(&bgen, None).expect_err("reserved compression should fail");

    assert_metadata_error_contains(error, "compression");
}
