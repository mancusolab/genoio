// pattern: Imperative Shell

use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

const FLAG_LAYOUT2: u32 = 2 << 2;
const FLAG_SAMPLE_IDENTIFIERS: u32 = 1 << 31;
const FLAG_ZLIB_COMPRESSION: u32 = 1;
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
    write_layout2_probability_block_header(
        writer,
        ProbabilityBlockHeader {
            n_samples,
            allele_count,
            min_ploidy: 2,
            max_ploidy: 2,
            sample_ploidies: &[2, 2],
            phased: 0,
            bit_depth: 0,
        },
    )
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
    let c = 10_u32
        .checked_add(
            u32::try_from(header.sample_ploidies.len())
                .expect("sample ploidy count should fit u32"),
        )
        .expect("probability block length should fit u32");
    writer.write_all(&c.to_le_bytes())?;
    writer.write_all(&header.n_samples.to_le_bytes())?;
    writer.write_all(&header.allele_count.to_le_bytes())?;
    writer.write_all(&header.min_ploidy.to_le_bytes())?;
    writer.write_all(&header.max_ploidy.to_le_bytes())?;
    writer.write_all(header.sample_ploidies)?;
    writer.write_all(&header.phased.to_le_bytes())?;
    writer.write_all(&header.bit_depth.to_le_bytes())?;
    Ok(())
}

fn write_placeholder_compressed_probability_block(writer: &mut impl Write) -> io::Result<()> {
    let compressed_payload = [0x11_u8, 0x22, 0x33];
    let c = 4_u32
        .checked_add(
            u32::try_from(compressed_payload.len())
                .expect("compressed payload length should fit u32"),
        )
        .expect("compressed block length should fit u32");
    let d = 10_u32;
    writer.write_all(&c.to_le_bytes())?;
    writer.write_all(&d.to_le_bytes())?;
    writer.write_all(&compressed_payload)
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
        write_placeholder_compressed_probability_block(writer)
            .expect("first compressed probability block should write");
        write_layout2_variant_identifying_data(writer, "var2", "rs2", "2", 20, &["C", "T"])
            .expect("second variant identifying data should write");
        write_placeholder_compressed_probability_block(writer)
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
