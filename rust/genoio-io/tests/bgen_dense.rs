// pattern: Imperative Shell

use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

const FLAG_LAYOUT2: u32 = 2 << 2;
const FLAG_SAMPLE_IDENTIFIERS: u32 = 1 << 31;

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
    writer.write_all(&10_u32.to_le_bytes())?;
    writer.write_all(&n_samples.to_le_bytes())?;
    writer.write_all(&allele_count.to_le_bytes())?;
    writer.write_all(&2_u8.to_le_bytes())?;
    writer.write_all(&2_u8.to_le_bytes())?;
    writer.write_all(&0_u8.to_le_bytes())?;
    writer.write_all(&0_u8.to_le_bytes())?;
    Ok(())
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
