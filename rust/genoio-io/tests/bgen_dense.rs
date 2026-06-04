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

    writer.write_all(&(block_len as u32).to_le_bytes())?;
    writer.write_all(&(sample_ids.len() as u32).to_le_bytes())?;
    for sample_id in sample_ids {
        writer.write_all(&(sample_id.len() as u16).to_le_bytes())?;
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
    writer.write_all(&(id.len() as u16).to_le_bytes())?;
    writer.write_all(id.as_bytes())?;
    writer.write_all(&(rsid.len() as u16).to_le_bytes())?;
    writer.write_all(rsid.as_bytes())?;
    writer.write_all(&(chrom.len() as u16).to_le_bytes())?;
    writer.write_all(chrom.as_bytes())?;
    writer.write_all(&pos.to_le_bytes())?;
    writer.write_all(&(alleles.len() as u16).to_le_bytes())?;
    for allele in alleles {
        writer.write_all(&(allele.len() as u32).to_le_bytes())?;
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
    writer.write_all(&(allele_count as u16).to_le_bytes())?;
    writer.write_all(&2_u8.to_le_bytes())?;
    writer.write_all(&2_u8.to_le_bytes())?;
    writer.write_all(&0_u8.to_le_bytes())?;
    writer.write_all(&0_u8.to_le_bytes())?;
    Ok(())
}

fn write_two_sample_two_variant_bgen(path: &Path) {
    let mut bgen = Vec::new();
    write_bgen_header(&mut bgen, 2, 2, FLAG_LAYOUT2, true).expect("header should write");
    write_sample_identifier_block(&mut bgen, &["sample_1", "sample_2"])
        .expect("sample block should write");
    let variant_offset = (bgen.len() - 4) as u32;
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
