use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use genoio_core::VariantWindow;

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

fn write_vcf(path: &Path) {
    write_text(
        path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t20\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1
2\t30\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/0
1\t40\trs4\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/0
",
    );
}

fn write_plink_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let bed = dir.join("tiny.bed");
    let bim = dir.join("tiny.bim");
    let fam = dir.join("tiny.fam");
    fs::write(&bed, [0x6c, 0x1b, 0x01, 0x04, 0x0d, 0x03, 0x00])
        .expect("bed fixture should be written");
    write_text(
        &bim,
        "\
1 rs1 0 10 G A
1 rs2 0 20 T C
2 rs3 0 30 A G
1 rs4 0 40 C T
",
    );
    write_text(
        &fam,
        "\
F1 S1 0 0 1 -9
F1 S2 0 0 2 -9
",
    );
    (bed, bim, fam)
}

#[test]
fn vcf_dense_window_uses_retained_variant_order_after_filters() {
    let dir = unique_dir("vcf-block-window");
    let path = dir.join("blocks.vcf");
    write_vcf(&path);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "1"}
    }))
    .expect("filter should parse");

    let block = genoio_io::read_vcf_dense_windowed(
        &path,
        None,
        Some(&filter),
        Some(VariantWindow { start: 1, len: 2 }),
    )
    .expect("windowed vcf should decode");

    assert_eq!(
        block
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs2", "rs4"]
    );
    assert_eq!(block.values, vec![1.0, 0.0, 2.0, 0.0]);
}

#[test]
fn plink1_dense_window_uses_retained_variant_order_after_filters() {
    let dir = unique_dir("plink-block-window");
    let (bed, bim, fam) = write_plink_fixture(&dir);
    let filter = genoio_core::VariantFilter::from_json_value(serde_json::json!({
        "op": "predicate",
        "name": "chrom",
        "params": {"value": "1"}
    }))
    .expect("filter should parse");

    let block = genoio_io::read_plink1_dense_windowed(
        &bed,
        &bim,
        &fam,
        None,
        Some(&filter),
        Some(VariantWindow { start: 1, len: 2 }),
    )
    .expect("windowed plink should decode");

    assert_eq!(
        block
            .variants
            .iter()
            .map(|variant| variant.id.as_str())
            .collect::<Vec<_>>(),
        vec!["rs2", "rs4"]
    );
    assert_eq!(block.values, vec![0.0, 2.0, 0.0, 2.0]);
}
