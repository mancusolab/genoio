use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock should be after unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("genoio-{name}-{nanos}"));
    fs::create_dir(&dir).expect("test temp dir should be created");
    dir
}

fn write_file(path: &Path, contents: &str) {
    fs::write(path, contents).expect("test fixture should be written");
}

#[test]
fn vcf_metadata_preserves_header_sample_and_variant_order() {
    let dir = unique_dir("vcf-metadata");
    let path = dir.join("tiny.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1
1\t20\trs2\tC\tT,A\t.\tPASS\t.\tGT\t0/1\t0/0\t1/2
2\t30\tindel1\tAT\tA\t.\tPASS\t.\tGT\t0/1\t./.\t0/0
",
    );

    let metadata = genoio_io::read_vcf_metadata(&path).expect("vcf metadata should parse");

    assert_eq!(
        metadata
            .samples
            .iter()
            .map(|sample| sample.iid.as_str())
            .collect::<Vec<_>>(),
        vec!["s1", "s2", "s3"]
    );
    assert_eq!(
        metadata
            .variants
            .iter()
            .map(|variant| (variant.chrom.as_str(), variant.pos, variant.id.as_str()))
            .collect::<Vec<_>>(),
        vec![("1", 10, "rs1"), ("1", 20, "rs2"), ("2", 30, "indel1")]
    );
    assert_eq!(metadata.variants[0].a0, "A");
    assert_eq!(metadata.variants[0].a1, "G");
    assert_eq!(metadata.variants[1].a1, "T");
    assert_eq!(metadata.variants[1].alt_allele.as_deref(), Some("T,A"));
    assert!(metadata.capabilities.supports_geno);
    assert!(!metadata.capabilities.supports_haplo);
}

#[test]
fn vcf_haplotype_capability_requires_phased_genotype_evidence() {
    let dir = unique_dir("vcf-phased");
    let path = dir.join("phased.vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1
",
    );

    let metadata = genoio_io::read_vcf_metadata(&path).expect("vcf metadata should parse");

    assert!(metadata.capabilities.supports_geno);
    assert!(metadata.capabilities.supports_haplo);
    assert!(metadata.capabilities.phased);
}

#[test]
fn vcf_haplotype_capability_is_detected_from_records_not_extension() {
    let dir = unique_dir("vcf-phased-extension");
    let path = dir.join("phased.not-vcf");
    write_file(
        &path,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1
",
    );

    let metadata = genoio_io::read_vcf_metadata(&path).expect("vcf metadata should parse");

    assert!(metadata.capabilities.supports_geno);
    assert!(metadata.capabilities.supports_haplo);
    assert!(metadata.capabilities.phased);
}

#[test]
fn plink1_metadata_normalizes_fam_and_bim_records_without_reading_bed_payload() {
    let dir = unique_dir("plink1-metadata");
    let bed = dir.join("tiny.bed");
    let bim = dir.join("tiny.bim");
    let fam = dir.join("tiny.fam");
    fs::write(&bed, [0x6c, 0x1b, 0x01, 0x00]).expect("bed header should be written");
    write_file(
        &bim,
        "\
1 rs1 0 10 G A
1 rs2 0 20 T C
2 indel1 0 30 A AT
",
    );
    write_file(
        &fam,
        "\
F1 S1 0 0 1 -9
F1 S2 S1 0 2 1.5
F2 S3 0 0 0 2.0
",
    );

    let metadata =
        genoio_io::read_plink1_metadata(&bed, &bim, &fam).expect("plink1 metadata should parse");

    assert_eq!(
        metadata
            .samples
            .iter()
            .map(|sample| (sample.fid.as_deref(), sample.iid.as_str()))
            .collect::<Vec<_>>(),
        vec![(Some("F1"), "S1"), (Some("F1"), "S2"), (Some("F2"), "S3")]
    );
    assert_eq!(metadata.samples[1].father.as_deref(), Some("S1"));
    assert_eq!(metadata.samples[1].mother, None);
    assert_eq!(metadata.variants[0].a0, "A");
    assert_eq!(metadata.variants[0].a1, "G");
    assert_eq!(metadata.variants[0].source_a0, "A");
    assert_eq!(metadata.variants[0].source_a1, "G");
    assert!(metadata.capabilities.supports_geno);
    assert!(!metadata.capabilities.supports_haplo);
    assert!(!metadata.capabilities.phased);
}

#[test]
fn plink1_metadata_rejects_malformed_fam_and_bim_lines() {
    let dir = unique_dir("plink1-malformed");
    let bed = dir.join("tiny.bed");
    let bim = dir.join("tiny.bim");
    let fam = dir.join("tiny.fam");
    fs::write(&bed, [0x6c, 0x1b, 0x01, 0x00]).expect("bed header should be written");

    write_file(&bim, "1 rs1 0 10 G A\n");
    write_file(&fam, "F1 S1 0 0 1\n");
    assert!(genoio_io::read_plink1_metadata(&bed, &bim, &fam).is_err());

    write_file(&bim, "1 rs1 0 10 G\n");
    write_file(&fam, "F1 S1 0 0 1 -9\n");
    assert!(genoio_io::read_plink1_metadata(&bed, &bim, &fam).is_err());
}
