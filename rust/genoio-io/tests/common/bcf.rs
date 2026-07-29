// pattern: Imperative Shell

use std::fs;
use std::path::Path;

use noodles_core::Position;
use noodles_vcf::{
    header::record::value::{
        map::{
            format::{Number, Type},
            Contig, Format,
        },
        Map,
    },
    variant::{
        io::Write as _,
        record::samples::keys::key,
        record_buf::{samples::sample::Value, samples::Keys, AlternateBases, Ids, Samples},
    },
};

pub(crate) fn write_genotype_dosage_fixture(path: &Path) {
    let file = fs::File::create(path).expect("test BCF should be created");
    let mut writer = noodles_bcf::io::Writer::new(file);
    let ds_format = Map::<Format>::builder()
        .set_number(Number::Count(1))
        .set_type(Type::Float)
        .set_description("Expected alternate allele dosage")
        .build()
        .expect("DS format should build");
    let header = noodles_vcf::Header::builder()
        .add_contig("1", Map::<Contig>::new())
        .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
        .add_format("DS", ds_format)
        .add_sample_name("s1")
        .add_sample_name("s2")
        .build();

    writer
        .write_header(&header)
        .expect("test BCF header should be written");
    for record in [
        record("rs1", 10, "A", "G", [("0/0", 0.1), ("0/1", 0.9)]),
        record("rs2", 20, "C", "T", [("0/1", 1.2), ("1/1", 1.8)]),
        record("rs3", 30, "G", "A", [("1/1", 1.9), ("0/0", 0.2)]),
    ] {
        writer
            .write_variant_record(&header, &record)
            .expect("test BCF record should be written");
    }
}

pub(crate) fn write_haplotype_fixture(path: &Path) {
    let file = fs::File::create(path).expect("test BCF should be created");
    let mut writer = noodles_bcf::io::Writer::new(file);
    let header = noodles_vcf::Header::builder()
        .add_contig("1", Map::<Contig>::new())
        .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
        .add_sample_name("s1")
        .add_sample_name("s2")
        .build();

    writer
        .write_header(&header)
        .expect("test BCF header should be written");
    for record in [
        gt_record("rs1", 10, "A", "G", ["0|1", "1|0"]),
        gt_record("rs2", 20, "C", "T", ["1|1", "0|0"]),
        gt_record("rs3", 30, "G", "A", ["0|0", "0|1"]),
    ] {
        writer
            .write_variant_record(&header, &record)
            .expect("test BCF record should be written");
    }
}

pub(crate) fn write_haplotype_filter_fixture(path: &Path) {
    let file = fs::File::create(path).expect("test BCF should be created");
    let mut writer = noodles_bcf::io::Writer::new(file);
    let header = noodles_vcf::Header::builder()
        .add_contig("1", Map::<Contig>::new())
        .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
        .add_sample_name("s1")
        .add_sample_name("s2")
        .build();

    writer
        .write_header(&header)
        .expect("test BCF header should be written");
    for record in [
        gt_record("retained", 10, "A", "G", ["0|1", "1|0"]),
        gt_record("unphased_drop", 20, "C", "T", ["0/0", "0/0"]),
    ] {
        writer
            .write_variant_record(&header, &record)
            .expect("test BCF record should be written");
    }
}

fn record(
    id: &str,
    pos: usize,
    reference_bases: &str,
    alternate_bases: &str,
    calls: [(&str, f32); 2],
) -> noodles_vcf::variant::RecordBuf {
    let ids: Ids = [id.to_owned()].into_iter().collect();
    let keys: Keys = [String::from(key::GENOTYPE), "DS".to_owned()]
        .into_iter()
        .collect();
    let samples = Samples::new(
        keys,
        calls
            .into_iter()
            .map(|(gt, ds)| vec![Some(Value::from(gt)), Some(Value::from(ds))])
            .collect(),
    );

    noodles_vcf::variant::RecordBuf::builder()
        .set_reference_sequence_name("1")
        .set_variant_start(Position::try_from(pos).expect("position should be valid"))
        .set_ids(ids)
        .set_reference_bases(reference_bases)
        .set_alternate_bases(AlternateBases::from(vec![alternate_bases.to_owned()]))
        .set_samples(samples)
        .build()
}

fn gt_record(
    id: &str,
    pos: usize,
    reference_bases: &str,
    alternate_bases: &str,
    calls: [&str; 2],
) -> noodles_vcf::variant::RecordBuf {
    let ids: Ids = [id.to_owned()].into_iter().collect();
    let keys: Keys = [String::from(key::GENOTYPE)].into_iter().collect();
    let samples = Samples::new(
        keys,
        calls
            .into_iter()
            .map(|gt| vec![Some(Value::from(gt))])
            .collect(),
    );

    noodles_vcf::variant::RecordBuf::builder()
        .set_reference_sequence_name("1")
        .set_variant_start(Position::try_from(pos).expect("position should be valid"))
        .set_ids(ids)
        .set_reference_bases(reference_bases)
        .set_alternate_bases(AlternateBases::from(vec![alternate_bases.to_owned()]))
        .set_samples(samples)
        .build()
}
