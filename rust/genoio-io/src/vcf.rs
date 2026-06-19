// pattern: Imperative Shell

use std::path::{Path, PathBuf};

use genoio_core::{
    DenseGenotypeMatrix, GenoioError, MetadataOutput, SampleRecord, SourceCapabilities,
    SparseGenotypeMatrix, VariantFilter, VariantWindow,
};
use noodles_bcf as bcf;
use noodles_vcf as noodles;
use noodles_vcf::variant::record::samples::{
    keys::key,
    series::{value::genotype::Phasing as NoodlesGenotypePhasing, Value as NoodlesSampleValue},
};

use self::fast::metadata_variant_record_from_variant_record;
use crate::error::Result;

mod bcf_fast;
mod fast;

/// Read VCF/BCF sample and variant metadata without returning genotypes.
pub fn read_vcf_metadata(path: &Path) -> Result<MetadataOutput> {
    if is_bcf_path(path) {
        return read_bcf_metadata(path);
    }
    fast::read_vcf_metadata(path)
}

fn read_bcf_metadata(path: &Path) -> Result<MetadataOutput> {
    let file = std::fs::File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf open error: {error}")))?;
    let mut reader = bcf::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|error| GenoioError::invalid_source(path, format!("bcf header error: {error}")))?;
    let samples = sample_records_from_noodles_header(&header);

    let mut variants = Vec::new();
    let mut has_phased_genotype_evidence = false;
    // Reuse noodles' lazy BCF record buffer so metadata scans do not allocate a
    // full RecordBuf for each variant before genotype decoding exists.
    let mut record = bcf::Record::default();
    loop {
        let n = reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("bcf record error: {error}"))
        })?;
        if n == 0 {
            break;
        }

        if !has_phased_genotype_evidence
            && noodles_record_has_phased_genotype(path, &header, &record)?
        {
            has_phased_genotype_evidence = true;
        }
        variants.push(metadata_variant_record_from_variant_record(
            path, &header, &record,
        )?);
    }

    let capabilities = if has_phased_genotype_evidence {
        SourceCapabilities::phased_genotypes()
    } else {
        SourceCapabilities::genotype_only()
    };

    Ok(MetadataOutput {
        samples,
        variants,
        capabilities,
    })
}

/// Read retained VCF/BCF diploid genotypes as a dense sample-by-variant matrix.
pub fn read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_dense_windowed(path, requested_samples, variant_filter, None, false)
}

/// Read retained VCF/BCF diploid genotypes as dense values over an optional block window.
pub fn read_vcf_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_dense_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
        None,
    )
}

/// Read retained VCF/BCF diploid genotypes with optional BGZF decompression threads.
pub fn read_vcf_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return read_bcf_dense_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            matrix_only,
        );
    }

    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) =
            fast::try_empty_vcf_dense(path, requested_samples, matrix_only, threads)?
        {
            return Ok(matrix);
        }
        return fast_path_declined(path, "empty dense read");
    }

    // Concrete region predicates can be pushed into the tabix/CSI index. More
    // complex region expressions are evaluated during a permissive full scan.
    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            if let Some(matrix) = fast::try_read_vcf_dense_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                matrix_only,
                threads,
            )? {
                return Ok(matrix);
            }
            return fast_path_declined(path, "indexed dense read");
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    if let Some(matrix) = fast::try_read_vcf_dense(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
        threads,
    )? {
        return Ok(matrix);
    }

    fast_path_declined(path, "dense read")
}

fn read_bcf_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    bcf_fast::read_dense_windowed(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
    )
}

/// Read retained VCF/BCF FORMAT/DS values as dense sample-by-variant dosages.
pub fn read_vcf_dosage_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_dosage_dense_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
        None,
    )
}

/// Read retained VCF/BCF FORMAT/DS values with optional BGZF decompression threads.
pub fn read_vcf_dosage_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return read_bcf_dosage_dense_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            matrix_only,
        );
    }

    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) =
            fast::try_empty_vcf_dense(path, requested_samples, matrix_only, threads)?
        {
            return Ok(matrix);
        }
        return fast_path_declined(path, "empty dosage dense read");
    }

    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            if let Some(matrix) = fast::try_read_vcf_dosage_dense_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                matrix_only,
                threads,
            )? {
                return Ok(matrix);
            }
            return fast_path_declined(path, "indexed dosage dense read");
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    if let Some(matrix) = fast::try_read_vcf_dosage_dense(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
        threads,
    )? {
        return Ok(matrix);
    }

    fast_path_declined(path, "dosage dense read")
}

fn read_bcf_dosage_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    bcf_fast::read_dosage_dense_windowed(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
    )
}

/// Read retained VCF/BCF diploid genotypes as sparse CSC.
pub fn read_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_vcf_sparse_windowed(path, requested_samples, variant_filter, None)
}

/// Read retained VCF/BCF diploid genotypes as sparse CSC over an optional block window.
pub fn read_vcf_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    read_vcf_sparse_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        None,
    )
}

/// Read retained VCF/BCF diploid genotypes as sparse CSC with optional BGZF decompression threads.
pub fn read_vcf_sparse_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return read_bcf_sparse_windowed(path, requested_samples, variant_filter, variant_window);
    }

    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) = fast::try_empty_vcf_sparse(path, requested_samples, threads)? {
            return Ok(matrix);
        }
        return fast_path_declined(path, "empty sparse read");
    }

    // Keep sparse and dense region behavior identical so both paths retain the
    // same variants and fail the same way for unindexed compressed inputs.
    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            if let Some(matrix) = fast::try_read_vcf_sparse_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                threads,
            )? {
                return Ok(matrix);
            }
            return fast_path_declined(path, "indexed sparse read");
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    if let Some(matrix) = fast::try_read_vcf_sparse(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        threads,
    )? {
        return Ok(matrix);
    }

    fast_path_declined(path, "sparse read")
}

fn read_bcf_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    bcf_fast::read_sparse_windowed(path, requested_samples, variant_filter, variant_window)
}

/// Read phased VCF/BCF diploid genotypes as dense haplotype rows.
pub fn read_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_haplotypes_dense_windowed(path, requested_samples, variant_filter, None, false)
}

/// Read phased VCF/BCF diploid genotypes as dense haplotype rows over a block window.
pub fn read_vcf_haplotypes_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    read_vcf_haplotypes_dense_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
        None,
    )
}

/// Read phased VCF/BCF diploid genotypes as dense haplotype rows with optional BGZF decompression threads.
pub fn read_vcf_haplotypes_dense_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return read_bcf_haplotypes_dense_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
            matrix_only,
        );
    }

    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) =
            fast::try_empty_vcf_haplotypes_dense(path, requested_samples, matrix_only, threads)?
        {
            return Ok(matrix);
        }
        return fast_path_declined(path, "empty haplotype dense read");
    }

    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            if let Some(matrix) = fast::try_read_vcf_haplotypes_dense_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                matrix_only,
                threads,
            )? {
                return Ok(matrix);
            }
            return fast_path_declined(path, "indexed haplotype dense read");
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    if let Some(matrix) = fast::try_read_vcf_haplotypes_dense(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
        threads,
    )? {
        return Ok(matrix);
    }

    fast_path_declined(path, "haplotype dense read")
}

fn read_bcf_haplotypes_dense_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    bcf_fast::read_haplotypes_dense_windowed(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        matrix_only,
    )
}

/// Read phased VCF/BCF diploid genotypes as sparse haplotype rows.
pub fn read_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_vcf_haplotypes_sparse_windowed(path, requested_samples, variant_filter, None)
}

/// Read phased VCF/BCF diploid genotypes as sparse haplotype rows over a block window.
pub fn read_vcf_haplotypes_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    read_vcf_haplotypes_sparse_windowed_with_threads(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        None,
    )
}

/// Read phased VCF/BCF diploid genotypes as sparse haplotype rows with optional BGZF decompression threads.
pub fn read_vcf_haplotypes_sparse_windowed_with_threads(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    validate_threaded_read_support(path, threads)?;
    if is_bcf_path(path) {
        return read_bcf_haplotypes_sparse_windowed(
            path,
            requested_samples,
            variant_filter,
            variant_window,
        );
    }

    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) =
            fast::try_empty_vcf_haplotypes_sparse(path, requested_samples, threads)?
        {
            return Ok(matrix);
        }
        return fast_path_declined(path, "empty haplotype sparse read");
    }

    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            if let Some(matrix) = fast::try_read_vcf_haplotypes_sparse_indexed(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                threads,
            )? {
                return Ok(matrix);
            }
            return fast_path_declined(path, "indexed haplotype sparse read");
        }
    }

    reject_unindexed_compressed_region(path, variant_filter)?;
    if let Some(matrix) = fast::try_read_vcf_haplotypes_sparse(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        threads,
    )? {
        return Ok(matrix);
    }

    fast_path_declined(path, "haplotype sparse read")
}

fn read_bcf_haplotypes_sparse_windowed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    bcf_fast::read_haplotypes_sparse_windowed(
        path,
        requested_samples,
        variant_filter,
        variant_window,
    )
}

fn fast_path_declined<T>(path: &Path, context: &str) -> Result<T> {
    Err(GenoioError::internal_contract(format!(
        "noodles VCF fast path declined supported {context} for {}",
        path.display()
    )))
}

fn validate_threaded_read_support(path: &Path, threads: Option<usize>) -> Result<()> {
    let Some(threads) = threads else {
        return Ok(());
    };
    if threads == 0 {
        return Err(GenoioError::invalid_source(
            path,
            "vcf thread count must be greater than zero",
        ));
    }
    if is_bcf_path(path) {
        return Err(GenoioError::invalid_source(
            path,
            "threaded BCF reads are not supported by the noodles backend",
        ));
    }
    if !is_compressed_vcf(path) {
        return Err(GenoioError::invalid_source(
            path,
            "threaded reads are only supported for compressed VCF sources",
        ));
    }
    Ok(())
}

fn sample_records_from_noodles_header(header: &noodles::Header) -> Vec<SampleRecord> {
    header
        .sample_names()
        .iter()
        .map(|sample| plain_sample_record(sample.to_string()))
        .collect()
}

fn plain_sample_record(iid: String) -> SampleRecord {
    SampleRecord {
        fid: None,
        iid,
        father: None,
        mother: None,
        sex: None,
        phenotype: None,
        source_sample_index: None,
        haplotype_index: None,
    }
}

fn noodles_record_has_phased_genotype<R>(
    path: &Path,
    header: &noodles::Header,
    record: &R,
) -> Result<bool>
where
    R: noodles::variant::Record + ?Sized,
{
    let samples = record.samples().map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf samples error: {error}"))
    })?;
    let Some(gt_series_result) = samples.select(header, key::GENOTYPE) else {
        return Ok(false);
    };
    let gt_series = gt_series_result.map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf genotype series error: {error}"))
    })?;

    for value_result in gt_series.iter(header) {
        let Some(NoodlesSampleValue::Genotype(genotype)) = value_result.map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf genotype value error: {error}"))
        })?
        else {
            continue;
        };

        for allele_result in genotype.iter() {
            let (_, phasing) = allele_result.map_err(|error| {
                GenoioError::invalid_source(path, format!("vcf genotype allele error: {error}"))
            })?;
            if phasing == NoodlesGenotypePhasing::Phased {
                return Ok(true);
            }
        }
    }

    Ok(false)
}

fn finite_qual(qual: f32) -> Option<f32> {
    qual.is_finite().then_some(qual)
}

fn reject_unindexed_compressed_region(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
) -> Result<()> {
    if variant_filter
        .and_then(VariantFilter::concrete_region_pushdown)
        .is_none()
        || !is_compressed_vcf(path)
    {
        return Ok(());
    }
    if has_vcf_index(path) {
        return Ok(());
    }
    Err(GenoioError::invalid_source(
        path,
        "region filter on compressed VCF requires an index",
    ))
}

fn has_vcf_index(path: &Path) -> bool {
    companion_index_path(path, "tbi").exists() || companion_index_path(path, "csi").exists()
}

fn companion_index_path(path: &Path, index_extension: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}", path.to_string_lossy(), index_extension))
}

fn is_compressed_vcf(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| matches!(extension, "gz" | "bgz"))
}

fn is_bcf_path(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bcf"))
}

fn haplotype_sample_records(
    samples: &[SampleRecord],
    source_indices: &[usize],
) -> Vec<SampleRecord> {
    let mut haplotype_samples = Vec::with_capacity(samples.len() * 2);
    for (sample, source_index) in samples.iter().zip(source_indices) {
        let original_source_index = sample.source_sample_index.unwrap_or(*source_index);
        for haplotype_index in 0..2 {
            let mut haplotype_sample = sample.clone();
            haplotype_sample.source_sample_index = Some(original_source_index);
            haplotype_sample.haplotype_index = Some(haplotype_index);
            haplotype_samples.push(haplotype_sample);
        }
    }
    haplotype_samples
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;
    use noodles_core::Position;
    use noodles_vcf::{
        self as noodles_vcf,
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

    fn write_test_bcf(path: &Path) {
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
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs1", 10, "A", &["G"], ["0/0", "0/1"]),
            )
            .expect("first BCF record should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs2", 20, "C", &["T"], ["1/1", "./."]),
            )
            .expect("second BCF record should be written");
    }

    fn write_test_bcf_with_ds(path: &Path) {
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
        writer
            .write_variant_record(
                &header,
                &bcf_ds_test_record(
                    "rs1",
                    10,
                    "A",
                    &["G"],
                    [("0/0", Some(0.2)), ("0/1", Some(1.4))],
                ),
            )
            .expect("first BCF record should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_ds_test_record("rs2", 20, "C", &["T"], [("1/1", Some(2.0)), ("./.", None)]),
            )
            .expect("second BCF record should be written");
    }

    fn write_test_bcf_phased(path: &Path) {
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
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs1", 10, "A", &["G"], ["0|1", "1|1"]),
            )
            .expect("first phased BCF record should be written");
        writer
            .write_variant_record(
                &header,
                &bcf_test_record("rs2", 20, "C", &["T"], ["1|0", "0|0"]),
            )
            .expect("second phased BCF record should be written");
    }

    fn bcf_test_record(
        id: &str,
        pos: usize,
        reference_bases: &str,
        alternate_bases: &[&str],
        genotypes: [&str; 2],
    ) -> noodles_vcf::variant::RecordBuf {
        let ids: Ids = [id.to_string()].into_iter().collect();
        let keys: Keys = [String::from(key::GENOTYPE)].into_iter().collect();
        let samples = Samples::new(
            keys,
            genotypes
                .into_iter()
                .map(|gt| vec![Some(Value::from(gt))])
                .collect(),
        );

        noodles_vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("1")
            .set_variant_start(Position::try_from(pos).expect("position should be valid"))
            .set_ids(ids)
            .set_reference_bases(reference_bases)
            .set_alternate_bases(AlternateBases::from(
                alternate_bases
                    .iter()
                    .copied()
                    .map(str::to_string)
                    .collect::<Vec<_>>(),
            ))
            .set_samples(samples)
            .build()
    }

    fn bcf_ds_test_record(
        id: &str,
        pos: usize,
        reference_bases: &str,
        alternate_bases: &[&str],
        calls: [(&str, Option<f32>); 2],
    ) -> noodles_vcf::variant::RecordBuf {
        let ids: Ids = [id.to_string()].into_iter().collect();
        let keys: Keys = [String::from(key::GENOTYPE), "DS".to_string()]
            .into_iter()
            .collect();
        let samples = Samples::new(
            keys,
            calls
                .into_iter()
                .map(|(gt, ds)| vec![Some(Value::from(gt)), ds.map(Value::from)])
                .collect(),
        );

        noodles_vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("1")
            .set_variant_start(Position::try_from(pos).expect("position should be valid"))
            .set_ids(ids)
            .set_reference_bases(reference_bases)
            .set_alternate_bases(AlternateBases::from(
                alternate_bases
                    .iter()
                    .copied()
                    .map(str::to_string)
                    .collect::<Vec<_>>(),
            ))
            .set_samples(samples)
            .build()
    }

    #[test]
    fn bcf_dense_gt_reads_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let dense = read_bcf_dense_windowed(file.path(), None, None, None, false)
            .expect("BCF dense GT should decode");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 2);
        assert_eq!(dense.values, vec![0.0, 2.0, 1.0, 0.0]);
        assert_eq!(dense.missing_mask, vec![false, false, false, true]);
        assert_eq!(
            dense
                .variants
                .iter()
                .map(|variant| variant.id.as_str())
                .collect::<Vec<_>>(),
            vec!["rs1", "rs2"]
        );
    }

    #[test]
    fn bcf_dense_gt_rejects_threaded_reads() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let error =
            read_vcf_dense_windowed_with_threads(file.path(), None, None, None, false, Some(2))
                .expect_err("BCF should reject explicit thread count");

        assert!(error
            .to_string()
            .contains("threaded BCF reads are not supported"));
    }

    #[test]
    fn bcf_dense_gt_applies_retained_windows() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());

        let dense = read_bcf_dense_windowed(
            file.path(),
            None,
            None,
            Some(VariantWindow { start: 1, len: 1 }),
            false,
        )
        .expect("BCF dense GT window should decode");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 1);
        assert_eq!(dense.values, vec![2.0, 0.0]);
        assert_eq!(dense.missing_mask, vec![false, true]);
        assert_eq!(dense.variants[0].id, "rs2");
    }

    #[test]
    fn bcf_dense_gt_filters_stats_after_sample_selection() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());
        let samples = vec!["s2".to_string()];
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "mac",
            "params": {"min": 1}
        }))
        .expect("filter should parse");

        let dense =
            read_bcf_dense_windowed(file.path(), Some(&samples), Some(&filter), None, false)
                .expect("BCF dense GT filter should decode");

        assert_eq!(dense.n_samples, 1);
        assert_eq!(dense.n_variants, 1);
        assert_eq!(dense.values, vec![1.0]);
        assert_eq!(dense.missing_mask, vec![false]);
        assert_eq!(dense.samples[0].iid, "s2");
        assert_eq!(dense.variants[0].id, "rs1");
        assert_eq!(dense.variants[0].mac, Some(1));
        assert_eq!(dense.variants[0].n_called, Some(1));
    }

    #[test]
    fn bcf_dense_ds_reads_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_with_ds(file.path());

        let dense = read_bcf_dosage_dense_windowed(file.path(), None, None, None, false)
            .expect("BCF dense DS should decode");

        assert_eq!(dense.n_samples, 2);
        assert_eq!(dense.n_variants, 2);
        assert_eq!(dense.values, vec![0.2, 2.0, 1.4, 0.0]);
        assert_eq!(dense.missing_mask, vec![false, false, false, true]);
        assert_eq!(
            dense
                .variants
                .iter()
                .map(|variant| variant.id.as_str())
                .collect::<Vec<_>>(),
            vec!["rs1", "rs2"]
        );
    }

    #[test]
    fn bcf_sparse_gt_reads_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf(file.path());
        let filter = VariantFilter::from_json_value(serde_json::json!({
            "op": "predicate",
            "name": "mac",
            "params": {"min": 1}
        }))
        .expect("filter should parse");

        let sparse = read_bcf_sparse_windowed(file.path(), None, Some(&filter), None)
            .expect("BCF sparse GT should decode");

        assert_eq!(sparse.n_rows, 2);
        assert_eq!(sparse.n_cols, 1);
        assert_eq!(sparse.indptr, vec![0, 1]);
        assert_eq!(sparse.indices, vec![1]);
        assert_eq!(sparse.data, vec![1.0]);
        assert_eq!(sparse.variants[0].id, "rs1");
    }

    #[test]
    fn bcf_dense_haplotypes_read_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_phased(file.path());

        let dense = read_bcf_haplotypes_dense_windowed(file.path(), None, None, None, false)
            .expect("BCF dense haplotypes should decode");

        assert_eq!(dense.n_samples, 4);
        assert_eq!(dense.n_variants, 2);
        assert_eq!(dense.values, vec![0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]);
        assert_eq!(dense.missing_mask, vec![false; 8]);
        assert_eq!(dense.samples[0].iid, "s1");
        assert_eq!(dense.samples[0].haplotype_index, Some(0));
        assert_eq!(dense.samples[1].haplotype_index, Some(1));
    }

    #[test]
    fn bcf_sparse_haplotypes_read_lazy_noodles_records() {
        let file = tempfile::Builder::new()
            .suffix(".bcf")
            .tempfile()
            .expect("temp BCF should be created");
        write_test_bcf_phased(file.path());

        let sparse = read_bcf_haplotypes_sparse_windowed(file.path(), None, None, None)
            .expect("BCF sparse haplotypes should decode");

        assert_eq!(sparse.n_rows, 4);
        assert_eq!(sparse.n_cols, 2);
        assert_eq!(sparse.indptr, vec![0, 1, 2]);
        assert_eq!(sparse.indices, vec![0, 0]);
        assert_eq!(sparse.data, vec![1.0, 1.0]);
        assert!(sparse.variants[0].flipped);
        assert!(!sparse.variants[1].flipped);
    }
}
