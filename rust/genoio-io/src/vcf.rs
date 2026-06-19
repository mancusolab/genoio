// pattern: Imperative Shell

use std::ffi::CString;
use std::path::{Path, PathBuf};

use genoio_core::{
    append_sparse_column, attach_variant_stats, compute_dosage_variant_stats,
    flip_haplotype_values_to_minor_allele, flip_values_to_minor_allele,
    reject_sparse_missing_values, select_samples_source_order, variant_stats_from_counts,
    DenseGenotypeMatrix, GenoioError, MetadataOutput, PartialFilterDecision, RegionPredicate,
    SampleRecord, SourceCapabilities, SparseGenotypeMatrix, VariantFilter, VariantRecord,
    VariantStats, VariantWindow,
};
use noodles_bcf as bcf;
use noodles_vcf as noodles;
use noodles_vcf::variant::record::samples::{
    keys::key,
    series::{value::genotype::Phasing as NoodlesGenotypePhasing, Value as NoodlesSampleValue},
};
use rust_htslib::bcf::{
    record::{GenotypeAllele, Numeric},
    IndexedReader, Read, Reader,
};
use rust_htslib::htslib;

use self::fast::metadata_variant_record_from_variant_record;
use crate::error::Result;
use crate::matrix::{
    empty_sparse_matrix, finish_dense_matrix, finish_variant_major_dense_matrix, DenseMatrixParts,
    VariantMajorDenseParts,
};
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

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
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) =
            fast::try_empty_vcf_dense(path, requested_samples, matrix_only, threads)?
        {
            return Ok(matrix);
        }
        let (reader, _original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
        return empty_vcf_dense(path, reader.header(), requested_samples, matrix_only);
    }

    // Region filters can be pushed into htslib only when the expression shape
    // is a concrete safe region and the compressed source has an index.
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
            return read_indexed_vcf_dense(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                matrix_only,
                threads,
            );
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

    let (mut reader, original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
    read_vcf_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        original_source_indices.as_deref(),
        matrix_only,
        &mut reader,
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
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) =
            fast::try_empty_vcf_dense(path, requested_samples, matrix_only, threads)?
        {
            return Ok(matrix);
        }
        let (reader, _original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
        return empty_vcf_dense(path, reader.header(), requested_samples, matrix_only);
    }

    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            return read_indexed_vcf_dosage_dense(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                matrix_only,
                threads,
            );
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

    let (mut reader, original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
    read_vcf_dosage_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        original_source_indices.as_deref(),
        matrix_only,
        &mut reader,
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
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) = fast::try_empty_vcf_sparse(path, requested_samples, threads)? {
            return Ok(matrix);
        }
        let (reader, _original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
        return empty_vcf_sparse(path, reader.header(), requested_samples);
    }

    // Keep sparse and dense region behavior identical so both paths retain the
    // same variants and fail the same way for unindexed compressed inputs.
    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            return read_indexed_vcf_sparse(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                threads,
            );
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

    let (mut reader, original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
    read_vcf_sparse_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        original_source_indices.as_deref(),
        &mut reader,
    )
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
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) =
            fast::try_empty_vcf_haplotypes_dense(path, requested_samples, matrix_only, threads)?
        {
            return Ok(matrix);
        }
        let (reader, _original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
        return empty_vcf_haplotypes_dense(path, reader.header(), requested_samples, matrix_only);
    }

    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            return read_indexed_vcf_haplotypes_dense(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                matrix_only,
                threads,
            );
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

    let (mut reader, original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
    read_vcf_haplotypes_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        original_source_indices.as_deref(),
        matrix_only,
        &mut reader,
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
    if variant_filter.is_some_and(VariantFilter::is_always_false) {
        if let Some(matrix) =
            fast::try_empty_vcf_haplotypes_sparse(path, requested_samples, threads)?
        {
            return Ok(matrix);
        }
        let (reader, _original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
        return empty_vcf_haplotypes_sparse(path, reader.header(), requested_samples);
    }

    if let Some(region) = variant_filter.and_then(VariantFilter::concrete_region_pushdown) {
        if has_vcf_index(path) {
            return read_indexed_vcf_haplotypes_sparse(
                path,
                requested_samples,
                variant_filter,
                variant_window,
                &region,
                threads,
            );
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

    let (mut reader, original_source_indices) = open_vcf_reader(path, threads, requested_samples)?;
    read_vcf_haplotypes_sparse_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        original_source_indices.as_deref(),
        &mut reader,
    )
}

fn open_vcf_reader(
    path: &Path,
    threads: Option<usize>,
    requested_samples: Option<&[String]>,
) -> Result<(Reader, Option<Vec<usize>>)> {
    let mut reader = Reader::from_path(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf reader error: {error}")))?;
    set_vcf_reader_threads(path, &mut reader, threads, "vcf reader")?;
    let original_source_indices = set_vcf_reader_samples(path, &mut reader, requested_samples)?;
    Ok((reader, original_source_indices))
}

fn set_vcf_reader_samples(
    path: &Path,
    reader: &mut Reader,
    requested_samples: Option<&[String]>,
) -> Result<Option<Vec<usize>>> {
    let Some(requested_samples) = requested_samples else {
        return Ok(None);
    };

    let source_samples = sample_records_from_header(reader.header());
    let selection = select_samples_source_order(&source_samples, Some(requested_samples), path)?;
    if selection.samples.len() == source_samples.len() {
        return Ok(None);
    }
    let sample_list = selection
        .samples
        .iter()
        .map(|sample| sample.iid.as_str())
        .collect::<Vec<_>>()
        .join(",");
    let sample_list = CString::new(sample_list).map_err(|_| {
        GenoioError::invalid_source(path, "vcf sample IDs must not contain NUL bytes")
    })?;
    let status =
        unsafe { htslib::bcf_hdr_set_samples(reader.header().as_ptr(), sample_list.as_ptr(), 0) };
    match status {
        0 => Ok(Some(selection.source_indices)),
        positive if positive > 0 => Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf sample subset rejected sample at list index {}",
                positive - 1
            ),
        )),
        _ => Err(GenoioError::invalid_source(
            path,
            "vcf sample subset setup error",
        )),
    }
}

fn open_indexed_vcf_reader(path: &Path, threads: Option<usize>) -> Result<IndexedReader> {
    let mut reader = IndexedReader::from_path(path).map_err(|error| {
        GenoioError::invalid_source(path, format!("indexed vcf reader error: {error}"))
    })?;
    set_vcf_reader_threads(path, &mut reader, threads, "indexed vcf reader")?;
    Ok(reader)
}

fn set_vcf_reader_threads<R: Read>(
    path: &Path,
    reader: &mut R,
    threads: Option<usize>,
    label: &str,
) -> Result<()> {
    let Some(threads) = threads else {
        return Ok(());
    };
    if threads == 0 {
        return Err(GenoioError::invalid_source(
            path,
            "vcf thread count must be greater than zero",
        ));
    }
    reader.set_threads(threads).map_err(|error| {
        GenoioError::invalid_source(path, format!("{label} thread setup error: {error}"))
    })
}

fn read_indexed_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    let mut reader = open_indexed_vcf_reader(path, threads)?;
    let header = reader.header().clone();
    let rid = match header.name2rid(region.chrom.as_bytes()) {
        Ok(rid) => rid,
        Err(_) => return empty_vcf_dense(path, &header, requested_samples, matrix_only),
    };
    reader
        .fetch(
            rid,
            u64::from(region.start - 1),
            Some(u64::from(region.end - 1)),
        )
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf region fetch error: {error}"))
        })?;

    read_vcf_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        None,
        matrix_only,
        &mut reader,
    )
}

fn read_indexed_vcf_dosage_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    let mut reader = open_indexed_vcf_reader(path, threads)?;
    let header = reader.header().clone();
    let rid = match header.name2rid(region.chrom.as_bytes()) {
        Ok(rid) => rid,
        Err(_) => return empty_vcf_dense(path, &header, requested_samples, matrix_only),
    };
    reader
        .fetch(
            rid,
            u64::from(region.start - 1),
            Some(u64::from(region.end - 1)),
        )
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf region fetch error: {error}"))
        })?;

    read_vcf_dosage_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        None,
        matrix_only,
        &mut reader,
    )
}

fn read_indexed_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    let mut reader = open_indexed_vcf_reader(path, threads)?;
    let header = reader.header().clone();
    let rid = match header.name2rid(region.chrom.as_bytes()) {
        Ok(rid) => rid,
        Err(_) => return empty_vcf_sparse(path, &header, requested_samples),
    };
    reader
        .fetch(
            rid,
            u64::from(region.start - 1),
            Some(u64::from(region.end - 1)),
        )
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf region fetch error: {error}"))
        })?;

    read_vcf_sparse_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        None,
        &mut reader,
    )
}

fn read_indexed_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    let mut reader = open_indexed_vcf_reader(path, threads)?;
    let header = reader.header().clone();
    let rid = match header.name2rid(region.chrom.as_bytes()) {
        Ok(rid) => rid,
        Err(_) => {
            return empty_vcf_haplotypes_dense(path, &header, requested_samples, matrix_only);
        }
    };
    reader
        .fetch(
            rid,
            u64::from(region.start - 1),
            Some(u64::from(region.end - 1)),
        )
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf region fetch error: {error}"))
        })?;

    read_vcf_haplotypes_dense_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        None,
        matrix_only,
        &mut reader,
    )
}

fn read_indexed_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    let mut reader = open_indexed_vcf_reader(path, threads)?;
    let header = reader.header().clone();
    let rid = match header.name2rid(region.chrom.as_bytes()) {
        Ok(rid) => rid,
        Err(_) => return empty_vcf_haplotypes_sparse(path, &header, requested_samples),
    };
    reader
        .fetch(
            rid,
            u64::from(region.start - 1),
            Some(u64::from(region.end - 1)),
        )
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf region fetch error: {error}"))
        })?;

    read_vcf_haplotypes_sparse_records(
        path,
        requested_samples,
        variant_filter,
        variant_window,
        None,
        &mut reader,
    )
}

fn read_vcf_dense_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    original_source_indices: Option<&[usize]>,
    matrix_only: bool,
    reader: &mut R,
) -> Result<DenseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let mut selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    apply_original_source_indices(&mut selection.samples, original_source_indices);
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
    let mut n_variants = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    for record_result in reader.records() {
        // Check before pulling another record; otherwise block reads still pay
        // to scan one extra variant after the requested retained window.
        if retention.window_is_satisfied() {
            break;
        }
        let record = record_result.map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf record error: {error}"))
        })?;
        let variant = if !matrix_only || variant_filter.is_some() {
            Some(variant_record_from_record(path, &header, &record)?)
        } else {
            None
        };
        let partial_decision = match (variant_filter, variant.as_ref()) {
            (Some(filter), Some(variant)) => filter.partial_decision(variant),
            (Some(_), None) => {
                return Err(GenoioError::internal_contract(
                    "vcf filter requires variant metadata",
                ));
            }
            (None, _) => PartialFilterDecision::Accept,
        };
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        validate_dense_biallelic_record(path, &record)?;
        let (current_values, current_missing, stats) = if needs_genotype_decision {
            let decoded =
                decode_diploid_genotype_with_stats(path, &record, &selection.source_indices)?;
            (decoded.values, decoded.missing, Some(decoded.stats))
        } else {
            let decoded = decode_diploid_genotype_values(path, &record, &selection.source_indices)?;
            (decoded.values, decoded.missing, None)
        };
        if needs_genotype_decision {
            let variant = variant.as_ref().ok_or_else(|| {
                GenoioError::internal_contract("vcf filter requires variant metadata")
            })?;
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }
        if !matrix_only {
            let mut variant = variant.ok_or_else(|| {
                GenoioError::internal_contract("vcf metadata output requires variant metadata")
            })?;
            if let Some(stats) = stats {
                attach_variant_stats(&mut variant, stats);
            }
            variants.push(variant);
        }
        n_variants += 1;
        variant_major_values.extend(current_values);
        variant_major_missing.extend(current_missing);
    }

    let n_samples = selection.samples.len();
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

fn read_vcf_dosage_dense_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    original_source_indices: Option<&[usize]>,
    matrix_only: bool,
    reader: &mut R,
) -> Result<DenseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let mut selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    apply_original_source_indices(&mut selection.samples, original_source_indices);
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
    let mut n_variants = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    for record_result in reader.records() {
        if retention.window_is_satisfied() {
            break;
        }
        let record = record_result.map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf record error: {error}"))
        })?;
        let mut variant = if !matrix_only || variant_filter.is_some() {
            Some(variant_record_from_record(path, &header, &record)?)
        } else {
            None
        };
        let partial_decision = match (variant_filter, variant.as_ref()) {
            (Some(filter), Some(variant)) => filter.partial_decision(variant),
            (Some(_), None) => {
                return Err(GenoioError::internal_contract(
                    "vcf filter requires variant metadata",
                ));
            }
            (None, _) => PartialFilterDecision::Accept,
        };
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        validate_dense_biallelic_record(path, &record)?;
        let decoded = decode_ds_dosage_record(path, &record, &selection.source_indices)?;
        if needs_genotype_decision {
            let stats = compute_dosage_variant_stats(&decoded.values, &decoded.missing)?;
            let filter_variant = variant.as_ref().ok_or_else(|| {
                GenoioError::internal_contract("vcf filter requires variant metadata")
            })?;
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(filter_variant, Some(&stats))),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if !matrix_only {
                let variant = variant.as_mut().ok_or_else(|| {
                    GenoioError::internal_contract("vcf metadata output requires variant metadata")
                })?;
                attach_variant_stats(variant, stats);
            }
        }
        if !matrix_only {
            let variant = variant.ok_or_else(|| {
                GenoioError::internal_contract("vcf metadata output requires variant metadata")
            })?;
            variants.push(variant);
        }
        n_variants += 1;
        variant_major_values.extend(decoded.values);
        variant_major_missing.extend(decoded.missing);
    }

    let n_samples = selection.samples.len();
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples: selection.samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

fn read_vcf_haplotypes_dense_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    original_source_indices: Option<&[usize]>,
    matrix_only: bool,
    reader: &mut R,
) -> Result<DenseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let mut selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    apply_original_source_indices(&mut selection.samples, original_source_indices);
    let mut diagnostics = selection.diagnostics;

    let mut variants = Vec::new();
    let mut variant_major_values = Vec::new();
    let mut variant_major_missing = Vec::new();
    let mut n_variants = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    for record_result in reader.records() {
        // Haplotype reads use the same retained-window semantics as genotype
        // reads, but each retained sample contributes two output rows.
        if retention.window_is_satisfied() {
            break;
        }
        let record = record_result.map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf record error: {error}"))
        })?;
        let variant = if !matrix_only || variant_filter.is_some() {
            Some(variant_record_from_record(path, &header, &record)?)
        } else {
            None
        };
        let partial_decision = match (variant_filter, variant.as_ref()) {
            (Some(filter), Some(variant)) => filter.partial_decision(variant),
            (Some(_), None) => {
                return Err(GenoioError::internal_contract(
                    "vcf filter requires variant metadata",
                ));
            }
            (None, _) => PartialFilterDecision::Accept,
        };
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        validate_dense_biallelic_record(path, &record)?;
        let stats = if needs_genotype_decision {
            let decoded =
                decode_diploid_genotype_with_stats(path, &record, &selection.source_indices)?;
            Some(decoded.stats)
        } else {
            None
        };
        if needs_genotype_decision {
            let variant = variant.as_ref().ok_or_else(|| {
                GenoioError::internal_contract("vcf filter requires variant metadata")
            })?;
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }
        if !matrix_only {
            let mut variant = variant.ok_or_else(|| {
                GenoioError::internal_contract("vcf metadata output requires variant metadata")
            })?;
            if let Some(stats) = stats {
                attach_variant_stats(&mut variant, stats);
            }
            variants.push(variant);
        }
        let decoded = decode_phased_haplotype_record(path, &record, &selection.source_indices)?;
        n_variants += 1;
        variant_major_values.extend(decoded.haplotype_values);
        variant_major_missing.extend(decoded.haplotype_missing);
    }

    let samples = if matrix_only {
        Vec::new()
    } else {
        haplotype_sample_records(&selection.samples, &selection.source_indices)
    };
    let n_samples = selection.samples.len() * 2;
    diagnostics.retained_variants = n_variants;
    finish_variant_major_dense_matrix(
        VariantMajorDenseParts {
            n_samples,
            n_variants,
            variant_major_values,
            variant_major_missing,
            samples,
            variants,
            diagnostics,
        },
        matrix_only,
    )
}

fn read_vcf_haplotypes_sparse_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    original_source_indices: Option<&[usize]>,
    reader: &mut R,
) -> Result<SparseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let mut selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    apply_original_source_indices(&mut selection.samples, original_source_indices);
    let mut diagnostics = selection.diagnostics;

    let samples = haplotype_sample_records(&selection.samples, &selection.source_indices);
    let n_samples = samples.len();
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();
    let mut retention = RetainedVariantState::new(variant_window);
    for record_result in reader.records() {
        // Stop before reading the next record so sparse blocks do not scan the
        // remainder of large VCF/BCF files after the block is filled.
        if retention.window_is_satisfied() {
            break;
        }
        let record = record_result.map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf record error: {error}"))
        })?;
        let mut variant = variant_record_from_record(path, &header, &record)?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        validate_dense_biallelic_record(path, &record)?;
        let stats = if needs_genotype_decision {
            let decoded =
                decode_diploid_genotype_with_stats(path, &record, &selection.source_indices)?;
            Some(decoded.stats)
        } else {
            None
        };
        if needs_genotype_decision {
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        let decoded = decode_phased_haplotype_record(path, &record, &selection.source_indices)?;
        reject_sparse_missing_values(&decoded.haplotype_missing)?;
        let mut haplotype_values = decoded.haplotype_values;
        flip_haplotype_values_to_minor_allele(&mut haplotype_values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &haplotype_values);
        variants.push(variant);
    }

    let n_variants = variants.len();
    diagnostics.retained_variants = n_variants;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        samples,
        variants,
        diagnostics,
    )
}

fn read_vcf_sparse_records<R: Read>(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    original_source_indices: Option<&[usize]>,
    reader: &mut R,
) -> Result<SparseGenotypeMatrix> {
    let header = reader.header().clone();
    let all_samples = sample_records_from_header(&header);
    let mut selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    apply_original_source_indices(&mut selection.samples, original_source_indices);
    let mut diagnostics = selection.diagnostics;

    let n_samples = selection.samples.len();
    let mut indptr = vec![0];
    let mut indices = Vec::new();
    let mut data = Vec::new();
    let mut variants = Vec::new();
    let mut retention = RetainedVariantState::new(variant_window);
    for record_result in reader.records() {
        // Stop before reading the next record so sparse blocks do not scan the
        // remainder of large VCF/BCF files after the block is filled.
        if retention.window_is_satisfied() {
            break;
        }
        let record = record_result.map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf record error: {error}"))
        })?;
        let mut variant = variant_record_from_record(path, &header, &record)?;
        let partial_decision = variant_filter.map_or(PartialFilterDecision::Accept, |filter| {
            filter.partial_decision(&variant)
        });
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);

        validate_dense_biallelic_record(path, &record)?;
        let (mut current_values, current_missing, stats) = if needs_genotype_decision {
            let decoded =
                decode_diploid_genotype_with_stats(path, &record, &selection.source_indices)?;
            (decoded.values, decoded.missing, Some(decoded.stats))
        } else {
            let decoded = decode_diploid_genotype_values(path, &record, &selection.source_indices)?;
            (decoded.values, decoded.missing, None)
        };
        if needs_genotype_decision {
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }
        if let Some(stats) = stats {
            attach_variant_stats(&mut variant, stats);
        }
        reject_sparse_missing_values(&current_missing)?;
        flip_values_to_minor_allele(&mut current_values, &mut variant);
        append_sparse_column(&mut indptr, &mut indices, &mut data, &current_values);
        variants.push(variant);
    }

    let n_variants = variants.len();
    diagnostics.retained_variants = n_variants;
    SparseGenotypeMatrix::new(
        n_samples,
        n_variants,
        indptr,
        indices,
        data,
        selection.samples,
        variants,
        diagnostics,
    )
}

fn empty_vcf_dense(
    path: &Path,
    header: &rust_htslib::bcf::header::HeaderView,
    requested_samples: Option<&[String]>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let all_samples = sample_records_from_header(header);
    let mut selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let n_samples = selection.samples.len();
    selection.diagnostics.retained_variants = 0;
    finish_dense_matrix(
        DenseMatrixParts {
            n_samples,
            n_variants: 0,
            values: Vec::new(),
            missing_mask: Vec::new(),
            samples: selection.samples,
            variants: Vec::new(),
            diagnostics: selection.diagnostics,
        },
        matrix_only,
    )
}

fn empty_vcf_sparse(
    path: &Path,
    header: &rust_htslib::bcf::header::HeaderView,
    requested_samples: Option<&[String]>,
) -> Result<SparseGenotypeMatrix> {
    let all_samples = sample_records_from_header(header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    empty_sparse_matrix(selection.samples, selection.diagnostics)
}

fn empty_vcf_haplotypes_dense(
    path: &Path,
    header: &rust_htslib::bcf::header::HeaderView,
    requested_samples: Option<&[String]>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let all_samples = sample_records_from_header(header);
    let mut selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let n_samples = selection.samples.len() * 2;
    let samples = if matrix_only {
        Vec::new()
    } else {
        haplotype_sample_records(&selection.samples, &selection.source_indices)
    };
    selection.diagnostics.retained_variants = 0;
    finish_dense_matrix(
        DenseMatrixParts {
            n_samples,
            n_variants: 0,
            values: Vec::new(),
            missing_mask: Vec::new(),
            samples,
            variants: Vec::new(),
            diagnostics: selection.diagnostics,
        },
        matrix_only,
    )
}

fn empty_vcf_haplotypes_sparse(
    path: &Path,
    header: &rust_htslib::bcf::header::HeaderView,
    requested_samples: Option<&[String]>,
) -> Result<SparseGenotypeMatrix> {
    let all_samples = sample_records_from_header(header);
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    let samples = haplotype_sample_records(&selection.samples, &selection.source_indices);
    empty_sparse_matrix(samples, selection.diagnostics)
}

fn validate_dense_biallelic_record(path: &Path, record: &rust_htslib::bcf::Record) -> Result<()> {
    let allele_count = record.allele_count();
    if allele_count == 2 {
        return Ok(());
    }
    let record_id = record.id();
    let id = String::from_utf8_lossy(&record_id);
    let reason = if allele_count > 2 {
        "multi-ALT records are not supported"
    } else {
        "records with fewer than two alleles are not supported"
    };
    Err(GenoioError::invalid_source(
        path,
        format!(
            "vcf dense reads require biallelic records; record {id} has {allele_count} alleles: {reason}"
        ),
    ))
}

fn sample_records_from_header(header: &rust_htslib::bcf::header::HeaderView) -> Vec<SampleRecord> {
    header
        .samples()
        .iter()
        .map(|sample| plain_sample_record(String::from_utf8_lossy(sample).into_owned()))
        .collect()
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

fn apply_original_source_indices(
    samples: &mut [SampleRecord],
    original_source_indices: Option<&[usize]>,
) {
    let Some(original_source_indices) = original_source_indices else {
        return;
    };

    debug_assert_eq!(samples.len(), original_source_indices.len());
    for (sample, source_index) in samples.iter_mut().zip(original_source_indices) {
        sample.source_sample_index = Some(*source_index);
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

fn variant_record_from_record(
    path: &Path,
    header: &rust_htslib::bcf::header::HeaderView,
    record: &rust_htslib::bcf::Record,
) -> Result<VariantRecord> {
    let rid = record.rid().ok_or_else(|| {
        GenoioError::invalid_source(path, "vcf record is missing a chromosome id")
    })?;
    let chrom =
        String::from_utf8_lossy(header.rid2name(rid).map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf rid error: {error}"))
        })?)
        .into_owned();
    let pos = u32::try_from(record.pos() + 1)
        .map_err(|_| GenoioError::invalid_source(path, "vcf record position is out of range"))?;
    let id = String::from_utf8_lossy(&record.id()).into_owned();
    let alleles = record.alleles();
    if alleles.len() < 2 {
        return Err(GenoioError::invalid_source(
            path,
            format!("vcf record {chrom}:{pos} has fewer than two alleles"),
        ));
    }
    let ref_allele = String::from_utf8_lossy(alleles[0]).into_owned();
    let alt_values = alleles[1..]
        .iter()
        .map(|allele| String::from_utf8_lossy(allele).into_owned())
        .collect::<Vec<_>>();
    let first_alt = alt_values[0].clone();
    let alt_allele = alt_values.join(",");

    Ok(VariantRecord {
        chrom,
        pos,
        id,
        a0: ref_allele.clone(),
        a1: first_alt.clone(),
        ref_allele: Some(ref_allele.clone()),
        alt_allele: Some(alt_allele),
        source_a0: ref_allele,
        source_a1: first_alt,
        flipped: false,
        qual: finite_qual(record.qual()),
        af: None,
        maf: None,
        mac: None,
        missing_rate: None,
        n_called: None,
    })
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DiploidGtClass {
    HomRef,
    Het,
    HomAlt,
    Missing,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct RawDiploidGtCall {
    value: f32,
    class: DiploidGtClass,
}

impl RawDiploidGtCall {
    fn is_missing(self) -> bool {
        self.class == DiploidGtClass::Missing
    }
}

fn decode_raw_diploid_gt_call(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    genotype: &[i32],
) -> Result<RawDiploidGtCall> {
    if genotype.len() != 2 {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf record {} has non-diploid GT with {} alleles",
                String::from_utf8_lossy(&record.id()),
                genotype.len()
            ),
        ));
    }

    let mut alt_count = 0_u8;
    for encoded in genotype {
        match decode_raw_gt_allele(*encoded) {
            RawGtAllele::Missing => {
                return Ok(RawDiploidGtCall {
                    value: 0.0,
                    class: DiploidGtClass::Missing,
                });
            }
            RawGtAllele::Reference => {}
            RawGtAllele::Alternate => alt_count += 1,
            RawGtAllele::Unsupported(other) => {
                return Err(GenoioError::invalid_source(
                    path,
                    format!(
                        "vcf record {} has multiallelic GT allele index {other}",
                        String::from_utf8_lossy(&record.id())
                    ),
                ));
            }
        }
    }

    let class = match alt_count {
        0 => DiploidGtClass::HomRef,
        1 => DiploidGtClass::Het,
        2 => DiploidGtClass::HomAlt,
        _ => unreachable!("two diploid GT alleles can only produce dosage 0, 1, or 2"),
    };
    Ok(RawDiploidGtCall {
        value: f32::from(alt_count),
        class,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RawGtAllele {
    Missing,
    Reference,
    Alternate,
    Unsupported(i32),
}

fn decode_raw_gt_allele(encoded: i32) -> RawGtAllele {
    match encoded {
        0 | 1 => RawGtAllele::Missing,
        value if value > 1 => match (value >> 1) - 1 {
            0 => RawGtAllele::Reference,
            1 => RawGtAllele::Alternate,
            allele => RawGtAllele::Unsupported(allele),
        },
        value => RawGtAllele::Unsupported(value),
    }
}

struct DecodedDiploidGenotypeWithStats {
    values: Vec<f32>,
    missing: Vec<bool>,
    stats: VariantStats,
}

struct DecodedDiploidGenotypeValues {
    values: Vec<f32>,
    missing: Vec<bool>,
}

fn decode_diploid_genotype_values(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    source_indices: &[usize],
) -> Result<DecodedDiploidGenotypeValues> {
    let genotypes = record.format(b"GT").integer().map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf genotype error: {error}"))
    })?;
    let mut values = Vec::with_capacity(source_indices.len());
    let mut missing = Vec::with_capacity(source_indices.len());

    for source_index in source_indices {
        let call = decode_raw_diploid_gt_call(path, record, genotypes[*source_index])?;
        values.push(call.value);
        missing.push(call.is_missing());
    }

    Ok(DecodedDiploidGenotypeValues { values, missing })
}

fn decode_diploid_genotype_with_stats(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    source_indices: &[usize],
) -> Result<DecodedDiploidGenotypeWithStats> {
    let genotypes = record.format(b"GT").integer().map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf genotype error: {error}"))
    })?;
    let mut values = Vec::with_capacity(source_indices.len());
    let mut missing = Vec::with_capacity(source_indices.len());
    let mut hom_ref_count = 0_u64;
    let mut het_count = 0_u64;
    let mut hom_alt_count = 0_u64;
    let mut missing_count = 0_u64;

    // Genotype filters need per-variant counts before matrix materialization, so
    // this path fuses GT decoding with hard-call counting to avoid a second pass.
    for source_index in source_indices {
        let call = decode_raw_diploid_gt_call(path, record, genotypes[*source_index])?;
        match call.class {
            DiploidGtClass::HomRef => hom_ref_count += 1,
            DiploidGtClass::Het => het_count += 1,
            DiploidGtClass::HomAlt => hom_alt_count += 1,
            DiploidGtClass::Missing => missing_count += 1,
        }
        values.push(call.value);
        missing.push(call.is_missing());
    }

    let stats = variant_stats_from_counts(hom_ref_count, het_count, hom_alt_count, missing_count)?;
    Ok(DecodedDiploidGenotypeWithStats {
        values,
        missing,
        stats,
    })
}

fn decode_ds_dosage_record(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    source_indices: &[usize],
) -> Result<DecodedDiploidGenotypeValues> {
    let dosages = record.format(b"DS").float().map_err(|error| {
        GenoioError::unsupported(format!(
            "vcf dosage reads require FORMAT/DS values: {error}"
        ))
    })?;
    let mut values = Vec::with_capacity(source_indices.len());
    let mut missing = Vec::with_capacity(source_indices.len());

    for source_index in source_indices {
        let dosage = dosages[*source_index];
        if dosage.len() != 1 {
            return Err(GenoioError::invalid_source(
                path,
                format!(
                    "vcf record {} has FORMAT/DS with {} values for a sample; expected one",
                    String::from_utf8_lossy(&record.id()),
                    dosage.len()
                ),
            ));
        }
        let value = dosage[0];
        if value.is_missing() {
            values.push(0.0);
            missing.push(true);
            continue;
        }
        if !value.is_finite() || !(0.0..=2.0).contains(&value) {
            return Err(GenoioError::invalid_source(
                path,
                format!(
                    "vcf record {} has invalid FORMAT/DS value {value}; expected finite value in [0, 2]",
                    String::from_utf8_lossy(&record.id())
                ),
            ));
        }
        values.push(value);
        missing.push(false);
    }

    Ok(DecodedDiploidGenotypeValues { values, missing })
}

struct PhasedHaplotypeRecord {
    haplotype_values: Vec<f32>,
    haplotype_missing: Vec<bool>,
}

fn decode_phased_haplotype_record(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    source_indices: &[usize],
) -> Result<PhasedHaplotypeRecord> {
    let genotypes = record.genotypes().map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf genotype error: {error}"))
    })?;
    let mut haplotype_values = Vec::with_capacity(source_indices.len() * 2);
    let mut haplotype_missing = Vec::with_capacity(source_indices.len() * 2);

    for source_index in source_indices {
        let genotype = genotypes.get(*source_index);
        let (values, missing) = decode_phased_diploid_gt(path, record, &genotype)?;
        haplotype_values.extend(values);
        haplotype_missing.extend(missing);
    }

    Ok(PhasedHaplotypeRecord {
        haplotype_values,
        haplotype_missing,
    })
}

fn decode_phased_diploid_gt(
    path: &Path,
    record: &rust_htslib::bcf::Record,
    genotype: &rust_htslib::bcf::record::Genotype,
) -> Result<([f32; 2], [bool; 2])> {
    if genotype.len() != 2 {
        return Err(GenoioError::invalid_source(
            path,
            format!(
                "vcf record {} has non-diploid GT with {} alleles",
                String::from_utf8_lossy(&record.id()),
                genotype.len()
            ),
        ));
    }

    let mut values = [0.0, 0.0];
    let mut missing = [false, false];
    for (allele_index, allele) in genotype.iter().enumerate() {
        if allele_index > 0
            && matches!(
                allele,
                GenotypeAllele::UnphasedMissing | GenotypeAllele::Unphased(_)
            )
        {
            return Err(GenoioError::unsupported(format!(
                "vcf haplotype read record {} contains an unphased GT separator in a retained haplotype variant",
                String::from_utf8_lossy(&record.id())
            )));
        }
        match allele {
            GenotypeAllele::PhasedMissing | GenotypeAllele::UnphasedMissing => {
                missing[allele_index] = true;
            }
            GenotypeAllele::Phased(index) | GenotypeAllele::Unphased(index) => match index {
                0 => {}
                1 => values[allele_index] = 1.0,
                other => {
                    return Err(GenoioError::invalid_source(
                        path,
                        format!(
                            "vcf record {} has multiallelic GT allele index {other}",
                            String::from_utf8_lossy(&record.id())
                        ),
                    ));
                }
            },
        }
    }

    Ok((values, missing))
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

    use rust_htslib::bcf::{Read, Reader};

    use super::*;

    fn write_test_vcf(contents: &str) -> tempfile::NamedTempFile {
        let file = tempfile::Builder::new()
            .suffix(".vcf")
            .tempfile()
            .expect("temp VCF should be created");
        fs::write(file.path(), contents).expect("test VCF should be written");
        file
    }

    #[test]
    fn diploid_gt_decode_accumulates_stats_in_same_pass() {
        let file = write_test_vcf(
            "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.
",
        );
        let path = file.path();
        let mut reader = Reader::from_path(path).expect("test VCF should open");
        let record = reader
            .records()
            .next()
            .expect("record should exist")
            .expect("record should decode");

        let decoded = decode_diploid_genotype_with_stats(path, &record, &[0, 2, 3])
            .expect("GT should decode");

        assert_eq!(decoded.values, vec![0.0, 2.0, 0.0]);
        assert_eq!(decoded.missing, vec![false, false, true]);
        assert_eq!(decoded.stats.n_called, 2);
        assert_eq!(decoded.stats.mac, Some(2.0));
        assert_eq!(decoded.stats.maf, Some(0.5));
        assert!((decoded.stats.missing_rate - (1.0 / 3.0)).abs() < f64::EPSILON);
        assert!(decoded.stats.polymorphic);
    }
}
