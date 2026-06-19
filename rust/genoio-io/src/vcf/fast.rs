//! Narrow VCF fast path built on noodles lazy records.
//!
//! This module is intentionally conservative: it accelerates full-scan metadata
//! and GT reads for common text VCFs, and lets the htslib path remain the
//! correctness path for BCF, indexed, threaded, and unsupported operations.

// pattern: Mixed (unavoidable)
// Reason: This performance path keeps reader setup close to decode routing so
// buffer ownership and htslib fallback boundaries stay explicit.

use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use genoio_core::{
    compute_dosage_variant_stats, select_samples_source_order, DenseGenotypeMatrix,
    DenseSampleSelection, GenoioError, MetadataOutput, PartialFilterDecision, RegionPredicate,
    SourceCapabilities, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};
use noodles_bgzf as bgzf;
use noodles_core::{Position, Region};
use noodles_csi::{self as csi, BinningIndex};
use noodles_tabix as tabix;
use noodles_vcf as noodles;

use crate::error::Result;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use self::ds::{decode_ds_record, DsDecodeBuffers};
use self::gt::{
    decode_gt_record, decode_phased_gt_dense_record, record_has_phased_gt_evidence,
    GtDecodeBuffers, GtStatsMode, HaplotypeDenseDecodeBuffers,
};
use self::header::read_sample_records_from_header;
use self::output::{can_write_sample_major_directly, FastDenseOutput};
use self::record::{metadata_variant_record_from_record, variant_record_from_record};
use self::sparse::{read_haplotype_sparse_records, read_sparse_records};

use super::{haplotype_sample_records, is_compressed_vcf};

mod ds;
mod format;
mod gt;
mod header;
mod output;
mod record;
mod sparse;

const VCF_FAST_BUFFER_SIZE: usize = 1 << 20;

type FastVcfReader = noodles::io::Reader<BufReader<MultiGzDecoder<File>>>;
type PlainVcfReader = noodles::io::Reader<BufReader<File>>;
type IndexChunk = csi::binning_index::index::reference_sequence::bin::Chunk;

struct FastVcfInput {
    reader: FastVcfReader,
    source_sample_count: usize,
    selection: DenseSampleSelection,
}

#[derive(Clone, Copy)]
struct DenseReadSource<'a> {
    sample_count: usize,
    region: Option<&'a RegionPredicate>,
}

pub(super) fn read_vcf_metadata(path: &Path) -> Result<MetadataOutput> {
    // Metadata reads can avoid strict full-header parsing. They only need the
    // #CHROM line plus record fields, which keeps real-world VCF headers from
    // forcing the slower htslib path.
    if is_compressed_vcf(path) {
        let mut reader = open_lazy_reader(path)?;
        read_metadata_records(path, &mut reader)
    } else {
        let mut reader = open_plain_reader(path)?;
        read_metadata_records(path, &mut reader)
    }
}

pub(super) fn try_read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<Option<DenseGenotypeMatrix>> {
    if !is_fast_vcf_supported(path, variant_filter, threads) {
        return Ok(None);
    }

    let FastVcfInput {
        mut reader,
        source_sample_count,
        selection,
    } = open_fast_vcf_input(path, requested_samples)?;

    read_dense_records(
        path,
        variant_filter,
        variant_window,
        matrix_only,
        DenseReadSource {
            sample_count: source_sample_count,
            region: None,
        },
        &selection,
        &mut reader,
    )
    .map(Some)
}

pub(super) fn try_read_vcf_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<Option<DenseGenotypeMatrix>> {
    if !is_fast_indexed_vcf_supported(path, threads) {
        return Ok(None);
    }

    let chunks = index_chunks_for_region(path, region)?;
    let mut bgzf_reader = open_bgzf_reader(path)?;
    let Ok(all_samples) = read_sample_records_from_header(path, &mut bgzf_reader) else {
        return Ok(None);
    };
    let source_sample_count = all_samples.len();
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;

    let Some(chunks) = chunks else {
        return empty_dense_from_selection(selection, matrix_only).map(Some);
    };

    // Query chunks directly rather than using noodles' strict IndexedReader.
    // The fast path only needs the #CHROM line, and real VCFs can have header
    // records that htslib accepts but noodles' full header parser rejects.
    let query = csi::io::Query::new(&mut bgzf_reader, chunks);
    let mut reader = noodles::io::Reader::new(query);
    read_dense_records(
        path,
        variant_filter,
        variant_window,
        matrix_only,
        DenseReadSource {
            sample_count: source_sample_count,
            region: Some(region),
        },
        &selection,
        &mut reader,
    )
    .map(Some)
}

pub(super) fn try_read_vcf_dosage_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<Option<DenseGenotypeMatrix>> {
    if !is_fast_vcf_supported(path, variant_filter, threads) {
        return Ok(None);
    }

    let FastVcfInput {
        mut reader,
        source_sample_count,
        selection,
    } = open_fast_vcf_input(path, requested_samples)?;

    read_dosage_dense_records(
        path,
        variant_filter,
        variant_window,
        matrix_only,
        source_sample_count,
        &selection,
        &mut reader,
    )
    .map(Some)
}

pub(super) fn try_read_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<Option<DenseGenotypeMatrix>> {
    if !is_fast_vcf_supported(path, variant_filter, threads) {
        return Ok(None);
    }

    let FastVcfInput {
        mut reader,
        selection,
        ..
    } = open_fast_vcf_input(path, requested_samples)?;

    read_haplotype_dense_records(
        path,
        variant_filter,
        variant_window,
        matrix_only,
        selection,
        &mut reader,
    )
    .map(Some)
}

pub(super) fn try_read_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    threads: Option<usize>,
) -> Result<Option<SparseGenotypeMatrix>> {
    if !is_fast_vcf_supported(path, variant_filter, threads) {
        return Ok(None);
    }

    let FastVcfInput {
        mut reader,
        selection,
        ..
    } = open_fast_vcf_input(path, requested_samples)?;

    read_sparse_records(path, variant_filter, variant_window, selection, &mut reader).map(Some)
}

pub(super) fn try_read_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    threads: Option<usize>,
) -> Result<Option<SparseGenotypeMatrix>> {
    if !is_fast_vcf_supported(path, variant_filter, threads) {
        return Ok(None);
    }

    let FastVcfInput {
        mut reader,
        selection,
        ..
    } = open_fast_vcf_input(path, requested_samples)?;

    read_haplotype_sparse_records(path, variant_filter, variant_window, selection, &mut reader)
        .map(Some)
}

fn is_fast_vcf_supported(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    threads: Option<usize>,
) -> bool {
    is_compressed_vcf(path)
        && threads.is_none()
        && !variant_filter.is_some_and(VariantFilter::has_region_predicate)
}

fn is_fast_indexed_vcf_supported(path: &Path, threads: Option<usize>) -> bool {
    is_compressed_vcf(path) && threads.is_none()
}

fn open_fast_vcf_input(path: &Path, requested_samples: Option<&[String]>) -> Result<FastVcfInput> {
    let mut reader = open_lazy_reader(path)?;
    let all_samples = read_sample_records_from_header(path, reader.get_mut())?;
    let source_sample_count = all_samples.len();
    let selection = select_samples_source_order(&all_samples, requested_samples, path)?;
    Ok(FastVcfInput {
        reader,
        source_sample_count,
        selection,
    })
}

fn noodles_region_from_predicate(path: &Path, region: &RegionPredicate) -> Result<Region> {
    let start = position_from_u32(path, region.start, "start")?;
    let end = position_from_u32(path, region.end, "end")?;
    Ok(Region::new(region.chrom.as_str(), start..=end))
}

fn position_from_u32(path: &Path, value: u32, label: &str) -> Result<Position> {
    let value = usize::try_from(value).map_err(|_| {
        GenoioError::invalid_source(
            path,
            format!("vcf region {label} coordinate is out of range"),
        )
    })?;
    Position::try_from(value).map_err(|error| {
        GenoioError::invalid_source(
            path,
            format!("vcf region {label} coordinate is invalid: {error}"),
        )
    })
}

fn empty_dense_from_selection(
    selection: DenseSampleSelection,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    FastDenseOutput::new(selection.samples.len(), 0, false).finish(
        0,
        selection.samples,
        Vec::new(),
        diagnostics,
        matrix_only,
    )
}

fn variant_in_region(variant: &genoio_core::VariantRecord, region: &RegionPredicate) -> bool {
    variant.chrom == region.chrom && variant.pos >= region.start && variant.pos <= region.end
}

fn open_bgzf_reader(path: &Path) -> Result<bgzf::io::Reader<File>> {
    File::open(path)
        .map(bgzf::io::Reader::new)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf open error: {error}")))
}

fn index_chunks_for_region(
    path: &Path,
    region: &RegionPredicate,
) -> Result<Option<Vec<IndexChunk>>> {
    let index = read_associated_index(path).map_err(|error| {
        GenoioError::invalid_source(path, format!("vcf index read error: {error}"))
    })?;
    let region = noodles_region_from_predicate(path, region)?;
    let Some(header) = index.header() else {
        return Ok(None);
    };
    let Some(reference_sequence_id) = header
        .reference_sequence_names()
        .get_index_of(region.name())
    else {
        return Ok(None);
    };
    index
        .query(reference_sequence_id, region.interval())
        .map(Some)
        .map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf index query error: {error}"))
        })
}

fn read_associated_index(path: &Path) -> std::io::Result<Box<dyn BinningIndex>> {
    match tabix::fs::read(companion_index_path(path, "tbi")) {
        Ok(index) => Ok(Box::new(index)),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            csi::fs::read(companion_index_path(path, "csi"))
                .map(|index| Box::new(index) as Box<dyn BinningIndex>)
        }
        Err(error) => Err(error),
    }
}

fn companion_index_path(path: &Path, extension: &str) -> PathBuf {
    let mut raw = OsString::from(path);
    raw.push(".");
    raw.push(extension);
    PathBuf::from(raw)
}

fn open_lazy_reader(path: &Path) -> Result<FastVcfReader> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf open error: {error}")))?;
    let reader = BufReader::with_capacity(VCF_FAST_BUFFER_SIZE, MultiGzDecoder::new(file));
    Ok(noodles::io::Reader::new(reader))
}

fn open_plain_reader(path: &Path) -> Result<PlainVcfReader> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf open error: {error}")))?;
    let reader = BufReader::with_capacity(VCF_FAST_BUFFER_SIZE, file);
    Ok(noodles::io::Reader::new(reader))
}

fn read_metadata_records<R: BufRead>(
    path: &Path,
    reader: &mut noodles::io::Reader<R>,
) -> Result<MetadataOutput> {
    let samples = read_sample_records_from_header(path, reader.get_mut())?;
    let mut variants = Vec::new();
    let mut has_phased_genotype_evidence = false;
    let mut record = noodles::Record::default();

    loop {
        if reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf fast record error: {error}"))
        })? == 0
        {
            break;
        }

        if !has_phased_genotype_evidence && record_has_phased_gt_evidence(&record) {
            has_phased_genotype_evidence = true;
        }
        variants.push(metadata_variant_record_from_record(path, &record)?);
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

fn read_dense_records<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics.clone();
    let n_samples = selection.samples.len();
    let output_variant_capacity = variant_window.map_or(0, |window| window.len);
    let direct_sample_major = matrix_only
        && variant_window.is_some()
        && can_write_sample_major_directly(selection, source.sample_count, variant_filter);
    let mut output = FastDenseOutput::new(n_samples, output_variant_capacity, direct_sample_major);
    let mut output_variant_count = 0;
    let mut variants = Vec::with_capacity(if matrix_only {
        0
    } else {
        output_variant_capacity
    });
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = noodles::Record::default();
    let mut decoded = GtDecodeBuffers::with_capacity(selection.source_indices.len());

    loop {
        if retention.window_is_satisfied() {
            break;
        }
        if reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf fast record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = variant_record_from_record(path, &record)?;
        // Tabix/CSI chunks can include neighboring records from the same BGZF
        // block. Keep the fast path's exact region contract independent of the
        // lower-level chunk boundaries.
        if source
            .region
            .is_some_and(|region| !variant_in_region(&variant, region))
        {
            continue;
        }
        let partial_decision = variant_filter
            .map(|filter| filter.partial_decision(&variant))
            .unwrap_or(PartialFilterDecision::Accept);
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        // `record.samples()` borrows noodles' reusable record buffer, so decode
        // selected GTs completely before the next `read_record` call.
        decode_gt_record(
            path,
            &record,
            &selection.source_indices,
            GtStatsMode::from_needed(needs_genotype_decision),
            &mut decoded,
        )?;
        let stats = if needs_genotype_decision {
            decoded.stats()
        } else {
            None
        };

        if let Some(stats) = stats {
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
        }
        if !matrix_only {
            if let Some(stats) = stats {
                genoio_core::attach_variant_stats(&mut variant, stats);
            }
            variants.push(variant);
        }

        output.write_variant(output_variant_count, decoded.values(), decoded.missing())?;
        output_variant_count += 1;
    }

    diagnostics.retained_variants = output_variant_count;
    output.finish(
        output_variant_count,
        selection.samples.clone(),
        variants,
        diagnostics,
        matrix_only,
    )
}

fn read_dosage_dense_records<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    source_sample_count: usize,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics.clone();
    let n_samples = selection.samples.len();
    let output_variant_capacity = variant_window.map_or(0, |window| window.len);
    let direct_sample_major = matrix_only
        && variant_window.is_some()
        && can_write_sample_major_directly(selection, source_sample_count, variant_filter);
    let mut output = FastDenseOutput::new(n_samples, output_variant_capacity, direct_sample_major);
    let mut output_variant_count = 0;
    let mut variants = Vec::with_capacity(if matrix_only {
        0
    } else {
        output_variant_capacity
    });
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = noodles::Record::default();
    let mut decoded = DsDecodeBuffers::with_capacity(selection.source_indices.len());

    loop {
        if retention.window_is_satisfied() {
            break;
        }
        if reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf fast record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = variant_record_from_record(path, &record)?;
        let partial_decision = variant_filter
            .map(|filter| filter.partial_decision(&variant))
            .unwrap_or(PartialFilterDecision::Accept);
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        decode_ds_record(path, &record, &selection.source_indices, &mut decoded)?;
        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        let stats = if needs_genotype_decision {
            Some(compute_dosage_variant_stats(
                decoded.values(),
                decoded.missing(),
            )?)
        } else {
            None
        };

        if let Some(stats) = stats {
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, Some(&stats))),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if !matrix_only {
                genoio_core::attach_variant_stats(&mut variant, stats);
            }
        }
        if !matrix_only {
            variants.push(variant);
        }

        output.write_variant(output_variant_count, decoded.values(), decoded.missing())?;
        output_variant_count += 1;
    }

    diagnostics.retained_variants = output_variant_count;
    output.finish(
        output_variant_count,
        selection.samples.clone(),
        variants,
        diagnostics,
        matrix_only,
    )
}

fn read_haplotype_dense_records<R: std::io::BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    selection: DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    let DenseSampleSelection {
        source_indices,
        samples,
        mut diagnostics,
    } = selection;
    let n_samples = samples.len() * 2;
    let output_variant_capacity = variant_window.map_or(0, |window| window.len);
    let mut output = FastDenseOutput::new(n_samples, output_variant_capacity, false);
    let mut variants = Vec::with_capacity(if matrix_only {
        0
    } else {
        output_variant_capacity
    });
    let mut output_variant_count = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = noodles::Record::default();
    let mut decoded = HaplotypeDenseDecodeBuffers::with_capacity(source_indices.len());

    loop {
        if retention.window_is_satisfied() {
            break;
        }
        if reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("vcf fast record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = variant_record_from_record(path, &record)?;
        let partial_decision = variant_filter
            .map(|filter| filter.partial_decision(&variant))
            .unwrap_or(PartialFilterDecision::Accept);
        match retention.metadata_decision(partial_decision, &mut diagnostics) {
            MetadataRetentionAction::Include | MetadataRetentionAction::DecodeGenotypes => {}
            MetadataRetentionAction::Skip => continue,
            MetadataRetentionAction::Stop => break,
        }

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        decode_phased_gt_dense_record(
            path,
            &record,
            &source_indices,
            GtStatsMode::from_needed(needs_genotype_decision),
            &mut decoded,
        )?;

        if needs_genotype_decision {
            let stats = decoded.stats();
            match retention.genotype_decision(
                variant_filter.is_none_or(|filter| filter.evaluate(&variant, stats.as_ref())),
                &mut diagnostics,
            ) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if !matrix_only {
                if let Some(stats) = stats {
                    genoio_core::attach_variant_stats(&mut variant, stats);
                }
            }
        }
        if !matrix_only {
            variants.push(variant);
        }

        // Dense haplotype reads expose source allele-1 indicators. Sparse
        // haplotype output flips columns to minor allele later to reduce nnz.
        output.write_variant(output_variant_count, decoded.values(), decoded.missing())?;
        output_variant_count += 1;
    }

    let samples = if matrix_only {
        Vec::new()
    } else {
        haplotype_sample_records(&samples, &source_indices)
    };
    diagnostics.retained_variants = output_variant_count;
    output.finish(
        output_variant_count,
        samples,
        variants,
        diagnostics,
        matrix_only,
    )
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    #[test]
    fn dense_fast_path_accepts_metadata_reads() {
        assert!(is_fast_vcf_supported(
            Path::new("example.vcf.gz"),
            None,
            None,
        ));
    }
}
