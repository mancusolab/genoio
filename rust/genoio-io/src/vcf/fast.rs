//! Narrow VCF GT fast path built on noodles lazy records.
//!
//! This module is intentionally conservative: it accelerates dense matrix-only
//! hardcall reads for common biallelic diploid VCFs, and lets the htslib path
//! remain the correctness path for unsupported operations.

// pattern: Mixed (unavoidable)
// Reason: This performance path keeps reader setup close to decode routing so
// buffer ownership and htslib fallback boundaries stay explicit.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use flate2::read::MultiGzDecoder;
use genoio_core::{
    select_samples_source_order, DenseGenotypeMatrix, DenseSampleSelection, GenoioError,
    PartialFilterDecision, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};
use noodles_vcf as noodles;

use crate::error::Result;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use self::gt::{decode_gt_record, GtDecodeBuffers, GtStatsMode};
use self::header::read_sample_records_from_header;
use self::output::{can_write_sample_major_directly, FastDenseOutput};
use self::record::variant_record_from_record;
use self::sparse::{read_haplotype_sparse_records, read_sparse_records};

use super::is_compressed_vcf;

mod gt;
mod header;
mod output;
mod record;
mod sparse;

const VCF_FAST_BUFFER_SIZE: usize = 1 << 20;

type FastVcfReader = noodles::io::Reader<BufReader<MultiGzDecoder<File>>>;

struct FastVcfInput {
    reader: FastVcfReader,
    source_sample_count: usize,
    selection: DenseSampleSelection,
}

pub(super) fn try_read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<Option<DenseGenotypeMatrix>> {
    if !is_fast_dense_supported(path, variant_filter, matrix_only, threads) {
        return Ok(None);
    }
    let Some(variant_window) = variant_window else {
        return Ok(None);
    };

    let FastVcfInput {
        mut reader,
        source_sample_count,
        selection,
    } = open_fast_vcf_input(path, requested_samples)?;

    read_dense_records(
        path,
        variant_filter,
        variant_window,
        source_sample_count,
        &selection,
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

fn is_fast_dense_supported(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    matrix_only: bool,
    threads: Option<usize>,
) -> bool {
    matrix_only && is_fast_vcf_supported(path, variant_filter, threads)
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

fn open_lazy_reader(path: &Path) -> Result<FastVcfReader> {
    let file = File::open(path)
        .map_err(|error| GenoioError::invalid_source(path, format!("vcf open error: {error}")))?;
    let reader = BufReader::with_capacity(VCF_FAST_BUFFER_SIZE, MultiGzDecoder::new(file));
    Ok(noodles::io::Reader::new(reader))
}

fn read_dense_records<R: std::io::BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: VariantWindow,
    source_sample_count: usize,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics.clone();
    let n_samples = selection.samples.len();
    let output_variant_capacity = variant_window.len;
    let direct_sample_major =
        can_write_sample_major_directly(selection, source_sample_count, variant_filter);
    let mut output = FastDenseOutput::new(n_samples, output_variant_capacity, direct_sample_major);
    let mut output_variant_count = 0;
    let mut retention = RetainedVariantState::new(Some(variant_window));
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

        let variant = variant_record_from_record(path, &record)?;
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

        output.write_variant(output_variant_count, &decoded)?;
        output_variant_count += 1;
    }

    diagnostics.retained_variants = output_variant_count;
    output.finish(output_variant_count, selection.samples.clone(), diagnostics)
}
