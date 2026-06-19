//! Narrow VCF GT fast path built on noodles lazy records.
//!
//! This module is intentionally conservative: it accelerates full-scan GT reads
//! for common biallelic diploid VCFs, and lets the htslib path remain the
//! correctness path for indexed, threaded, and unsupported operations.

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

use self::gt::{
    decode_gt_record, decode_phased_gt_dense_record, GtDecodeBuffers, GtStatsMode,
    HaplotypeDenseDecodeBuffers,
};
use self::header::read_sample_records_from_header;
use self::output::{can_write_sample_major_directly, FastDenseOutput};
use self::record::variant_record_from_record;
use self::sparse::{read_haplotype_sparse_records, read_sparse_records};

use super::{haplotype_sample_records, is_compressed_vcf};

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
