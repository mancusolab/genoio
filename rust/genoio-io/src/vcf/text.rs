//! Text VCF backend built on noodles lazy records.
//!
//! This module is intentionally conservative: it accelerates common text VCF
//! scans while BCF uses a separate typed backend. Threaded compressed VCF reads
//! use noodles' BGZF block decompression; record parsing remains ordered and
//! single-consumer.

// pattern: Mixed (unavoidable)
// Reason: Hot record loops interleave lazy record iteration, filtering, and
// output staging so buffers can be reused without extra abstraction overhead.

use std::io::BufRead;
use std::path::Path;

use genoio_core::{
    DenseGenotypeMatrix, DenseSampleSelection, GenoioError, MetadataOutput, PartialFilterDecision,
    RegionPredicate, SourceCapabilities, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};
use noodles_vcf as noodles;

use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::hardcall::evaluate_hardcall_counts_filter;
use crate::matrix::empty_sparse_matrix;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use self::ds::{decode_ds_record, DsDecodeBuffers};
use self::gt::{
    decode_gt_record, decode_phased_gt_dense_record, text_record_has_phased_genotype,
    GtDecodeBuffers, GtStatsMode, HaplotypeDenseDecodeBuffers,
};
use self::header::read_sample_records_from_header;
use self::output::{can_write_sample_major_directly, TextDenseOutput};
use self::record::{
    skip_variant_for_region, validate_biallelic_variant, variant_record_from_text_record,
};
use self::source::{
    ensure_text_indexed_vcf_supported, ensure_text_vcf_supported, open_compressed_reader,
    open_plain_reader, open_text_sample_selection, open_text_vcf_input,
    with_indexed_text_vcf_input, with_threaded_indexed_text_vcf_input, DenseReadSource,
    TextVcfInput, TextVcfSource,
};
use self::sparse::{read_haplotype_sparse_records, read_sparse_records};

use super::{haplotype_sample_records, is_compressed_vcf};

mod ds;
mod format;
mod gt;
mod header;
mod output;
mod record;
mod source;
mod sparse;

pub(in crate::vcf) use self::record::variant_record_from_noodles_variant_record;

pub(super) const VCF_TEXT_BUFFER_SIZE: usize = 1 << 20;

impl TextVcfSource {
    fn read_dense(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => {
                read_dense_from_input(path, variant_filter, variant_window, matrix_only, input)
            }
            Self::ThreadedCompressed(input) => {
                read_dense_from_input(path, variant_filter, variant_window, matrix_only, input)
            }
            Self::Plain(input) => {
                read_dense_from_input(path, variant_filter, variant_window, matrix_only, input)
            }
        }
    }

    fn read_dosage_dense(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => read_dosage_dense_from_input(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                input,
            ),
            Self::ThreadedCompressed(input) => read_dosage_dense_from_input(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                input,
            ),
            Self::Plain(input) => read_dosage_dense_from_input(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                input,
            ),
        }
    }

    fn read_haplotype_dense(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => read_haplotype_dense_from_input(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                input,
            ),
            Self::ThreadedCompressed(input) => read_haplotype_dense_from_input(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                input,
            ),
            Self::Plain(input) => read_haplotype_dense_from_input(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                input,
            ),
        }
    }

    fn read_sparse(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
    ) -> Result<SparseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => {
                read_sparse_from_input(path, variant_filter, variant_window, input)
            }
            Self::ThreadedCompressed(input) => {
                read_sparse_from_input(path, variant_filter, variant_window, input)
            }
            Self::Plain(input) => {
                read_sparse_from_input(path, variant_filter, variant_window, input)
            }
        }
    }

    fn read_haplotype_sparse(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
    ) -> Result<SparseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => {
                read_haplotype_sparse_from_input(path, variant_filter, variant_window, input)
            }
            Self::ThreadedCompressed(input) => {
                read_haplotype_sparse_from_input(path, variant_filter, variant_window, input)
            }
            Self::Plain(input) => {
                read_haplotype_sparse_from_input(path, variant_filter, variant_window, input)
            }
        }
    }

    fn into_selection(self) -> DenseSampleSelection {
        match self {
            Self::Compressed(input) => input.selection,
            Self::ThreadedCompressed(input) => input.selection,
            Self::Plain(input) => input.selection,
        }
    }
}

// Keep the threaded/unthreaded indexed reader choice at setup time. A small
// macro avoids per-record dynamic dispatch without repeating the same match in
// every output-mode entry point.
macro_rules! with_indexed_text_vcf_input_for_threads {
    (
        $path:expr,
        $requested_samples:expr,
        $region:expr,
        $threads:expr,
        |$input:ident, $reader:ident| $read:block,
        $empty:expr $(,)?
    ) => {{
        match $threads {
            Some(threads) => with_threaded_indexed_text_vcf_input(
                $path,
                $requested_samples,
                $region,
                threads,
                |$input, $reader| $read,
                $empty,
            ),
            None => with_indexed_text_vcf_input(
                $path,
                $requested_samples,
                $region,
                |$input, $reader| $read,
                $empty,
            ),
        }
    }};
}

pub(super) fn read_vcf_metadata(path: &Path) -> Result<MetadataOutput> {
    // Metadata reads can avoid strict full-header parsing. They only need the
    // #CHROM line plus record fields, which keeps real-world VCF headers from
    // forcing a strict parser onto the hot path.
    if is_compressed_vcf(path) {
        let mut reader = open_compressed_reader(path)?;
        read_metadata_records(path, &mut reader)
    } else {
        let mut reader = open_plain_reader(path)?;
        read_metadata_records(path, &mut reader)
    }
}

pub(super) fn read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_dense(
        path,
        variant_filter,
        variant_window,
        matrix_only,
    )
}

pub(super) fn empty_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_dense_from_selection(selection, matrix_only)
}

pub(super) fn empty_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_sparse_from_selection(selection)
}

pub(super) fn empty_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_haplotype_dense_from_selection(selection, matrix_only)
}

pub(super) fn empty_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_haplotype_sparse_from_selection(selection)
}

fn read_dense_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    input: TextVcfInput<R>,
) -> Result<DenseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        source_sample_count,
        selection,
    } = input;
    read_dense_records(
        path,
        variant_filter,
        variant_window,
        matrix_only,
        DenseReadSource::full_scan(source_sample_count),
        &selection,
        &mut reader,
    )
}

pub(super) fn read_vcf_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    with_indexed_text_vcf_input_for_threads!(
        path,
        requested_samples,
        region,
        threads,
        |input, reader| {
            read_dense_records(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                input.dense_source(),
                &input.selection,
                reader,
            )
        },
        |selection| empty_dense_from_selection(selection, matrix_only),
    )
}

pub(super) fn read_vcf_dosage_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_dosage_dense(
        path,
        variant_filter,
        variant_window,
        matrix_only,
    )
}

pub(super) fn read_vcf_dosage_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    with_indexed_text_vcf_input_for_threads!(
        path,
        requested_samples,
        region,
        threads,
        |input, reader| {
            read_dosage_dense_records(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                input.dense_source(),
                &input.selection,
                reader,
            )
        },
        |selection| empty_dense_from_selection(selection, matrix_only),
    )
}

fn read_dosage_dense_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    input: TextVcfInput<R>,
) -> Result<DenseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        source_sample_count,
        selection,
    } = input;
    read_dosage_dense_records(
        path,
        variant_filter,
        variant_window,
        matrix_only,
        DenseReadSource::full_scan(source_sample_count),
        &selection,
        &mut reader,
    )
}

pub(super) fn read_vcf_haplotypes_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    with_indexed_text_vcf_input_for_threads!(
        path,
        requested_samples,
        region,
        threads,
        |input, reader| {
            read_haplotype_dense_records(
                path,
                variant_filter,
                variant_window,
                matrix_only,
                Some(input.region),
                input.selection,
                reader,
            )
        },
        |selection| empty_haplotype_dense_from_selection(selection, matrix_only),
    )
}

pub(super) fn read_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_haplotype_dense(
        path,
        variant_filter,
        variant_window,
        matrix_only,
    )
}

fn read_haplotype_dense_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
    input: TextVcfInput<R>,
) -> Result<DenseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        selection,
        ..
    } = input;
    read_haplotype_dense_records(
        path,
        variant_filter,
        variant_window,
        matrix_only,
        None,
        selection,
        &mut reader,
    )
}

pub(super) fn read_vcf_sparse_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    with_indexed_text_vcf_input_for_threads!(
        path,
        requested_samples,
        region,
        threads,
        |input, reader| {
            read_sparse_records(
                path,
                variant_filter,
                variant_window,
                Some(input.region),
                input.selection,
                reader,
            )
        },
        empty_sparse_from_selection,
    )
}

pub(super) fn read_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_sparse(
        path,
        variant_filter,
        variant_window,
    )
}

fn read_sparse_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    input: TextVcfInput<R>,
) -> Result<SparseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        selection,
        ..
    } = input;
    read_sparse_records(
        path,
        variant_filter,
        variant_window,
        None,
        selection,
        &mut reader,
    )
}

pub(super) fn read_vcf_haplotypes_sparse_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    with_indexed_text_vcf_input_for_threads!(
        path,
        requested_samples,
        region,
        threads,
        |input, reader| {
            read_haplotype_sparse_records(
                path,
                variant_filter,
                variant_window,
                Some(input.region),
                input.selection,
                reader,
            )
        },
        empty_haplotype_sparse_from_selection,
    )
}

pub(super) fn read_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_haplotype_sparse(
        path,
        variant_filter,
        variant_window,
    )
}

fn read_haplotype_sparse_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    input: TextVcfInput<R>,
) -> Result<SparseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        selection,
        ..
    } = input;
    read_haplotype_sparse_records(
        path,
        variant_filter,
        variant_window,
        None,
        selection,
        &mut reader,
    )
}

fn empty_dense_from_selection(
    selection: DenseSampleSelection,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    TextDenseOutput::new(selection.samples.len(), 0, false).finish(
        0,
        selection.samples,
        Vec::new(),
        diagnostics,
        matrix_only,
    )
}

fn empty_sparse_from_selection(selection: DenseSampleSelection) -> Result<SparseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    empty_sparse_matrix(selection.samples, diagnostics)
}

fn empty_haplotype_dense_from_selection(
    selection: DenseSampleSelection,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len() * 2;
    let samples = if matrix_only {
        Vec::new()
    } else {
        haplotype_sample_records(&selection.samples, &selection.source_indices)
    };
    TextDenseOutput::new(n_samples, 0, false).finish(
        0,
        samples,
        Vec::new(),
        diagnostics,
        matrix_only,
    )
}

fn empty_haplotype_sparse_from_selection(
    selection: DenseSampleSelection,
) -> Result<SparseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let samples = haplotype_sample_records(&selection.samples, &selection.source_indices);
    empty_sparse_matrix(samples, diagnostics)
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
            GenoioError::invalid_source(path, format!("text VCF record error: {error}"))
        })? == 0
        {
            break;
        }

        if !has_phased_genotype_evidence && text_record_has_phased_genotype(&record) {
            has_phased_genotype_evidence = true;
        }
        variants.push(variant_record_from_text_record(path, &record)?);
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
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity, direct_sample_major);
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
            GenoioError::invalid_source(path, format!("text VCF record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = variant_record_from_text_record(path, &record)?;
        // Tabix/CSI chunks can include neighboring records from the same BGZF
        // block. Keep the text backend's exact region contract independent of the
        // lower-level chunk boundaries.
        if skip_variant_for_region(&variant, source.region) {
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
        validate_biallelic_variant(path, &variant)?;

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        // `record.samples()` borrows noodles' reusable record buffer, so decode
        // selected GTs completely before the next `read_record` call.
        let stats_mode = match (needs_genotype_decision, matrix_only) {
            (true, true) => GtStatsMode::Counts,
            (true, false) => GtStatsMode::Compute,
            (false, _) => GtStatsMode::Skip,
        };
        decode_gt_record(
            path,
            &record,
            &selection.source_indices,
            stats_mode,
            &mut decoded,
        )?;

        if needs_genotype_decision {
            let filter = variant_filter.ok_or_else(|| {
                GenoioError::internal_contract("genotype decision requires a variant filter")
            })?;
            let (retain_variant, stats) = if matrix_only {
                let counts = decoded.counts().ok_or_else(|| {
                    GenoioError::internal_contract("matrix-only vcf GT filter missing counts")
                })?;
                evaluate_hardcall_counts_filter(
                    counts,
                    filter,
                    filter.genotype_filter_plan(),
                    Some(&variant),
                    false,
                )?
            } else {
                let stats = decoded
                    .stats()
                    .ok_or_else(|| GenoioError::internal_contract("vcf GT filter missing stats"))?;
                (filter.evaluate(&variant, Some(&stats)), Some(stats))
            };
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if let Some(stats) = stats {
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

fn read_dosage_dense_records<R: BufRead>(
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
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity, direct_sample_major);
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
            GenoioError::invalid_source(path, format!("text VCF record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = variant_record_from_text_record(path, &record)?;
        if skip_variant_for_region(&variant, source.region) {
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
        validate_biallelic_variant(path, &variant)?;

        decode_ds_record(path, &record, &selection.source_indices, &mut decoded)?;
        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        if needs_genotype_decision {
            let filter = variant_filter.ok_or_else(|| {
                GenoioError::internal_contract("genotype decision requires a variant filter")
            })?;
            let (retain_variant, stats) = evaluate_dosage_filter(
                decoded.values(),
                decoded.missing(),
                filter,
                &variant,
                !matrix_only,
            )?;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
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
    source_region: Option<&RegionPredicate>,
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
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity, false);
    let mut variants = Vec::with_capacity(if matrix_only {
        0
    } else {
        output_variant_capacity
    });
    let mut output_variant_count = 0;
    let mut retention = RetainedVariantState::new(variant_window);
    let mut record = noodles::Record::default();
    let mut decoded = HaplotypeDenseDecodeBuffers::with_capacity(source_indices.len());
    let mut stats_decoded = GtDecodeBuffers::with_capacity(source_indices.len());

    loop {
        if retention.window_is_satisfied() {
            break;
        }
        if reader.read_record(&mut record).map_err(|error| {
            GenoioError::invalid_source(path, format!("text VCF record error: {error}"))
        })? == 0
        {
            break;
        }

        let mut variant = variant_record_from_text_record(path, &record)?;
        if skip_variant_for_region(&variant, source_region) {
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
        validate_biallelic_variant(path, &variant)?;

        let needs_genotype_decision =
            matches!(partial_decision, PartialFilterDecision::NeedGenotypes);
        if needs_genotype_decision {
            // Genotype-stat filters are evaluated on diploid dosage before
            // enforcing phased output. This lets filters drop unphased records
            // without surfacing a haplotype decode error.
            decode_gt_record(
                path,
                &record,
                &source_indices,
                GtStatsMode::Compute,
                &mut stats_decoded,
            )?;
            let stats = stats_decoded.stats();
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
        decode_phased_gt_dense_record(
            path,
            &record,
            &source_indices,
            GtStatsMode::Skip,
            &mut decoded,
        )?;
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
    fn dense_text_backend_accepts_metadata_reads() {
        assert!(ensure_text_vcf_supported(Path::new("example.vcf.gz"), None).is_ok());
    }
}
