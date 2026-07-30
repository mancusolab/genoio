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
    DenseDiagnostics, DenseGenotypeMatrix, DenseMissingPolicy, DenseSampleSelection, GenoioError,
    MetadataOutput, RegionPredicate, SampleMetadataBuffers, SourceCapabilities,
    SparseGenotypeMatrix, VariantFilter, VariantMetadataBuffers, VariantMetadataView, VariantStats,
    VariantWindow,
};
use noodles_vcf as noodles;

use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::hardcall::evaluate_hardcall_counts_filter;
use crate::retention::{RetainedVariantState, RetentionAction};

use self::ds::{decode_ds_record, DsDecodeBuffers};
use self::gt::{
    decode_gt_record, decode_phased_gt_dense_record, text_record_has_phased_genotype,
    GtDecodeBuffers, GtStatsMode, HaplotypeDenseDecodeBuffers,
};
use self::header::read_sample_records_from_header;
use self::output::TextDenseOutput;
use self::record::{
    append_public_variant_metadata_from_text_record, prepare_text_candidate, PreparedTextCandidate,
    TextCandidateAction, TextVariantView,
};
use self::source::{
    ensure_text_indexed_vcf_supported, ensure_text_vcf_supported, open_compressed_reader,
    open_plain_reader, open_text_sample_selection, open_text_vcf_input,
    with_indexed_text_vcf_input, with_threaded_indexed_text_vcf_input, DenseReadSource,
    TextVcfInput, TextVcfSource,
};
use self::sparse::{
    read_haplotype_sparse_records_with_metadata, read_sparse_records_with_metadata,
};

use super::{haplotype_sample_records, is_compressed_vcf};

mod ds;
mod format;
mod gt;
mod header;
mod output;
mod record;
mod session;
mod source;
mod sparse;

pub(in crate::vcf) use self::record::append_public_variant_metadata_from_noodles_variant_record;
pub(crate) use self::session::TextVcfBlockSession;

pub(super) const VCF_TEXT_BUFFER_SIZE: usize = 1 << 20;
const VCF_METADATA_INITIAL_VARIANT_CAPACITY: usize = 4096;
const VCF_TEXT_INITIAL_MATRIX_VARIANT_CAPACITY: usize = 4096;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct VcfMetadataReturn {
    samples: bool,
    variants: bool,
}

impl VcfMetadataReturn {
    fn matrix_only(self) -> bool {
        !self.samples && !self.variants
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantMetadataSinkKind {
    Output,
    None,
}

impl VariantMetadataSinkKind {
    const fn for_output(metadata_return: VcfMetadataReturn) -> Self {
        if metadata_return.variants {
            Self::Output
        } else {
            Self::None
        }
    }
}

enum VariantMetadataSink {
    Output(Box<VariantMetadataBuffers>),
    None,
}

impl VariantMetadataSink {
    fn new(kind: VariantMetadataSinkKind, capacity: usize) -> Self {
        match kind {
            VariantMetadataSinkKind::Output => {
                Self::Output(Box::new(VariantMetadataBuffers::with_capacity(capacity)))
            }
            VariantMetadataSinkKind::None => Self::None,
        }
    }

    fn push_view<V: VariantMetadataView + ?Sized>(&mut self, variant: &V) -> Result<()> {
        match self {
            Self::Output(variants) => variants.push_view(variant)?,
            Self::None => {}
        }
        Ok(())
    }

    fn push_view_with_optional_stats_and_orientation<V: VariantMetadataView + ?Sized>(
        &mut self,
        variant: &V,
        stats: Option<VariantStats>,
        flipped: bool,
    ) -> Result<()> {
        match self {
            Self::Output(variants) => {
                let row_index = variants.len();
                if flipped {
                    variants.push_flipped_view(variant)?;
                } else {
                    variants.push_view(variant)?;
                }
                if let Some(stats) = stats {
                    variants.attach_stats(row_index, stats)?;
                }
            }
            Self::None => {}
        }
        Ok(())
    }

    fn push_view_with_stats<V: VariantMetadataView + ?Sized>(
        &mut self,
        variant: &V,
        stats: VariantStats,
    ) -> Result<()> {
        self.push_view_with_optional_stats_and_orientation(variant, Some(stats), false)
    }

    fn into_output(self) -> Result<Option<VariantMetadataBuffers>> {
        match self {
            Self::Output(variants) => Ok(Some(*variants)),
            Self::None => Ok(None),
        }
    }
}

/// Choose the initial retained-variant capacity for text matrix output.
///
/// Indexed/windowed reads know the exact variant bound; full scans start with a
/// small nonzero capacity so the first retained variants do not repeatedly
/// grow empty vectors.
fn dense_output_variant_capacity(variant_window: Option<VariantWindow>) -> usize {
    variant_window.map_or(VCF_TEXT_INITIAL_MATRIX_VARIANT_CAPACITY, |window| {
        window.len
    })
}

/// Append one decoded dense text variant, using the no-missing fast path when valid.
fn write_dense_text_variant(
    output: &mut TextDenseOutput,
    values: &[f32],
    missing_indices: &[usize],
    missing_policy: DenseMissingPolicy,
) -> Result<()> {
    if missing_indices.is_empty() {
        return output.write_variant_no_missing_direct(values);
    }
    output.write_variant(values, missing_indices, missing_policy)
}

impl TextVcfSource {
    fn read_dense(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
    ) -> Result<DenseGenotypeMatrix> {
        self.read_dense_with_metadata(
            path,
            variant_filter,
            variant_window,
            missing_policy,
            metadata_return,
            VariantMetadataSinkKind::for_output(metadata_return),
        )
    }

    fn read_dense_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => read_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::ThreadedCompressed(input) => read_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::Plain(input) => read_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
        }
    }

    fn read_dosage_dense(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
    ) -> Result<DenseGenotypeMatrix> {
        self.read_dosage_dense_with_metadata(
            path,
            variant_filter,
            variant_window,
            missing_policy,
            metadata_return,
            VariantMetadataSinkKind::for_output(metadata_return),
        )
    }

    fn read_dosage_dense_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => read_dosage_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::ThreadedCompressed(input) => read_dosage_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::Plain(input) => read_dosage_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
        }
    }

    fn read_haplotype_dense(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
    ) -> Result<DenseGenotypeMatrix> {
        self.read_haplotype_dense_with_metadata(
            path,
            variant_filter,
            variant_window,
            missing_policy,
            metadata_return,
            VariantMetadataSinkKind::for_output(metadata_return),
        )
    }

    fn read_haplotype_dense_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => read_haplotype_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::ThreadedCompressed(input) => read_haplotype_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::Plain(input) => read_haplotype_dense_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                variant_sink_kind,
                input,
            ),
        }
    }

    fn read_sparse(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        metadata_return: VcfMetadataReturn,
    ) -> Result<SparseGenotypeMatrix> {
        self.read_sparse_with_metadata(
            path,
            variant_filter,
            variant_window,
            metadata_return,
            VariantMetadataSinkKind::for_output(metadata_return),
        )
    }

    fn read_sparse_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<SparseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => read_sparse_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::ThreadedCompressed(input) => read_sparse_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::Plain(input) => read_sparse_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                metadata_return,
                variant_sink_kind,
                input,
            ),
        }
    }

    fn read_haplotype_sparse(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        metadata_return: VcfMetadataReturn,
    ) -> Result<SparseGenotypeMatrix> {
        self.read_haplotype_sparse_with_metadata(
            path,
            variant_filter,
            variant_window,
            metadata_return,
            VariantMetadataSinkKind::for_output(metadata_return),
        )
    }

    fn read_haplotype_sparse_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<SparseGenotypeMatrix> {
        match self {
            Self::Compressed(input) => read_haplotype_sparse_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::ThreadedCompressed(input) => read_haplotype_sparse_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                metadata_return,
                variant_sink_kind,
                input,
            ),
            Self::Plain(input) => read_haplotype_sparse_with_metadata_from_input(
                path,
                variant_filter,
                variant_window,
                metadata_return,
                variant_sink_kind,
                input,
            ),
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

pub(super) fn read_vcf_public_metadata(path: &Path) -> Result<MetadataOutput> {
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

#[expect(
    clippy::too_many_arguments,
    reason = "VCF Output metadata path keeps sample and variant return choices explicit"
)]
pub(super) fn read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_dense(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

pub(super) fn empty_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_dense_output_from_selection(
        selection,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

pub(super) fn empty_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_sparse_output_from_selection(
        selection,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

pub(super) fn empty_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_haplotype_dense_output_from_selection(
        selection,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

pub(super) fn empty_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_haplotype_sparse_output_from_selection(
        selection,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

fn read_dense_with_metadata_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    input: TextVcfInput<R>,
) -> Result<DenseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        selection,
    } = input;
    read_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        variant_sink_kind,
        DenseReadSource::full_scan(),
        &selection,
        &mut reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Output boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    let metadata_return = VcfMetadataReturn {
        samples: return_samples,
        variants: return_variants,
    };
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
                missing_policy,
                metadata_return,
                input.dense_source(),
                &input.selection,
                reader,
            )
        },
        |selection| empty_dense_output_from_selection(selection, metadata_return),
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "VCF Output metadata path keeps sample and variant return choices explicit"
)]
pub(super) fn read_vcf_dosage_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_dosage_dense(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Output boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_dosage_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    let metadata_return = VcfMetadataReturn {
        samples: return_samples,
        variants: return_variants,
    };
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
                missing_policy,
                metadata_return,
                input.dense_source(),
                &input.selection,
                reader,
            )
        },
        |selection| empty_dense_output_from_selection(selection, metadata_return),
    )
}

fn read_dosage_dense_with_metadata_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    input: TextVcfInput<R>,
) -> Result<DenseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        selection,
    } = input;
    read_dosage_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        variant_sink_kind,
        DenseReadSource::full_scan(),
        &selection,
        &mut reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Output boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_haplotypes_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    let metadata_return = VcfMetadataReturn {
        samples: return_samples,
        variants: return_variants,
    };
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
                missing_policy,
                metadata_return,
                Some(input.region),
                input.selection,
                reader,
            )
        },
        |selection| empty_haplotype_dense_output_from_selection(selection, metadata_return),
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "VCF Output metadata path keeps sample and variant return choices explicit"
)]
pub(super) fn read_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_haplotype_dense(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

fn read_haplotype_dense_with_metadata_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    input: TextVcfInput<R>,
) -> Result<DenseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        selection,
        ..
    } = input;
    read_haplotype_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        variant_sink_kind,
        None,
        selection,
        &mut reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Output boundary carries region and threading explicitly"
)]
pub(super) fn read_vcf_sparse_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    let metadata_return = VcfMetadataReturn {
        samples: return_samples,
        variants: return_variants,
    };
    with_indexed_text_vcf_input_for_threads!(
        path,
        requested_samples,
        region,
        threads,
        |input, reader| {
            read_sparse_records_with_metadata(
                path,
                variant_filter,
                variant_window,
                Some(input.region),
                metadata_return,
                VariantMetadataSinkKind::for_output(metadata_return),
                input.selection,
                reader,
            )
        },
        |selection| empty_sparse_output_from_selection(selection, metadata_return),
    )
}

pub(super) fn read_vcf_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_sparse(
        path,
        variant_filter,
        variant_window,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

fn read_sparse_with_metadata_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    input: TextVcfInput<R>,
) -> Result<SparseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        selection,
        ..
    } = input;
    read_sparse_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        None,
        metadata_return,
        variant_sink_kind,
        selection,
        &mut reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Output boundary carries region and threading explicitly"
)]
pub(super) fn read_vcf_haplotypes_sparse_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_indexed_vcf_supported(path)?;
    let metadata_return = VcfMetadataReturn {
        samples: return_samples,
        variants: return_variants,
    };
    with_indexed_text_vcf_input_for_threads!(
        path,
        requested_samples,
        region,
        threads,
        |input, reader| {
            read_haplotype_sparse_records_with_metadata(
                path,
                variant_filter,
                variant_window,
                Some(input.region),
                metadata_return,
                VariantMetadataSinkKind::for_output(metadata_return),
                input.selection,
                reader,
            )
        },
        |selection| empty_haplotype_sparse_output_from_selection(selection, metadata_return),
    )
}

pub(super) fn read_vcf_haplotypes_sparse(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_haplotype_sparse(
        path,
        variant_filter,
        variant_window,
        VcfMetadataReturn {
            samples: return_samples,
            variants: return_variants,
        },
    )
}

fn read_haplotype_sparse_with_metadata_from_input<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    input: TextVcfInput<R>,
) -> Result<SparseGenotypeMatrix> {
    let TextVcfInput {
        mut reader,
        selection,
        ..
    } = input;
    read_haplotype_sparse_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        None,
        metadata_return,
        variant_sink_kind,
        selection,
        &mut reader,
    )
}

fn empty_dense_output_from_selection(
    selection: DenseSampleSelection,
    metadata_return: VcfMetadataReturn,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len();
    let samples = SampleMetadataBuffers::optional_from_records(
        &selection.samples,
        metadata_return.samples,
        false,
    )?;
    let variants = metadata_return
        .variants
        .then(|| VariantMetadataBuffers::with_capacity(0));
    TextDenseOutput::new(n_samples, 0).finish(0, samples, variants, diagnostics)
}

fn empty_sparse_output_from_selection(
    selection: DenseSampleSelection,
    metadata_return: VcfMetadataReturn,
) -> Result<SparseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len();
    let samples = SampleMetadataBuffers::optional_from_records(
        &selection.samples,
        metadata_return.samples,
        false,
    )?;
    let variants = metadata_return
        .variants
        .then(|| VariantMetadataBuffers::with_capacity(0));
    SparseGenotypeMatrix::new(
        n_samples,
        0,
        vec![0],
        Vec::new(),
        Vec::new(),
        samples,
        variants,
        diagnostics,
    )
}

fn empty_haplotype_dense_output_from_selection(
    selection: DenseSampleSelection,
    metadata_return: VcfMetadataReturn,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len() * 2;
    let haplotype_samples = haplotype_sample_records(&selection.samples, &selection.source_indices);
    let samples = SampleMetadataBuffers::optional_from_records(
        &haplotype_samples,
        metadata_return.samples,
        true,
    )?;
    let variants = metadata_return
        .variants
        .then(|| VariantMetadataBuffers::with_capacity(0));
    TextDenseOutput::new(n_samples, 0).finish(0, samples, variants, diagnostics)
}

fn empty_haplotype_sparse_output_from_selection(
    selection: DenseSampleSelection,
    metadata_return: VcfMetadataReturn,
) -> Result<SparseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len() * 2;
    let haplotype_samples = haplotype_sample_records(&selection.samples, &selection.source_indices);
    let samples = SampleMetadataBuffers::optional_from_records(
        &haplotype_samples,
        metadata_return.samples,
        true,
    )?;
    let variants = metadata_return
        .variants
        .then(|| VariantMetadataBuffers::with_capacity(0));
    SparseGenotypeMatrix::new(
        n_samples,
        0,
        vec![0],
        Vec::new(),
        Vec::new(),
        samples,
        variants,
        diagnostics,
    )
}

fn read_metadata_records<R: BufRead>(
    path: &Path,
    reader: &mut noodles::io::Reader<R>,
) -> Result<MetadataOutput> {
    let samples = read_sample_records_from_header(path, reader.get_mut())?;
    let mut variants = VariantMetadataBuffers::with_capacity(VCF_METADATA_INITIAL_VARIANT_CAPACITY);
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
        append_public_variant_metadata_from_text_record(path, &record, &mut variants)?;
    }

    let capabilities = if has_phased_genotype_evidence {
        SourceCapabilities::phased_genotypes()
    } else {
        SourceCapabilities::genotype_only()
    };

    Ok(MetadataOutput {
        samples: SampleMetadataBuffers::from_records(&samples, false)?,
        variants,
        capabilities,
    })
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_dense_records<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    read_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        VariantMetadataSinkKind::for_output(metadata_return),
        source,
        selection,
        reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_dense_records_with_metadata<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics.clone();
    let n_samples = selection.samples.len();
    let output_variant_capacity = dense_output_variant_capacity(variant_window);
    // Variant-major appends have better locality for text VCF's decoded
    // per-record values. Python exposes the public sample-by-variant shape via
    // layout-aware assembly, so avoid strided sample-major writes here.
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity);
    let mut output_variant_count = 0;
    let mut variants = VariantMetadataSink::new(variant_sink_kind, output_variant_capacity);
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

        // Indexed chunks can include neighboring records; shared candidate
        // preparation retains exact region post-filtering before any decode.
        let prepared = match prepare_text_candidate(
            path,
            &record,
            source.region,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            TextCandidateAction::Skip => continue,
            TextCandidateAction::Stop => break,
            TextCandidateAction::Decode(prepared) => prepared,
        };
        // `record.samples()` borrows noodles' reusable record buffer, so the
        // shared transition decodes selected GTs before the next read.
        let (variant, stats) = match process_text_gt_candidate(
            path,
            &record,
            &selection.source_indices,
            prepared,
            variant_filter,
            &mut retention,
            &mut diagnostics,
            metadata_return.matrix_only(),
            true,
            "GT",
            &mut decoded,
        )? {
            DecodedTextCandidate::Include { variant, stats } => (variant, stats),
            DecodedTextCandidate::Skip => continue,
            DecodedTextCandidate::Stop => break,
        };
        if let Some(stats) = stats {
            variants.push_view_with_stats(&variant, stats)?;
        } else {
            variants.push_view(&variant)?;
        }

        write_dense_text_variant(
            &mut output,
            decoded.values(),
            decoded.missing_indices(),
            missing_policy,
        )?;
        output_variant_count += 1;
    }

    diagnostics.retained_variants = output_variant_count;
    let samples = SampleMetadataBuffers::optional_from_records(
        &selection.samples,
        metadata_return.samples,
        false,
    )?;
    output.finish(
        output_variant_count,
        samples,
        variants.into_output()?,
        diagnostics,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_dosage_dense_records<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    read_dosage_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        VariantMetadataSinkKind::for_output(metadata_return),
        source,
        selection,
        reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_dosage_dense_records_with_metadata<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    let mut diagnostics = selection.diagnostics.clone();
    let n_samples = selection.samples.len();
    let output_variant_capacity = dense_output_variant_capacity(variant_window);
    // Keep text VCF dosage output aligned with the genotype path: append
    // retained variants contiguously and let the adapter expose public shape.
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity);
    let mut output_variant_count = 0;
    let mut variants = VariantMetadataSink::new(variant_sink_kind, output_variant_capacity);
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

        let prepared = match prepare_text_candidate(
            path,
            &record,
            source.region,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            TextCandidateAction::Skip => continue,
            TextCandidateAction::Stop => break,
            TextCandidateAction::Decode(prepared) => prepared,
        };
        let (variant, stats) = match process_text_ds_candidate(
            path,
            &record,
            &selection.source_indices,
            prepared,
            variant_filter,
            &mut retention,
            &mut diagnostics,
            metadata_return.matrix_only(),
            &mut decoded,
        )? {
            DecodedTextCandidate::Include { variant, stats } => (variant, stats),
            DecodedTextCandidate::Skip => continue,
            DecodedTextCandidate::Stop => break,
        };
        if let Some(stats) = stats {
            variants.push_view_with_stats(&variant, stats)?;
        } else {
            variants.push_view(&variant)?;
        }

        write_dense_text_variant(
            &mut output,
            decoded.values(),
            decoded.missing_indices(),
            missing_policy,
        )?;
        output_variant_count += 1;
    }

    diagnostics.retained_variants = output_variant_count;
    let samples = SampleMetadataBuffers::optional_from_records(
        &selection.samples,
        metadata_return.samples,
        false,
    )?;
    output.finish(
        output_variant_count,
        samples,
        variants.into_output()?,
        diagnostics,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_haplotype_dense_records<R: std::io::BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    source_region: Option<&RegionPredicate>,
    selection: DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    read_haplotype_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        VariantMetadataSinkKind::for_output(metadata_return),
        source_region,
        selection,
        reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_haplotype_dense_records_with_metadata<R: std::io::BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    variant_sink_kind: VariantMetadataSinkKind,
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
    let output_variant_capacity = dense_output_variant_capacity(variant_window);
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity);
    let mut variants = VariantMetadataSink::new(variant_sink_kind, output_variant_capacity);
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

        let prepared = match prepare_text_candidate(
            path,
            &record,
            source_region,
            variant_filter,
            &mut retention,
            &mut diagnostics,
        )? {
            TextCandidateAction::Skip => continue,
            TextCandidateAction::Stop => break,
            TextCandidateAction::Decode(prepared) => prepared,
        };
        // Genotype-stat filters are evaluated on diploid dosage before
        // enforcing phased output. Rejected unphased records therefore do not
        // surface a haplotype decode error.
        let (variant, stats) = match process_text_gt_candidate(
            path,
            &record,
            &source_indices,
            prepared,
            variant_filter,
            &mut retention,
            &mut diagnostics,
            metadata_return.matrix_only(),
            false,
            "haplotype",
            &mut stats_decoded,
        )? {
            DecodedTextCandidate::Include { variant, stats } => (variant, stats),
            DecodedTextCandidate::Skip => continue,
            DecodedTextCandidate::Stop => break,
        };
        if let Some(stats) = stats {
            variants.push_view_with_stats(&variant, stats)?;
        } else {
            variants.push_view(&variant)?;
        }
        decode_phased_gt_dense_record(
            path,
            &record,
            &source_indices,
            GtStatsMode::Skip,
            &mut decoded,
        )?;

        // Dense haplotype reads expose source allele-1 indicators. Sparse
        // haplotype output flips columns to minor allele later to reduce nnz.
        write_dense_text_variant(
            &mut output,
            decoded.values(),
            decoded.missing_indices(),
            missing_policy,
        )?;
        output_variant_count += 1;
    }

    let haplotype_samples = haplotype_sample_records(&samples, &source_indices);
    let samples = SampleMetadataBuffers::optional_from_records(
        &haplotype_samples,
        metadata_return.samples,
        true,
    )?;
    diagnostics.retained_variants = output_variant_count;
    output.finish(
        output_variant_count,
        samples,
        variants.into_output()?,
        diagnostics,
    )
}

fn evaluate_text_gt_filter<V: VariantMetadataView + ?Sized>(
    decoded: &GtDecodeBuffers,
    filter: &VariantFilter,
    variant: &V,
    matrix_only: bool,
    context: &str,
) -> Result<(bool, Option<VariantStats>)> {
    if matrix_only {
        let counts = decoded.counts().ok_or_else(|| {
            GenoioError::internal_contract(format!(
                "matrix-only vcf {context} filter missing counts"
            ))
        })?;
        return evaluate_hardcall_counts_filter(
            counts,
            filter,
            filter.genotype_filter_plan(),
            Some(variant),
            false,
        );
    }

    let stats = decoded.stats().ok_or_else(|| {
        GenoioError::internal_contract(format!("vcf {context} filter missing stats"))
    })?;
    Ok((filter.evaluate_view(variant, Some(&stats)), Some(stats)))
}

enum DecodedTextCandidate<'a> {
    Include {
        variant: TextVariantView<'a>,
        stats: Option<VariantStats>,
    },
    Skip,
    Stop,
}

#[expect(
    clippy::too_many_arguments,
    reason = "shared text-VCF transition receives borrowed record state plus reusable decode buffers"
)]
fn process_text_gt_candidate<'a>(
    path: &Path,
    record: &'a noodles::Record,
    source_indices: &[usize],
    prepared: PreparedTextCandidate<'a>,
    variant_filter: Option<&VariantFilter>,
    retention: &mut RetainedVariantState,
    diagnostics: &mut DenseDiagnostics,
    matrix_only: bool,
    decode_when_unfiltered: bool,
    context: &str,
    decoded: &mut GtDecodeBuffers,
) -> Result<DecodedTextCandidate<'a>> {
    let variant = prepared.variant;
    let needs_genotype_decision = prepared.needs_genotype_decision;
    if !needs_genotype_decision && !decode_when_unfiltered {
        return Ok(DecodedTextCandidate::Include {
            variant,
            stats: None,
        });
    }
    let stats_mode = match (needs_genotype_decision, matrix_only) {
        (true, true) => GtStatsMode::Counts,
        (true, false) => GtStatsMode::Compute,
        (false, _) => GtStatsMode::Skip,
    };
    decode_gt_record(path, record, source_indices, stats_mode, decoded)?;

    let stats = if needs_genotype_decision {
        let filter = variant_filter.ok_or_else(|| {
            GenoioError::internal_contract("genotype decision requires a variant filter")
        })?;
        let (retain_variant, stats) =
            evaluate_text_gt_filter(decoded, filter, &variant, matrix_only, context)?;
        match retention.genotype_decision(retain_variant, diagnostics) {
            RetentionAction::Include => stats,
            RetentionAction::Skip => return Ok(DecodedTextCandidate::Skip),
            RetentionAction::Stop => return Ok(DecodedTextCandidate::Stop),
        }
    } else {
        None
    };
    Ok(DecodedTextCandidate::Include { variant, stats })
}

#[expect(
    clippy::too_many_arguments,
    reason = "shared text-VCF dosage transition receives borrowed record state plus reusable decode buffers"
)]
fn process_text_ds_candidate<'a>(
    path: &Path,
    record: &'a noodles::Record,
    source_indices: &[usize],
    prepared: PreparedTextCandidate<'a>,
    variant_filter: Option<&VariantFilter>,
    retention: &mut RetainedVariantState,
    diagnostics: &mut DenseDiagnostics,
    matrix_only: bool,
    decoded: &mut DsDecodeBuffers,
) -> Result<DecodedTextCandidate<'a>> {
    let variant = prepared.variant;
    decode_ds_record(path, record, source_indices, decoded)?;

    let stats = if prepared.needs_genotype_decision {
        let filter = variant_filter.ok_or_else(|| {
            GenoioError::internal_contract("genotype decision requires a variant filter")
        })?;
        let (retain_variant, stats) = evaluate_dosage_filter(
            decoded.values(),
            decoded.missing_indices(),
            filter,
            &variant,
            !matrix_only,
        )?;
        match retention.genotype_decision(retain_variant, diagnostics) {
            RetentionAction::Include => stats,
            RetentionAction::Skip => return Ok(DecodedTextCandidate::Skip),
            RetentionAction::Stop => return Ok(DecodedTextCandidate::Stop),
        }
    } else {
        None
    };
    Ok(DecodedTextCandidate::Include { variant, stats })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::path::Path;

    use super::*;

    #[test]
    fn dense_text_backend_accepts_metadata_reads() {
        assert!(ensure_text_vcf_supported(Path::new("example.vcf.gz"), None).is_ok());
    }

    #[test]
    fn pbr_rust_textvcf_003_stateless_and_persistent_paths_share_gt_transition() {
        let path = Path::new("fixture.vcf");
        let mut reader =
            noodles::io::Reader::new(Cursor::new(b"1\t42\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n"));
        let mut record = noodles::Record::default();
        reader
            .read_record(&mut record)
            .expect("text VCF record should parse");
        let mut retention = RetainedVariantState::new(None);
        let mut diagnostics = genoio_core::DenseDiagnostics::default();
        let prepared = match prepare_text_candidate(
            path,
            &record,
            None,
            None,
            &mut retention,
            &mut diagnostics,
        )
        .expect("candidate preparation should succeed")
        {
            TextCandidateAction::Decode(prepared) => prepared,
            TextCandidateAction::Skip | TextCandidateAction::Stop => {
                panic!("unfiltered candidate should decode")
            }
        };
        let mut decoded = GtDecodeBuffers::with_capacity(1);

        let action = process_text_gt_candidate(
            path,
            &record,
            &[0],
            prepared,
            None,
            &mut retention,
            &mut diagnostics,
            false,
            true,
            "GT",
            &mut decoded,
        )
        .expect("shared GT transition should decode");

        assert!(matches!(action, DecodedTextCandidate::Include { .. }));
        assert_eq!(decoded.values(), &[1.0]);
    }
}
