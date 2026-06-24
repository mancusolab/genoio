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
    DenseGenotypeMatrix, DenseGenotypeMatrixArrowVariants, DenseMissingPolicy,
    DenseSampleSelection, GenoioError, MetadataArrowOutput, MetadataOutput, PartialFilterDecision,
    RegionPredicate, SourceCapabilities, SparseGenotypeMatrix, SparseGenotypeMatrixArrowVariants,
    VariantFilter, VariantMetadataArrowBuffers, VariantRecord, VariantStats, VariantWindow,
};
use noodles_vcf as noodles;

use crate::dosage_filter::evaluate_dosage_filter;
use crate::error::Result;
use crate::hardcall::evaluate_hardcall_counts_filter;
use crate::retention::{MetadataRetentionAction, RetainedVariantState, RetentionAction};

use self::ds::{decode_ds_record, DsDecodeBuffers};
use self::gt::{
    decode_gt_record, decode_phased_gt_dense_record, text_record_has_phased_genotype,
    GtDecodeBuffers, GtStatsMode, HaplotypeDenseDecodeBuffers,
};
use self::header::read_sample_records_from_header;
use self::output::{can_write_sample_major_directly, TextDenseOutput};
use self::record::{
    append_public_variant_metadata_from_text_record, skip_variant_for_region,
    validate_biallelic_variant, variant_record_from_text_record,
};
use self::source::{
    ensure_text_indexed_vcf_supported, ensure_text_vcf_supported, open_compressed_reader,
    open_plain_reader, open_text_sample_selection, open_text_vcf_input,
    with_indexed_text_vcf_input, with_threaded_indexed_text_vcf_input, DenseReadSource,
    TextVcfInput, TextVcfSource,
};
use self::sparse::{
    read_haplotype_sparse_records, read_haplotype_sparse_records_with_metadata,
    read_sparse_records, read_sparse_records_with_metadata, TextSparseReadOutput,
};

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
const VCF_METADATA_INITIAL_VARIANT_CAPACITY: usize = 4096;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct VcfMetadataReturn {
    samples: bool,
    variants: bool,
}

impl VcfMetadataReturn {
    const fn from_matrix_only(matrix_only: bool) -> Self {
        Self {
            samples: !matrix_only,
            variants: !matrix_only,
        }
    }

    fn matrix_only(self) -> bool {
        !self.samples && !self.variants
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantMetadataSinkKind {
    Arrow,
    None,
}

impl VariantMetadataSinkKind {
    const fn for_arrow_output(metadata_return: VcfMetadataReturn) -> Self {
        if metadata_return.variants {
            Self::Arrow
        } else {
            Self::None
        }
    }
}

enum VariantMetadataSink {
    Arrow(Box<VariantMetadataArrowBuffers>),
    None,
}

impl VariantMetadataSink {
    fn new(kind: VariantMetadataSinkKind, capacity: usize) -> Self {
        match kind {
            VariantMetadataSinkKind::Arrow => Self::Arrow(Box::new(
                VariantMetadataArrowBuffers::with_capacity(capacity),
            )),
            VariantMetadataSinkKind::None => Self::None,
        }
    }

    fn push_variant(&mut self, variant: VariantRecord) -> Result<()> {
        match self {
            Self::Arrow(variants) => variants.push_record(&variant)?,
            Self::None => {}
        }
        Ok(())
    }

    fn into_arrow(self) -> Result<Option<VariantMetadataArrowBuffers>> {
        match self {
            Self::Arrow(variants) => Ok(Some(*variants)),
            Self::None => Ok(None),
        }
    }
}

enum TextDenseReadOutput {
    Arrow(DenseGenotypeMatrixArrowVariants),
}

impl TextDenseReadOutput {
    fn into_records_or_matrix_only(
        self,
        _matrix_only: bool,
        context: &'static str,
    ) -> Result<DenseGenotypeMatrix> {
        match self {
            Self::Arrow(output) => output.into_matrix().map_err(|error| {
                GenoioError::internal_contract(format!(
                    "text VCF {context} Arrow-to-row compatibility conversion failed: {error}"
                ))
            }),
        }
    }

    fn into_arrow(self) -> DenseGenotypeMatrixArrowVariants {
        match self {
            Self::Arrow(output) => output,
        }
    }
}

impl TextVcfSource {
    fn read_dense(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        let metadata_return = VcfMetadataReturn::from_matrix_only(matrix_only);
        self.read_dense_with_metadata(
            path,
            variant_filter,
            variant_window,
            missing_policy,
            metadata_return,
            VariantMetadataSinkKind::for_arrow_output(metadata_return),
        )?
        .into_records_or_matrix_only(matrix_only, "dense")
    }

    fn read_dense_arrow_variants(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
    ) -> Result<DenseGenotypeMatrixArrowVariants> {
        Ok(self
            .read_dense_with_metadata(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                VariantMetadataSinkKind::for_arrow_output(metadata_return),
            )?
            .into_arrow())
    }

    fn read_dense_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<TextDenseReadOutput> {
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
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        let metadata_return = VcfMetadataReturn::from_matrix_only(matrix_only);
        self.read_dosage_dense_with_metadata(
            path,
            variant_filter,
            variant_window,
            missing_policy,
            metadata_return,
            VariantMetadataSinkKind::for_arrow_output(metadata_return),
        )?
        .into_records_or_matrix_only(matrix_only, "dosage dense")
    }

    fn read_dosage_dense_arrow_variants(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
    ) -> Result<DenseGenotypeMatrixArrowVariants> {
        Ok(self
            .read_dosage_dense_with_metadata(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                VariantMetadataSinkKind::for_arrow_output(metadata_return),
            )?
            .into_arrow())
    }

    fn read_dosage_dense_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<TextDenseReadOutput> {
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
        matrix_only: bool,
    ) -> Result<DenseGenotypeMatrix> {
        let metadata_return = VcfMetadataReturn::from_matrix_only(matrix_only);
        self.read_haplotype_dense_with_metadata(
            path,
            variant_filter,
            variant_window,
            missing_policy,
            metadata_return,
            VariantMetadataSinkKind::for_arrow_output(metadata_return),
        )?
        .into_records_or_matrix_only(matrix_only, "haplotype dense")
    }

    fn read_haplotype_dense_arrow_variants(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
    ) -> Result<DenseGenotypeMatrixArrowVariants> {
        Ok(self
            .read_haplotype_dense_with_metadata(
                path,
                variant_filter,
                variant_window,
                missing_policy,
                metadata_return,
                VariantMetadataSinkKind::for_arrow_output(metadata_return),
            )?
            .into_arrow())
    }

    fn read_haplotype_dense_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        missing_policy: DenseMissingPolicy,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<TextDenseReadOutput> {
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
    ) -> Result<SparseGenotypeMatrix> {
        self.read_sparse_with_metadata(
            path,
            variant_filter,
            variant_window,
            VcfMetadataReturn::from_matrix_only(false),
            VariantMetadataSinkKind::Arrow,
        )?
        .into_legacy_matrix("sparse")
    }

    fn read_sparse_arrow_variants(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        metadata_return: VcfMetadataReturn,
    ) -> Result<SparseGenotypeMatrixArrowVariants> {
        Ok(self
            .read_sparse_with_metadata(
                path,
                variant_filter,
                variant_window,
                metadata_return,
                VariantMetadataSinkKind::for_arrow_output(metadata_return),
            )?
            .into_arrow())
    }

    fn read_sparse_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<TextSparseReadOutput> {
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
    ) -> Result<SparseGenotypeMatrix> {
        self.read_haplotype_sparse_with_metadata(
            path,
            variant_filter,
            variant_window,
            VcfMetadataReturn::from_matrix_only(false),
            VariantMetadataSinkKind::Arrow,
        )?
        .into_legacy_matrix("haplotype sparse")
    }

    fn read_haplotype_sparse_arrow_variants(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        metadata_return: VcfMetadataReturn,
    ) -> Result<SparseGenotypeMatrixArrowVariants> {
        Ok(self
            .read_haplotype_sparse_with_metadata(
                path,
                variant_filter,
                variant_window,
                metadata_return,
                VariantMetadataSinkKind::for_arrow_output(metadata_return),
            )?
            .into_arrow())
    }

    fn read_haplotype_sparse_with_metadata(
        self,
        path: &Path,
        variant_filter: Option<&VariantFilter>,
        variant_window: Option<VariantWindow>,
        metadata_return: VcfMetadataReturn,
        variant_sink_kind: VariantMetadataSinkKind,
    ) -> Result<TextSparseReadOutput> {
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

pub(super) fn read_vcf_public_metadata_arrow(path: &Path) -> Result<MetadataArrowOutput> {
    // Metadata reads can avoid strict full-header parsing. They only need the
    // #CHROM line plus record fields, which keeps real-world VCF headers from
    // forcing a strict parser onto the hot path.
    if is_compressed_vcf(path) {
        let mut reader = open_compressed_reader(path)?;
        read_metadata_arrow_records(path, &mut reader)
    } else {
        let mut reader = open_plain_reader(path)?;
        read_metadata_arrow_records(path, &mut reader)
    }
}

pub(super) fn read_vcf_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_dense(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        matrix_only,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "VCF Arrow metadata path keeps sample and variant return choices explicit"
)]
pub(super) fn read_vcf_dense_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_dense_arrow_variants(
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
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_dense_arrow_from_selection(selection, VcfMetadataReturn::from_matrix_only(matrix_only))?
        .into_matrix()
}

pub(super) fn empty_vcf_dense_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_dense_arrow_from_selection(
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
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_sparse_arrow_from_selection(selection, VcfMetadataReturn::from_matrix_only(false))?
        .into_matrix()
}

pub(super) fn empty_vcf_sparse_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_sparse_arrow_from_selection(
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
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_haplotype_dense_arrow_from_selection(
        selection,
        VcfMetadataReturn::from_matrix_only(matrix_only),
    )?
    .into_matrix()
}

pub(super) fn empty_vcf_haplotypes_dense_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_haplotype_dense_arrow_from_selection(
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
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_haplotype_sparse_arrow_from_selection(
        selection,
        VcfMetadataReturn::from_matrix_only(false),
    )?
    .into_matrix()
}

pub(super) fn empty_vcf_haplotypes_sparse_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    let selection = open_text_sample_selection(path, requested_samples, threads)?;
    empty_haplotype_sparse_arrow_from_selection(
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
) -> Result<TextDenseReadOutput> {
    let TextVcfInput {
        mut reader,
        source_sample_count,
        selection,
    } = input;
    read_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        variant_sink_kind,
        DenseReadSource::full_scan(source_sample_count),
        &selection,
        &mut reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
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
                missing_policy,
                matrix_only,
                input.dense_source(),
                &input.selection,
                reader,
            )
        },
        |selection| {
            empty_dense_arrow_from_selection(
                selection,
                VcfMetadataReturn::from_matrix_only(matrix_only),
            )?
            .into_matrix()
        },
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Arrow boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_dense_arrow_variants_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
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
            read_dense_records_arrow_variants(
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
        |selection| empty_dense_arrow_from_selection(selection, metadata_return),
    )
}

pub(super) fn read_vcf_dosage_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_dosage_dense(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        matrix_only,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "VCF Arrow metadata path keeps sample and variant return choices explicit"
)]
pub(super) fn read_vcf_dosage_dense_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_dosage_dense_arrow_variants(
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
    reason = "indexed VCF boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_dosage_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
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
                missing_policy,
                matrix_only,
                input.dense_source(),
                &input.selection,
                reader,
            )
        },
        |selection| {
            empty_dense_arrow_from_selection(
                selection,
                VcfMetadataReturn::from_matrix_only(matrix_only),
            )?
            .into_matrix()
        },
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Arrow boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_dosage_dense_arrow_variants_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
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
            read_dosage_dense_records_arrow_variants(
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
        |selection| empty_dense_arrow_from_selection(selection, metadata_return),
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
) -> Result<TextDenseReadOutput> {
    let TextVcfInput {
        mut reader,
        source_sample_count,
        selection,
    } = input;
    read_dosage_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        variant_sink_kind,
        DenseReadSource::full_scan(source_sample_count),
        &selection,
        &mut reader,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_haplotypes_dense_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
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
                missing_policy,
                matrix_only,
                Some(input.region),
                input.selection,
                reader,
            )
        },
        |selection| {
            empty_haplotype_dense_arrow_from_selection(
                selection,
                VcfMetadataReturn::from_matrix_only(matrix_only),
            )?
            .into_matrix()
        },
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Arrow boundary carries region, threading, and dense missing policy explicitly"
)]
pub(super) fn read_vcf_haplotypes_dense_arrow_variants_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
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
            read_haplotype_dense_records_arrow_variants(
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
        |selection| empty_haplotype_dense_arrow_from_selection(selection, metadata_return),
    )
}

pub(super) fn read_vcf_haplotypes_dense(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    matrix_only: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrix> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_haplotype_dense(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        matrix_only,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "VCF Arrow metadata path keeps sample and variant return choices explicit"
)]
pub(super) fn read_vcf_haplotypes_dense_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_haplotype_dense_arrow_variants(
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
) -> Result<TextDenseReadOutput> {
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
        |selection| {
            empty_sparse_arrow_from_selection(
                selection,
                VcfMetadataReturn::from_matrix_only(false),
            )?
            .into_matrix()
        },
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Arrow boundary carries region and threading explicitly"
)]
pub(super) fn read_vcf_sparse_arrow_variants_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
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
                VariantMetadataSinkKind::for_arrow_output(metadata_return),
                input.selection,
                reader,
            )
            .map(TextSparseReadOutput::into_arrow)
        },
        |selection| empty_sparse_arrow_from_selection(selection, metadata_return),
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

pub(super) fn read_vcf_sparse_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_sparse_arrow_variants(
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
) -> Result<TextSparseReadOutput> {
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
        |selection| {
            empty_haplotype_sparse_arrow_from_selection(
                selection,
                VcfMetadataReturn::from_matrix_only(false),
            )?
            .into_matrix()
        },
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "indexed VCF Arrow boundary carries region and threading explicitly"
)]
pub(super) fn read_vcf_haplotypes_sparse_arrow_variants_indexed(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    region: &RegionPredicate,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
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
                VariantMetadataSinkKind::for_arrow_output(metadata_return),
                input.selection,
                reader,
            )
            .map(TextSparseReadOutput::into_arrow)
        },
        |selection| empty_haplotype_sparse_arrow_from_selection(selection, metadata_return),
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

pub(super) fn read_vcf_haplotypes_sparse_arrow_variants(
    path: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
    threads: Option<usize>,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    ensure_text_vcf_supported(path, threads)?;
    open_text_vcf_input(path, requested_samples, threads)?.read_haplotype_sparse_arrow_variants(
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
) -> Result<TextSparseReadOutput> {
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

fn empty_dense_arrow_from_selection(
    selection: DenseSampleSelection,
    metadata_return: VcfMetadataReturn,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len();
    let samples = if metadata_return.samples {
        selection.samples
    } else {
        Vec::new()
    };
    let variants = metadata_return
        .variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(0));
    TextDenseOutput::new(n_samples, 0, false).finish_arrow_variants(
        0,
        samples,
        variants,
        diagnostics,
    )
}

fn empty_sparse_arrow_from_selection(
    selection: DenseSampleSelection,
    metadata_return: VcfMetadataReturn,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len();
    let samples = if metadata_return.samples {
        selection.samples
    } else {
        Vec::new()
    };
    let variants = metadata_return
        .variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(0));
    SparseGenotypeMatrixArrowVariants::new(
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

fn empty_haplotype_dense_arrow_from_selection(
    selection: DenseSampleSelection,
    metadata_return: VcfMetadataReturn,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len() * 2;
    let samples = if metadata_return.samples {
        haplotype_sample_records(&selection.samples, &selection.source_indices)
    } else {
        Vec::new()
    };
    let variants = metadata_return
        .variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(0));
    TextDenseOutput::new(n_samples, 0, false).finish_arrow_variants(
        0,
        samples,
        variants,
        diagnostics,
    )
}

fn empty_haplotype_sparse_arrow_from_selection(
    selection: DenseSampleSelection,
    metadata_return: VcfMetadataReturn,
) -> Result<SparseGenotypeMatrixArrowVariants> {
    let mut diagnostics = selection.diagnostics;
    diagnostics.retained_variants = 0;
    let n_samples = selection.samples.len() * 2;
    let samples = if metadata_return.samples {
        haplotype_sample_records(&selection.samples, &selection.source_indices)
    } else {
        Vec::new()
    };
    let variants = metadata_return
        .variants
        .then(|| VariantMetadataArrowBuffers::with_capacity(0));
    SparseGenotypeMatrixArrowVariants::new(
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

fn read_metadata_arrow_records<R: BufRead>(
    path: &Path,
    reader: &mut noodles::io::Reader<R>,
) -> Result<MetadataArrowOutput> {
    let samples = read_sample_records_from_header(path, reader.get_mut())?;
    let mut variants =
        VariantMetadataArrowBuffers::with_capacity(VCF_METADATA_INITIAL_VARIANT_CAPACITY);
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

    Ok(MetadataArrowOutput {
        samples,
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
    matrix_only: bool,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    read_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        VcfMetadataReturn::from_matrix_only(matrix_only),
        VariantMetadataSinkKind::for_arrow_output(VcfMetadataReturn::from_matrix_only(matrix_only)),
        source,
        selection,
        reader,
    )?
    .into_records_or_matrix_only(matrix_only, "dense")
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_dense_records_arrow_variants<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    Ok(read_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        VariantMetadataSinkKind::for_arrow_output(metadata_return),
        source,
        selection,
        reader,
    )?
    .into_arrow())
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
) -> Result<TextDenseReadOutput> {
    let mut diagnostics = selection.diagnostics.clone();
    let n_samples = selection.samples.len();
    let output_variant_capacity = variant_window.map_or(0, |window| window.len);
    let direct_sample_major = metadata_return.matrix_only()
        && variant_window.is_some()
        && can_write_sample_major_directly(selection, source.sample_count, variant_filter);
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity, direct_sample_major);
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
        let stats_mode = match (needs_genotype_decision, metadata_return.matrix_only()) {
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
            let (retain_variant, stats) = evaluate_text_gt_filter(
                &decoded,
                filter,
                &variant,
                metadata_return.matrix_only(),
                "GT",
            )?;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if let Some(stats) = stats {
                genoio_core::attach_variant_stats(&mut variant, stats);
            }
        }
        variants.push_variant(variant)?;

        output.write_variant(
            output_variant_count,
            decoded.values(),
            decoded.missing_indices(),
            missing_policy,
        )?;
        output_variant_count += 1;
    }

    diagnostics.retained_variants = output_variant_count;
    let samples = if metadata_return.samples {
        selection.samples.clone()
    } else {
        Vec::new()
    };
    output
        .finish_arrow_variants(
            output_variant_count,
            samples,
            variants.into_arrow()?,
            diagnostics,
        )
        .map(TextDenseReadOutput::Arrow)
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
    matrix_only: bool,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    read_dosage_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        VcfMetadataReturn::from_matrix_only(matrix_only),
        VariantMetadataSinkKind::for_arrow_output(VcfMetadataReturn::from_matrix_only(matrix_only)),
        source,
        selection,
        reader,
    )?
    .into_records_or_matrix_only(matrix_only, "dosage dense")
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_dosage_dense_records_arrow_variants<R: BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    source: DenseReadSource<'_>,
    selection: &genoio_core::DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    Ok(read_dosage_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        VariantMetadataSinkKind::for_arrow_output(metadata_return),
        source,
        selection,
        reader,
    )?
    .into_arrow())
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
) -> Result<TextDenseReadOutput> {
    let mut diagnostics = selection.diagnostics.clone();
    let n_samples = selection.samples.len();
    let output_variant_capacity = variant_window.map_or(0, |window| window.len);
    let direct_sample_major = metadata_return.matrix_only()
        && variant_window.is_some()
        && can_write_sample_major_directly(selection, source.sample_count, variant_filter);
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity, direct_sample_major);
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
                decoded.missing_indices(),
                filter,
                &variant,
                !metadata_return.matrix_only(),
            )?;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if !metadata_return.matrix_only() {
                if let Some(stats) = stats {
                    genoio_core::attach_variant_stats(&mut variant, stats);
                }
            }
        }
        variants.push_variant(variant)?;

        output.write_variant(
            output_variant_count,
            decoded.values(),
            decoded.missing_indices(),
            missing_policy,
        )?;
        output_variant_count += 1;
    }

    diagnostics.retained_variants = output_variant_count;
    let samples = if metadata_return.samples {
        selection.samples.clone()
    } else {
        Vec::new()
    };
    output
        .finish_arrow_variants(
            output_variant_count,
            samples,
            variants.into_arrow()?,
            diagnostics,
        )
        .map(TextDenseReadOutput::Arrow)
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
    matrix_only: bool,
    source_region: Option<&RegionPredicate>,
    selection: DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrix> {
    read_haplotype_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        VcfMetadataReturn::from_matrix_only(matrix_only),
        VariantMetadataSinkKind::for_arrow_output(VcfMetadataReturn::from_matrix_only(matrix_only)),
        source_region,
        selection,
        reader,
    )?
    .into_records_or_matrix_only(matrix_only, "haplotype dense")
}

#[expect(
    clippy::too_many_arguments,
    reason = "record loop receives prevalidated output mode, sample selection, and reader state"
)]
fn read_haplotype_dense_records_arrow_variants<R: std::io::BufRead>(
    path: &Path,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    metadata_return: VcfMetadataReturn,
    source_region: Option<&RegionPredicate>,
    selection: DenseSampleSelection,
    reader: &mut noodles::io::Reader<R>,
) -> Result<DenseGenotypeMatrixArrowVariants> {
    Ok(read_haplotype_dense_records_with_metadata(
        path,
        variant_filter,
        variant_window,
        missing_policy,
        metadata_return,
        VariantMetadataSinkKind::for_arrow_output(metadata_return),
        source_region,
        selection,
        reader,
    )?
    .into_arrow())
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
) -> Result<TextDenseReadOutput> {
    let DenseSampleSelection {
        source_indices,
        samples,
        mut diagnostics,
    } = selection;
    let n_samples = samples.len() * 2;
    let output_variant_capacity = variant_window.map_or(0, |window| window.len);
    let mut output = TextDenseOutput::new(n_samples, output_variant_capacity, false);
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
            let stats_mode = if metadata_return.matrix_only() {
                GtStatsMode::Counts
            } else {
                GtStatsMode::Compute
            };
            decode_gt_record(
                path,
                &record,
                &source_indices,
                stats_mode,
                &mut stats_decoded,
            )?;
            let filter = variant_filter.ok_or_else(|| {
                GenoioError::internal_contract("genotype decision requires a variant filter")
            })?;
            let (retain_variant, stats) = evaluate_text_gt_filter(
                &stats_decoded,
                filter,
                &variant,
                metadata_return.matrix_only(),
                "haplotype",
            )?;
            match retention.genotype_decision(retain_variant, &mut diagnostics) {
                RetentionAction::Include => {}
                RetentionAction::Skip => continue,
                RetentionAction::Stop => break,
            }
            if !metadata_return.matrix_only() {
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
        variants.push_variant(variant)?;

        // Dense haplotype reads expose source allele-1 indicators. Sparse
        // haplotype output flips columns to minor allele later to reduce nnz.
        output.write_variant(
            output_variant_count,
            decoded.values(),
            decoded.missing_indices(),
            missing_policy,
        )?;
        output_variant_count += 1;
    }

    let samples = if metadata_return.samples {
        haplotype_sample_records(&samples, &source_indices)
    } else {
        Vec::new()
    };
    diagnostics.retained_variants = output_variant_count;
    output
        .finish_arrow_variants(
            output_variant_count,
            samples,
            variants.into_arrow()?,
            diagnostics,
        )
        .map(TextDenseReadOutput::Arrow)
}

fn evaluate_text_gt_filter(
    decoded: &GtDecodeBuffers,
    filter: &VariantFilter,
    variant: &VariantRecord,
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
    Ok((filter.evaluate(variant, Some(&stats)), Some(stats)))
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
