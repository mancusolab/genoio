// pattern: Imperative Shell

use std::fmt;
use std::path::PathBuf;

use genoio_core::{
    DenseDiagnostics, DenseGenotypeMatrix, DenseMissingPolicy, GenoioError, SparseGenotypeMatrix,
    VariantFilter,
};

use crate::bgen::BgenBlockSession;
use crate::error::Result;
use crate::plink::{Plink1BlockSession, Plink2BlockSession};
use crate::vcf::TextVcfBlockSession;

/// Owned source paths for one persistent block-reader session.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BlockSource {
    /// Variant Call Format input.
    Vcf { vcf: PathBuf },
    /// Binary Call Format input.
    Bcf { bcf: PathBuf },
    /// PLINK 1 binary genotype input and its metadata companions.
    Plink1 {
        bed: PathBuf,
        bim: PathBuf,
        fam: PathBuf,
    },
    /// PLINK 2 genotype input and its metadata companions.
    Plink2 {
        pgen: PathBuf,
        pvar: PathBuf,
        psam: PathBuf,
    },
    /// BGEN input and its optional sample metadata companion.
    Bgen {
        bgen: PathBuf,
        sample: Option<PathBuf>,
    },
}

/// Biological row representation requested from a block reader.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixKind {
    /// One row per selected sample.
    Genotype,
    /// One row per selected sample haplotype.
    Haplotype,
}

/// Source field used to populate matrix values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DosageSource {
    /// Decode discrete hard-call values.
    Hardcall,
    /// Decode dosage values supplied by the source.
    Dosage,
}

/// Backend-neutral options owned for the lifetime of a block-reader session.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockReadOptions {
    pub matrix_kind: MatrixKind,
    pub sparse: bool,
    pub requested_samples: Option<Vec<String>>,
    pub variant_filter: Option<VariantFilter>,
    pub dosage_source: DosageSource,
    pub missing_policy: DenseMissingPolicy,
    pub return_samples: bool,
    pub return_variants: bool,
}

/// One owned dense or sparse block returned by a persistent reader.
#[derive(Debug, Clone, PartialEq)]
pub enum BlockOutput {
    /// Dense matrix block with optional metadata and diagnostics.
    Dense(DenseGenotypeMatrix),
    /// Sparse CSC matrix block with optional metadata and diagnostics.
    Sparse(SparseGenotypeMatrix),
}

enum BlockBackend {
    Bgen(BgenBlockSession),
    Plink1(Plink1BlockSession),
    Plink2(Plink2BlockSession),
    TextVcf(TextVcfBlockSession),
}

/// Backend-neutral persistent reader that yields bounded genotype blocks.
pub struct BlockReader {
    backend: BlockBackend,
    lifecycle: BlockLifecycle,
}

impl fmt::Debug for BlockReader {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let backend = match &self.backend {
            BlockBackend::Bgen(_) => "bgen",
            BlockBackend::Plink1(_) => "plink1",
            BlockBackend::Plink2(_) => "plink2",
            BlockBackend::TextVcf(_) => "text-vcf",
        };
        let lifecycle = if self.lifecycle.eof { "eof" } else { "active" };
        formatter
            .debug_struct("BlockReader")
            .field("backend", &backend)
            .field("lifecycle", &lifecycle)
            .field("block_size", &self.lifecycle.block_size)
            .finish()
    }
}

impl BlockReader {
    /// Open one persistent source session.
    pub fn open(source: BlockSource, options: BlockReadOptions, block_size: usize) -> Result<Self> {
        let lifecycle = BlockLifecycle::new(block_size)?;
        let backend = match source {
            BlockSource::Bgen { bgen, sample } => {
                BlockBackend::Bgen(BgenBlockSession::open(bgen, sample, options)?)
            }
            BlockSource::Plink1 { bed, bim, fam } => {
                BlockBackend::Plink1(Plink1BlockSession::open(bed, bim, fam, options)?)
            }
            BlockSource::Plink2 { pgen, pvar, psam } => {
                BlockBackend::Plink2(Plink2BlockSession::open(pgen, pvar, psam, options)?)
            }
            BlockSource::Vcf { vcf } => {
                BlockBackend::TextVcf(TextVcfBlockSession::open(vcf, options)?)
            }
            BlockSource::Bcf { .. } => {
                return Err(GenoioError::unsupported(
                    "persistent block reads are not implemented for this source yet",
                ));
            }
        };
        Ok(Self { backend, lifecycle })
    }

    /// Decode the next retained block, or return `None` after terminal EOF.
    pub fn next_block(&mut self) -> Result<Option<BlockOutput>> {
        let Self { backend, lifecycle } = self;
        lifecycle.next_block(|block_size| match backend {
            BlockBackend::Bgen(session) => session.next_block(block_size),
            BlockBackend::Plink1(session) => session.next_block(block_size),
            BlockBackend::Plink2(session) => session.next_block(block_size),
            BlockBackend::TextVcf(session) => session.next_block(block_size),
        })
    }
}

impl BlockOutput {
    fn width(&self) -> usize {
        match self {
            Self::Dense(matrix) => matrix.n_variants,
            Self::Sparse(matrix) => matrix.n_cols,
        }
    }

    #[expect(
        dead_code,
        reason = "used by concrete block sessions introduced in later phases"
    )]
    pub(crate) fn diagnostics_mut(&mut self) -> &mut DenseDiagnostics {
        match self {
            Self::Dense(matrix) => &mut matrix.diagnostics,
            Self::Sparse(matrix) => &mut matrix.diagnostics,
        }
    }
}

fn validate_block_output(output: Option<&BlockOutput>, block_size: usize) -> Result<()> {
    let Some(output) = output else {
        return Ok(());
    };
    let width = output.width();
    if width == 0 || width > block_size {
        return Err(GenoioError::internal_contract(format!(
            "block output width {width} must be between 1 and block size {block_size}",
        )));
    }
    Ok(())
}

#[derive(Debug)]
pub(crate) struct BlockLifecycle {
    block_size: usize,
    eof: bool,
}

impl BlockLifecycle {
    pub(crate) fn new(block_size: usize) -> Result<Self> {
        if block_size == 0 {
            return Err(GenoioError::internal_contract(
                "block size must be greater than zero",
            ));
        }
        Ok(Self {
            block_size,
            eof: false,
        })
    }

    pub(crate) fn next_block<F>(&mut self, read_next: F) -> Result<Option<BlockOutput>>
    where
        F: FnOnce(usize) -> Result<Option<BlockOutput>>,
    {
        if self.eof {
            return Ok(None);
        }

        let output = read_next(self.block_size)?;
        validate_block_output(output.as_ref(), self.block_size)?;
        if output.is_none() {
            self.eof = true;
        }
        Ok(output)
    }
}

pub(crate) fn checked_dense_block_len(n_rows: usize, block_size: usize) -> Result<usize> {
    n_rows.checked_mul(block_size).ok_or_else(|| {
        GenoioError::internal_contract(format!(
            "dense block shape {n_rows} x {block_size} is out of range",
        ))
    })
}

pub(crate) fn checked_sparse_indptr_len(block_size: usize) -> Result<usize> {
    block_size
        .checked_add(1)
        .ok_or_else(|| GenoioError::internal_contract("sparse block column count is out of range"))
}

pub(crate) fn block_diagnostics_snapshot(
    cumulative: &DenseDiagnostics,
    block_width: usize,
) -> DenseDiagnostics {
    DenseDiagnostics {
        retained_variants: block_width,
        ..cumulative.clone()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;

    use genoio_core::{
        DenseDiagnostics, DenseGenotypeMatrix, DenseLayout, GenoioError, SampleMetadataBuffers,
        SampleRecord, SparseGenotypeMatrix, VariantMetadataBuffers,
    };

    use super::*;

    fn sample_metadata() -> SampleMetadataBuffers {
        SampleMetadataBuffers::from_records(
            &[
                SampleRecord {
                    fid: None,
                    iid: "sample-1".to_owned(),
                    father: None,
                    mother: None,
                    sex: None,
                    phenotype: None,
                    source_sample_index: None,
                    haplotype_index: None,
                },
                SampleRecord {
                    fid: None,
                    iid: "sample-2".to_owned(),
                    father: None,
                    mother: None,
                    sex: None,
                    phenotype: None,
                    source_sample_index: None,
                    haplotype_index: None,
                },
            ],
            false,
        )
        .expect("sample metadata should be valid")
    }

    fn variant_metadata(width: usize) -> VariantMetadataBuffers {
        let mut variants = VariantMetadataBuffers::with_capacity(width);
        for index in 0..width {
            variants
                .push(
                    "1",
                    i64::try_from(index + 1).expect("test index should fit in i64"),
                    &format!("variant-{}", index + 1),
                    "A",
                    "C",
                )
                .expect("variant metadata should be valid");
        }
        variants
    }

    fn diagnostics(width: usize) -> DenseDiagnostics {
        DenseDiagnostics {
            requested_samples: 2,
            retained_samples: 2,
            candidate_variants: width + 3,
            retained_variants: width,
            dropped_metadata_variants: 1,
            dropped_genotype_variants: 2,
            ..DenseDiagnostics::default()
        }
    }

    fn dense_output(width: usize) -> BlockOutput {
        BlockOutput::Dense(
            DenseGenotypeMatrix::new_with_layout(
                2,
                width,
                vec![0.0; 2 * width],
                DenseLayout::SampleMajor,
                Some(sample_metadata()),
                Some(variant_metadata(width)),
                diagnostics(width),
            )
            .expect("dense matrix should be valid"),
        )
    }

    fn sparse_output(width: usize) -> BlockOutput {
        BlockOutput::Sparse(
            SparseGenotypeMatrix::new(
                2,
                width,
                vec![0; width + 1],
                Vec::new(),
                Vec::new(),
                Some(sample_metadata()),
                Some(variant_metadata(width)),
                diagnostics(width),
            )
            .expect("sparse matrix should be valid"),
        )
    }

    #[test]
    fn pbr_rust_block_001_block_output_width_uses_established_matrix_dimensions() {
        assert_eq!(dense_output(3).width(), 3);
        assert_eq!(sparse_output(4).width(), 4);
    }

    #[test]
    fn pbr_rust_block_001_block_output_validation_accepts_partial_and_exact_widths() {
        validate_block_output(Some(&dense_output(2)), 4)
            .expect("partial dense block should be valid");
        validate_block_output(Some(&sparse_output(4)), 4)
            .expect("exact-width sparse block should be valid");
        validate_block_output(None, 4).expect("EOF should be valid");
    }

    #[test]
    fn pbr_rust_block_001_block_output_validation_rejects_zero_and_over_widths() {
        let zero_error = validate_block_output(Some(&dense_output(0)), 4)
            .expect_err("zero-width output should fail");
        let over_error = validate_block_output(Some(&sparse_output(5)), 4)
            .expect_err("over-width output should fail");

        assert!(matches!(zero_error, GenoioError::InternalContract { .. }));
        assert!(matches!(over_error, GenoioError::InternalContract { .. }));
    }

    #[test]
    fn pbr_rust_block_001_dense_output_retains_owned_matrix_values_metadata_and_diagnostics() {
        let output = dense_output(2);

        let BlockOutput::Dense(matrix) = output else {
            panic!("expected dense output");
        };
        assert_eq!(matrix.values, vec![0.0; 4]);
        assert_eq!(matrix.samples, Some(sample_metadata()));
        assert_eq!(matrix.variants, Some(variant_metadata(2)));
        assert_eq!(matrix.diagnostics, diagnostics(2));
    }

    #[test]
    fn pbr_rust_block_001_sparse_output_retains_owned_matrix_values_metadata_and_diagnostics() {
        let output = sparse_output(2);

        let BlockOutput::Sparse(matrix) = output else {
            panic!("expected sparse output");
        };
        assert_eq!(matrix.indptr, vec![0, 0, 0]);
        assert!(matrix.indices.is_empty());
        assert!(matrix.data.is_empty());
        assert_eq!(matrix.samples, Some(sample_metadata()));
        assert_eq!(matrix.variants, Some(variant_metadata(2)));
        assert_eq!(matrix.diagnostics, diagnostics(2));
    }

    #[test]
    fn pbr_rust_alloc_001_checked_block_capacities_follow_requested_dimensions() {
        assert_eq!(
            checked_dense_block_len(3, 4).expect("dense dimensions should fit"),
            12
        );
        assert_eq!(
            checked_sparse_indptr_len(4).expect("sparse dimensions should fit"),
            5
        );
    }

    #[test]
    fn pbr_rust_alloc_001_checked_block_capacities_reject_arithmetic_overflow() {
        let dense_error = checked_dense_block_len(usize::MAX, 2)
            .expect_err("dense output length overflow should fail");
        let sparse_error = checked_sparse_indptr_len(usize::MAX)
            .expect_err("sparse indptr length overflow should fail");

        assert!(matches!(dense_error, GenoioError::InternalContract { .. }));
        assert!(matches!(sparse_error, GenoioError::InternalContract { .. }));
    }

    #[test]
    fn pbr_rust_diag_001_diagnostics_snapshot_preserves_cumulative_counts_and_uses_block_width() {
        let cumulative = DenseDiagnostics {
            requested_samples: 11,
            retained_samples: 7,
            missing_samples: 4,
            candidate_variants: 23,
            retained_variants: 19,
            dropped_metadata_variants: 2,
            dropped_genotype_variants: 2,
        };

        let snapshot = block_diagnostics_snapshot(&cumulative, 3);

        assert_eq!(
            snapshot,
            DenseDiagnostics {
                retained_variants: 3,
                ..cumulative.clone()
            }
        );
        assert_eq!(cumulative.retained_variants, 19);
    }

    #[test]
    fn pbr_rust_block_001_block_lifecycle_rejects_zero_block_size() {
        let error = BlockLifecycle::new(0).expect_err("zero block size should fail");

        assert!(matches!(error, GenoioError::InternalContract { .. }));
    }

    #[test]
    fn pbr_rust_eof_001_block_lifecycle_makes_immediate_eof_sticky() {
        let mut lifecycle = BlockLifecycle::new(2).expect("positive block size should be valid");
        let mut backend_invocations = 0;

        let first = lifecycle
            .next_block(|_| {
                backend_invocations += 1;
                Ok(None)
            })
            .expect("immediate EOF should be valid");
        let second = lifecycle
            .next_block(|_| {
                backend_invocations += 1;
                Ok(Some(dense_output(1)))
            })
            .expect("sticky EOF should be valid");

        assert!(first.is_none());
        assert!(second.is_none());
        assert_eq!(backend_invocations, 1);
    }

    #[test]
    fn pbr_rust_block_001_block_lifecycle_accepts_exact_then_partial_final_output() {
        let mut lifecycle = BlockLifecycle::new(3).expect("positive block size should be valid");
        let mut backend_invocations = 0;
        let mut script = VecDeque::from([
            Ok(Some(dense_output(3))),
            Ok(Some(sparse_output(2))),
            Ok(None),
        ]);

        let mut next = || {
            lifecycle.next_block(|block_size| {
                backend_invocations += 1;
                assert_eq!(block_size, 3);
                script.pop_front().expect("script should contain a result")
            })
        };

        assert_eq!(
            next()
                .expect("exact-width output should be valid")
                .map(|output| output.width()),
            Some(3)
        );
        assert_eq!(
            next()
                .expect("partial-final output should be valid")
                .map(|output| output.width()),
            Some(2)
        );
        assert!(next().expect("EOF should be valid").is_none());
        assert_eq!(backend_invocations, 3);
    }

    #[test]
    fn pbr_rust_block_001_block_lifecycle_rejects_invalid_nonterminal_widths() {
        let mut lifecycle = BlockLifecycle::new(2).expect("positive block size should be valid");
        let mut backend_invocations = 0;

        let zero_error = lifecycle
            .next_block(|_| {
                backend_invocations += 1;
                Ok(Some(dense_output(0)))
            })
            .expect_err("zero-width output should fail");
        let over_error = lifecycle
            .next_block(|_| {
                backend_invocations += 1;
                Ok(Some(sparse_output(3)))
            })
            .expect_err("over-width output should fail");

        assert!(matches!(zero_error, GenoioError::InternalContract { .. }));
        assert!(matches!(over_error, GenoioError::InternalContract { .. }));
        assert_eq!(backend_invocations, 2);
    }

    #[test]
    fn pbr_rust_eof_001_block_lifecycle_error_does_not_prefetch_or_consume_later_result() {
        let mut lifecycle = BlockLifecycle::new(2).expect("positive block size should be valid");
        let mut backend_invocations = 0;
        let mut script = VecDeque::from([
            Err(GenoioError::invalid_source(
                "<script>",
                "scripted backend failure",
            )),
            Ok(Some(dense_output(1))),
            Ok(None),
        ]);

        let error = lifecycle
            .next_block(|_| {
                backend_invocations += 1;
                script.pop_front().expect("script should contain a result")
            })
            .expect_err("backend error should be returned");

        assert!(matches!(error, GenoioError::InvalidSource { .. }));
        assert_eq!(backend_invocations, 1);
        assert_eq!(script.len(), 2);

        let output = lifecycle
            .next_block(|_| {
                backend_invocations += 1;
                script.pop_front().expect("script should contain a result")
            })
            .expect("later scripted result should remain available")
            .expect("later scripted result should be a block");

        assert_eq!(output.width(), 1);
        assert_eq!(backend_invocations, 2);
        assert_eq!(script.len(), 1);
    }
}
