// pattern: Imperative Shell

use std::path::PathBuf;

use genoio_core::{
    DenseDiagnostics, DenseGenotypeMatrix, DenseMissingPolicy, GenoioError, SparseGenotypeMatrix,
    VariantFilter,
};

use crate::error::Result;

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

#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "used by the lifecycle state added later in this phase"
    )
)]
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

#[cfg(test)]
mod tests {
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
    fn block_output_width_uses_established_matrix_dimensions() {
        assert_eq!(dense_output(3).width(), 3);
        assert_eq!(sparse_output(4).width(), 4);
    }

    #[test]
    fn block_output_validation_accepts_partial_and_exact_widths() {
        validate_block_output(Some(&dense_output(2)), 4)
            .expect("partial dense block should be valid");
        validate_block_output(Some(&sparse_output(4)), 4)
            .expect("exact-width sparse block should be valid");
        validate_block_output(None, 4).expect("EOF should be valid");
    }

    #[test]
    fn block_output_validation_rejects_zero_and_over_widths() {
        let zero_error = validate_block_output(Some(&dense_output(0)), 4)
            .expect_err("zero-width output should fail");
        let over_error = validate_block_output(Some(&sparse_output(5)), 4)
            .expect_err("over-width output should fail");

        assert!(matches!(zero_error, GenoioError::InternalContract { .. }));
        assert!(matches!(over_error, GenoioError::InternalContract { .. }));
    }

    #[test]
    fn dense_output_retains_owned_matrix_values_metadata_and_diagnostics() {
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
    fn sparse_output_retains_owned_matrix_values_metadata_and_diagnostics() {
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
}
