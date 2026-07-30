// pattern: Imperative Shell

use std::path::PathBuf;

use genoio_core::{DenseGenotypeMatrix, GenoioError, MetadataOutput, SparseGenotypeMatrix};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::options::{BlockOptions, DosageSource, ReadOptions};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MatrixKind {
    Genotype,
    Haplotype,
}

impl MatrixKind {
    pub(crate) fn from_str(value: &str) -> Result<Self, GenoioError> {
        match value {
            "geno" => Ok(Self::Genotype),
            "haplo" => Ok(Self::Haplotype),
            other => Err(GenoioError::unsupported(format!(
                "unsupported genotype kind: {other}"
            ))),
        }
    }
}

impl From<MatrixKind> for genoio_io::MatrixKind {
    fn from(value: MatrixKind) -> Self {
        match value {
            MatrixKind::Genotype => Self::Genotype,
            MatrixKind::Haplotype => Self::Haplotype,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SourceFormat {
    Vcf,
    Bcf,
    Plink1,
    Plink2,
    Bgen,
}

impl SourceFormat {
    pub(crate) fn from_str(format: &str) -> PyResult<Self> {
        match format {
            "vcf" => Ok(Self::Vcf),
            "bcf" => Ok(Self::Bcf),
            "plink1" => Ok(Self::Plink1),
            "plink2" => Ok(Self::Plink2),
            "bgen" => Ok(Self::Bgen),
            other => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unsupported source format: {other}"
            ))),
        }
    }
}

impl From<SourceFormat> for genoio_io::BlockFormat {
    fn from(value: SourceFormat) -> Self {
        match value {
            SourceFormat::Vcf => Self::Vcf,
            SourceFormat::Bcf => Self::Bcf,
            SourceFormat::Plink1 => Self::Plink1,
            SourceFormat::Plink2 => Self::Plink2,
            SourceFormat::Bgen => Self::Bgen,
        }
    }
}

pub(crate) enum SourceMembers {
    Vcf {
        path: PathBuf,
    },
    Bcf {
        path: PathBuf,
    },
    Plink1 {
        bed: PathBuf,
        bim: PathBuf,
        fam: PathBuf,
    },
    Plink2 {
        pgen: PathBuf,
        pvar: PathBuf,
        psam: PathBuf,
    },
    Bgen {
        bgen: PathBuf,
        sample: Option<PathBuf>,
    },
}

impl SourceMembers {
    fn format(&self) -> SourceFormat {
        match self {
            Self::Vcf { .. } => SourceFormat::Vcf,
            Self::Bcf { .. } => SourceFormat::Bcf,
            Self::Plink1 { .. } => SourceFormat::Plink1,
            Self::Plink2 { .. } => SourceFormat::Plink2,
            Self::Bgen { .. } => SourceFormat::Bgen,
        }
    }

    pub(crate) fn into_block_source(self) -> genoio_io::BlockSource {
        match self {
            Self::Vcf { path } => genoio_io::BlockSource::Vcf { vcf: path },
            Self::Bcf { path } => genoio_io::BlockSource::Bcf { bcf: path },
            Self::Plink1 { bed, bim, fam } => genoio_io::BlockSource::Plink1 { bed, bim, fam },
            Self::Plink2 { pgen, pvar, psam } => {
                genoio_io::BlockSource::Plink2 { pgen, pvar, psam }
            }
            Self::Bgen { bgen, sample } => genoio_io::BlockSource::Bgen { bgen, sample },
        }
    }
}

pub(crate) fn source_members(format: &str, members: &Bound<'_, PyDict>) -> PyResult<SourceMembers> {
    let source_format = SourceFormat::from_str(format)?;
    match source_format {
        SourceFormat::Vcf => Ok(SourceMembers::Vcf {
            path: member_path(members, "vcf")?,
        }),
        SourceFormat::Bcf => Ok(SourceMembers::Bcf {
            path: member_path(members, "bcf")?,
        }),
        SourceFormat::Plink1 => Ok(SourceMembers::Plink1 {
            bed: member_path(members, "bed")?,
            bim: member_path(members, "bim")?,
            fam: member_path(members, "fam")?,
        }),
        SourceFormat::Plink2 => Ok(SourceMembers::Plink2 {
            pgen: member_path(members, "pgen")?,
            pvar: member_path(members, "pvar")?,
            psam: member_path(members, "psam")?,
        }),
        SourceFormat::Bgen => Ok(SourceMembers::Bgen {
            bgen: member_path(members, "bgen")?,
            sample: optional_member_path(members, "sample")?,
        }),
    }
}

pub(crate) fn validate_read_support_impl(
    format: SourceFormat,
    kind: MatrixKind,
    dosage: DosageSource,
    sparse: bool,
) -> Result<(), GenoioError> {
    genoio_io::validate_block_support(format.into(), kind.into(), dosage.into(), sparse)
}

pub(crate) fn block_read_options(
    kind: MatrixKind,
    sparse: bool,
    options: BlockOptions,
) -> genoio_io::BlockReadOptions {
    genoio_io::BlockReadOptions {
        matrix_kind: kind.into(),
        sparse,
        requested_samples: options.requested_samples,
        variant_filter: options.variant_filter,
        dosage_source: options.dosage.into(),
        missing_policy: options.missing,
        return_samples: options.return_samples,
        return_variants: options.return_variants,
    }
}

pub(crate) fn read_source_metadata(source: &SourceMembers) -> Result<MetadataOutput, GenoioError> {
    match source {
        SourceMembers::Vcf { path } | SourceMembers::Bcf { path } => {
            genoio_io::read_vcf_public_metadata(path)
        }
        SourceMembers::Plink1 { bed, bim, fam } => genoio_io::read_plink1_metadata(bed, bim, fam),
        SourceMembers::Plink2 { pgen, pvar, psam } => {
            genoio_io::read_plink2_metadata(pgen, pvar, psam)
        }
        SourceMembers::Bgen { bgen, sample } => {
            genoio_io::read_bgen_metadata(bgen, sample.as_deref())
        }
    }
}

pub(crate) fn read_dense_matrix_for_py(
    source: &SourceMembers,
    kind: MatrixKind,
    options: &ReadOptions,
) -> Result<DenseGenotypeMatrix, GenoioError> {
    validate_read_support_impl(source.format(), kind, options.dosage, false)?;
    match (source, kind, options.dosage) {
        (
            SourceMembers::Vcf { path } | SourceMembers::Bcf { path },
            MatrixKind::Genotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_vcf_dense_windowed(
            path,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.return_samples,
            options.return_variants,
        ),
        (
            SourceMembers::Vcf { path } | SourceMembers::Bcf { path },
            MatrixKind::Genotype,
            DosageSource::Dosage,
        ) => genoio_io::read_vcf_dosage_dense_windowed(
            path,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.return_samples,
            options.return_variants,
        ),
        (
            SourceMembers::Vcf { path } | SourceMembers::Bcf { path },
            MatrixKind::Haplotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_vcf_haplotypes_dense_windowed(
            path,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.return_samples,
            options.return_variants,
        ),
        (SourceMembers::Plink1 { bed, bim, fam }, MatrixKind::Genotype, DosageSource::Hardcall) => {
            genoio_io::read_plink1_dense_windowed(
                bed,
                bim,
                fam,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.return_samples,
                options.return_variants,
            )
        }
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Genotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_plink2_dense_windowed(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.return_samples,
            options.return_variants,
        ),
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Genotype,
            DosageSource::Dosage,
        ) => genoio_io::read_plink2_dosage_dense_windowed(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.return_samples,
            options.return_variants,
        ),
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Haplotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_plink2_haplotypes_dense_windowed(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.return_samples,
            options.return_variants,
        ),
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Haplotype,
            DosageSource::Dosage,
        ) => genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.return_samples,
            options.return_variants,
        ),
        (SourceMembers::Bgen { bgen, sample }, MatrixKind::Genotype, DosageSource::Dosage) => {
            genoio_io::read_bgen_dosage_dense_windowed(
                bgen,
                sample.as_deref(),
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.return_samples,
                options.return_variants,
            )
        }
        (SourceMembers::Bgen { bgen, sample }, MatrixKind::Haplotype, DosageSource::Dosage) => {
            genoio_io::read_bgen_haplotypes_dosage_dense_windowed(
                bgen,
                sample.as_deref(),
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.return_samples,
                options.return_variants,
            )
        }
        _ => Err(GenoioError::internal_contract(
            "read support validation accepted unsupported dense dispatch",
        )),
    }
}

pub(crate) fn read_sparse_matrix_for_py(
    source: &SourceMembers,
    kind: MatrixKind,
    options: &ReadOptions,
) -> Result<SparseGenotypeMatrix, GenoioError> {
    validate_read_support_impl(source.format(), kind, options.dosage, true)?;
    match (source, kind, options.dosage) {
        (
            SourceMembers::Vcf { path } | SourceMembers::Bcf { path },
            MatrixKind::Genotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_vcf_sparse_windowed(
            path,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.return_samples,
            options.return_variants,
        ),
        (
            SourceMembers::Vcf { path } | SourceMembers::Bcf { path },
            MatrixKind::Haplotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_vcf_haplotypes_sparse_windowed(
            path,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.return_samples,
            options.return_variants,
        ),
        (SourceMembers::Plink1 { bed, bim, fam }, MatrixKind::Genotype, DosageSource::Hardcall) => {
            genoio_io::read_plink1_sparse_windowed(
                bed,
                bim,
                fam,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.return_samples,
                options.return_variants,
            )
        }
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Genotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_plink2_sparse_windowed(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.return_samples,
            options.return_variants,
        ),
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Haplotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_plink2_haplotypes_sparse_windowed(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.return_samples,
            options.return_variants,
        ),
        _ => Err(GenoioError::internal_contract(
            "read support validation accepted unsupported sparse dispatch",
        )),
    }
}

fn member_path(members: &Bound<'_, PyDict>, key: &str) -> PyResult<PathBuf> {
    let value = members.get_item(key)?.ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err(format!("missing source member: {key}"))
    })?;
    Ok(PathBuf::from(value.extract::<String>()?))
}

fn optional_member_path(members: &Bound<'_, PyDict>, key: &str) -> PyResult<Option<PathBuf>> {
    let Some(value) = members.get_item(key)? else {
        return Ok(None);
    };
    Ok(Some(PathBuf::from(value.extract::<String>()?)))
}
