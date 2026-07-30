// pattern: Imperative Shell
//! PLINK1 BED/BIM/FAM readers.
//!
//! This module coordinates companion-file validation, sample selection, variant
//! filtering, and dense or sparse output assembly. Binary BED decoding and text
//! metadata parsing live in submodules.

use std::fs;
use std::path::Path;

use genoio_core::{
    DenseGenotypeMatrix, DenseMissingPolicy, GenoioError, MetadataOutput, SampleMetadataBuffers,
    SourceCapabilities, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};

use crate::error::Result;

mod bed;
mod metadata;
mod session;

pub(crate) use session::Plink1BlockSession;

use metadata::{parse_bim_metadata, parse_fam};

/// Read PLINK1 metadata with variant metadata staged as columnar buffers.
pub fn read_plink1_metadata(bed: &Path, bim: &Path, fam: &Path) -> Result<MetadataOutput> {
    fs::metadata(bed).map_err(|source| GenoioError::Io {
        path: bed.to_path_buf(),
        source,
    })?;
    let samples = parse_fam(fam)?;
    let variants = parse_bim_metadata(bim)?;

    Ok(MetadataOutput {
        samples: SampleMetadataBuffers::from_records(&samples, false)?,
        variants,
        capabilities: SourceCapabilities::genotype_only(),
    })
}

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors dense read options plus metadata return choices"
)]
pub fn read_plink1_dense_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    missing_policy: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
) -> Result<DenseGenotypeMatrix> {
    let matrix_only = !return_samples && !return_variants;
    let retained_skip = variant_window.map_or(0, |window| window.start);
    let read_bim_records = !(matrix_only
        && variant_window.is_some()
        && variant_filter.is_none_or(|filter| {
            filter.requires_genotype_stats() && filter.is_genotype_stats_only()
        }));
    let options = crate::blocks::BlockReadOptions {
        matrix_kind: crate::blocks::MatrixKind::Genotype,
        sparse: false,
        requested_samples: requested_samples.map(<[String]>::to_vec),
        variant_filter: variant_filter.cloned(),
        dosage_source: crate::blocks::DosageSource::Hardcall,
        missing_policy,
        return_samples,
        return_variants,
    };
    let mut session = Plink1BlockSession::open_windowed(
        bed.to_path_buf(),
        bim.to_path_buf(),
        fam.to_path_buf(),
        options,
        retained_skip,
        read_bim_records,
    )?;
    let block_size = match variant_window {
        Some(window) => window.len,
        None => session.source_record_capacity(),
    };
    match session.next_dense_block(block_size)? {
        Some(matrix) => Ok(matrix),
        None => session.empty_dense_output(),
    }
}

#[expect(
    clippy::too_many_arguments,
    reason = "output facade mirrors sparse read options plus metadata return choices"
)]
pub fn read_plink1_sparse_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    return_samples: bool,
    return_variants: bool,
) -> Result<SparseGenotypeMatrix> {
    let retained_skip = variant_window.map_or(0, |window| window.start);
    let options = crate::blocks::BlockReadOptions {
        matrix_kind: crate::blocks::MatrixKind::Genotype,
        sparse: true,
        requested_samples: requested_samples.map(<[String]>::to_vec),
        variant_filter: variant_filter.cloned(),
        dosage_source: crate::blocks::DosageSource::Hardcall,
        missing_policy: DenseMissingPolicy::Raise,
        return_samples,
        return_variants,
    };
    let mut session = Plink1BlockSession::open_windowed(
        bed.to_path_buf(),
        bim.to_path_buf(),
        fam.to_path_buf(),
        options,
        retained_skip,
        true,
    )?;
    let block_size = match variant_window {
        Some(window) => window.len,
        None => session.source_record_capacity(),
    };
    match session.next_sparse_block(block_size)? {
        Some(matrix) => Ok(matrix),
        None => session.empty_sparse_output(),
    }
}
