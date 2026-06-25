// pattern: Imperative Shell
#![allow(dead_code, unused_imports)]

use std::path::Path;

use genoio_core::{
    DenseGenotypeMatrix, DenseMissingPolicy, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};

pub(crate) use super::output::{
    dense_missing_sample_major_output, dense_values_sample_major_output,
    sparse_values_dense_output, string_at, strings, variant_alt_allele, variant_id, variant_ids,
    variant_ref_allele, variants,
};
pub(crate) use ::genoio_io::{read_plink1_metadata, read_plink2_metadata, Result};

pub(crate) fn read_plink1_dense(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_plink1_dense_windowed(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        None,
        false,
    )
}

pub(crate) fn read_plink1_dense_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_plink1_dense_windowed(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_plink1_sparse(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_plink1_sparse_windowed(bed, bim, fam, requested_samples, variant_filter, None)
}

pub(crate) fn read_plink1_sparse_windowed(
    bed: &Path,
    bim: &Path,
    fam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    ::genoio_io::read_plink1_sparse_windowed(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )
}

pub(crate) fn read_plink2_dense(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<DenseGenotypeMatrix> {
    read_plink2_dense_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        None,
        false,
    )
}

pub(crate) fn read_plink2_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_plink2_dense_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_plink2_dosage_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_plink2_dosage_dense_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_plink2_sparse(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
) -> Result<SparseGenotypeMatrix> {
    read_plink2_sparse_windowed(pgen, pvar, psam, requested_samples, variant_filter, None)
}

pub(crate) fn read_plink2_sparse_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    ::genoio_io::read_plink2_sparse_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )
}

pub(crate) fn read_plink2_haplotypes_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_plink2_haplotypes_dense_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_plink2_haplotypes_dosage_dense_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_plink2_haplotypes_dosage_dense_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_plink2_haplotypes_sparse_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    ::genoio_io::read_plink2_haplotypes_sparse_windowed(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )
}
