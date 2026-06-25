// pattern: Imperative Shell
#![allow(dead_code, unused_imports)]

use std::path::Path;

use genoio_core::{
    DenseGenotypeMatrix, DenseMissingPolicy, SparseGenotypeMatrix, VariantFilter, VariantWindow,
};

pub(crate) use ::genoio_io::{
    read_bgen_metadata, read_bgen_metadata_arrow, read_plink1_metadata, read_plink1_metadata_arrow,
    read_plink2_metadata, read_plink2_metadata_arrow, Result,
};

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
    ::genoio_io::read_plink1_dense_windowed_with_arrow_variants(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )?
    .into_matrix()
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
    ::genoio_io::read_plink1_sparse_windowed_with_arrow_variants(
        bed,
        bim,
        fam,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )?
    .into_matrix()
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
    ::genoio_io::read_plink2_dense_windowed_with_arrow_variants(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )?
    .into_matrix()
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
    ::genoio_io::read_plink2_dosage_dense_windowed_with_arrow_variants(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )?
    .into_matrix()
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
    ::genoio_io::read_plink2_sparse_windowed_with_arrow_variants(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )?
    .into_matrix()
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
    ::genoio_io::read_plink2_haplotypes_dense_windowed_with_arrow_variants(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )?
    .into_matrix()
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
    ::genoio_io::read_plink2_haplotypes_dosage_dense_windowed_with_arrow_variants(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )?
    .into_matrix()
}

pub(crate) fn read_plink2_haplotypes_sparse_windowed(
    pgen: &Path,
    pvar: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
) -> Result<SparseGenotypeMatrix> {
    ::genoio_io::read_plink2_haplotypes_sparse_windowed_with_arrow_variants(
        pgen,
        pvar,
        psam,
        requested_samples,
        variant_filter,
        variant_window,
        true,
        true,
    )?
    .into_matrix()
}

pub(crate) fn read_bgen_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_bgen_dosage_dense_windowed_with_arrow_variants(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )?
    .into_matrix()
}

pub(crate) fn read_bgen_haplotypes_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_bgen_haplotypes_dosage_dense_windowed_with_arrow_variants(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )?
    .into_matrix()
}
