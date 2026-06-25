// pattern: Imperative Shell
#![allow(dead_code, unused_imports)]

use std::path::Path;

use genoio_core::{DenseGenotypeMatrix, DenseMissingPolicy, VariantFilter, VariantWindow};

pub(crate) use super::output::{
    dense_missing_sample_major_output, string_at, variant_a0, variant_a1, variant_alt_allele,
    variant_chrom, variant_id, variant_ids, variant_ref_allele, variants,
};
pub(crate) use ::genoio_io::{read_bgen_metadata, Result};

pub(crate) fn read_bgen_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_bgen_dosage_dense_windowed(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}

pub(crate) fn read_bgen_haplotypes_dosage_dense_windowed(
    bgen: &Path,
    sample: Option<&Path>,
    requested_samples: Option<&[String]>,
    variant_filter: Option<&VariantFilter>,
    variant_window: Option<VariantWindow>,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    ::genoio_io::read_bgen_haplotypes_dosage_dense_windowed(
        bgen,
        sample,
        requested_samples,
        variant_filter,
        variant_window,
        DenseMissingPolicy::Nan,
        !matrix_only,
        !matrix_only,
    )
}
