// pattern: Imperative Shell
//! PLINK2 read setup and empty-output helpers.
//!
//! Read contexts pair a PGEN header with PSAM sample selection before the dense,
//! dosage, haplotype, or sparse loops run. Empty-output helpers centralize the
//! matrix-only metadata elision rules.

use std::fs;
use std::path::Path;

use genoio_core::{
    select_samples_source_order, DenseDiagnostics, DenseGenotypeMatrix, DenseSampleSelection,
    GenoioError, SampleRecord, SparseGenotypeMatrix, VariantWindow,
};

use crate::error::Result;
use crate::matrix::{finish_dense_matrix, DenseMatrixParts};

use super::metadata::parse_psam;
use super::pgen::{
    read_supported_pgen_header, read_supported_pgen_header_prefix, validate_plink2_sample_count,
    PgenHeader,
};

pub(super) struct Plink2ReadContext {
    pub(super) header: PgenHeader,
    pub(super) selection: DenseSampleSelection,
    pub(super) all_samples_selected: bool,
}

impl Plink2ReadContext {
    pub(super) fn new(
        pgen: &Path,
        psam: &Path,
        requested_samples: Option<&[String]>,
    ) -> Result<Self> {
        let header = read_supported_pgen_header(pgen)?;
        Self::from_header(pgen, psam, requested_samples, header)
    }

    pub(super) fn new_prefix(
        pgen: &Path,
        psam: &Path,
        requested_samples: Option<&[String]>,
        decode_variant_ct: usize,
    ) -> Result<Self> {
        let header = read_supported_pgen_header_prefix(pgen, decode_variant_ct)?;
        Self::from_header(pgen, psam, requested_samples, header)
    }

    fn from_header(
        pgen: &Path,
        psam: &Path,
        requested_samples: Option<&[String]>,
        header: PgenHeader,
    ) -> Result<Self> {
        let selection = select_samples_for_header(pgen, psam, requested_samples, &header)?;
        Ok(Self {
            header,
            selection,
            all_samples_selected: requested_samples.is_none(),
        })
    }
}

pub(super) fn select_samples_for_header(
    pgen: &Path,
    psam: &Path,
    requested_samples: Option<&[String]>,
    header: &PgenHeader,
) -> Result<DenseSampleSelection> {
    let all_samples = parse_psam(psam)?;
    validate_plink2_sample_count(pgen, header, all_samples.len())?;
    select_samples_source_order(&all_samples, requested_samples, pgen)
}

pub(super) fn require_pvar(path: &Path) -> Result<()> {
    fs::metadata(path).map_err(|source| GenoioError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

pub(super) fn variant_output_capacity(
    header: &PgenHeader,
    variant_window: Option<VariantWindow>,
) -> usize {
    variant_window.map_or(header.variant_ct, |window| {
        window.len.min(header.variant_ct)
    })
}

pub(super) fn empty_dense_for_samples(
    samples: Vec<SampleRecord>,
    mut diagnostics: DenseDiagnostics,
    matrix_only: bool,
) -> Result<DenseGenotypeMatrix> {
    diagnostics.retained_variants = 0;
    finish_dense_matrix(
        DenseMatrixParts {
            n_samples: samples.len(),
            n_variants: 0,
            values: Vec::new(),
            samples,
            variants: Vec::new(),
            diagnostics,
        },
        matrix_only,
    )
}

pub(super) fn empty_sparse_for_selection(
    selection: DenseSampleSelection,
) -> Result<SparseGenotypeMatrix> {
    empty_sparse_for_samples(selection.samples, selection.diagnostics)
}

pub(super) fn empty_sparse_for_samples(
    samples: Vec<SampleRecord>,
    mut diagnostics: DenseDiagnostics,
) -> Result<SparseGenotypeMatrix> {
    diagnostics.retained_variants = 0;
    SparseGenotypeMatrix::new(
        samples.len(),
        0,
        vec![0],
        Vec::new(),
        Vec::new(),
        samples,
        Vec::new(),
        diagnostics,
    )
}

pub(super) fn expand_selected_samples_to_haplotypes(
    selection: &DenseSampleSelection,
) -> Vec<SampleRecord> {
    let mut haplotype_samples = Vec::with_capacity(selection.samples.len() * 2);
    for (sample, &source_index) in selection.samples.iter().zip(&selection.source_indices) {
        for haplotype_index in 0..2 {
            let mut haplotype_sample = sample.clone();
            haplotype_sample.source_sample_index = Some(source_index);
            haplotype_sample.haplotype_index = Some(haplotype_index);
            haplotype_samples.push(haplotype_sample);
        }
    }
    haplotype_samples
}

pub(super) fn matrix_only_source_window_diagnostics(
    n_samples: usize,
    n_variants: usize,
) -> DenseDiagnostics {
    DenseDiagnostics {
        requested_samples: n_samples,
        retained_samples: n_samples,
        candidate_variants: n_variants,
        retained_variants: n_variants,
        ..DenseDiagnostics::default()
    }
}
