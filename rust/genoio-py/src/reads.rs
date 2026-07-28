// pattern: Imperative Shell

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::errors::{catch_internal_panic, genoio_error_to_py, run_without_gil};
use crate::options::{read_options, DosageSource};
use crate::output::{dense_matrix_to_py, metadata_to_py, sparse_matrix_to_py};
use crate::source::{
    read_dense_matrix_for_py, read_source_metadata, read_sparse_matrix_for_py, source_members,
    validate_read_support_impl, MatrixKind, SourceFormat,
};

#[pyfunction]
/// Return the Rust IO backend name for Python diagnostics.
pub(crate) fn backend_name() -> &'static str {
    genoio_io::backend_name()
}

#[pyfunction]
/// Validate whether a source/read representation is supported by the Rust backend.
pub(crate) fn validate_read_support(
    format: &str,
    kind: &str,
    dosage: &str,
    sparse: bool,
) -> PyResult<()> {
    catch_internal_panic(|| {
        let format = SourceFormat::from_str(format)?;
        let kind = MatrixKind::from_str(kind).map_err(genoio_error_to_py)?;
        let dosage = DosageSource::from_str(dosage).map_err(genoio_error_to_py)?;
        validate_read_support_impl(format, kind, dosage, sparse).map_err(genoio_error_to_py)
    })
}

#[pyfunction]
/// Read source metadata and convert it into Python-owned Arrow stream wrappers.
pub(crate) fn read_metadata(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    catch_internal_panic(|| read_metadata_impl(py, format, members))
}

fn read_metadata_impl(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let source = source_members(format, members)?;
    let output = run_without_gil(py, || read_source_metadata(&source))?;

    metadata_to_py(py, output)
}

#[pyfunction]
/// Read dense genotypes from a resolved source.
pub(crate) fn read_dense(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    catch_internal_panic(|| read_dense_impl(py, format, members, options))
}

fn read_dense_impl(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let read_options = read_options(options)?;
    let source = source_members(format, members)?;
    let output = run_without_gil(py, || {
        read_dense_matrix_for_py(&source, MatrixKind::Genotype, &read_options)
    })?;

    dense_matrix_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
    )
}

#[pyfunction]
/// Read sparse genotypes from a resolved source.
pub(crate) fn read_sparse(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    catch_internal_panic(|| read_sparse_impl(py, format, members, options))
}

fn read_sparse_impl(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let read_options = read_options(options)?;
    let source = source_members(format, members)?;
    let output = run_without_gil(py, || {
        read_sparse_matrix_for_py(&source, MatrixKind::Genotype, &read_options)
    })?;

    sparse_matrix_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
    )
}

#[pyfunction]
/// Read dense haplotype rows from a resolved source.
pub(crate) fn read_haplotypes_dense(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    catch_internal_panic(|| read_haplotypes_dense_impl(py, format, members, options))
}

fn read_haplotypes_dense_impl(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let read_options = read_options(options)?;
    let source = source_members(format, members)?;
    let output = run_without_gil(py, || {
        read_dense_matrix_for_py(&source, MatrixKind::Haplotype, &read_options)
    })?;

    dense_matrix_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
    )
}

#[pyfunction]
/// Read sparse haplotype rows from a resolved source.
pub(crate) fn read_haplotypes_sparse(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    catch_internal_panic(|| read_haplotypes_sparse_impl(py, format, members, options))
}

fn read_haplotypes_sparse_impl(
    py: Python<'_>,
    format: &str,
    members: &Bound<'_, PyDict>,
    options: &Bound<'_, PyDict>,
) -> PyResult<Py<PyDict>> {
    let read_options = read_options(options)?;
    let source = source_members(format, members)?;
    let output = run_without_gil(py, || {
        read_sparse_matrix_for_py(&source, MatrixKind::Haplotype, &read_options)
    })?;

    sparse_matrix_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
    )
}
