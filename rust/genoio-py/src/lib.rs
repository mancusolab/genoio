// pattern: Imperative Shell

#![cfg_attr(
    not(test),
    deny(clippy::expect_used, clippy::panic, clippy::unwrap_used)
)]

//! PyO3 adapter for the private `genoio._rust` extension module.
//!
//! This crate is the boundary between the public Python API and the Rust
//! backend crates. Python owns user-facing ergonomics, overloads, docstrings,
//! and public exceptions. This adapter validates Python-owned dictionaries,
//! translates them into Rust-owned paths/options/filters, calls `genoio-io`,
//! and converts validated core structs back into Python dictionaries.
//!
//! Long-running reader calls release the GIL only after argument extraction has
//! produced Rust-owned values. No borrowed Python object crosses
//! `Python::allow_threads`; result construction resumes under the GIL.
//!
//! Expected backend failures are represented as `genoio_core::GenoioError` and
//! mapped to private Python exception classes. The pure-Python layer then maps
//! those private classes to public `genoio` exceptions. Rust panics are treated
//! as internal bugs and are contained at each `#[pyfunction]` entry point.
//!
//! Numeric arrays returned from this crate transfer Rust vector ownership to
//! NumPy. Dense missing policies are resolved while the IO layer builds matrix
//! values, so dense reads do not materialize Python missing masks. `usize`
//! indices are checked and converted to NumPy-compatible `int64` values before
//! ownership transfer.

use std::any::Any as PanicPayload;
use std::panic::{self, AssertUnwindSafe};
use std::path::PathBuf;
use std::sync::Arc;

use arrow_array::{ArrayRef, Int64Array, RecordBatch, StringArray, StructArray};
use arrow_schema::{ArrowError, DataType, Field, Schema};
use genoio_core::{
    DenseLayout, DenseMissingPolicy, GenoioError, SampleMetadataColumns, VariantMetadataColumns,
};
use numpy::{Element, PyArray1};
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{
    PyAny, PyBool, PyCapsule, PyDict, PyFloat, PyInt, PyList, PyModule, PyString, PyTuple,
};
use pyo3_arrow::ffi::{to_stream_pycapsule, ArrayIterator};

pyo3::create_exception!(genoio_py, RustInvalidSourceError, PyException);
pyo3::create_exception!(genoio_py, RustUnsupportedRepresentationError, PyException);
pyo3::create_exception!(genoio_py, RustInvalidOptionError, PyException);
pyo3::create_exception!(genoio_py, RustMissingDataError, PyException);
pyo3::create_exception!(genoio_py, RustSampleFilterError, PyException);
pyo3::create_exception!(genoio_py, RustInternalError, PyException);

const SPARSE_DOSAGE_BACKED_GENOTYPE_UNSUPPORTED: &str =
    "sparse dosage-backed genotype reads are intentionally unsupported";
const PLINK2_SPARSE_DOSAGE_BACKED_HAPLOTYPE_UNSUPPORTED: &str =
    "plink2 sparse haplotype reads are intentionally unsupported for dosage-backed sources; use dense haplotype reads with sparse=False";

#[pyclass(module = "genoio._rust", skip_from_py_object)]
#[derive(Clone)]
struct ArrowMetadataFrame {
    batch: RecordBatch,
}

#[pymethods]
impl ArrowMetadataFrame {
    #[pyo3(signature = (requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        requested_schema: Option<Bound<'py, PyCapsule>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let schema = self.batch.schema();
        let fields = schema.fields().clone();
        let metadata = schema.metadata.clone();
        let arrays = vec![Ok(
            Arc::new(StructArray::from(self.batch.clone())) as ArrayRef
        )];
        let array_reader = Box::new(ArrayIterator::new(
            arrays,
            Field::new_struct("", fields, false)
                .with_metadata(metadata)
                .into(),
        ));

        to_stream_pycapsule(py, array_reader, requested_schema).map_err(PyErr::from)
    }
}

#[pyfunction]
/// Return the Rust IO backend name for Python diagnostics.
fn backend_name() -> &'static str {
    genoio_io::backend_name()
}

#[pyfunction]
/// Validate whether a source/read representation is supported by the Rust backend.
fn validate_read_support(format: &str, kind: &str, dosage: &str, sparse: bool) -> PyResult<()> {
    catch_internal_panic(|| {
        let format = SourceFormat::from_str(format)?;
        let kind = MatrixKind::from_str(kind).map_err(genoio_error_to_py)?;
        let dosage = DosageSource::from_str(dosage).map_err(genoio_error_to_py)?;
        validate_read_support_impl(format, kind, dosage, sparse).map_err(genoio_error_to_py)
    })
}

#[pyfunction]
/// Read source metadata and convert it into Python dictionaries.
fn read_metadata(
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
fn read_dense(
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
        read_dense_matrix(&source, MatrixKind::Genotype, &read_options)
    })?;

    dense_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
        false,
    )
}

#[pyfunction]
/// Read sparse genotypes from a resolved source.
fn read_sparse(
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
        read_sparse_matrix(&source, MatrixKind::Genotype, &read_options)
    })?;

    sparse_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
        false,
    )
}

#[pyfunction]
/// Read dense haplotype rows from a resolved source.
fn read_haplotypes_dense(
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
        read_dense_matrix(&source, MatrixKind::Haplotype, &read_options)
    })?;

    dense_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
        true,
    )
}

#[pyfunction]
/// Read sparse haplotype rows from a resolved source.
fn read_haplotypes_sparse(
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
        read_sparse_matrix(&source, MatrixKind::Haplotype, &read_options)
    })?;

    sparse_to_py(
        py,
        output,
        read_options.return_samples,
        read_options.return_variants,
        true,
    )
}

fn catch_internal_panic<T>(f: impl FnOnce() -> PyResult<T>) -> PyResult<T> {
    match panic::catch_unwind(AssertUnwindSafe(f)) {
        // Expected errors are normal Result/PyErr values and pass through
        // unchanged. Only a Rust panic becomes RustInternalError.
        Ok(result) => result,
        Err(payload) => Err(RustInternalError::new_err(format!(
            "internal Rust backend panic: {}",
            panic_message(payload.as_ref())
        ))),
    }
}

fn run_without_gil<T: Send>(
    py: Python<'_>,
    f: impl FnOnce() -> Result<T, GenoioError> + Send,
) -> PyResult<T> {
    // Only Rust-owned source paths, filters, options, and output structs may
    // cross this boundary. Python argument parsing and result construction stay
    // under the GIL.
    py.detach(f).map_err(genoio_error_to_py)
}

fn panic_message(payload: &(dyn PanicPayload + Send)) -> &str {
    if let Some(message) = payload.downcast_ref::<&'static str>() {
        message
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.as_str()
    } else {
        "unknown panic payload"
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MatrixKind {
    Genotype,
    Haplotype,
}

impl MatrixKind {
    fn from_str(value: &str) -> Result<Self, GenoioError> {
        match value {
            "geno" => Ok(Self::Genotype),
            "haplo" => Ok(Self::Haplotype),
            other => Err(GenoioError::unsupported(format!(
                "unsupported genotype kind: {other}"
            ))),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SourceFormat {
    Vcf,
    Bcf,
    Plink1,
    Plink2,
    Bgen,
}

impl SourceFormat {
    fn from_str(format: &str) -> PyResult<Self> {
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

enum SourceMembers {
    Vcf {
        format: SourceFormat,
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
            Self::Vcf { format, .. } => *format,
            Self::Plink1 { .. } => SourceFormat::Plink1,
            Self::Plink2 { .. } => SourceFormat::Plink2,
            Self::Bgen { .. } => SourceFormat::Bgen,
        }
    }
}

fn source_members(format: &str, members: &Bound<'_, PyDict>) -> PyResult<SourceMembers> {
    let source_format = SourceFormat::from_str(format)?;
    match source_format {
        SourceFormat::Vcf | SourceFormat::Bcf => {
            let key = if format == "vcf" { "vcf" } else { "bcf" };
            Ok(SourceMembers::Vcf {
                format: source_format,
                path: member_path(members, key)?,
            })
        }
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

fn validate_read_support_impl(
    format: SourceFormat,
    kind: MatrixKind,
    dosage: DosageSource,
    sparse: bool,
) -> Result<(), GenoioError> {
    match (format, kind, dosage, sparse) {
        (
            SourceFormat::Vcf | SourceFormat::Bcf,
            MatrixKind::Genotype,
            DosageSource::Hardcall,
            _,
        )
        | (
            SourceFormat::Vcf | SourceFormat::Bcf,
            MatrixKind::Haplotype,
            DosageSource::Hardcall,
            _,
        )
        | (
            SourceFormat::Vcf | SourceFormat::Bcf,
            MatrixKind::Genotype,
            DosageSource::Dosage,
            false,
        )
        | (
            SourceFormat::Plink1,
            MatrixKind::Genotype,
            DosageSource::Hardcall,
            _,
        )
        | (
            SourceFormat::Plink2,
            MatrixKind::Genotype,
            DosageSource::Hardcall,
            _,
        )
        | (
            SourceFormat::Plink2,
            MatrixKind::Haplotype,
            DosageSource::Hardcall,
            _,
        )
        | (
            SourceFormat::Plink2,
            MatrixKind::Genotype,
            DosageSource::Dosage,
            false,
        )
        | (
            SourceFormat::Plink2,
            MatrixKind::Haplotype,
            DosageSource::Dosage,
            false,
        )
        | (
            SourceFormat::Bgen,
            MatrixKind::Genotype,
            DosageSource::Dosage,
            false,
        )
        | (
            SourceFormat::Bgen,
            MatrixKind::Haplotype,
            DosageSource::Dosage,
            false,
        ) => Ok(()),
        (
            SourceFormat::Vcf | SourceFormat::Bcf,
            MatrixKind::Haplotype,
            DosageSource::Dosage,
            false,
        ) => Err(GenoioError::unsupported(
            "VCF haplotype dosage reads are unsupported because VCF haplotype support is hardcall GT-based",
        )),
        (
            SourceFormat::Vcf | SourceFormat::Bcf | SourceFormat::Plink1 | SourceFormat::Plink2,
            MatrixKind::Genotype,
            DosageSource::Dosage,
            true,
        ) => Err(GenoioError::unsupported(
            SPARSE_DOSAGE_BACKED_GENOTYPE_UNSUPPORTED,
        )),
        (
            SourceFormat::Vcf | SourceFormat::Bcf,
            MatrixKind::Haplotype,
            DosageSource::Dosage,
            true,
        ) => Err(GenoioError::unsupported(
            "sparse haplotype reads are intentionally unsupported for dosage-backed sources; use dense haplotype reads with sparse=False",
        )),
        (
            SourceFormat::Plink2,
            MatrixKind::Haplotype,
            DosageSource::Dosage,
            true,
        ) => Err(GenoioError::unsupported(
            PLINK2_SPARSE_DOSAGE_BACKED_HAPLOTYPE_UNSUPPORTED,
        )),
        (SourceFormat::Plink1, MatrixKind::Genotype, DosageSource::Dosage, false) => {
            Err(GenoioError::unsupported(
                "plink1 does not support dosage-backed genotype reads",
            ))
        }
        (SourceFormat::Plink1, MatrixKind::Haplotype, _, _) => Err(GenoioError::unsupported(
            "unsupported haplotype format: plink1",
        )),
        (SourceFormat::Bgen, MatrixKind::Genotype, _, true) => Err(GenoioError::unsupported(
            "bgen sparse genotype reads are not implemented",
        )),
        (SourceFormat::Bgen, MatrixKind::Haplotype, _, true) => Err(GenoioError::unsupported(
            "bgen sparse haplotype reads are not implemented; use dense haplotype reads with sparse=False",
        )),
        (SourceFormat::Bgen, MatrixKind::Genotype, DosageSource::Hardcall, false) => {
            Err(GenoioError::unsupported(
                "bgen hardcall genotype reads are not implemented; use dosage=\"dosage\"",
            ))
        }
        (SourceFormat::Bgen, MatrixKind::Haplotype, DosageSource::Hardcall, false) => {
            Err(GenoioError::unsupported(
                "bgen hardcall haplotype reads are not implemented; use dosage=\"dosage\" for source-encoded phased haplotype dosage",
            ))
        }
    }
}

fn read_source_metadata(
    source: &SourceMembers,
) -> Result<genoio_core::MetadataOutput, GenoioError> {
    match source {
        SourceMembers::Vcf { path, .. } => genoio_io::read_vcf_metadata(path),
        SourceMembers::Plink1 { bed, bim, fam } => genoio_io::read_plink1_metadata(bed, bim, fam),
        SourceMembers::Plink2 { pgen, pvar, psam } => {
            genoio_io::read_plink2_metadata(pgen, pvar, psam)
        }
        SourceMembers::Bgen { bgen, sample } => {
            genoio_io::read_bgen_metadata(bgen, sample.as_deref())
        }
    }
}

fn read_dense_matrix(
    source: &SourceMembers,
    kind: MatrixKind,
    options: &ReadOptions,
) -> Result<genoio_core::DenseGenotypeMatrix, GenoioError> {
    validate_read_support_impl(source.format(), kind, options.dosage, false)?;
    match (source, kind, options.dosage) {
        (SourceMembers::Vcf { path, .. }, MatrixKind::Genotype, DosageSource::Hardcall) => {
            genoio_io::read_vcf_dense_windowed_with_missing_policy(
                path,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.matrix_only,
            )
        }
        (SourceMembers::Vcf { path, .. }, MatrixKind::Genotype, DosageSource::Dosage) => {
            genoio_io::read_vcf_dosage_dense_windowed_with_missing_policy(
                path,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.matrix_only,
            )
        }
        (SourceMembers::Vcf { path, .. }, MatrixKind::Haplotype, DosageSource::Hardcall) => {
            genoio_io::read_vcf_haplotypes_dense_windowed_with_missing_policy(
                path,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.matrix_only,
            )
        }
        (SourceMembers::Plink1 { bed, bim, fam }, MatrixKind::Genotype, DosageSource::Hardcall) => {
            genoio_io::read_plink1_dense_windowed_with_missing_policy(
                bed,
                bim,
                fam,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.matrix_only,
            )
        }
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Genotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_plink2_dense_windowed_with_missing_policy(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.matrix_only,
        ),
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Genotype,
            DosageSource::Dosage,
        ) => genoio_io::read_plink2_dosage_dense_windowed_with_missing_policy(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.matrix_only,
        ),
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Haplotype,
            DosageSource::Hardcall,
        ) => genoio_io::read_plink2_haplotypes_dense_windowed_with_missing_policy(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.matrix_only,
        ),
        (
            SourceMembers::Plink2 { pgen, pvar, psam },
            MatrixKind::Haplotype,
            DosageSource::Dosage,
        ) => genoio_io::read_plink2_haplotypes_dosage_dense_windowed_with_missing_policy(
            pgen,
            pvar,
            psam,
            options.requested_samples.as_deref(),
            options.variant_filter.as_ref(),
            options.variant_window,
            options.missing,
            options.matrix_only,
        ),
        (SourceMembers::Bgen { bgen, sample }, MatrixKind::Genotype, DosageSource::Dosage) => {
            genoio_io::read_bgen_dosage_dense_windowed_with_missing_policy(
                bgen,
                sample.as_deref(),
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.matrix_only,
            )
        }
        (SourceMembers::Bgen { bgen, sample }, MatrixKind::Haplotype, DosageSource::Dosage) => {
            genoio_io::read_bgen_haplotypes_dosage_dense_windowed_with_missing_policy(
                bgen,
                sample.as_deref(),
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
                options.missing,
                options.matrix_only,
            )
        }
        _ => Err(GenoioError::internal_contract(
            "read support validation accepted unsupported dense dispatch",
        )),
    }
}

fn read_sparse_matrix(
    source: &SourceMembers,
    kind: MatrixKind,
    options: &ReadOptions,
) -> Result<genoio_core::SparseGenotypeMatrix, GenoioError> {
    validate_read_support_impl(source.format(), kind, options.dosage, true)?;
    match (source, kind, options.dosage) {
        (SourceMembers::Vcf { path, .. }, MatrixKind::Genotype, DosageSource::Hardcall) => {
            genoio_io::read_vcf_sparse_windowed(
                path,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
            )
        }
        (SourceMembers::Vcf { path, .. }, MatrixKind::Haplotype, DosageSource::Hardcall) => {
            genoio_io::read_vcf_haplotypes_sparse_windowed(
                path,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
            )
        }
        (SourceMembers::Plink1 { bed, bim, fam }, MatrixKind::Genotype, DosageSource::Hardcall) => {
            genoio_io::read_plink1_sparse_windowed(
                bed,
                bim,
                fam,
                options.requested_samples.as_deref(),
                options.variant_filter.as_ref(),
                options.variant_window,
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
        ),
        _ => Err(GenoioError::internal_contract(
            "read support validation accepted unsupported sparse dispatch",
        )),
    }
}

#[pymodule]
/// PyO3 extension module used by the pure-Python public API.
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "RustInvalidSourceError",
        module.py().get_type::<RustInvalidSourceError>(),
    )?;
    module.add(
        "RustUnsupportedRepresentationError",
        module.py().get_type::<RustUnsupportedRepresentationError>(),
    )?;
    module.add(
        "RustInvalidOptionError",
        module.py().get_type::<RustInvalidOptionError>(),
    )?;
    module.add(
        "RustMissingDataError",
        module.py().get_type::<RustMissingDataError>(),
    )?;
    module.add(
        "RustSampleFilterError",
        module.py().get_type::<RustSampleFilterError>(),
    )?;
    module.add(
        "RustInternalError",
        module.py().get_type::<RustInternalError>(),
    )?;
    module.add_function(wrap_pyfunction!(backend_name, module)?)?;
    module.add_function(wrap_pyfunction!(validate_read_support, module)?)?;
    module.add_function(wrap_pyfunction!(read_metadata, module)?)?;
    module.add_function(wrap_pyfunction!(read_dense, module)?)?;
    module.add_function(wrap_pyfunction!(read_sparse, module)?)?;
    module.add_function(wrap_pyfunction!(read_haplotypes_dense, module)?)?;
    module.add_function(wrap_pyfunction!(read_haplotypes_sparse, module)?)?;
    module.add_class::<ArrowMetadataFrame>()?;
    Ok(())
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

struct ReadOptions {
    requested_samples: Option<Vec<String>>,
    variant_filter: Option<genoio_core::VariantFilter>,
    variant_window: Option<genoio_core::VariantWindow>,
    dosage: DosageSource,
    missing: DenseMissingPolicy,
    return_samples: bool,
    return_variants: bool,
    matrix_only: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DosageSource {
    Hardcall,
    Dosage,
}

impl DosageSource {
    fn from_str(value: &str) -> Result<Self, GenoioError> {
        match value {
            "hardcall" => Ok(Self::Hardcall),
            "dosage" => Ok(Self::Dosage),
            other => Err(GenoioError::invalid_filter(format!(
                "unsupported dosage source: {other}"
            ))),
        }
    }
}

fn read_options(options: &Bound<'_, PyDict>) -> PyResult<ReadOptions> {
    Ok(ReadOptions {
        requested_samples: samples_option(options)?,
        variant_filter: variants_option(options)?,
        variant_window: variant_window_option(options)?,
        dosage: dosage_option(options)?,
        missing: missing_policy_option(options)?,
        return_samples: bool_option(options, "return_samples")?,
        return_variants: bool_option(options, "return_variants")?,
        matrix_only: required_bool_option(options, "matrix_only")?,
    })
}

fn samples_option(options: &Bound<'_, PyDict>) -> PyResult<Option<Vec<String>>> {
    let Some(value) = options.get_item("samples")? else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    value.extract::<Vec<String>>().map(Some)
}

fn variants_option(options: &Bound<'_, PyDict>) -> PyResult<Option<genoio_core::VariantFilter>> {
    let Some(value) = options.get_item("variants")? else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    // Python validates and serializes the public FilterExpr tree; Rust parses
    // the JSON-compatible IR again so malformed direct extension calls fail at
    // the same boundary as normal reads.
    let json = py_to_json_value(&value)?;
    genoio_core::VariantFilter::from_json_value(json)
        .map(Some)
        .map_err(genoio_error_to_py)
}

fn variant_window_option(
    options: &Bound<'_, PyDict>,
) -> PyResult<Option<genoio_core::VariantWindow>> {
    let Some(value) = options.get_item("variant_window")? else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    let dict = value.cast::<PyDict>()?;
    let start = dict
        .get_item("start")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing variant_window.start"))?
        .extract::<usize>()?;
    let len = dict
        .get_item("len")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing variant_window.len"))?
        .extract::<usize>()?;
    Ok(Some(genoio_core::VariantWindow { start, len }))
}

fn dosage_option(options: &Bound<'_, PyDict>) -> PyResult<DosageSource> {
    let value = options
        .get_item("dosage")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("missing option: dosage"))?;
    DosageSource::from_str(value.extract::<String>()?.as_str()).map_err(genoio_error_to_py)
}

fn missing_policy_option(options: &Bound<'_, PyDict>) -> PyResult<DenseMissingPolicy> {
    let Some(value) = options.get_item("missing")? else {
        return Ok(DenseMissingPolicy::Raise);
    };
    if value.is_none() {
        return Ok(DenseMissingPolicy::Raise);
    }
    dense_missing_policy_from_str(value.extract::<String>()?.as_str()).map_err(genoio_error_to_py)
}

fn dense_missing_policy_from_str(value: &str) -> Result<DenseMissingPolicy, GenoioError> {
    match value {
        "raise" => Ok(DenseMissingPolicy::Raise),
        "nan" => Ok(DenseMissingPolicy::Nan),
        "impute" => Ok(DenseMissingPolicy::Impute),
        other => Err(GenoioError::invalid_filter(format!(
            "unsupported missing-data policy: {other}"
        ))),
    }
}

fn bool_option(options: &Bound<'_, PyDict>, key: &str) -> PyResult<bool> {
    let Some(value) = options.get_item(key)? else {
        return Ok(false);
    };
    if value.is_none() {
        return Ok(false);
    }
    value.extract::<bool>()
}

fn required_bool_option(options: &Bound<'_, PyDict>, key: &str) -> PyResult<bool> {
    let value = options
        .get_item(key)?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("missing option: {key}")))?;
    if value.is_none() {
        return Err(pyo3::exceptions::PyTypeError::new_err(format!(
            "{key} must be a bool"
        )));
    }
    value.extract::<bool>()
}

fn metadata_to_py(py: Python<'_>, output: genoio_core::MetadataOutput) -> PyResult<Py<PyDict>> {
    let genoio_core::MetadataOutput {
        samples,
        variants,
        capabilities: source_capabilities,
    } = output;
    let dict = PyDict::new(py);

    let capabilities = PyDict::new(py);
    capabilities.set_item("supports_geno", source_capabilities.supports_geno)?;
    capabilities.set_item("supports_haplo", source_capabilities.supports_haplo)?;
    capabilities.set_item("phased", source_capabilities.phased)?;

    dict.set_item("samples", sample_records_to_py(py, samples, false)?)?;
    dict.set_item("variants", variant_records_to_py(py, variants)?)?;
    dict.set_item("capabilities", capabilities)?;
    Ok(dict.unbind())
}

fn dense_to_py(
    py: Python<'_>,
    output: genoio_core::DenseGenotypeMatrix,
    return_samples: bool,
    return_variants: bool,
    include_haplotype_sample_columns: bool,
) -> PyResult<Py<PyDict>> {
    let genoio_core::DenseGenotypeMatrix {
        n_samples,
        n_variants,
        values,
        layout,
        samples,
        variants,
        diagnostics: dense_diagnostics,
    } = output;

    let dict = PyDict::new(py);
    dict.set_item("values", f32_vec_to_numpy(py, values)?)?;
    dict.set_item("shape", (n_samples, n_variants))?;
    // Python assembly uses this tag to reshape variant-major buffers as views
    // instead of requiring Rust to allocate a physically transposed matrix.
    dict.set_item("values_layout", dense_layout_to_py(layout))?;
    if return_samples {
        dict.set_item(
            "samples",
            sample_records_to_py(py, samples, include_haplotype_sample_columns)?,
        )?;
    }
    if return_variants {
        dict.set_item("variants", variant_records_to_py(py, variants)?)?;
    }

    let diagnostics = PyDict::new(py);
    diagnostics.set_item("requested_samples", dense_diagnostics.requested_samples)?;
    diagnostics.set_item("retained_samples", dense_diagnostics.retained_samples)?;
    diagnostics.set_item("missing_samples", dense_diagnostics.missing_samples)?;
    diagnostics.set_item("candidate_variants", dense_diagnostics.candidate_variants)?;
    diagnostics.set_item("retained_variants", dense_diagnostics.retained_variants)?;
    diagnostics.set_item(
        "dropped_metadata_variants",
        dense_diagnostics.dropped_metadata_variants,
    )?;
    diagnostics.set_item(
        "dropped_genotype_variants",
        dense_diagnostics.dropped_genotype_variants,
    )?;
    dict.set_item("diagnostics", diagnostics)?;

    Ok(dict.unbind())
}

fn dense_layout_to_py(layout: DenseLayout) -> &'static str {
    match layout {
        DenseLayout::SampleMajor => "sample_major",
        DenseLayout::VariantMajor => "variant_major",
    }
}

fn sparse_to_py(
    py: Python<'_>,
    output: genoio_core::SparseGenotypeMatrix,
    return_samples: bool,
    return_variants: bool,
    include_haplotype_sample_columns: bool,
) -> PyResult<Py<PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("indptr", usize_vec_to_numpy_i64(py, output.indptr)?)?;
    dict.set_item("indices", usize_vec_to_numpy_i64(py, output.indices)?)?;
    dict.set_item("data", f32_vec_to_numpy(py, output.data)?)?;
    dict.set_item("shape", (output.n_rows, output.n_cols))?;
    if return_samples {
        dict.set_item(
            "samples",
            sample_records_to_py(py, output.samples, include_haplotype_sample_columns)?,
        )?;
    }
    if return_variants {
        dict.set_item("variants", variant_records_to_py(py, output.variants)?)?;
    }

    let diagnostics = PyDict::new(py);
    diagnostics.set_item("requested_samples", output.diagnostics.requested_samples)?;
    diagnostics.set_item("retained_samples", output.diagnostics.retained_samples)?;
    diagnostics.set_item("missing_samples", output.diagnostics.missing_samples)?;
    diagnostics.set_item("candidate_variants", output.diagnostics.candidate_variants)?;
    diagnostics.set_item("retained_variants", output.diagnostics.retained_variants)?;
    diagnostics.set_item(
        "dropped_metadata_variants",
        output.diagnostics.dropped_metadata_variants,
    )?;
    diagnostics.set_item(
        "dropped_genotype_variants",
        output.diagnostics.dropped_genotype_variants,
    )?;
    dict.set_item("diagnostics", diagnostics)?;

    Ok(dict.unbind())
}

fn f32_vec_to_numpy(py: Python<'_>, values: Vec<f32>) -> PyResult<Bound<'_, PyAny>> {
    Ok(vec_to_numpy(py, values))
}

fn usize_vec_to_numpy_i64(py: Python<'_>, values: Vec<usize>) -> PyResult<Bound<'_, PyAny>> {
    let values = values
        .into_iter()
        .map(|value| {
            i64::try_from(value).map_err(|_| {
                pyo3::exceptions::PyOverflowError::new_err(
                    "array index exceeds supported NumPy int64 range",
                )
            })
        })
        .collect::<PyResult<Vec<i64>>>()?;
    Ok(vec_to_numpy(py, values))
}

fn vec_to_numpy<'py, T>(py: Python<'py>, values: Vec<T>) -> Bound<'py, PyAny>
where
    T: Element,
{
    PyArray1::from_vec(py, values).into_any()
}

fn sample_records_to_py(
    py: Python<'_>,
    samples: Vec<genoio_core::SampleRecord>,
    include_haplotype_columns: bool,
) -> PyResult<Py<ArrowMetadataFrame>> {
    let columns = SampleMetadataColumns::from_records(samples, include_haplotype_columns);
    Py::new(
        py,
        ArrowMetadataFrame {
            batch: sample_columns_to_arrow_batch(columns)?,
        },
    )
}

fn variant_records_to_py(
    py: Python<'_>,
    variants: Vec<genoio_core::VariantRecord>,
) -> PyResult<Py<ArrowMetadataFrame>> {
    let columns = VariantMetadataColumns::from_records(variants);
    Py::new(
        py,
        ArrowMetadataFrame {
            batch: variant_columns_to_arrow_batch(columns)?,
        },
    )
}

fn sample_columns_to_arrow_batch(columns: SampleMetadataColumns) -> PyResult<RecordBatch> {
    let mut fields = vec![
        Field::new("fid", DataType::Utf8, true),
        Field::new("iid", DataType::Utf8, false),
        Field::new("father", DataType::Utf8, true),
        Field::new("mother", DataType::Utf8, true),
        Field::new("sex", DataType::Utf8, true),
        Field::new("phenotype", DataType::Utf8, true),
    ];
    let mut arrays: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from_iter(columns.fids)) as ArrayRef,
        Arc::new(StringArray::from_iter_values(columns.iids)) as ArrayRef,
        Arc::new(StringArray::from_iter(columns.fathers)) as ArrayRef,
        Arc::new(StringArray::from_iter(columns.mothers)) as ArrayRef,
        Arc::new(StringArray::from_iter(columns.sexes)) as ArrayRef,
        Arc::new(StringArray::from_iter(columns.phenotypes)) as ArrayRef,
    ];

    if let Some(source_sample_indices) = columns.source_sample_indices {
        fields.push(Field::new("source_sample_index", DataType::Int64, true));
        arrays.push(Arc::new(optional_usize_to_i64_array(source_sample_indices)?) as ArrayRef);
    }
    if let Some(haplotype_indices) = columns.haplotype_indices {
        fields.push(Field::new("haplotype_index", DataType::Int64, true));
        arrays.push(Arc::new(optional_usize_to_i64_array(haplotype_indices)?) as ArrayRef);
    }

    RecordBatch::try_new(Arc::new(Schema::new(fields)), arrays).map_err(arrow_error_to_py)
}

fn variant_columns_to_arrow_batch(columns: VariantMetadataColumns) -> PyResult<RecordBatch> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::Int64, false),
        Field::new("id", DataType::Utf8, false),
        Field::new("a0", DataType::Utf8, false),
        Field::new("a1", DataType::Utf8, false),
    ]));
    let positions = columns
        .positions
        .into_iter()
        .map(i64::from)
        .collect::<Vec<_>>();
    let arrays: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from_iter_values(columns.chroms)) as ArrayRef,
        Arc::new(Int64Array::from(positions)) as ArrayRef,
        Arc::new(StringArray::from_iter_values(columns.ids)) as ArrayRef,
        Arc::new(StringArray::from_iter_values(columns.a0s)) as ArrayRef,
        Arc::new(StringArray::from_iter_values(columns.a1s)) as ArrayRef,
    ];

    RecordBatch::try_new(schema, arrays).map_err(arrow_error_to_py)
}

fn optional_usize_to_i64_array(values: Vec<Option<usize>>) -> PyResult<Int64Array> {
    let values = values
        .into_iter()
        .map(|value| {
            value
                .map(|value| {
                    i64::try_from(value).map_err(|_| {
                        pyo3::exceptions::PyOverflowError::new_err(
                            "array index exceeds supported Arrow int64 range",
                        )
                    })
                })
                .transpose()
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(Int64Array::from(values))
}

fn arrow_error_to_py(error: ArrowError) -> PyErr {
    RustInternalError::new_err(format!("failed to build Arrow metadata payload: {error}"))
}

fn py_to_json_value(value: &Bound<'_, PyAny>) -> PyResult<serde_json::Value> {
    // PyO3 exposes bool as an int subclass, so check bool before PyInt to keep
    // filter IR values semantically stable.
    if value.is_none() {
        return Ok(serde_json::Value::Null);
    }
    if value.cast::<PyBool>().is_ok() {
        return Ok(serde_json::Value::Bool(value.extract::<bool>()?));
    }
    if value.cast::<PyString>().is_ok() {
        return Ok(serde_json::Value::String(value.extract::<String>()?));
    }
    if value.cast::<PyInt>().is_ok() {
        return Ok(serde_json::Value::Number(value.extract::<i64>()?.into()));
    }
    if value.cast::<PyFloat>().is_ok() {
        let number = serde_json::Number::from_f64(value.extract::<f64>()?).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("filter IR contains a non-finite float")
        })?;
        return Ok(serde_json::Value::Number(number));
    }
    if let Ok(dict) = value.cast::<PyDict>() {
        let mut object = serde_json::Map::new();
        for (key, item) in dict.iter() {
            object.insert(key.extract::<String>()?, py_to_json_value(&item)?);
        }
        return Ok(serde_json::Value::Object(object));
    }
    if let Ok(list) = value.cast::<PyList>() {
        return list
            .iter()
            .map(|item| py_to_json_value(&item))
            .collect::<PyResult<Vec<_>>>()
            .map(serde_json::Value::Array);
    }
    if let Ok(tuple) = value.cast::<PyTuple>() {
        return tuple
            .iter()
            .map(|item| py_to_json_value(&item))
            .collect::<PyResult<Vec<_>>>()
            .map(serde_json::Value::Array);
    }
    Err(pyo3::exceptions::PyValueError::new_err(
        "filter IR must contain only JSON-compatible values",
    ))
}

fn genoio_error_to_py(error: GenoioError) -> PyErr {
    match error {
        GenoioError::Io { .. } | GenoioError::InvalidSource { .. } => {
            RustInvalidSourceError::new_err(error.to_string())
        }
        GenoioError::UnsupportedRepresentation { .. } => {
            RustUnsupportedRepresentationError::new_err(error.to_string())
        }
        GenoioError::SampleFilter { .. } => RustSampleFilterError::new_err(error.to_string()),
        GenoioError::MissingData { .. } => RustMissingDataError::new_err(error.to_string()),
        GenoioError::InvalidFilter { .. } => RustInvalidOptionError::new_err(error.to_string()),
        GenoioError::InternalContract { .. } => RustInternalError::new_err(error.to_string()),
    }
}
