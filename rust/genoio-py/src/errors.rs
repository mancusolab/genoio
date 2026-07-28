// pattern: Imperative Shell

use std::any::Any as PanicPayload;
use std::panic::{self, AssertUnwindSafe};

use genoio_core::GenoioError;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

pyo3::create_exception!(genoio_py, RustInvalidSourceError, PyException);
pyo3::create_exception!(genoio_py, RustUnsupportedRepresentationError, PyException);
pyo3::create_exception!(genoio_py, RustInvalidOptionError, PyException);
pyo3::create_exception!(genoio_py, RustMissingDataError, PyException);
pyo3::create_exception!(genoio_py, RustSampleFilterError, PyException);
pyo3::create_exception!(genoio_py, RustInternalError, PyException);

pub(crate) fn catch_internal_panic<T>(f: impl FnOnce() -> PyResult<T>) -> PyResult<T> {
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

pub(crate) fn run_without_gil<T: Send>(
    py: Python<'_>,
    f: impl FnOnce() -> Result<T, GenoioError> + Send,
) -> PyResult<T> {
    // Only Rust-owned source paths, filters, options, and output structs may
    // cross this boundary. Python argument parsing and result construction stay
    // under the GIL.
    py.detach(f).map_err(genoio_error_to_py)
}

pub(crate) fn genoio_error_to_py(error: GenoioError) -> PyErr {
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

fn panic_message(payload: &(dyn PanicPayload + Send)) -> &str {
    if let Some(message) = payload.downcast_ref::<&'static str>() {
        message
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.as_str()
    } else {
        "unknown panic payload"
    }
}
