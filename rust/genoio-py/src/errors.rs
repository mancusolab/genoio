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

#[cfg(test)]
mod tests {
    use std::sync::mpsc;
    use std::thread;
    use std::time::Duration;

    use genoio_core::GenoioError;
    use pyo3::Python;

    use super::{catch_internal_panic, run_without_gil, RustInternalError};

    #[test]
    fn pbr_rust_panic_001_catch_internal_panic_maps_to_rust_internal_error() {
        Python::attach(|py| {
            let error = catch_internal_panic::<()>(|| panic!("intentional boundary panic"))
                .expect_err("intentional panic should be translated");

            assert!(error.is_instance_of::<RustInternalError>(py));
            assert!(error.to_string().contains("intentional boundary panic"));
        });
    }

    #[test]
    fn pbr_gil_001_run_without_gil_allows_another_python_thread_to_attach() {
        let (result, worker) = Python::attach(|py| {
            let (attached_tx, attached_rx) = mpsc::channel();
            let worker = thread::spawn(move || {
                Python::attach(|_| {
                    let _ = attached_tx.send(());
                });
            });

            let result = run_without_gil(py, move || {
                attached_rx
                    .recv_timeout(Duration::from_secs(2))
                    .map_err(|_| {
                        GenoioError::internal_contract(
                            "another Python thread could not attach while Rust was waiting",
                        )
                    })?;
                Ok(())
            });
            (result, worker)
        });

        worker
            .join()
            .expect("Python attachment worker should not panic");
        result.expect("run_without_gil should release the GIL while Rust is waiting");
    }
}
