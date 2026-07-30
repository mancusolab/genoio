// pattern: Imperative Shell

use std::sync::{Mutex, MutexGuard};

use genoio_core::GenoioError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::errors::{catch_internal_panic, run_without_gil};
use crate::options::block_options;
use crate::output::block_output_to_py;
use crate::source::{block_read_options, source_members, MatrixKind};

#[pyclass(frozen, module = "genoio._rust", name = "_BlockReader")]
pub(crate) struct PyBlockReader {
    reader: Mutex<Option<genoio_io::BlockReader>>,
    return_samples: bool,
    return_variants: bool,
}

#[pymethods]
impl PyBlockReader {
    #[new]
    #[pyo3(signature = (format, members, kind, sparse, options, block_size))]
    fn new(
        py: Python<'_>,
        format: &str,
        members: &Bound<'_, PyDict>,
        kind: &str,
        sparse: bool,
        options: &Bound<'_, PyDict>,
        block_size: usize,
    ) -> PyResult<Self> {
        catch_internal_panic(|| {
            let source = source_members(format, members)?.into_block_source();
            let kind = MatrixKind::from_str(kind).map_err(crate::errors::genoio_error_to_py)?;
            let options = block_options(options)?;
            let return_samples = options.return_samples;
            let return_variants = options.return_variants;
            let options = block_read_options(kind, sparse, options);
            let reader = run_without_gil(py, move || {
                genoio_io::BlockReader::open(source, options, block_size)
            })?;

            Ok(Self {
                reader: Mutex::new(Some(reader)),
                return_samples,
                return_variants,
            })
        })
    }

    fn next_block(&self, py: Python<'_>) -> PyResult<Option<Py<PyDict>>> {
        catch_internal_panic(|| {
            let output = run_without_gil(py, || {
                let mut guard = self.lock_reader()?;
                let output = match guard.as_mut() {
                    Some(reader) => reader.next_block(),
                    None => Ok(None),
                };
                drop(guard);
                output
            })?;

            output
                .map(|output| {
                    block_output_to_py(py, output, self.return_samples, self.return_variants)
                })
                .transpose()
        })
    }

    fn close(&self, py: Python<'_>) -> PyResult<()> {
        catch_internal_panic(|| {
            run_without_gil(py, || {
                let mut guard = self.lock_reader()?;
                let reader = guard.take();
                drop(guard);
                drop(reader);
                Ok(())
            })
        })
    }
}

impl PyBlockReader {
    fn lock_reader(&self) -> Result<MutexGuard<'_, Option<genoio_io::BlockReader>>, GenoioError> {
        self.reader.lock().map_err(|_| {
            GenoioError::internal_contract(
                "native block reader lock is poisoned; the session was closed fail-closed",
            )
        })
    }
}

#[expect(
    dead_code,
    reason = "compile-time bridge ownership assertions are intentionally not called"
)]
fn assert_bridge_type_contracts() {
    fn assert_send<T: Send>() {}
    fn assert_send_sync<T: Send + Sync>() {}

    assert_send::<genoio_io::BlockReader>();
    assert_send_sync::<PyBlockReader>();
}

#[cfg(test)]
mod tests {
    use std::panic::{self, AssertUnwindSafe};

    use pyo3::Python;

    use super::PyBlockReader;
    use crate::errors::RustInternalError;

    #[test]
    fn pbr_py_error_001_poisoned_reader_fails_closed_as_internal_error() {
        let reader = PyBlockReader {
            reader: std::sync::Mutex::new(None),
            return_samples: false,
            return_variants: false,
        };
        let panic_result = panic::catch_unwind(AssertUnwindSafe(|| {
            let _guard = reader
                .reader
                .lock()
                .expect("test reader mutex should initially lock");
            panic!("poison test reader");
        }));
        assert!(panic_result.is_err());

        Python::attach(|py| {
            let error = reader
                .next_block(py)
                .expect_err("poisoned reader should fail closed");
            assert!(error.is_instance_of::<RustInternalError>(py));
        });
    }
}
