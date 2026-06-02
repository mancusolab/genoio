// pattern: Imperative Shell

use pyo3::prelude::*;

#[pyfunction]
fn backend_name() -> &'static str {
    genoio_io::backend_name()
}

#[pymodule]
fn _rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(backend_name, module)?)?;
    Ok(())
}
