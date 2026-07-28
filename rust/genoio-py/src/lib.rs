// pattern: Imperative Shell

#![cfg_attr(
    not(test),
    deny(clippy::expect_used, clippy::panic, clippy::unwrap_used)
)]

//! PyO3 adapter for the private `genoio._rust` extension module.
//!
//! Python-facing registration stays here. Argument parsing, source dispatch,
//! reader entry points, error translation, and output conversion are private
//! modules so the extension boundary remains explicit.

mod errors;
mod options;
mod output;
mod reads;
mod source;

use pyo3::prelude::*;

use errors::{
    RustInternalError, RustInvalidOptionError, RustInvalidSourceError, RustMissingDataError,
    RustSampleFilterError, RustUnsupportedRepresentationError,
};
use output::ArrowMetadataFrame;
use reads::{
    backend_name, read_dense, read_haplotypes_dense, read_haplotypes_sparse, read_metadata,
    read_sparse, validate_read_support,
};

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
