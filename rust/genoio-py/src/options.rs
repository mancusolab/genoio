// pattern: Imperative Shell

use genoio_core::{DenseMissingPolicy, GenoioError, VariantFilter, VariantWindow};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBool, PyDict, PyFloat, PyInt, PyList, PyString, PyTuple};

use crate::errors::genoio_error_to_py;

pub(crate) struct ReadOptions {
    pub(crate) requested_samples: Option<Vec<String>>,
    pub(crate) variant_filter: Option<VariantFilter>,
    pub(crate) variant_window: Option<VariantWindow>,
    pub(crate) dosage: DosageSource,
    pub(crate) missing: DenseMissingPolicy,
    pub(crate) return_samples: bool,
    pub(crate) return_variants: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DosageSource {
    Hardcall,
    Dosage,
}

impl DosageSource {
    pub(crate) fn from_str(value: &str) -> Result<Self, GenoioError> {
        match value {
            "hardcall" => Ok(Self::Hardcall),
            "dosage" => Ok(Self::Dosage),
            other => Err(GenoioError::invalid_filter(format!(
                "unsupported dosage source: {other}"
            ))),
        }
    }
}

pub(crate) fn read_options(options: &Bound<'_, PyDict>) -> PyResult<ReadOptions> {
    Ok(ReadOptions {
        requested_samples: samples_option(options)?,
        variant_filter: variants_option(options)?,
        variant_window: variant_window_option(options)?,
        dosage: dosage_option(options)?,
        missing: missing_policy_option(options)?,
        return_samples: bool_option(options, "return_samples")?,
        return_variants: bool_option(options, "return_variants")?,
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

fn variants_option(options: &Bound<'_, PyDict>) -> PyResult<Option<VariantFilter>> {
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
    VariantFilter::from_json_value(json)
        .map(Some)
        .map_err(genoio_error_to_py)
}

fn variant_window_option(options: &Bound<'_, PyDict>) -> PyResult<Option<VariantWindow>> {
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
    Ok(Some(VariantWindow { start, len }))
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
