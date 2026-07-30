// pattern: Functional Core

use std::collections::BTreeSet;

use serde::Deserialize;
use serde_json::Value;

use crate::GenoioError;

#[derive(Debug, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
pub(super) enum RawExpr {
    Predicate {
        name: String,
        #[serde(default)]
        params: Value,
    },
    And {
        left: Box<RawExpr>,
        right: Box<RawExpr>,
    },
    Or {
        left: Box<RawExpr>,
        right: Box<RawExpr>,
    },
    Not {
        expr: Box<RawExpr>,
    },
}

pub(super) fn expect_no_params(params: &Value) -> Result<(), GenoioError> {
    let object = params_object(params)?;
    if object.is_empty() {
        Ok(())
    } else {
        Err(GenoioError::invalid_source(
            "<filter>",
            "predicate does not accept parameters",
        ))
    }
}

pub(super) fn required_string(params: &Value, key: &str) -> Result<String, GenoioError> {
    match params_object(params)?.get(key) {
        Some(Value::String(value)) if !value.is_empty() => Ok(value.clone()),
        _ => Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} must be a non-empty string"),
        )),
    }
}

pub(super) fn required_string_set(
    params: &Value,
    key: &str,
) -> Result<BTreeSet<String>, GenoioError> {
    match params_object(params)?.get(key) {
        Some(Value::Array(values)) => {
            let mut set = BTreeSet::new();
            for value in values {
                let Value::String(text) = value else {
                    return Err(GenoioError::invalid_source(
                        "<filter>",
                        format!("predicate parameter {key:?} must contain only strings"),
                    ));
                };
                if !set.insert(text.clone()) {
                    return Err(GenoioError::invalid_source(
                        "<filter>",
                        format!("predicate parameter {key:?} must not contain duplicates"),
                    ));
                }
            }
            Ok(set)
        }
        _ => Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} must be a string array"),
        )),
    }
}

pub(super) fn optional_rate(params: &Value, key: &str) -> Result<Option<f32>, GenoioError> {
    match params_object(params)?.get(key) {
        Some(value) => Ok(Some(value_to_rate(key, value)?)),
        None => Ok(None),
    }
}

pub(super) fn required_rate(params: &Value, key: &str) -> Result<f32, GenoioError> {
    match optional_rate(params, key)? {
        Some(value) => Ok(value),
        None => Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} is required"),
        )),
    }
}

pub(super) fn optional_nonnegative_f32(
    params: &Value,
    key: &str,
) -> Result<Option<f32>, GenoioError> {
    let Some(value) = params.get(key) else {
        return Ok(None);
    };
    let number = value.as_f64().ok_or_else(|| {
        GenoioError::invalid_source(
            "<filter>",
            format!("{key} must be a non-negative finite number"),
        )
    })?;
    if !number.is_finite() || number < 0.0 {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("{key} must be a non-negative finite number"),
        ));
    }
    Ok(Some(number as f32))
}

pub(super) fn optional_u32(params: &Value, key: &str) -> Result<Option<u32>, GenoioError> {
    match params_object(params)?.get(key) {
        Some(value) => {
            let Some(number) = value.as_u64() else {
                return Err(GenoioError::invalid_source(
                    "<filter>",
                    format!("predicate parameter {key:?} must be a non-negative integer"),
                ));
            };
            Ok(Some(u32::try_from(number).map_err(|_| {
                GenoioError::invalid_source(
                    "<filter>",
                    format!("predicate parameter {key:?} is out of range"),
                )
            })?))
        }
        None => Ok(None),
    }
}

pub(super) fn validate_range<T: PartialOrd>(
    name: &str,
    min: Option<T>,
    max: Option<T>,
) -> Result<(), GenoioError> {
    if min.is_none() && max.is_none() {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("{name} predicate requires at least one threshold"),
        ));
    }
    if min.zip(max).is_some_and(|(min, max)| min > max) {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("{name} predicate min must be <= max"),
        ));
    }
    Ok(())
}

pub(super) fn parse_region(value: &str) -> Result<(String, u32, u32), GenoioError> {
    let Some((chrom, coordinates)) = value.split_once(':') else {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "invalid region syntax; expected chrom:start-end",
        ));
    };
    let Some((start_text, end_text)) = coordinates.split_once('-') else {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "invalid region syntax; expected chrom:start-end",
        ));
    };
    let start = start_text.parse::<u32>().map_err(|error| {
        GenoioError::invalid_source(
            "<filter>",
            format!("invalid region start coordinate: {error}"),
        )
    })?;
    let end = end_text.parse::<u32>().map_err(|error| {
        GenoioError::invalid_source(
            "<filter>",
            format!("invalid region end coordinate: {error}"),
        )
    })?;
    if chrom.is_empty() || start == 0 || end < start {
        return Err(GenoioError::invalid_source(
            "<filter>",
            "invalid region coordinates; expected 1-based start <= end",
        ));
    }
    Ok((chrom.to_string(), start, end))
}

fn params_object(params: &Value) -> Result<&serde_json::Map<String, Value>, GenoioError> {
    match params {
        Value::Object(object) => Ok(object),
        _ => Err(GenoioError::invalid_source(
            "<filter>",
            "predicate params must be a JSON object",
        )),
    }
}

fn value_to_rate(key: &str, value: &Value) -> Result<f32, GenoioError> {
    let Some(number) = value.as_f64() else {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} must be numeric"),
        ));
    };
    if !(0.0..=1.0).contains(&number) {
        return Err(GenoioError::invalid_source(
            "<filter>",
            format!("predicate parameter {key:?} must be between 0 and 1"),
        ));
    }
    Ok(number as f32)
}
