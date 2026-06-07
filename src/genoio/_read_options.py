# pattern: Functional Core

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import DTypeLike

from ._errors import InvalidOptionError, UnsupportedRepresentation
from ._filters import FilterExpr, id_in

_SUPPORTED_KINDS = {"geno", "haplo"}
_SUPPORTED_DOSAGE_SOURCES = {"hardcall", "dosage"}
_SUPPORTED_MISSING_POLICIES = {"nan", "raise", "impute"}

_READ_OPTION_DEFAULTS: dict[str, Any] = {
    "kind": "geno",
    "dosage": "hardcall",
    "sparse": False,
    "variants": None,
    "samples": None,
    "missing": None,
    "dtype": "float32",
    "return_samples": False,
    "return_variants": False,
}


@dataclass(frozen=True, slots=True)
class _ReadOptions:
    kind: str
    dosage: str
    sparse: object
    variants: Any
    samples: list[str] | tuple[str, ...] | set[str] | None
    missing: object
    dtype: DTypeLike
    return_samples: bool
    return_variants: bool


@dataclass(frozen=True, slots=True)
class _ValidatedReadOptions:
    # Values that have crossed this boundary are safe to send to Rust or use
    # for Python-side matrix assembly without rechecking user input types.
    dtype: np.dtype[Any]
    missing: str
    sparse_format: str | None
    variant_filter_ir: dict[str, Any] | None


def _validate_read_options(options: _ReadOptions) -> _ValidatedReadOptions:
    _validate_kind(options.kind)
    _validate_dosage_source(options.dosage)
    sparse_format = _validate_sparse(options.sparse)
    normalized_missing = _normalize_missing(options.missing, sparse_format=sparse_format)
    variant_filter_ir = _variant_filter_ir(options.variants)
    normalized_dtype = _normalize_dtype(options.dtype)
    _validate_sparse_missing_compatibility(sparse_format, normalized_missing)
    _validate_missing_dtype_compatibility(normalized_missing, normalized_dtype)
    _validate_sample_filter(options.samples)
    _validate_bool_option("return_samples", options.return_samples)
    _validate_bool_option("return_variants", options.return_variants)
    return _ValidatedReadOptions(
        dtype=normalized_dtype,
        missing=normalized_missing,
        sparse_format=sparse_format,
        variant_filter_ir=variant_filter_ir,
    )


def _validate_kind(kind: str) -> None:
    if not isinstance(kind, str) or kind not in _SUPPORTED_KINDS:
        raise UnsupportedRepresentation(f"unsupported genotype kind: {kind}")


def _validate_dosage_source(dosage: str) -> None:
    if not isinstance(dosage, str) or dosage not in _SUPPORTED_DOSAGE_SOURCES:
        raise InvalidOptionError(f"unsupported dosage source: {dosage!r}")


def _validate_sparse(sparse: object) -> str | None:
    match sparse:
        case False:
            return None
        case True:
            return "csc"
        case "csc":
            return "csc"
        case "csr":
            return "csr"
        case _:
            raise InvalidOptionError(f"unsupported sparse option: {sparse!r}")


def _variant_filter_ir(variants: Any) -> dict[str, Any] | None:
    if variants is None:
        return None
    if isinstance(variants, FilterExpr):
        return variants.to_ir()
    if callable(variants):
        raise InvalidOptionError("variants must be a serializable filter expression or variant ID iterable")
    if isinstance(variants, str) or not isinstance(variants, Iterable):
        raise InvalidOptionError("variants must be a serializable filter expression or variant ID iterable")
    variant_ids: list[str] = []
    for variant_id in variants:
        if not isinstance(variant_id, str):
            raise InvalidOptionError("variant ID filters must contain only strings")
        variant_ids.append(variant_id)
    # ID lists share the same Rust filter path as declarative predicates.
    # This keeps variant selection order source-defined for both forms.
    return id_in(variant_ids).to_ir()


def _validate_sample_filter(samples: list[str] | tuple[str, ...] | set[str] | None) -> None:
    if samples is None:
        return
    if not isinstance(samples, list | tuple | set):
        raise InvalidOptionError("samples must be a list, tuple, or set of sample IDs")
    if any(not isinstance(sample, str) for sample in samples):
        raise InvalidOptionError("samples must contain only sample ID strings")
    if len(samples) != len(set(samples)):
        raise InvalidOptionError("samples must not contain duplicate sample IDs")


def _validate_missing(missing: object) -> str:
    if not isinstance(missing, str) or missing not in _SUPPORTED_MISSING_POLICIES:
        raise InvalidOptionError(f"unsupported missing-data policy: {missing}")
    return missing


def _normalize_missing(missing: object, *, sparse_format: str | None) -> str:
    if missing is None:
        return "raise" if sparse_format is not None else "nan"
    return _validate_missing(missing)


def _normalize_dtype(dtype: DTypeLike) -> np.dtype[Any]:
    try:
        return np.dtype(dtype)
    except TypeError as error:
        raise InvalidOptionError(f"invalid dtype: {dtype!r}") from error


def _validate_missing_dtype_compatibility(missing: str, dtype: np.dtype[Any]) -> None:
    if missing in {"nan", "impute"} and not np.issubdtype(dtype, np.floating):
        raise InvalidOptionError(f'missing="{missing}" requires a floating dtype')


def _validate_sparse_missing_compatibility(sparse_format: str | None, missing: str) -> None:
    if sparse_format is not None and missing in {"nan", "impute"}:
        raise InvalidOptionError("this release does not store sparse missing values; use missing='raise'")


def _validate_bool_option(name: str, value: bool) -> None:
    if not isinstance(value, bool):
        raise InvalidOptionError(f"{name} must be a bool")


def _read_options_with_defaults(read_options: Mapping[str, object]) -> _ReadOptions:
    unknown = set(read_options) - set(_READ_OPTION_DEFAULTS)
    if unknown:
        keys = ", ".join(sorted(unknown))
        raise InvalidOptionError(f"unsupported option(s): {keys}")
    # Iterator methods accept **read_options, so they reconstruct the same
    # option object that read() would have built from explicit parameters.
    merged: dict[str, Any] = {**_READ_OPTION_DEFAULTS, **read_options}
    return _ReadOptions(
        kind=merged["kind"],
        dosage=merged["dosage"],
        sparse=merged["sparse"],
        variants=merged["variants"],
        samples=merged["samples"],
        missing=merged["missing"],
        dtype=merged["dtype"],
        return_samples=merged["return_samples"],
        return_variants=merged["return_variants"],
    )
