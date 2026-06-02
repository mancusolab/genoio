# pattern: Mixed (needs refactoring)
# Reason: Public entrypoints call resolve_source(), which performs filesystem checks;
# Dataset validation helpers remain pure.

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ._errors import InvalidOptionError, UnsupportedRepresentation
from ._source import ResolvedSource, resolve_source

_SUPPORTED_KINDS = {"geno", "haplo"}
_SUPPORTED_SPARSE = {False, True, "csc", "csr"}
_SUPPORTED_MISSING_POLICIES = {"nan", "raise", "impute"}


@dataclass(frozen=True)
class Dataset:
    source: ResolvedSource

    def read(
        self,
        *,
        kind: str = "geno",
        sparse: bool | str = False,
        variants: Any = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: str = "nan",
        dtype: Any = "float32",
        return_samples: bool = False,
        return_variants: bool = False,
    ) -> Any:
        _validate_read_options(
            kind=kind,
            sparse=sparse,
            variants=variants,
            samples=samples,
            missing=missing,
            dtype=dtype,
            return_samples=return_samples,
            return_variants=return_variants,
        )
        raise NotImplementedError("genotype reading is implemented in a later phase")

    def samples(self, **options: Any) -> Any:
        _reject_options(options)
        raise NotImplementedError("sample metadata reading is implemented in a later phase")

    def variants(self, *, stats: Any = None, **options: Any) -> Any:
        _validate_variant_stats(stats)
        _reject_options(options)
        raise NotImplementedError("variant metadata reading is implemented in a later phase")

    def blocks(self, size: int, **read_options: Any) -> Any:
        _validate_block_size(size)
        _validate_read_options_from_mapping(read_options)
        raise NotImplementedError("block iteration is implemented in a later phase")


def open(path: str | Path, format: str | None = None) -> Dataset:
    return Dataset(source=resolve_source(path, format=format))


def read(
    path: str | Path,
    *,
    format: str | None = None,
    kind: str = "geno",
    sparse: bool | str = False,
    variants: Any = None,
    samples: list[str] | tuple[str, ...] | set[str] | None = None,
    missing: str = "nan",
    dtype: Any = "float32",
    return_samples: bool = False,
    return_variants: bool = False,
) -> Any:
    return open(path, format=format).read(
        kind=kind,
        sparse=sparse,
        variants=variants,
        samples=samples,
        missing=missing,
        dtype=dtype,
        return_samples=return_samples,
        return_variants=return_variants,
    )


def samples(path: str | Path, *, format: str | None = None, **options: Any) -> Any:
    return open(path, format=format).samples(**options)


def variants(path: str | Path, *, format: str | None = None, stats: Any = None, **options: Any) -> Any:
    return open(path, format=format).variants(stats=stats, **options)


def _validate_read_options(
    *,
    kind: str,
    sparse: bool | str,
    variants: Any,
    samples: list[str] | tuple[str, ...] | set[str] | None,
    missing: str,
    dtype: Any,
    return_samples: bool,
    return_variants: bool,
) -> None:
    _validate_kind(kind)
    _validate_sparse(sparse)
    _validate_variant_filter(variants)
    _validate_sample_filter(samples)
    _validate_missing(missing)
    normalized_dtype = _normalize_dtype(dtype)
    _validate_missing_dtype_compatibility(missing, normalized_dtype)
    _validate_bool_option("return_samples", return_samples)
    _validate_bool_option("return_variants", return_variants)


def _reject_options(options: dict[str, Any]) -> None:
    if options:
        keys = ", ".join(sorted(options))
        raise InvalidOptionError(f"unsupported option(s): {keys}")


def _validate_variant_stats(stats: Any) -> None:
    if stats is not None:
        raise InvalidOptionError("variant stats are not implemented until a later phase")


def _validate_kind(kind: str) -> None:
    if not isinstance(kind, str) or kind not in _SUPPORTED_KINDS:
        raise UnsupportedRepresentation(f"unsupported genotype kind: {kind}")


def _validate_sparse(sparse: bool | str) -> None:
    if not isinstance(sparse, bool | str) or sparse not in _SUPPORTED_SPARSE:
        raise UnsupportedRepresentation(f"unsupported sparse representation: {sparse}")


def _validate_variant_filter(variants: Any) -> None:
    if variants is None:
        return
    to_ir = getattr(variants, "to_ir", None)
    if to_ir is None:
        return
    ir = to_ir()
    if not isinstance(ir, dict):
        raise InvalidOptionError("variant filter must serialize to a dictionary IR")


def _validate_sample_filter(samples: list[str] | tuple[str, ...] | set[str] | None) -> None:
    if samples is None:
        return
    if not isinstance(samples, list | tuple | set):
        raise InvalidOptionError("samples must be a list, tuple, or set of sample IDs")
    if any(not isinstance(sample, str) for sample in samples):
        raise InvalidOptionError("samples must contain only sample ID strings")
    if len(samples) != len(set(samples)):
        raise InvalidOptionError("samples must not contain duplicate sample IDs")


def _validate_missing(missing: str) -> None:
    if not isinstance(missing, str) or missing not in _SUPPORTED_MISSING_POLICIES:
        raise InvalidOptionError(f"unsupported missing-data policy: {missing}")


def _normalize_dtype(dtype: Any) -> np.dtype[Any]:
    try:
        return np.dtype(dtype)
    except TypeError as error:
        raise InvalidOptionError(f"invalid dtype: {dtype!r}") from error


def _validate_missing_dtype_compatibility(missing: str, dtype: np.dtype[Any]) -> None:
    if missing in {"nan", "impute"} and not np.issubdtype(dtype, np.floating):
        raise InvalidOptionError(f'missing="{missing}" requires a floating dtype')


def _validate_bool_option(name: str, value: bool) -> None:
    if not isinstance(value, bool):
        raise InvalidOptionError(f"{name} must be a bool")


def _validate_block_size(size: int) -> None:
    if not isinstance(size, int) or isinstance(size, bool) or size < 1:
        raise InvalidOptionError("block size must be a positive integer")


def _validate_read_options_from_mapping(read_options: dict[str, Any]) -> None:
    defaults = {
        "kind": "geno",
        "sparse": False,
        "variants": None,
        "samples": None,
        "missing": "nan",
        "dtype": "float32",
        "return_samples": False,
        "return_variants": False,
    }
    unknown = set(read_options) - set(defaults)
    if unknown:
        keys = ", ".join(sorted(unknown))
        raise InvalidOptionError(f"unsupported option(s): {keys}")
    _validate_read_options(**(defaults | read_options))
