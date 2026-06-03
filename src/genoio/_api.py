# pattern: Mixed (needs refactoring)
# Reason: Public entrypoints call resolve_source(), which performs filesystem checks;
# Dataset validation helpers remain pure.

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from . import _rust
from ._assembly import (
    dense_array_from_rust,
    read_result_tuple,
    samples_frame,
    sparse_matrix_from_rust,
    variants_frame,
)
from ._errors import (
    InvalidOptionError,
    InvalidSourceError,
    MissingDataError,
    SampleFilterError,
    UnsupportedFormatError,
    UnsupportedRepresentation,
)
from ._filters import FilterExpr, id_in
from ._source import ResolvedSource, resolve_source

_SUPPORTED_KINDS = {"geno", "haplo"}
_SUPPORTED_MISSING_POLICIES = {"nan", "raise", "impute"}


class _DefaultMissing:
    def __repr__(self) -> str:
        return "DEFAULT_MISSING"


_DEFAULT_MISSING = _DefaultMissing()


@dataclass(frozen=True)
class _ValidatedReadOptions:
    dtype: np.dtype[Any]
    missing: str
    sparse_format: str | None
    variant_filter_ir: dict[str, Any] | None


@dataclass(frozen=True)
class Dataset:
    source: ResolvedSource
    _metadata_cache: dict[str, Any] | None = field(default=None, init=False, compare=False, repr=False)

    def read(
        self,
        *,
        kind: str = "geno",
        sparse: bool | str = False,
        variants: Any = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Any = _DEFAULT_MISSING,
        dtype: Any = "float32",
        return_samples: bool = False,
        return_variants: bool = False,
    ) -> Any:
        validated_options = _validate_read_options(
            kind=kind,
            sparse=sparse,
            variants=variants,
            samples=samples,
            missing=missing,
            dtype=dtype,
            return_samples=return_samples,
            return_variants=return_variants,
        )
        _reject_deferred_plink2_decode(self.source.format.value)
        if kind != "geno":
            capabilities = self._metadata()["capabilities"]
            if not capabilities["supports_haplo"]:
                raise _unsupported_haplotype_source(self.source.format.value)

        return self._read_validated(
            kind=kind,
            samples=samples,
            return_samples=return_samples,
            return_variants=return_variants,
            validated_options=validated_options,
            variant_window=None,
        )

    def samples(self, **options: Any) -> Any:
        _reject_options(options)
        return samples_frame(self._metadata()["samples"])

    def variants(self, *, stats: Any = None, **options: Any) -> Any:
        _validate_variant_stats(stats)
        _reject_options(options)
        return variants_frame(self._metadata()["variants"])

    def blocks(self, size: int, **read_options: Any) -> Any:
        _validate_block_size(size)
        normalized_options = _read_options_with_defaults(read_options)
        validated_options = _validate_read_options(**normalized_options)
        _reject_deferred_plink2_decode(self.source.format.value)
        if normalized_options["kind"] != "geno":
            capabilities = self._metadata()["capabilities"]
            if not capabilities["supports_haplo"]:
                raise _unsupported_haplotype_source(self.source.format.value)
        return self._block_iterator(size, normalized_options, validated_options)

    def _block_iterator(
        self,
        size: int,
        read_options: dict[str, Any],
        validated_options: _ValidatedReadOptions,
    ) -> Any:
        start = 0
        while True:
            block = self._read_block_from_rust(
                kind=read_options["kind"],
                samples=read_options["samples"],
                return_samples=read_options["return_samples"],
                return_variants=read_options["return_variants"],
                validated_options=validated_options,
                variant_window={"start": start, "len": size},
            )
            genotype_matrix = block[0] if isinstance(block, tuple) else block
            if genotype_matrix.shape[1] == 0:
                break
            yield block
            if genotype_matrix.shape[1] < size:
                break
            start += size

    def _read_block_from_rust(
        self,
        *,
        kind: str,
        samples: list[str] | tuple[str, ...] | set[str] | None,
        return_samples: bool,
        return_variants: bool,
        validated_options: _ValidatedReadOptions,
        variant_window: dict[str, int],
    ) -> Any:
        return self._read_validated(
            kind=kind,
            samples=samples,
            return_samples=return_samples,
            return_variants=return_variants,
            validated_options=validated_options,
            variant_window=variant_window,
        )

    def _read_validated(
        self,
        *,
        kind: str,
        samples: list[str] | tuple[str, ...] | set[str] | None,
        return_samples: bool,
        return_variants: bool,
        validated_options: _ValidatedReadOptions,
        variant_window: dict[str, int] | None,
    ) -> Any:
        members = {key: str(path) for key, path in self.source.members.items()}
        options = {
            "samples": None if samples is None else list(samples),
            "variants": validated_options.variant_filter_ir,
            "variant_window": variant_window,
            "return_samples": return_samples,
            "return_variants": return_variants,
        }
        if validated_options.sparse_format is None:
            rust_result = (
                self._read_dense_from_rust(members, options)
                if kind == "geno"
                else self._read_haplotypes_dense_from_rust(members, options)
            )
            genotype_matrix = dense_array_from_rust(
                values=rust_result["values"],
                shape=tuple(rust_result["shape"]),
                missing_mask=rust_result["missing_mask"],
                missing=validated_options.missing,
                dtype=validated_options.dtype,
            )
        else:
            rust_result = (
                self._read_sparse_from_rust(members, options)
                if kind == "geno"
                else self._read_haplotypes_sparse_from_rust(members, options)
            )
            genotype_matrix = sparse_matrix_from_rust(
                indptr=rust_result["indptr"],
                indices=rust_result["indices"],
                data=rust_result["data"],
                shape=tuple(rust_result["shape"]),
                dtype=validated_options.dtype,
                sparse_format=validated_options.sparse_format,
            )
        sample_metadata = samples_frame(rust_result["samples"]) if return_samples else None
        variant_metadata = variants_frame(rust_result["variants"]) if return_variants else None
        return read_result_tuple(
            genotype_matrix,
            sample_metadata,
            variant_metadata,
            return_samples=return_samples,
            return_variants=return_variants,
        )

    def _read_dense_from_rust(self, members: dict[str, str], options: dict[str, Any]) -> dict[str, Any]:
        try:
            return _rust.read_dense(self.source.format.value, members, options)
        except ValueError as error:
            raise _public_read_error(error) from error

    def _read_haplotypes_dense_from_rust(
        self,
        members: dict[str, str],
        options: dict[str, Any],
    ) -> dict[str, Any]:
        try:
            return _rust.read_haplotypes_dense(self.source.format.value, members, options)
        except ValueError as error:
            raise _public_haplotype_read_error(error) from error

    def _read_haplotypes_sparse_from_rust(
        self,
        members: dict[str, str],
        options: dict[str, Any],
    ) -> dict[str, Any]:
        try:
            return _rust.read_haplotypes_sparse(self.source.format.value, members, options)
        except ValueError as error:
            raise _public_haplotype_read_error(error) from error

    def _read_sparse_from_rust(self, members: dict[str, str], options: dict[str, Any]) -> dict[str, Any]:
        try:
            return _rust.read_sparse(self.source.format.value, members, options)
        except ValueError as error:
            raise _public_read_error(error) from error

    def _metadata(self) -> dict[str, Any]:
        if self._metadata_cache is None:
            members = {key: str(path) for key, path in self.source.members.items()}
            try:
                metadata = _rust.read_metadata(self.source.format.value, members)
            except ValueError as error:
                from ._errors import InvalidSourceError

                raise InvalidSourceError(str(error)) from error
            object.__setattr__(self, "_metadata_cache", metadata)
        return self._metadata_cache


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
    missing: Any = _DEFAULT_MISSING,
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


def blocks(path: str | Path, size: int, *, format: str | None = None, **read_options: Any) -> Any:
    return open(path, format=format).blocks(size, **read_options)


def _validate_read_options(
    *,
    kind: str,
    sparse: bool | str,
    variants: Any,
    samples: list[str] | tuple[str, ...] | set[str] | None,
    missing: Any,
    dtype: Any,
    return_samples: bool,
    return_variants: bool,
) -> _ValidatedReadOptions:
    _validate_kind(kind)
    sparse_format = _validate_sparse(sparse)
    normalized_missing = _normalize_missing(missing, sparse_format=sparse_format)
    variant_filter_ir = _validate_variant_filter(variants)
    normalized_dtype = _normalize_dtype(dtype)
    _validate_sparse_missing_compatibility(sparse_format, normalized_missing)
    _validate_missing_dtype_compatibility(normalized_missing, normalized_dtype)
    _validate_sample_filter(samples)
    _validate_bool_option("return_samples", return_samples)
    _validate_bool_option("return_variants", return_variants)
    return _ValidatedReadOptions(
        dtype=normalized_dtype,
        missing=normalized_missing,
        sparse_format=sparse_format,
        variant_filter_ir=variant_filter_ir,
    )


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


def _reject_deferred_plink2_decode(source_format: str) -> None:
    if source_format == "plink2":
        raise UnsupportedFormatError("PLINK2 decode is deferred pending docs/plink2-parser-strategy.md")


def _unsupported_haplotype_source(source_format: str) -> UnsupportedRepresentation:
    if source_format == "vcf":
        return UnsupportedRepresentation(
            'VCF source has no phased GT evidence; kind="haplo" requires phased retained variants'
        )
    return UnsupportedRepresentation(f"{source_format} does not support haplo reads")


def _validate_sparse(sparse: bool | str) -> str | None:
    if sparse is False:
        return None
    if sparse is True:
        return "csc"
    if isinstance(sparse, str) and sparse in {"csc", "csr"}:
        return sparse
    raise InvalidOptionError(f"unsupported sparse option: {sparse!r}")


def _validate_variant_filter(variants: Any) -> dict[str, Any] | None:
    if variants is None:
        return None
    return _variant_filter_ir(variants)


def _variant_filter_ir(variants: Any) -> dict[str, Any] | None:
    if variants is None:
        return None
    if isinstance(variants, FilterExpr):
        return variants.to_ir()
    if callable(variants):
        raise InvalidOptionError("variants must be a serializable filter expression or variant ID iterable")
    if isinstance(variants, str) or not isinstance(variants, Iterable):
        raise InvalidOptionError("variants must be a serializable filter expression or variant ID iterable")
    try:
        return id_in(list(variants)).to_ir()
    except InvalidOptionError:
        raise
    except TypeError as error:
        raise InvalidOptionError("variants must be a serializable filter expression or variant ID iterable") from error


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


def _normalize_missing(missing: Any, *, sparse_format: str | None) -> str:
    if missing is _DEFAULT_MISSING:
        return "raise" if sparse_format is not None else "nan"
    _validate_missing(missing)
    return missing


def _normalize_dtype(dtype: Any) -> np.dtype[Any]:
    try:
        return np.dtype(dtype)
    except TypeError as error:
        raise InvalidOptionError(f"invalid dtype: {dtype!r}") from error


def _validate_missing_dtype_compatibility(missing: str, dtype: np.dtype[Any]) -> None:
    if isinstance(missing, str) and missing in {"nan", "impute"} and not np.issubdtype(dtype, np.floating):
        raise InvalidOptionError(f'missing="{missing}" requires a floating dtype')


def _validate_sparse_missing_compatibility(sparse_format: str | None, missing: str) -> None:
    if sparse_format is not None and missing in {"nan", "impute"}:
        raise InvalidOptionError("this release does not store sparse missing values; use missing='raise'")


def _public_read_error(error: ValueError) -> Exception:
    message = str(error)
    if "missing requested sample" in message:
        return SampleFilterError(message)
    if "sparse missing values" in message:
        return MissingDataError(message)
    return InvalidSourceError(message)


def _public_haplotype_read_error(error: ValueError) -> Exception:
    message = str(error)
    if "unphased" in message or "unsupported haplotype format" in message:
        return UnsupportedRepresentation(message)
    if "sparse missing values" in message:
        return MissingDataError(message)
    return _public_read_error(error)


def _validate_bool_option(name: str, value: bool) -> None:
    if not isinstance(value, bool):
        raise InvalidOptionError(f"{name} must be a bool")


def _validate_block_size(size: int) -> None:
    if not isinstance(size, int) or isinstance(size, bool) or size < 1:
        raise InvalidOptionError("block size must be a positive integer")


def _read_options_with_defaults(read_options: dict[str, Any]) -> dict[str, Any]:
    defaults = {
        "kind": "geno",
        "sparse": False,
        "variants": None,
        "samples": None,
        "missing": _DEFAULT_MISSING,
        "dtype": "float32",
        "return_samples": False,
        "return_variants": False,
    }
    unknown = set(read_options) - set(defaults)
    if unknown:
        keys = ", ".join(sorted(unknown))
        raise InvalidOptionError(f"unsupported option(s): {keys}")
    return defaults | read_options
