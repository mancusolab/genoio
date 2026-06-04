# pattern: Functional Core

from __future__ import annotations

import math
import re
from collections.abc import Iterable
from dataclasses import dataclass

from ._errors import InvalidOptionError

_REGION_PATTERN = re.compile(r"^[^:\s]+:[0-9]+-[0-9]+$")
_GENOTYPE_RATE_RANGE = (0.0, 1.0)

JsonScalar = str | int | float | bool | None
JsonValue = JsonScalar | list["JsonValue"] | dict[str, "JsonValue"]
ParamValue = JsonScalar | tuple[JsonScalar, ...]


class FilterExpr:
    r"""Serializable variant filter expression.

    Filter expressions can be combined with ``&`` and ``|`` and negated with
    ``~``. They are converted to a JSON-compatible intermediate representation
    before being evaluated in Rust. Expressions are declarative: Python
    callbacks are not supported, and a valid expression may retain zero
    variants.
    """

    def __and__(self, other: FilterExpr) -> FilterExpr:
        return AndExpr(self, _ensure_expression(other))

    def __or__(self, other: FilterExpr) -> FilterExpr:
        return OrExpr(self, _ensure_expression(other))

    def __invert__(self) -> FilterExpr:
        return NotExpr(self)

    def to_ir(self) -> dict[str, JsonValue]:
        r"""Return the JSON-compatible filter representation consumed by Rust.

        Treat this representation as an inspection aid. Rust validates and may
        simplify the expression before readers evaluate it.

        **Returns:**

        Dictionary containing only JSON-compatible values.
        """
        raise NotImplementedError


@dataclass(frozen=True)
class PredicateExpr(FilterExpr):
    r"""Leaf filter predicate with validated parameter values.

    **Arguments:**

    - `name`: predicate name understood by the Rust filter evaluator.
    - `params`: validated parameter key/value pairs.
    """

    name: str
    params: tuple[tuple[str, ParamValue], ...] = ()

    def to_ir(self) -> dict[str, JsonValue]:
        r"""Return this predicate as filter IR."""
        return {
            "op": "predicate",
            "name": self.name,
            "params": {key: _json_param(value) for key, value in self.params},
        }


@dataclass(frozen=True)
class AndExpr(FilterExpr):
    r"""Boolean conjunction of two filter expressions."""

    left: FilterExpr
    right: FilterExpr

    def to_ir(self) -> dict[str, JsonValue]:
        r"""Return this conjunction as filter IR."""
        return {"op": "and", "left": self.left.to_ir(), "right": self.right.to_ir()}


@dataclass(frozen=True)
class OrExpr(FilterExpr):
    r"""Boolean disjunction of two filter expressions."""

    left: FilterExpr
    right: FilterExpr

    def to_ir(self) -> dict[str, JsonValue]:
        r"""Return this disjunction as filter IR."""
        return {"op": "or", "left": self.left.to_ir(), "right": self.right.to_ir()}


@dataclass(frozen=True)
class NotExpr(FilterExpr):
    r"""Boolean negation of one filter expression."""

    expr: FilterExpr

    def to_ir(self) -> dict[str, JsonValue]:
        r"""Return this negation as filter IR."""
        return {"op": "not", "expr": self.expr.to_ir()}


def chrom(value: str) -> FilterExpr:
    r"""Keep variants on chromosome `value`.

    **Arguments:**

    - `value`: non-empty chromosome label.

    **Returns:**

    `FilterExpr` predicate.
    """
    if not isinstance(value, str) or not value:
        raise InvalidOptionError("chrom filter requires a non-empty chromosome string")
    return PredicateExpr("chrom", (("value", value),))


def region(value: str) -> FilterExpr:
    r"""Keep variants in a 1-based inclusive `"chrom:start-end"` region.

    **Arguments:**

    - `value`: region string with 1-based inclusive coordinates.

    **Returns:**

    `FilterExpr` predicate.
    """
    _validate_region(value)
    return PredicateExpr("region", (("value", value),))


def snp() -> FilterExpr:
    r"""Keep single-nucleotide variants with one REF and one ALT base.

    **Returns:**

    `FilterExpr` predicate.
    """
    return PredicateExpr("snp")


def biallelic() -> FilterExpr:
    r"""Keep variants with exactly one ALT allele.

    **Returns:**

    `FilterExpr` predicate.
    """
    return PredicateExpr("biallelic")


def maf(*, min: float | None = None, max: float | None = None) -> FilterExpr:
    r"""Keep variants by minor allele frequency over called genotypes.

    **Arguments:**

    - `min`: optional inclusive lower MAF threshold.
    - `max`: optional inclusive upper MAF threshold.

    **Returns:**

    `FilterExpr` predicate.
    """
    return PredicateExpr("maf", _validate_float_range("maf", min=min, max=max))


def qual(*, min: float | None = None, max: float | None = None) -> FilterExpr:
    r"""Keep variants by source QUAL value.

    **Arguments:**

    - `min`: optional inclusive lower QUAL threshold.
    - `max`: optional inclusive upper QUAL threshold.

    **Returns:**

    `FilterExpr` predicate.
    """
    return PredicateExpr("qual", _validate_nonnegative_float_range("qual", min=min, max=max))


def mac(*, min: int | None = None, max: int | None = None) -> FilterExpr:
    r"""Keep variants by minor allele count over called genotypes.

    **Arguments:**

    - `min`: optional inclusive lower MAC threshold.
    - `max`: optional inclusive upper MAC threshold.

    **Returns:**

    `FilterExpr` predicate.
    """
    return PredicateExpr("mac", _validate_int_range("mac", min=min, max=max))


def missing_rate(max: float) -> FilterExpr:
    r"""Keep variants with missing-call rate at or below `max`.

    **Arguments:**

    - `max`: inclusive upper missing-call-rate threshold.

    **Returns:**

    `FilterExpr` predicate.
    """
    value = _validate_rate("missing_rate max", max)
    return PredicateExpr("missing_rate", (("max", value),))


def polymorphic() -> FilterExpr:
    r"""Keep variants with nonzero minor allele count.

    **Returns:**

    `FilterExpr` predicate.
    """
    return PredicateExpr("polymorphic")


def id_in(values: Iterable[str]) -> FilterExpr:
    r"""Keep variants whose source ID is in `values`.

    Sets are accepted for convenience and are sorted before serialization so
    the generated filter IR is deterministic.

    **Arguments:**

    - `values`: iterable of unique variant ID strings.

    **Returns:**

    `FilterExpr` predicate.
    """
    if isinstance(values, str) or not isinstance(values, Iterable):
        raise InvalidOptionError("id_in requires an iterable of variant ID strings")
    normalized = tuple(sorted(values) if isinstance(values, set) else values)
    if any(not isinstance(value, str) for value in normalized):
        raise InvalidOptionError("id_in values must contain only variant ID strings")
    if len(normalized) != len(set(normalized)):
        raise InvalidOptionError("id_in values must not contain duplicate variant IDs")
    return PredicateExpr("id_in", (("values", normalized),))


def _validate_region(value: str) -> None:
    if not isinstance(value, str) or not _REGION_PATTERN.fullmatch(value):
        raise InvalidOptionError(f"invalid region syntax: {value!r}; expected 'chrom:start-end'")

    _, coordinates = value.split(":", 1)
    start_text, end_text = coordinates.split("-", 1)
    start = int(start_text)
    end = int(end_text)
    if start < 1 or end < start:
        raise InvalidOptionError(f"invalid region coordinates: {value!r}; expected 1-based start <= end")


def _validate_float_range(name: str, *, min: float | None, max: float | None) -> tuple[tuple[str, ParamValue], ...]:
    if min is None and max is None:
        raise InvalidOptionError(f"{name} requires at least one threshold")
    min_value = None if min is None else _validate_rate(f"{name} min", min)
    max_value = None if max is None else _validate_rate(f"{name} max", max)
    if min_value is not None and max_value is not None and min_value > max_value:
        raise InvalidOptionError(f"{name} min must be <= max")
    return tuple((key, value) for key, value in (("min", min_value), ("max", max_value)) if value is not None)


def _validate_int_range(name: str, *, min: int | None, max: int | None) -> tuple[tuple[str, ParamValue], ...]:
    if min is None and max is None:
        raise InvalidOptionError(f"{name} requires at least one threshold")
    min_value = None if min is None else _validate_nonnegative_int(f"{name} min", min)
    max_value = None if max is None else _validate_nonnegative_int(f"{name} max", max)
    if min_value is not None and max_value is not None and min_value > max_value:
        raise InvalidOptionError(f"{name} min must be <= max")
    return tuple((key, value) for key, value in (("min", min_value), ("max", max_value)) if value is not None)


def _validate_nonnegative_float_range(
    name: str, *, min: float | None, max: float | None
) -> tuple[tuple[str, ParamValue], ...]:
    if min is None and max is None:
        raise InvalidOptionError(f"{name} requires at least one threshold")
    min_value = None if min is None else _validate_nonnegative_float(f"{name} min", min)
    max_value = None if max is None else _validate_nonnegative_float(f"{name} max", max)
    if min_value is not None and max_value is not None and min_value > max_value:
        raise InvalidOptionError(f"{name} min must be <= max")
    return tuple((key, value) for key, value in (("min", min_value), ("max", max_value)) if value is not None)


def _validate_rate(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise InvalidOptionError(f"{name} must be a number between 0 and 1")
    normalized = float(value)
    lower, upper = _GENOTYPE_RATE_RANGE
    if not math.isfinite(normalized) or normalized < lower or normalized > upper:
        raise InvalidOptionError(f"{name} must be between 0 and 1")
    return normalized


def _validate_nonnegative_float(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise InvalidOptionError(f"{name} must be a non-negative number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0.0:
        raise InvalidOptionError(f"{name} must be a non-negative number")
    return normalized


def _validate_nonnegative_int(name: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidOptionError(f"{name} must be a non-negative integer")
    if value < 0:
        raise InvalidOptionError(f"{name} must be a non-negative integer")
    return value


def _ensure_expression(value: FilterExpr) -> FilterExpr:
    if not isinstance(value, FilterExpr):
        raise TypeError(f"expected FilterExpr, got {type(value).__name__}")
    return value


def _json_param(value: ParamValue) -> JsonValue:
    if isinstance(value, tuple):
        return list(value)
    return value
