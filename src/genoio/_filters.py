# pattern: Functional Core

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from ._errors import InvalidOptionError

_REGION_PATTERN = re.compile(r"^[^:\s]+:[0-9]+-[0-9]+$")


@dataclass(frozen=True)
class Expression:
    op: str
    args: tuple[Expression, ...] = ()
    value: Any = None
    options: tuple[tuple[str, Any], ...] = ()

    def __and__(self, other: Expression) -> Expression:
        return Expression("and", args=(self, _ensure_expression(other)))

    def __or__(self, other: Expression) -> Expression:
        return Expression("or", args=(self, _ensure_expression(other)))

    def __invert__(self) -> Expression:
        return Expression("not", args=(self,))

    def to_ir(self) -> dict[str, Any]:
        if self.op == "not":
            return {"op": self.op, "arg": self.args[0].to_ir()}
        if self.args:
            return {"op": self.op, "args": [arg.to_ir() for arg in self.args]}
        if self.options:
            return {"op": self.op, **dict(self.options)}
        return {"op": self.op, "value": self.value}


def chrom(value: str) -> Expression:
    return Expression("chrom", value=value)


def region(value: str) -> Expression:
    _validate_region(value)
    return Expression("region", value=value)


def snp() -> Expression:
    return Expression("snp")


def biallelic() -> Expression:
    return Expression("biallelic")


def maf(*, min: float | None = None, max: float | None = None) -> Expression:
    return _range_expression("maf", min=min, max=max)


def mac(*, min: int | None = None, max: int | None = None) -> Expression:
    return _range_expression("mac", min=min, max=max)


def missing_rate(*, min: float | None = None, max: float | None = None) -> Expression:
    return _range_expression("missing_rate", min=min, max=max)


def polymorphic() -> Expression:
    return Expression("polymorphic")


def id_in(values: list[str] | tuple[str, ...] | set[str]) -> Expression:
    return Expression("id_in", value=tuple(values))


def _range_expression(op: str, *, min: float | int | None, max: float | int | None) -> Expression:
    options = tuple((key, value) for key, value in (("min", min), ("max", max)) if value is not None)
    return Expression(op, options=options)


def _validate_region(value: str) -> None:
    if not isinstance(value, str) or not _REGION_PATTERN.fullmatch(value):
        raise InvalidOptionError(f"invalid region syntax: {value!r}; expected 'chrom:start-end'")

    _, coordinates = value.split(":", 1)
    start_text, end_text = coordinates.split("-", 1)
    start = int(start_text)
    end = int(end_text)
    if start < 1 or end < start:
        raise InvalidOptionError(f"invalid region coordinates: {value!r}; expected 1-based start <= end")


def _ensure_expression(value: Expression) -> Expression:
    if not isinstance(value, Expression):
        raise TypeError(f"expected Expression, got {type(value).__name__}")
    return value
