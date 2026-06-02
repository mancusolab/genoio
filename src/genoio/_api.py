# pattern: Mixed (needs refactoring)
# Reason: Public entrypoints call resolve_source(), which performs filesystem checks;
# Dataset validation helpers remain pure.

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ._errors import InvalidOptionError, UnsupportedRepresentation
from ._source import ResolvedSource, resolve_source

_SUPPORTED_REPRESENTATIONS = {None, "numpy", "polars", "scipy", "csr", "csc"}


@dataclass(frozen=True)
class Dataset:
    source: ResolvedSource

    def read(self, *, representation: str | None = None, **options: Any) -> Any:
        _validate_representation(representation)
        _reject_options(options)
        raise NotImplementedError("genotype reading is implemented in a later phase")

    def samples(self, **options: Any) -> Any:
        _reject_options(options)
        raise NotImplementedError("sample metadata reading is implemented in a later phase")

    def variants(self, **options: Any) -> Any:
        _reject_options(options)
        raise NotImplementedError("variant metadata reading is implemented in a later phase")

    def blocks(self, *, representation: str | None = None, **options: Any) -> Any:
        _validate_representation(representation)
        _reject_options(options)
        raise NotImplementedError("block iteration is implemented in a later phase")


def open(path: str | Path, format: str | None = None) -> Dataset:
    return Dataset(source=resolve_source(path, format=format))


def read(path: str | Path, *, format: str | None = None, representation: str | None = None, **options: Any) -> Any:
    return open(path, format=format).read(representation=representation, **options)


def samples(path: str | Path, *, format: str | None = None, **options: Any) -> Any:
    return open(path, format=format).samples(**options)


def variants(path: str | Path, *, format: str | None = None, **options: Any) -> Any:
    return open(path, format=format).variants(**options)


def _validate_representation(representation: str | None) -> None:
    if representation not in _SUPPORTED_REPRESENTATIONS:
        raise UnsupportedRepresentation(f"unsupported representation: {representation}")


def _reject_options(options: dict[str, Any]) -> None:
    if options:
        keys = ", ".join(sorted(options))
        raise InvalidOptionError(f"unsupported option(s): {keys}")
