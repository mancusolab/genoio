# pattern: Mixed (unavoidable)
# Reason: Phase 1 explicitly defines source resolution as one cohesive filesystem shell module.

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from ._errors import (
    AmbiguousSourceError,
    InvalidSourceError,
    MissingCompanionFileError,
    UnsupportedFormatError,
)


class SourceFormat(Enum):
    VCF = "vcf"
    BCF = "bcf"
    PLINK1 = "plink1"
    PLINK2 = "plink2"


@dataclass(frozen=True)
class ResolvedSource:
    format: SourceFormat
    path: Path
    members: Mapping[str, Path]
    prefix: Path | None = None


_PLINK1_SUFFIXES = {"bed": ".bed", "bim": ".bim", "fam": ".fam"}
_PLINK2_SUFFIXES = {"pgen": ".pgen", "pvar": ".pvar", "psam": ".psam"}
_MEMBER_SUFFIX_TO_FORMAT = {
    **{suffix: SourceFormat.PLINK1 for suffix in _PLINK1_SUFFIXES.values()},
    **{suffix: SourceFormat.PLINK2 for suffix in _PLINK2_SUFFIXES.values()},
}


def resolve_source(path: str | Path, format: str | SourceFormat | None = None) -> ResolvedSource:
    source_path = Path(path)
    requested_format = _normalize_format(format)

    if requested_format in {SourceFormat.PLINK1, SourceFormat.PLINK2}:
        return _resolve_plink(source_path, requested_format)

    if requested_format in {SourceFormat.VCF, SourceFormat.BCF}:
        return _resolve_single_file(source_path, requested_format)

    detected_format = _detect_single_file_format(source_path)
    if detected_format is not None:
        return _resolve_single_file(source_path, detected_format)

    member_format = _MEMBER_SUFFIX_TO_FORMAT.get(source_path.suffix)
    if member_format is not None:
        return _resolve_plink(source_path, member_format)

    if source_path.suffix:
        raise UnsupportedFormatError(f"unsupported source extension: {source_path}")

    return _resolve_prefix(source_path)


def _normalize_format(format: str | SourceFormat | None) -> SourceFormat | None:
    if format is None or isinstance(format, SourceFormat):
        return format
    normalized = format.lower()
    aliases = {
        "vcf": SourceFormat.VCF,
        "bcf": SourceFormat.BCF,
        "plink1": SourceFormat.PLINK1,
        "bed": SourceFormat.PLINK1,
        "plink2": SourceFormat.PLINK2,
        "pgen": SourceFormat.PLINK2,
    }
    try:
        return aliases[normalized]
    except KeyError as error:
        raise UnsupportedFormatError(f"unsupported format: {format}") from error


def _detect_single_file_format(path: Path) -> SourceFormat | None:
    if path.name.endswith(".vcf.gz") or path.suffix == ".vcf":
        return SourceFormat.VCF
    if path.suffix == ".bcf":
        return SourceFormat.BCF
    return None


def _resolve_single_file(path: Path, format: SourceFormat) -> ResolvedSource:
    if not path.exists():
        raise InvalidSourceError(f"source path does not exist: {path}")
    if not path.is_file():
        raise InvalidSourceError(f"source path is not a file: {path}")

    expected_format = _detect_single_file_format(path)
    if expected_format is None:
        raise UnsupportedFormatError(f"unsupported source extension: {path}")
    if expected_format is not format:
        raise UnsupportedFormatError(f"source path {path} is not {format.value}")

    member_key = "vcf" if format is SourceFormat.VCF else "bcf"
    return ResolvedSource(format=format, path=path, members={member_key: path}, prefix=None)


def _resolve_prefix(path: Path) -> ResolvedSource:
    plink1_members = _existing_members(path, _PLINK1_SUFFIXES)
    plink2_members = _existing_members(path, _PLINK2_SUFFIXES)
    plink1_complete = set(plink1_members) == set(_PLINK1_SUFFIXES)
    plink2_complete = set(plink2_members) == set(_PLINK2_SUFFIXES)

    if plink1_complete and plink2_complete:
        raise AmbiguousSourceError(f"source prefix matches multiple formats: {path}")
    if plink1_members:
        return _resolve_plink_prefix(path, SourceFormat.PLINK1)
    if plink2_members:
        return _resolve_plink_prefix(path, SourceFormat.PLINK2)
    if path.exists():
        raise UnsupportedFormatError(f"unsupported source extension: {path}")
    raise InvalidSourceError(f"source path does not exist: {path}")


def _resolve_plink(path: Path, format: SourceFormat) -> ResolvedSource:
    prefix = _prefix_for_plink_path(path, format)
    return _resolve_plink_prefix(prefix, format)


def _prefix_for_plink_path(path: Path, format: SourceFormat) -> Path:
    suffixes = _PLINK1_SUFFIXES if format is SourceFormat.PLINK1 else _PLINK2_SUFFIXES
    if path.suffix in suffixes.values():
        if not path.exists():
            raise InvalidSourceError(f"source path does not exist: {path}")
        return path.with_suffix("")
    return path


def _resolve_plink_prefix(prefix: Path, format: SourceFormat) -> ResolvedSource:
    suffixes = _PLINK1_SUFFIXES if format is SourceFormat.PLINK1 else _PLINK2_SUFFIXES
    members = _existing_members(prefix, suffixes)
    missing = sorted(set(suffixes) - set(members))
    if missing:
        missing_paths = ", ".join(str(prefix.with_suffix(suffixes[key])) for key in missing)
        raise MissingCompanionFileError(f"missing companion file(s): {missing_paths}")
    primary_key = "bed" if format is SourceFormat.PLINK1 else "pgen"
    return ResolvedSource(format=format, path=members[primary_key], members=members, prefix=prefix)


def _existing_members(prefix: Path, suffixes: Mapping[str, str]) -> dict[str, Path]:
    members: dict[str, Path] = {}
    for key, suffix in suffixes.items():
        member_path = prefix.with_suffix(suffix)
        if member_path.exists():
            if not member_path.is_file():
                raise InvalidSourceError(f"source member is not a file: {member_path}")
            members[key] = member_path
    return members
