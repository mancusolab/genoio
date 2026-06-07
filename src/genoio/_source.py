# pattern: Mixed (unavoidable)
# Reason: Phase 1 explicitly defines source resolution as one cohesive filesystem shell module.

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from types import MappingProxyType

from ._errors import InvalidSourceError, MissingCompanionFileError, UnsupportedFormatError


class SourceFormat(Enum):
    r"""Supported on-disk genotype source formats."""

    VCF = "vcf"
    BCF = "bcf"
    BGEN = "bgen"
    PLINK1 = "plink1"
    PLINK2 = "plink2"


@dataclass(frozen=True, slots=True)
class ResolvedSource:
    r"""Resolved source path and required member files.

    `members` maps logical member names such as `"vcf"`, `"bed"`, or `"pgen"`
    to concrete paths. PLINK sources also carry the shared prefix.

    **Attributes:**

    - `format`: detected or requested source format.
    - `path`: primary input path used by the backend reader.
    - `members`: required source members keyed by logical role.
    - `prefix`: PLINK prefix, or `None` for single-file sources.
    """

    format: SourceFormat
    path: Path
    members: Mapping[str, Path]
    prefix: Path | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "members", MappingProxyType(dict(self.members)))


_PLINK1_SUFFIXES = {"bed": ".bed", "bim": ".bim", "fam": ".fam"}
_PLINK2_SUFFIXES = {"pgen": ".pgen", "pvar": ".pvar", "psam": ".psam"}
_PLINK_SUFFIXES = {
    SourceFormat.PLINK1: _PLINK1_SUFFIXES,
    SourceFormat.PLINK2: _PLINK2_SUFFIXES,
}
_PLINK2_COMPRESSED_PVAR_SUFFIX = ".pvar.zst"
_BGEN_SUFFIX = ".bgen"
_BGEN_SAMPLE_SUFFIX = ".sample"


def resolve_vcf(path: str | Path) -> ResolvedSource:
    r"""Resolve a VCF or BCF source file.

    **Arguments:**

    - `path`: `.vcf`, `.vcf.gz`, or `.bcf` path.

    **Returns:**

    `ResolvedSource` with one VCF/BCF member.

    **Raises:**

    - `genoio.SourceResolutionError`: if the path is missing, is not a file, or
      does not have a supported VCF/BCF extension.
    """
    source_path = Path(path)
    detected_format = _detect_single_file_format(source_path)
    if detected_format is None:
        raise UnsupportedFormatError(f"unsupported source extension: {source_path}")
    return _resolve_single_file(source_path, detected_format)


def resolve_bfile(path: str | Path) -> ResolvedSource:
    r"""Resolve a PLINK1 BED/BIM/FAM file set from a prefix or member path.

    **Arguments:**

    - `path`: PLINK1 prefix or one `.bed`, `.bim`, or `.fam` member path.

    **Returns:**

    `ResolvedSource` with `.bed`, `.bim`, and `.fam` members.
    """
    return _resolve_plink(Path(path), SourceFormat.PLINK1)


def resolve_bgen(path: str | Path) -> ResolvedSource:
    r"""Resolve a BGEN source file from a prefix or `.bgen` member path.

    **Arguments:**

    - `path`: BGEN prefix or `.bgen` member path.

    **Returns:**

    `ResolvedSource` with the required `.bgen` member and optional same-prefix
    `.sample` member.
    """
    source_path = Path(path)
    if source_path.suffix and source_path.suffix != _BGEN_SUFFIX and source_path.exists():
        raise UnsupportedFormatError(f"source path {source_path} is not bgen")
    prefix = source_path.with_suffix("") if source_path.suffix == _BGEN_SUFFIX else source_path
    return _resolve_bgen_prefix(prefix)


def resolve_pfile(path: str | Path) -> ResolvedSource:
    r"""Resolve a PLINK2 PGEN/PVAR/PSAM file set from a prefix or member path.

    **Arguments:**

    - `path`: PLINK2 prefix or one `.pgen`, `.pvar`, `.pvar.zst`, or `.psam`
      member path.

    **Returns:**

    `ResolvedSource` with `.pgen`, `.psam`, and either `.pvar` or `.pvar.zst`
    as the logical `"pvar"` member.
    """
    return _resolve_plink(Path(path), SourceFormat.PLINK2)


def _resolve_bgen_prefix(prefix: Path) -> ResolvedSource:
    bgen_path = _append_suffix(prefix, _BGEN_SUFFIX)
    if not bgen_path.exists():
        raise InvalidSourceError(f"source path does not exist: {bgen_path}")
    if not bgen_path.is_file():
        raise InvalidSourceError(f"source path is not a file: {bgen_path}")

    members = {"bgen": bgen_path}
    sample_path = _append_suffix(prefix, _BGEN_SAMPLE_SUFFIX)
    if sample_path.exists():
        if not sample_path.is_file():
            raise InvalidSourceError(f"source member is not a file: {sample_path}")
        members["sample"] = sample_path
    return ResolvedSource(format=SourceFormat.BGEN, path=bgen_path, members=members, prefix=prefix)


def _append_suffix(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


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


def _resolve_plink(path: Path, format: SourceFormat) -> ResolvedSource:
    suffixes = _plink_suffixes(format)
    if path.suffix and path.suffix not in suffixes.values() and not _is_compressed_pvar(path, format):
        raise UnsupportedFormatError(f"source path {path} is not {format.value}")
    prefix = _prefix_for_plink_path(path, format)
    return _resolve_plink_prefix(prefix, format)


def _prefix_for_plink_path(path: Path, format: SourceFormat) -> Path:
    suffixes = _plink_suffixes(format)
    if _is_compressed_pvar(path, format):
        if not path.exists():
            raise InvalidSourceError(f"source path does not exist: {path}")
        return path.with_suffix("").with_suffix("")
    if path.suffix in suffixes.values():
        if not path.exists():
            raise InvalidSourceError(f"source path does not exist: {path}")
        return path.with_suffix("")
    return path


def _resolve_plink_prefix(prefix: Path, format: SourceFormat) -> ResolvedSource:
    suffixes = _plink_suffixes(format)
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
        if key == "pvar" and not member_path.exists():
            compressed_member = Path(f"{member_path}.zst")
            if compressed_member.exists():
                member_path = compressed_member
        if member_path.exists():
            if not member_path.is_file():
                raise InvalidSourceError(f"source member is not a file: {member_path}")
            members[key] = member_path
    return members


def _plink_suffixes(format: SourceFormat) -> Mapping[str, str]:
    return _PLINK_SUFFIXES[format]


def _is_compressed_pvar(path: Path, format: SourceFormat) -> bool:
    return format is SourceFormat.PLINK2 and path.name.endswith(_PLINK2_COMPRESSED_PVAR_SUFFIX)
