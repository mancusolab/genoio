# pattern: Mixed
# Reason: Public dataset constructors perform filesystem source resolution;
# Dataset validation helpers remain pure.

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import DTypeLike

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
    UnsupportedRepresentation,
)
from ._filters import FilterExpr, id_in
from ._source import ResolvedSource, resolve_bfile, resolve_bgen, resolve_pfile, resolve_vcf

_SUPPORTED_KINDS = {"geno", "haplo"}
_SUPPORTED_DOSAGE_SOURCES = {"hardcall", "dosage"}
_SUPPORTED_MISSING_POLICIES = {"nan", "raise", "impute"}


class _DefaultMissing:
    def __repr__(self) -> str:
        return "DEFAULT_MISSING"


_DEFAULT_MISSING = _DefaultMissing()


@dataclass(frozen=True)
class _ReadOptions:
    kind: str
    dosage: str
    sparse: bool | str
    variants: Any
    samples: list[str] | tuple[str, ...] | set[str] | None
    missing: Any
    dtype: DTypeLike
    return_samples: bool
    return_variants: bool


@dataclass(frozen=True)
class _ValidatedReadOptions:
    dtype: np.dtype[Any]
    missing: str
    sparse_format: str | None
    variant_filter_ir: dict[str, Any] | None


@dataclass(frozen=True)
class Dataset:
    r"""Resolved genotype dataset with metadata, whole-read, and block-read methods.

    Constructed by [`genoio.vcf`][], [`genoio.bfile`][], [`genoio.bgen`][], or
    [`genoio.pfile`][]. The object caches source metadata after the first
    metadata-dependent operation, but matrix reads are executed on each call.

    **Attributes:**

    - `source`: resolved source format, primary path, companion files, and
      optional PLINK prefix.
    """

    source: ResolvedSource
    _metadata_cache: dict[str, Any] | None = field(default=None, init=False, compare=False, repr=False)
    _samples_frame_cache: Any | None = field(default=None, init=False, compare=False, repr=False)
    _variants_frame_cache: Any | None = field(default=None, init=False, compare=False, repr=False)

    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: bool | str = False,
        variants: Any = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Any = _DEFAULT_MISSING,
        dtype: DTypeLike = "float32",
        return_samples: bool = False,
        return_variants: bool = False,
    ) -> Any:
        r"""Read a genotype or haplotype matrix from this dataset.

        Dense reads return a NumPy array with shape `(samples, variants)`.
        Sparse reads return SciPy CSC by default or CSR when `sparse="csr"`.
        Set `return_samples` or `return_variants` to return Polars metadata
        frames with the matrix.

        **Arguments:**

        - `kind`: `"geno"` for diploid sample-by-variant genotype values or
          `"haplo"` for phased haplotype rows. Haplotype reads currently require
          phased VCF.
        - `dosage`: `"hardcall"` reads allele counts from hard calls.
          `"dosage"` reads dosage-backed genotype values when the source
          supports them. This release supports dense VCF `FORMAT/DS` and
          PLINK2 unphased biallelic dosage reads, and dense BGEN Layout 2
          biallelic diploid dosage reads. Phased BGEN records are collapsed to
          expected diploid A1 dosage. BGEN sample IDs must be embedded in the
          `.bgen` file or supplied by a companion `.sample` file. Concrete
          BGEN region filters use a same-path `.bgen.bgi` index when present.
          Haplotype and sparse reads only support `"hardcall"`.
        - `sparse`: `False` for dense NumPy, `True` or `"csc"` for CSC,
          `"csr"` for CSR.
        - `variants`: filter expression from `genoio` or iterable of variant
          IDs. `None` keeps all variants.
        - `samples`: optional sample ID keep list. Retained rows stay in source
          order.
        - `missing`: `"nan"`, `"raise"`, or `"impute"`. The default is `"nan"`
          for dense reads and `"raise"` for sparse reads.
        - `dtype`: NumPy dtype for the returned matrix.
        - `return_samples`: include a sample metadata frame.
        - `return_variants`: include a variant metadata frame.

        **Returns:**

        Matrix alone, or a tuple containing the matrix and requested metadata
        frames.

        **Raises:**

        - `genoio.InvalidOptionError`: if read options are invalid.
        - `genoio.UnsupportedRepresentation`: if the requested representation
          is unavailable for the source.
        - `genoio.InvalidSourceError`: if the source cannot be decoded.
        - `genoio.MissingDataError`: if retained missing calls conflict with
          the requested missing-data policy.
        """
        read_options = _ReadOptions(
            kind=kind,
            dosage=dosage,
            sparse=sparse,
            variants=variants,
            samples=samples,
            missing=missing,
            dtype=dtype,
            return_samples=return_samples,
            return_variants=return_variants,
        )
        validated_options = _validate_read_options(read_options)
        self._validate_source_supports_read(
            read_options.kind,
            read_options.dosage,
            validated_options.sparse_format,
        )

        return self._read_validated(
            read_options=read_options,
            validated_options=validated_options,
            variant_window=None,
        )

    def samples(self, **options: Any) -> Any:
        r"""Return sample metadata as a Polars DataFrame.

        Columns are `fid`, `iid`, `father`, `mother`, `sex`, and `phenotype`.
        Rows are ordered as they appear in the source. Haplotype reads that
        return sample metadata add `source_sample_index` and `haplotype_index`
        columns to map haplotype rows back to source samples.

        **Returns:**

        Polars DataFrame with source sample metadata in source order.
        """
        _reject_options(options)
        if self._samples_frame_cache is None:
            object.__setattr__(self, "_samples_frame_cache", samples_frame(self._metadata()["samples"]))
        return self._samples_frame_cache

    def variants(self, *, stats: Any = None, **options: Any) -> Any:
        r"""Return variant metadata as a Polars DataFrame.

        Columns are `chrom`, `pos`, `id`, `a0`, and `a1`. Rows are ordered as
        they appear in the source; variant frames returned by matrix reads are
        ordered to match matrix columns after filtering. The `a1` allele is the
        allele counted by returned genotype values.

        The `stats` argument is reserved for future metadata-stat controls.
        Passing it currently raises `genoio.InvalidOptionError`.

        **Arguments:**

        - `stats`: reserved; must be `None`.

        **Returns:**

        Polars DataFrame with source variant metadata in source order.
        """
        _validate_variant_stats(stats)
        _reject_options(options)
        if self._variants_frame_cache is None:
            object.__setattr__(self, "_variants_frame_cache", variants_frame(self._metadata()["variants"]))
        return self._variants_frame_cache

    def blocks(self, size: int, **read_options: Any) -> Any:
        r"""Yield consecutive variant blocks from this dataset.

        Each yielded block has at most `size` variants and follows the same
        return contract as [`genoio.Dataset.read`][]. Blocks are ordered by
        source variant order after any retained-variant filtering. BGEN dosage
        blocks with a concrete region filter use a same-path `.bgen.bgi` index
        when present.

        **Arguments:**

        - `size`: maximum number of variants per yielded block.
        - `read_options`: forwarded to [`genoio.Dataset.read`][].

        **Returns:**

        Iterator yielding matrices or matrix/metadata tuples.
        """
        _validate_block_size(size)
        normalized_options = _read_options_with_defaults(read_options)
        validated_options = _validate_read_options(normalized_options)
        self._validate_source_supports_read(
            normalized_options.kind,
            normalized_options.dosage,
            validated_options.sparse_format,
        )
        return self._block_iterator(size, normalized_options, validated_options)

    def _block_iterator(
        self,
        size: int,
        read_options: _ReadOptions,
        validated_options: _ValidatedReadOptions,
    ) -> Any:
        start = 0
        while True:
            # Variant windows are expressed in retained-variant coordinates.
            # Rust applies metadata/genotype filters before deciding whether a
            # retained variant falls into this block.
            block = self._read_validated(
                read_options=read_options,
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

    def _read_validated(
        self,
        *,
        read_options: _ReadOptions,
        validated_options: _ValidatedReadOptions,
        variant_window: dict[str, int] | None,
    ) -> Any:
        members = {key: str(path) for key, path in self.source.members.items()}
        options = {
            "samples": None if read_options.samples is None else list(read_options.samples),
            "variants": validated_options.variant_filter_ir,
            "variant_window": variant_window,
            "dosage": read_options.dosage,
            "return_samples": read_options.return_samples,
            "return_variants": read_options.return_variants,
            "matrix_only": (
                not read_options.return_samples
                and not read_options.return_variants
                and read_options.samples is None
                and validated_options.variant_filter_ir is None
            ),
        }
        if validated_options.sparse_format is None:
            rust_result = (
                self._read_dense_from_rust(members, options)
                if read_options.kind == "geno"
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
                if read_options.kind == "geno"
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
        # Metadata frames are assembled only when requested. Large PLINK2
        # block reads can otherwise avoid parsing full variant metadata.
        sample_metadata = samples_frame(rust_result["samples"]) if read_options.return_samples else None
        variant_metadata = variants_frame(rust_result["variants"]) if read_options.return_variants else None
        return read_result_tuple(
            genotype_matrix,
            sample_metadata,
            variant_metadata,
            return_samples=read_options.return_samples,
            return_variants=read_options.return_variants,
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
                if _is_unsupported_bgen_representation_error(str(error)):
                    raise UnsupportedRepresentation(str(error)) from error
                raise InvalidSourceError(str(error)) from error
            object.__setattr__(self, "_metadata_cache", metadata)
        assert self._metadata_cache is not None
        return self._metadata_cache

    def _validate_source_supports_kind(self, kind: str) -> None:
        if kind == "geno":
            return
        if self.source.format.value == "bgen":
            raise _unsupported_haplotype_source(self.source.format.value)
        capabilities = self._metadata()["capabilities"]
        if not capabilities["supports_haplo"]:
            raise _unsupported_haplotype_source(self.source.format.value)

    def _validate_source_supports_read(
        self,
        kind: str,
        dosage: str,
        sparse_format: str | None,
    ) -> None:
        if self.source.format.value == "bgen":
            if kind == "haplo":
                raise UnsupportedRepresentation(
                    'bgen haplo reads are not implemented; use kind="geno" with dosage="dosage"'
                )
            if dosage == "hardcall":
                if sparse_format is not None:
                    raise UnsupportedRepresentation(
                        'bgen sparse genotype reads are not implemented; use sparse=False with dosage="dosage"'
                    )
                raise UnsupportedRepresentation('bgen hardcall genotype reads are not implemented; use dosage="dosage"')
            return
        self._validate_source_supports_kind(kind)
        self._validate_source_supports_dosage(kind, dosage)

    def _validate_source_supports_dosage(self, kind: str, dosage: str) -> None:
        if dosage == "hardcall":
            return
        if kind == "haplo":
            raise UnsupportedRepresentation('kind="haplo" does not support dosage-backed reads')
        if self.source.format.value in {"vcf", "plink2"}:
            return
        raise UnsupportedRepresentation(f"{self.source.format.value} does not support dosage-backed genotype reads")


def vcf(path: str | Path) -> Dataset:
    r"""Resolve a VCF/BCF file and return a reusable dataset.

    **Arguments:**

    - `path`: `.vcf`, `.vcf.gz`, or `.bcf` path.

    **Returns:**

    [`genoio.Dataset`][] backed by the VCF/BCF source.

    **Raises:**

    - `genoio.SourceResolutionError`: if the path cannot be used as VCF/BCF.
    """
    return Dataset(source=resolve_vcf(path))


def bfile(path: str | Path) -> Dataset:
    r"""Resolve a PLINK1 BED/BIM/FAM file set and return a reusable dataset.

    `path` may be the shared prefix or one `.bed`, `.bim`, or `.fam` member.

    **Arguments:**

    - `path`: PLINK1 prefix or member path.

    **Returns:**

    [`genoio.Dataset`][] backed by the PLINK1 source.
    """
    return Dataset(source=resolve_bfile(path))


def bgen(path: str | Path) -> Dataset:
    r"""Resolve a BGEN source and return a reusable dataset.

    `path` may be the shared prefix or the `.bgen` member. If a same-prefix
    `.sample` file exists, it is recorded as an optional companion. Concrete
    region filters look for a same-path bgenix SQLite index beside the BGEN
    member, for example `cohort.bgen.bgi`.

    **Arguments:**

    - `path`: BGEN prefix or `.bgen` member path.

    **Returns:**

    [`genoio.Dataset`][] backed by the BGEN source.
    """
    return Dataset(source=resolve_bgen(path))


def pfile(path: str | Path) -> Dataset:
    r"""Resolve a PLINK2 PGEN/PVAR/PSAM file set and return a reusable dataset.

    `path` may be the shared prefix or one `.pgen`, `.pvar`, `.pvar.zst`, or
    `.psam` member. If both `.pvar` and `.pvar.zst` exist for a prefix,
    uncompressed `.pvar` is preferred.

    **Arguments:**

    - `path`: PLINK2 prefix or member path.

    **Returns:**

    [`genoio.Dataset`][] backed by the PLINK2 source.
    """
    return Dataset(source=resolve_pfile(path))


def _validate_read_options(options: _ReadOptions) -> _ValidatedReadOptions:
    _validate_kind(options.kind)
    _validate_dosage_source(options.dosage)
    sparse_format = _validate_sparse(options.sparse)
    normalized_missing = _normalize_missing(options.missing, sparse_format=sparse_format)
    variant_filter_ir = _variant_filter_ir(options.variants)
    normalized_dtype = _normalize_dtype(options.dtype)
    _validate_dosage_compatibility(options.dosage, sparse_format)
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


def _validate_dosage_source(dosage: str) -> None:
    if not isinstance(dosage, str) or dosage not in _SUPPORTED_DOSAGE_SOURCES:
        raise InvalidOptionError(f"unsupported dosage source: {dosage!r}")


def _validate_dosage_compatibility(
    dosage: str,
    sparse_format: str | None,
) -> None:
    if dosage == "hardcall":
        return
    if sparse_format is not None:
        raise UnsupportedRepresentation("sparse dosage-backed genotype reads are not implemented")


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


def _validate_missing(missing: str) -> None:
    if not isinstance(missing, str) or missing not in _SUPPORTED_MISSING_POLICIES:
        raise InvalidOptionError(f"unsupported missing-data policy: {missing}")


def _normalize_missing(missing: Any, *, sparse_format: str | None) -> str:
    if missing is _DEFAULT_MISSING:
        return "raise" if sparse_format is not None else "nan"
    _validate_missing(missing)
    return missing


def _normalize_dtype(dtype: DTypeLike) -> np.dtype[Any]:
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
    if "bgen" in message and "not implemented" in message and "dosage reads" not in message:
        return UnsupportedRepresentation(message)
    if _is_unsupported_bgen_representation_error(message):
        return UnsupportedRepresentation(message)
    if (
        "FORMAT/DS" in message
        or "pgen does not contain dosage values" in message
        or "pgen record does not contain dosage values" in message
        or ("unsupported pgen" in message and "dosage" in message)
    ):
        return UnsupportedRepresentation(message)
    return InvalidSourceError(message)


def _public_haplotype_read_error(error: ValueError) -> Exception:
    message = str(error)
    if (
        "unphased" in message
        or "unsupported haplotype format" in message
        or ("bgen" in message and "not implemented" in message)
    ):
        return UnsupportedRepresentation(message)
    if "sparse missing values" in message:
        return MissingDataError(message)
    return _public_read_error(error)


def _is_unsupported_bgen_representation_error(message: str) -> bool:
    unsupported_markers = (
        "unsupported bgen",
        "bgen layout",
        "bgen metadata parsing requires layout",
        "bgen compression value is reserved",
    )
    representation_markers = (
        "multiallelic",
        "phased",
        "variable-ploidy",
        "bit depth",
        "layout",
        "compression",
    )
    return any(marker in message for marker in unsupported_markers) and any(
        marker in message for marker in representation_markers
    )


def _validate_bool_option(name: str, value: bool) -> None:
    if not isinstance(value, bool):
        raise InvalidOptionError(f"{name} must be a bool")


def _validate_block_size(size: int) -> None:
    if not isinstance(size, int) or isinstance(size, bool) or size < 1:
        raise InvalidOptionError("block size must be a positive integer")


def _read_options_with_defaults(read_options: dict[str, Any]) -> _ReadOptions:
    defaults = {
        "kind": "geno",
        "dosage": "hardcall",
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
    merged: dict[str, Any] = defaults | read_options
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
