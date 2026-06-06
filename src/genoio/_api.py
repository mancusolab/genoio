# pattern: Mixed
# Reason: Public dataset constructors perform filesystem source resolution;
# Dataset validation helpers remain pure.

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, TypeVar, cast, overload

import numpy as np
import polars as pl
from numpy.typing import DTypeLike

from . import _rust
from ._assembly import (
    MatrixResult,
    ReadResult,
    SparseMatrixResult,
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
_Region = TypeVar("_Region")
_SPARSE_DOSAGE_BACKED_GENOTYPE_UNSUPPORTED = "sparse dosage-backed genotype reads are intentionally unsupported"
_SPARSE_DOSAGE_BACKED_HAPLOTYPE_UNSUPPORTED = (
    "sparse haplotype reads are intentionally unsupported for dosage-backed sources; "
    "use dense haplotype reads with sparse=False"
)
_PLINK2_SPARSE_DOSAGE_BACKED_HAPLOTYPE_UNSUPPORTED = (
    "plink2 sparse haplotype reads are intentionally unsupported for dosage-backed sources; "
    "use dense haplotype reads with sparse=False"
)
_RUST_ERROR_MAP = (
    (_rust.RustSampleFilterError, SampleFilterError),
    (_rust.RustMissingDataError, MissingDataError),
    (_rust.RustUnsupportedRepresentationError, UnsupportedRepresentation),
    (_rust.RustInvalidOptionError, InvalidOptionError),
    (_rust.RustInvalidSourceError, InvalidSourceError),
)
_RUST_PUBLIC_ERROR_TYPES = tuple(error_type for error_type, _ in _RUST_ERROR_MAP)


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
    dtype: np.dtype[Any]
    missing: str
    sparse_format: str | None
    variant_filter_ir: dict[str, Any] | None


@dataclass(frozen=True, slots=True)
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
    _samples_frame_cache: pl.DataFrame | None = field(default=None, init=False, compare=False, repr=False)
    _variants_frame_cache: pl.DataFrame | None = field(default=None, init=False, compare=False, repr=False)

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: Literal[False] = False,
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: Literal[False] = False,
        return_variants: Literal[False] = False,
    ) -> np.ndarray: ...

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: Literal[True, "csc", "csr"],
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: Literal[False] = False,
        return_variants: Literal[False] = False,
    ) -> SparseMatrixResult: ...

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: Literal[False] = False,
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: Literal[True],
        return_variants: Literal[False] = False,
    ) -> tuple[np.ndarray, pl.DataFrame]: ...

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: Literal[False] = False,
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: Literal[False] = False,
        return_variants: Literal[True],
    ) -> tuple[np.ndarray, pl.DataFrame]: ...

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: Literal[False] = False,
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: Literal[True],
        return_variants: Literal[True],
    ) -> tuple[np.ndarray, pl.DataFrame, pl.DataFrame]: ...

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: Literal[True, "csc", "csr"],
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: Literal[True],
        return_variants: Literal[False] = False,
    ) -> tuple[SparseMatrixResult, pl.DataFrame]: ...

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: Literal[True, "csc", "csr"],
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: Literal[False] = False,
        return_variants: Literal[True],
    ) -> tuple[SparseMatrixResult, pl.DataFrame]: ...

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: Literal[True, "csc", "csr"],
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: Literal[True],
        return_variants: Literal[True],
    ) -> tuple[SparseMatrixResult, pl.DataFrame, pl.DataFrame]: ...

    @overload
    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: bool | str = False,
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: bool = False,
        return_variants: bool = False,
    ) -> ReadResult: ...

    def read(
        self,
        *,
        kind: str = "geno",
        dosage: str = "hardcall",
        sparse: bool | str = False,
        variants: FilterExpr | Iterable[str] | None = None,
        samples: list[str] | tuple[str, ...] | set[str] | None = None,
        missing: Literal["nan", "raise", "impute"] | None = None,
        dtype: DTypeLike = "float32",
        return_samples: bool = False,
        return_variants: bool = False,
    ) -> ReadResult:
        r"""Read a genotype or haplotype matrix from this dataset.

        Dense reads return a NumPy array with shape `(samples, variants)`.
        Sparse reads return SciPy CSC by default or CSR when `sparse="csr"`.
        Set `return_samples` or `return_variants` to return Polars metadata
        frames with the matrix.

        **Arguments:**

        - `kind`: Matrix row layout. `"geno"` returns one row per retained
          sample, with diploid genotype values in each cell. `"haplo"` returns
          one row per source haplotype, so diploid samples contribute two rows.
          Haplotype reads require phased records in the source.
        - `dosage`: Value source for each matrix cell. `"hardcall"` reads allele
          counts from called genotypes. `"dosage"` reads expected allele counts
          from dosage/probability fields when the source format supports them.
          `genoio` does not convert dosages into hard calls.
        - `sparse`: Output storage. `False` returns a dense NumPy array.
          `True` and `"csc"` return a SciPy CSC matrix; `"csr"` returns a SciPy
          CSR matrix. Sparse reads require `missing="raise"` because this
          release does not store sparse missing-value masks.
        - `variants`: Variants to keep. Pass a `genoio` filter expression to
          filter by metadata or genotype predicates, or pass an iterable of
          variant IDs to keep matching IDs. `None` keeps all variants. Retained
          columns stay in source order, not request order.
        - `samples`: Sample IDs to keep. Pass a list, tuple, or set of sample
          IDs. `None` keeps all samples. Retained rows stay in source order;
          duplicate requested IDs are rejected.
        - `missing`: Missing-call policy. `None` uses `"nan"` for dense reads
          and `"raise"` for sparse reads. `"nan"` stores missing calls as
          `np.nan`, `"raise"` fails if retained calls are missing, and
          `"impute"` fills missing calls with the retained variant mean.
        - `dtype`: NumPy dtype for returned matrix values. Missing policies
          that write `np.nan` or imputed means require a floating dtype.
        - `return_samples`: When `True`, return a sample metadata frame with
          the matrix. Haplotype reads include columns that map haplotype rows
          back to source samples.
        - `return_variants`: When `True`, return a variant metadata frame for
          the retained matrix columns.

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

    def samples(self, **options: object) -> pl.DataFrame:
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
        samples = self._samples_frame_cache
        assert samples is not None
        return samples

    def variants(self, *, stats: object = None, **options: object) -> pl.DataFrame:
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
        variants = self._variants_frame_cache
        assert variants is not None
        return variants

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[False] = False,
        return_variants: Literal[False] = False,
        **read_options: object,
    ) -> Iterator[SparseMatrixResult]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[True],
        return_variants: Literal[True],
        **read_options: object,
    ) -> Iterator[tuple[SparseMatrixResult, pl.DataFrame, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[True],
        return_variants: Literal[False] = False,
        **read_options: object,
    ) -> Iterator[tuple[SparseMatrixResult, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[False] = False,
        return_variants: Literal[True],
        **read_options: object,
    ) -> Iterator[tuple[SparseMatrixResult, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        return_samples: Literal[True],
        return_variants: Literal[True],
        **read_options: object,
    ) -> Iterator[tuple[np.ndarray, pl.DataFrame, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        return_samples: Literal[True],
        return_variants: Literal[False] = False,
        **read_options: object,
    ) -> Iterator[tuple[np.ndarray, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        return_samples: Literal[False] = False,
        return_variants: Literal[True],
        **read_options: object,
    ) -> Iterator[tuple[np.ndarray, pl.DataFrame]]: ...

    @overload
    def iter_blocks(self, size: int, **read_options: object) -> Iterator[np.ndarray]: ...

    def iter_blocks(self, size: int, **read_options: object) -> Iterator[ReadResult]:
        r"""Yield consecutive variant blocks from this dataset.

        Each yielded block has at most `size` variants and follows the same
        return contract as [`genoio.Dataset.read`][]. Blocks are fixed-width
        retained-variant chunks ordered by source variant order after any
        filtering. BGEN dosage blocks with a concrete region filter use a
        same-path `.bgen.bgi` index when present. Haplotype blocks follow the
        same source-encoded representation rules as [`genoio.Dataset.read`][].

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

    @overload
    def iter_regions(
        self,
        regions: Iterable[_Region],
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[False] = False,
        return_variants: Literal[False] = False,
        **read_options: object,
    ) -> Iterator[tuple[_Region, SparseMatrixResult]]: ...

    @overload
    def iter_regions(
        self,
        regions: Iterable[_Region],
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[True],
        return_variants: Literal[True],
        **read_options: object,
    ) -> Iterator[tuple[_Region, tuple[SparseMatrixResult, pl.DataFrame, pl.DataFrame]]]: ...

    @overload
    def iter_regions(
        self,
        regions: Iterable[_Region],
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[True],
        return_variants: Literal[False] = False,
        **read_options: object,
    ) -> Iterator[tuple[_Region, tuple[SparseMatrixResult, pl.DataFrame]]]: ...

    @overload
    def iter_regions(
        self,
        regions: Iterable[_Region],
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[False] = False,
        return_variants: Literal[True],
        **read_options: object,
    ) -> Iterator[tuple[_Region, tuple[SparseMatrixResult, pl.DataFrame]]]: ...

    @overload
    def iter_regions(
        self,
        regions: Iterable[_Region],
        *,
        return_samples: Literal[True],
        return_variants: Literal[True],
        **read_options: object,
    ) -> Iterator[tuple[_Region, tuple[np.ndarray, pl.DataFrame, pl.DataFrame]]]: ...

    @overload
    def iter_regions(
        self,
        regions: Iterable[_Region],
        *,
        return_samples: Literal[True],
        return_variants: Literal[False] = False,
        **read_options: object,
    ) -> Iterator[tuple[_Region, tuple[np.ndarray, pl.DataFrame]]]: ...

    @overload
    def iter_regions(
        self,
        regions: Iterable[_Region],
        *,
        return_samples: Literal[False] = False,
        return_variants: Literal[True],
        **read_options: object,
    ) -> Iterator[tuple[_Region, tuple[np.ndarray, pl.DataFrame]]]: ...

    @overload
    def iter_regions(
        self,
        regions: Iterable[_Region],
        **read_options: object,
    ) -> Iterator[tuple[_Region, np.ndarray]]: ...

    def iter_regions(
        self,
        regions: Iterable[_Region],
        **read_options: object,
    ) -> Iterator[tuple[_Region, ReadResult]]:
        r"""Yield one read result per requested region filter.

        Each yielded item is `(region, result)`, where `region` is the original
        object from `regions` and `result` follows the same return contract as
        [`genoio.Dataset.read`][]. Concrete VCF/BCF and BGEN region filters use
        the same indexed pushdown paths as normal reads when an index is
        present. Haplotype region reads follow the same source-encoded
        representation rules as [`genoio.Dataset.read`][].

        **Arguments:**

        - `regions`: iterable of region filter expressions.
        - `read_options`: forwarded to [`genoio.Dataset.read`][], except
          `variants`, which is supplied by each region.

        **Returns:**

        Iterator yielding `(region, matrix_or_tuple)` pairs.
        """
        if "variants" in read_options:
            raise InvalidOptionError("iter_regions supplies variants from the regions argument")
        return self._region_iterator(regions, read_options)

    def _region_iterator(
        self,
        regions: Iterable[_Region],
        read_options: Mapping[str, object],
    ) -> Iterator[tuple[_Region, ReadResult]]:
        read = cast(Any, self.read)
        for region in regions:
            yield region, cast(ReadResult, read(variants=region, **read_options))

    def _block_iterator(
        self,
        size: int,
        read_options: _ReadOptions,
        validated_options: _ValidatedReadOptions,
    ) -> Iterator[ReadResult]:
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
            genotype_matrix = cast(MatrixResult, block[0] if isinstance(block, tuple) else block)
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
    ) -> ReadResult:
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
            rust_result = self._read_from_rust(read_options.kind, False, members, options)
            genotype_matrix = dense_array_from_rust(
                values=rust_result["values"],
                shape=tuple(rust_result["shape"]),
                missing_mask=rust_result["missing_mask"],
                missing=validated_options.missing,
                dtype=validated_options.dtype,
            )
        else:
            rust_result = self._read_from_rust(read_options.kind, True, members, options)
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
        sample_metadata = (
            samples_frame(
                rust_result["samples"],
                include_haplotype_columns=read_options.kind == "haplo",
            )
            if read_options.return_samples
            else None
        )
        variant_metadata = variants_frame(rust_result["variants"]) if read_options.return_variants else None
        return read_result_tuple(
            genotype_matrix,
            sample_metadata,
            variant_metadata,
            return_samples=read_options.return_samples,
            return_variants=read_options.return_variants,
        )

    def _read_from_rust(
        self,
        kind: str,
        sparse: bool,
        members: dict[str, str],
        options: dict[str, Any],
    ) -> dict[str, Any]:
        if kind == "haplo":
            read = _rust.read_haplotypes_sparse if sparse else _rust.read_haplotypes_dense
        else:
            read = _rust.read_sparse if sparse else _rust.read_dense
        try:
            return read(self.source.format.value, members, options)
        except _RUST_PUBLIC_ERROR_TYPES as error:
            raise _public_rust_error(error) from error
        except ValueError as error:
            raise _public_read_error(error) from error

    def _metadata(self) -> dict[str, Any]:
        if self._metadata_cache is None:
            members = {key: str(path) for key, path in self.source.members.items()}
            try:
                metadata = _rust.read_metadata(self.source.format.value, members)
            except _RUST_PUBLIC_ERROR_TYPES as error:
                raise _public_rust_error(error) from error
            except ValueError as error:
                raise _public_read_error(error) from error
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
        if kind == "haplo" and sparse_format is not None and self.source.format.value == "bgen":
            raise UnsupportedRepresentation(
                f"{self.source.format.value} sparse haplotype reads are not implemented; "
                "use dense haplotype reads with sparse=False"
            )
        if self.source.format.value == "bgen":
            if kind == "haplo":
                if dosage == "dosage" and sparse_format is None:
                    return
                raise UnsupportedRepresentation(
                    'bgen hardcall haplotype reads are not implemented; use dosage="dosage" for '
                    "source-encoded phased haplotype dosage"
                )
            if dosage == "hardcall":
                if sparse_format is not None:
                    raise UnsupportedRepresentation(
                        'bgen sparse genotype reads are not implemented; use sparse=False with dosage="dosage"'
                    )
                raise UnsupportedRepresentation('bgen hardcall genotype reads are not implemented; use dosage="dosage"')
            return
        if kind == "haplo" and dosage == "dosage" and self.source.format.value == "vcf":
            raise UnsupportedRepresentation(
                "VCF haplotype dosage reads are unsupported because VCF haplotype support is hardcall GT-based"
            )
        if kind == "haplo" and self.source.format.value == "plink2":
            self._validate_source_supports_dosage(kind, dosage)
            return
        self._validate_source_supports_kind(kind)
        self._validate_source_supports_dosage(kind, dosage)

    def _validate_source_supports_dosage(self, kind: str, dosage: str) -> None:
        if dosage == "hardcall":
            return
        if kind == "haplo":
            if self.source.format.value in {"plink2", "bgen"}:
                return
            raise UnsupportedRepresentation(
                f"{self.source.format.value} does not support dosage-backed haplotype reads"
            )
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
    _validate_dosage_compatibility(options.kind, options.dosage, sparse_format)
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


def _reject_options(options: Mapping[str, object]) -> None:
    if options:
        keys = ", ".join(sorted(options))
        raise InvalidOptionError(f"unsupported option(s): {keys}")


def _validate_variant_stats(stats: object) -> None:
    if stats is not None:
        raise InvalidOptionError("variant stats are not implemented until a later phase")


def _validate_kind(kind: str) -> None:
    if not isinstance(kind, str) or kind not in _SUPPORTED_KINDS:
        raise UnsupportedRepresentation(f"unsupported genotype kind: {kind}")


def _validate_dosage_source(dosage: str) -> None:
    if not isinstance(dosage, str) or dosage not in _SUPPORTED_DOSAGE_SOURCES:
        raise InvalidOptionError(f"unsupported dosage source: {dosage!r}")


def _validate_dosage_compatibility(
    kind: str,
    dosage: str,
    sparse_format: str | None,
) -> None:
    if dosage == "hardcall":
        return
    if sparse_format is not None:
        if kind == "haplo":
            raise UnsupportedRepresentation(_SPARSE_DOSAGE_BACKED_HAPLOTYPE_UNSUPPORTED)
        raise UnsupportedRepresentation(_SPARSE_DOSAGE_BACKED_GENOTYPE_UNSUPPORTED)


def _unsupported_haplotype_source(source_format: str) -> UnsupportedRepresentation:
    if source_format == "vcf":
        return UnsupportedRepresentation(
            'VCF source has no phased GT evidence; kind="haplo" requires phased retained variants'
        )
    return UnsupportedRepresentation(f"{source_format} does not support haplo reads")


def _validate_sparse(sparse: object) -> str | None:
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
    if isinstance(missing, str) and missing in {"nan", "impute"} and not np.issubdtype(dtype, np.floating):
        raise InvalidOptionError(f'missing="{missing}" requires a floating dtype')


def _validate_sparse_missing_compatibility(sparse_format: str | None, missing: str) -> None:
    if sparse_format is not None and missing in {"nan", "impute"}:
        raise InvalidOptionError("this release does not store sparse missing values; use missing='raise'")


def _public_rust_error(error: Exception) -> Exception:
    for rust_error_type, public_error_type in _RUST_ERROR_MAP:
        if isinstance(error, rust_error_type):
            return public_error_type(str(error))
    return InvalidSourceError(str(error))


def _public_read_error(error: ValueError) -> Exception:
    return InvalidSourceError(str(error))


def _validate_bool_option(name: str, value: bool) -> None:
    if not isinstance(value, bool):
        raise InvalidOptionError(f"{name} must be a bool")


def _validate_block_size(size: int) -> None:
    if not isinstance(size, int) or isinstance(size, bool) or size < 1:
        raise InvalidOptionError("block size must be a positive integer")


def _read_options_with_defaults(read_options: Mapping[str, object]) -> _ReadOptions:
    unknown = set(read_options) - set(_READ_OPTION_DEFAULTS)
    if unknown:
        keys = ", ".join(sorted(unknown))
        raise InvalidOptionError(f"unsupported option(s): {keys}")
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
