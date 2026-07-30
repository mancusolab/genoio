# pattern: Mixed
# Reason: Public dataset constructors perform filesystem source resolution;
# Dataset validation helpers remain pure.

from __future__ import annotations

import sys
from collections.abc import Callable, Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import TracebackType
from typing import Any, Literal, Self, TypeVar, cast, overload

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
    InternalError,
    InvalidOptionError,
    InvalidSourceError,
    MissingDataError,
    SampleFilterError,
    UnsupportedRepresentation,
)
from ._filters import FilterExpr
from ._read_options import (
    _read_options_with_defaults,
    _ReadOptions,
    _validate_read_options,
    _ValidatedReadOptions,
)
from ._source import ResolvedSource, resolve_bfile, resolve_bgen, resolve_pfile, resolve_vcf

_Region = TypeVar("_Region")
_NativeResult = TypeVar("_NativeResult")
_BlockResult = TypeVar("_BlockResult", covariant=True)
_RUST_ERROR_MAP = (
    (_rust.RustInternalError, InternalError),
    (_rust.RustSampleFilterError, SampleFilterError),
    (_rust.RustMissingDataError, MissingDataError),
    (_rust.RustUnsupportedRepresentationError, UnsupportedRepresentation),
    (_rust.RustInvalidOptionError, InvalidOptionError),
    (_rust.RustInvalidSourceError, InvalidSourceError),
)
# Rust exposes private exception classes so Python can preserve the public
# genoio error hierarchy without parsing backend message text.
_RUST_PUBLIC_ERROR_TYPES = tuple(error_type for error_type, _ in _RUST_ERROR_MAP)


@dataclass(frozen=True, slots=True)
class _ReadPayload:
    # Internal reads keep matrix and metadata separate. Public methods decide
    # later whether to expose metadata by returning a tuple.
    matrix: MatrixResult
    samples: pl.DataFrame | None
    variants: pl.DataFrame | None

    def to_result(self, options: _ReadOptions) -> ReadResult:
        return read_result_tuple(
            self.matrix,
            self.samples,
            self.variants,
            return_samples=options.return_samples,
            return_variants=options.return_variants,
        )


class BlockIterator(Iterator[_BlockResult]):
    r"""Iterate over blocks while owning one persistent native reader.

    `Dataset.iter_blocks()` constructs this object. Plain iteration closes the
    reader on exhaustion. Use it as a context manager when control flow may
    stop early so the native reader closes when the `with` block exits.

    `close()` is idempotent. A closed iterator remains exhausted.
    """

    __slots__ = (
        "_dataset",
        "_size",
        "_read_options",
        "_validated_options",
        "_members",
        "_options",
        "_native",
        "_closed",
        "__weakref__",
    )

    def __init__(
        self,
        dataset: Dataset,
        size: int,
        read_options: _ReadOptions,
        validated_options: _ValidatedReadOptions,
        members: dict[str, str],
        options: dict[str, Any],
    ) -> None:
        self._dataset = dataset
        self._size = size
        self._read_options = read_options
        self._validated_options = validated_options
        self._members = members
        self._options = options
        self._native: _rust._BlockReader | None = None
        self._closed = False

    def __iter__(self) -> Self:
        return self

    def __next__(self) -> _BlockResult:
        if self._closed:
            raise StopIteration
        try:
            native = self._native
            if native is None:
                native = _call_native(
                    _rust._BlockReader,
                    self._dataset.source.format.value,
                    self._members,
                    self._read_options.kind,
                    self._validated_options.sparse_format is not None,
                    self._options,
                    self._size,
                )
                self._native = native
            rust_result = _call_native(native.next_block)
            if rust_result is None:
                self.close()
                raise StopIteration
            payload = self._dataset._assemble_read_payload(
                rust_result,
                read_options=self._read_options,
                validated_options=self._validated_options,
            )
            return cast(_BlockResult, payload.to_result(self._read_options))
        except BaseException as primary:
            self._close_after_failure(primary)
            raise

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> Literal[False]:
        del exc_type, traceback
        try:
            self.close()
        except BaseException as close_error:
            if exc_value is None:
                raise
            exc_value.add_note(_close_failure_note(close_error))
        return False

    def close(self) -> None:
        """Close the native reader without opening an unadvanced iterator."""
        if self._closed:
            return
        self._closed = True
        native, self._native = self._native, None
        if native is not None:
            _call_native(native.close)

    def _close_after_failure(self, primary: BaseException) -> None:
        try:
            self.close()
        except BaseException as close_error:
            primary.add_note(_close_failure_note(close_error))

    def __del__(self) -> None:
        try:
            self.close()
        except BaseException:
            pass


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
        - `genoio.InternalError`: if the compiled backend reports an internal
          invariant failure.
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
    ) -> BlockIterator[SparseMatrixResult]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[True],
        return_variants: Literal[True],
        **read_options: object,
    ) -> BlockIterator[tuple[SparseMatrixResult, pl.DataFrame, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[True],
        return_variants: Literal[False] = False,
        **read_options: object,
    ) -> BlockIterator[tuple[SparseMatrixResult, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        sparse: Literal[True, "csc", "csr"],
        return_samples: Literal[False] = False,
        return_variants: Literal[True],
        **read_options: object,
    ) -> BlockIterator[tuple[SparseMatrixResult, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        return_samples: Literal[True],
        return_variants: Literal[True],
        **read_options: object,
    ) -> BlockIterator[tuple[np.ndarray, pl.DataFrame, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        return_samples: Literal[True],
        return_variants: Literal[False] = False,
        **read_options: object,
    ) -> BlockIterator[tuple[np.ndarray, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        *,
        return_samples: Literal[False] = False,
        return_variants: Literal[True],
        **read_options: object,
    ) -> BlockIterator[tuple[np.ndarray, pl.DataFrame]]: ...

    @overload
    def iter_blocks(
        self,
        size: int,
        **read_options: object,
    ) -> BlockIterator[np.ndarray]: ...

    def iter_blocks(
        self,
        size: int,
        **read_options: object,
    ) -> BlockIterator[ReadResult]:
        r"""Yield consecutive variant blocks from this dataset.

        Each yielded block has at most `size` variants and follows the same
        return contract as [`genoio.Dataset.read`][]. Blocks are fixed-width
        retained-variant chunks ordered by source variant order after any
        filtering. BGEN dosage blocks with a concrete region filter use a
        same-path `.bgen.bgi` index when present. Haplotype blocks follow the
        same source-encoded representation rules as [`genoio.Dataset.read`][].

        Source opening and record decoding are lazy: source and record errors
        can be raised when the iterator is first advanced or between yielded
        blocks. Exhaustion and read failure close the native reader. When a loop
        may stop before exhaustion, use the returned iterator as a context
        manager:

        ```python
        with dataset.iter_blocks(10_000, return_variants=True) as blocks:
            for matrix, variants in blocks:
                if analysis_is_complete(matrix, variants):
                    break
        ```

        **Arguments:**

        - `size`: maximum number of variants per yielded block.
        - `read_options`: forwarded to [`genoio.Dataset.read`][].

        **Returns:**

        [`genoio.BlockIterator`][] yielding matrices or matrix/metadata tuples.

        **Raises:**

        - `genoio.InvalidOptionError`: if `size` or a read option is invalid.
        - `genoio.InvalidSourceError`: if lazy source opening or decoding fails.
        - `genoio.MissingDataError`: if a yielded block violates the requested
          missing-data policy.
        """
        _validate_block_size(size)
        normalized_options = _read_options_with_defaults(read_options)
        validated_options = _validate_read_options(normalized_options)
        self._validate_source_supports_read(
            normalized_options.kind,
            normalized_options.dosage,
            validated_options.sparse_format,
        )
        members, options = self._native_read_arguments(
            read_options=normalized_options,
            validated_options=validated_options,
            variant_window=None,
        )
        return BlockIterator(
            self,
            size,
            normalized_options,
            validated_options,
            members,
            options,
        )

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

    def _read_validated(
        self,
        *,
        read_options: _ReadOptions,
        validated_options: _ValidatedReadOptions,
        variant_window: dict[str, int] | None,
    ) -> ReadResult:
        # Whole reads expose the public return contract after assembling the
        # stateless native payload.
        return self._read_payload(
            read_options=read_options,
            validated_options=validated_options,
            variant_window=variant_window,
        ).to_result(read_options)

    def _native_read_arguments(
        self,
        *,
        read_options: _ReadOptions,
        validated_options: _ValidatedReadOptions,
        variant_window: dict[str, int] | None,
    ) -> tuple[dict[str, str], dict[str, Any]]:
        members = {key: str(path) for key, path in self.source.members.items()}
        options = {
            "samples": None if read_options.samples is None else list(read_options.samples),
            "variants": validated_options.variant_filter_ir,
            "variant_window": variant_window,
            "dosage": read_options.dosage,
            "missing": validated_options.missing,
            "return_samples": read_options.return_samples,
            "return_variants": read_options.return_variants,
            "matrix_only": not read_options.return_samples and not read_options.return_variants,
        }
        return members, options

    def _assemble_read_payload(
        self,
        rust_result: dict[str, Any],
        *,
        read_options: _ReadOptions,
        validated_options: _ValidatedReadOptions,
    ) -> _ReadPayload:
        if validated_options.sparse_format is None:
            genotype_matrix = dense_array_from_rust(
                values=rust_result["values"],
                shape=tuple(rust_result["shape"]),
                dtype=validated_options.dtype,
                values_layout=str(rust_result.get("values_layout", "sample_major")),
            )
        else:
            genotype_matrix = sparse_matrix_from_rust(
                indptr=rust_result["indptr"],
                indices=rust_result["indices"],
                data=rust_result["data"],
                shape=tuple(rust_result["shape"]),
                dtype=validated_options.dtype,
                sparse_format=validated_options.sparse_format,
            )
        sample_metadata = (
            samples_frame(
                rust_result["samples"],
                include_haplotype_columns=read_options.kind == "haplo",
            )
            if read_options.return_samples
            else None
        )
        variant_metadata = variants_frame(rust_result["variants"]) if read_options.return_variants else None
        return _ReadPayload(genotype_matrix, sample_metadata, variant_metadata)

    def _read_payload(
        self,
        *,
        read_options: _ReadOptions,
        validated_options: _ValidatedReadOptions,
        variant_window: dict[str, int] | None,
    ) -> _ReadPayload:
        members, options = self._native_read_arguments(
            read_options=read_options,
            validated_options=validated_options,
            variant_window=variant_window,
        )
        rust_result = self._read_from_rust(
            read_options.kind,
            validated_options.sparse_format is not None,
            members,
            options,
        )
        return self._assemble_read_payload(
            rust_result,
            read_options=read_options,
            validated_options=validated_options,
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
        return _call_native(read, self.source.format.value, members, options)

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

    def _validate_source_supports_read(
        self,
        kind: str,
        dosage: str,
        sparse_format: str | None,
    ) -> None:
        try:
            _rust.validate_read_support(
                self.source.format.value,
                kind,
                dosage,
                sparse_format is not None,
            )
        except _RUST_PUBLIC_ERROR_TYPES as error:
            raise _public_rust_error(error) from error
        except ValueError as error:
            raise _public_read_error(error) from error


def vcf(path: str | Path) -> Dataset:
    r"""Resolve a VCF/BCF file and return a reusable dataset.

    **Arguments:**

    - `path`: `.vcf`, `.vcf.gz`, `.vcf.bgz`, or `.bcf` path.

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


def _reject_options(options: Mapping[str, object]) -> None:
    if options:
        keys = ", ".join(sorted(options))
        raise InvalidOptionError(f"unsupported option(s): {keys}")


def _validate_variant_stats(stats: object) -> None:
    if stats is not None:
        raise InvalidOptionError("variant stats are not implemented until a later phase")


def _call_native(
    native_call: Callable[..., _NativeResult],
    /,
    *args: Any,
    **kwargs: Any,
) -> _NativeResult:
    try:
        return native_call(*args, **kwargs)
    except _RUST_PUBLIC_ERROR_TYPES as error:
        raise _public_rust_error(error) from error
    except ValueError as error:
        raise _public_read_error(error) from error


def _close_failure_note(error: BaseException) -> str:
    return f"native block reader close failed: {type(error).__name__}: {error}"


def _public_rust_error(error: Exception) -> Exception:
    for rust_error_type, public_error_type in _RUST_ERROR_MAP:
        if isinstance(error, rust_error_type):
            return public_error_type(str(error))
    return InvalidSourceError(str(error))


def _public_read_error(error: ValueError) -> Exception:
    return InvalidSourceError(str(error))


def _validate_block_size(size: int) -> None:
    if not isinstance(size, int) or isinstance(size, bool) or size < 1:
        raise InvalidOptionError("block size must be a positive integer")
    if size > sys.maxsize:
        raise InvalidOptionError(f"block size exceeds this platform's limit of {sys.maxsize}")
