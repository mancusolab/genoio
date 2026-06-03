# pattern: Functional Core


class GenoioError(Exception):
    r"""Base class for public genoio errors."""


class SourceResolutionError(GenoioError):
    r"""Base class for source resolution errors."""


class MissingCompanionFileError(SourceResolutionError):
    r"""Raised when a multi-file source is missing required companion files."""


class UnsupportedFormatError(SourceResolutionError):
    r"""Raised when a source format is not supported."""


class InvalidSourceError(SourceResolutionError):
    r"""Raised when a source path cannot be used."""


class UnsupportedRepresentation(GenoioError):
    r"""Raised when a requested output representation is not supported."""


class InvalidOptionError(GenoioError):
    r"""Raised when a public API option is invalid."""


class MissingDataError(GenoioError):
    r"""Raised when requested data is unavailable."""


class SampleFilterError(GenoioError):
    r"""Raised when a sample keep list cannot be satisfied."""
