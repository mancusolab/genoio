# pattern: Functional Core


class GenoioError(Exception):
    """Base class for public genoio errors."""


class SourceResolutionError(GenoioError):
    """Base class for source resolution errors."""


class AmbiguousSourceError(SourceResolutionError):
    """Raised when a path can resolve to multiple source formats."""


class MissingCompanionFileError(SourceResolutionError):
    """Raised when a multi-file source is missing required companion files."""


class UnsupportedFormatError(SourceResolutionError):
    """Raised when a source format is not supported."""


class InvalidSourceError(SourceResolutionError):
    """Raised when a source path cannot be used."""


class UnsupportedRepresentation(GenoioError):
    """Raised when a requested output representation is not supported."""


class InvalidOptionError(GenoioError):
    """Raised when a public API option is invalid."""


class MissingDataError(GenoioError):
    """Raised when requested data is unavailable."""


class SampleFilterError(GenoioError):
    """Raised when a sample keep list cannot be satisfied."""
