# Errors

`genoio` reports expected failures through public Python exception classes.
Source discovery problems, unsupported formats or representations, invalid
options, missing data policies, and sample-filter failures are separate classes
so callers can decide whether to retry, change options, or reject an input file.
Lower-level Rust parser details are preserved in the exception message when
malformed companion files or unsupported genotype records are encountered.

Unexpected compiled-backend invariant failures are reported as
`genoio.InternalError`. Treat these as bugs rather than user-input errors.

::: genoio.GenoioError

::: genoio.SourceResolutionError

::: genoio.MissingCompanionFileError

::: genoio.UnsupportedFormatError

::: genoio.InvalidSourceError

::: genoio.UnsupportedRepresentation

::: genoio.InvalidOptionError

::: genoio.MissingDataError

::: genoio.SampleFilterError

::: genoio.InternalError
