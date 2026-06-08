# Errors

`genoio` reports expected failures through public Python exception classes. Use
these classes when callers need to distinguish source-resolution failures,
invalid options, unsupported representations, missing data, or sample-filter
errors.

`genoio.InternalError` reports an unexpected backend invariant failure. Treat it
as a bug.

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
