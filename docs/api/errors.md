# Errors

These exceptions describe failures at the Python API boundary. They separate
source discovery problems, unsupported formats or representations, invalid
options, missing data policies, and sample-filter failures so callers can decide
whether to retry, change options, or reject an input file. Lower-level parser
details are preserved in the exception message when malformed companion files or
unsupported genotype records are encountered.

::: genoio.GenoioError

::: genoio.SourceResolutionError

::: genoio.MissingCompanionFileError

::: genoio.UnsupportedFormatError

::: genoio.InvalidSourceError

::: genoio.UnsupportedRepresentation

::: genoio.InvalidOptionError

::: genoio.MissingDataError

::: genoio.SampleFilterError
