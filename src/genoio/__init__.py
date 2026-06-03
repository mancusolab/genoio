# pattern: Functional Core

"""Read VCF, PLINK1, and PLINK2 genotype matrices from Python."""

from ._api import Dataset, blocks, open, read, samples, variants
from ._errors import (
    AmbiguousSourceError,
    GenoioError,
    InvalidOptionError,
    InvalidSourceError,
    MissingCompanionFileError,
    MissingDataError,
    SampleFilterError,
    UnsupportedFormatError,
    UnsupportedRepresentation,
)
from ._filters import biallelic, chrom, id_in, mac, maf, missing_rate, polymorphic, qual, region, snp

__all__ = [
    "Dataset",
    "open",
    "read",
    "samples",
    "variants",
    "blocks",
    "chrom",
    "region",
    "snp",
    "biallelic",
    "maf",
    "qual",
    "mac",
    "missing_rate",
    "polymorphic",
    "id_in",
    "GenoioError",
    "AmbiguousSourceError",
    "MissingCompanionFileError",
    "UnsupportedFormatError",
    "InvalidSourceError",
    "UnsupportedRepresentation",
    "InvalidOptionError",
    "MissingDataError",
    "SampleFilterError",
]
