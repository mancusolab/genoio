# pattern: Functional Core

"""Read VCF, PLINK1, and PLINK2 genotype matrices from Python."""

from ._api import Dataset, bfile, pfile, vcf
from ._errors import (
    GenoioError,
    InvalidOptionError,
    InvalidSourceError,
    MissingCompanionFileError,
    MissingDataError,
    SampleFilterError,
    SourceResolutionError,
    UnsupportedFormatError,
    UnsupportedRepresentation,
)
from ._filters import FilterExpr, biallelic, chrom, id_in, mac, maf, missing_rate, polymorphic, qual, region, snp

__all__ = [
    "Dataset",
    "vcf",
    "bfile",
    "pfile",
    "FilterExpr",
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
    "SourceResolutionError",
    "MissingCompanionFileError",
    "UnsupportedFormatError",
    "InvalidSourceError",
    "UnsupportedRepresentation",
    "InvalidOptionError",
    "MissingDataError",
    "SampleFilterError",
]
