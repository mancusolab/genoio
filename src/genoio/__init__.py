# pattern: Functional Core

"""Read VCF, BGEN, PLINK1, and PLINK2 genotype matrices from Python."""

from importlib.metadata import PackageNotFoundError, version

from ._api import BlockIterator, Dataset, bfile, bgen, pfile, vcf
from ._errors import (
    GenoioError,
    InternalError,
    InvalidOptionError,
    InvalidSourceError,
    MissingCompanionFileError,
    MissingDataError,
    SampleFilterError,
    SourceResolutionError,
    UnsupportedFormatError,
    UnsupportedRepresentation,
)
from ._filters import (
    FilterExpr,
    biallelic,
    chrom,
    id_in,
    mac,
    maf,
    missing_rate,
    polymorphic,
    qual,
    region,
    snp,
)

try:
    __version__ = version("genoio")
except PackageNotFoundError:
    __version__ = "0.4.1"

__all__ = [
    "__version__",
    "BlockIterator",
    "Dataset",
    "vcf",
    "bfile",
    "bgen",
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
    "InternalError",
]
