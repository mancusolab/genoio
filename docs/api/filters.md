# Filters

Filters are composable predicates over variant metadata and, when requested,
genotype-derived statistics. Metadata-only filters can usually be applied before
genotype decoding, while MAF, MAC, missing-rate, and polymorphism filters require
called genotype statistics for the selected samples. Use these constructors with
`Dataset.read` or `Dataset.blocks` when you need the returned matrix and variant
metadata to reflect the same retained variant set.

::: genoio.FilterExpr

::: genoio.chrom

::: genoio.region

::: genoio.snp

::: genoio.biallelic

::: genoio.qual

::: genoio.maf

::: genoio.mac

::: genoio.missing_rate

::: genoio.polymorphic

::: genoio.id_in
