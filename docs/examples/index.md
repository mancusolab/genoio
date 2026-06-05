# Examples

These examples show the two main iteration patterns in `genoio`.

The main difference from many Python genetics workflows is where filtering
happens. Instead of loading a full matrix and then slicing it in Python, you
build filter objects such as [`genoio.maf(max=0.05)`](../api/filters.md#genoio.maf)
or [`genoio.region("22:20000000-21000000")`](../api/filters.md#genoio.region).
`genoio` sends those filters to Rust with the read request, so the reader can
skip records before constructing the matrix when the source format allows it.

- [GWAS](gwas.md): one phenotype vector and one covariate matrix across
  fixed-width genotype blocks.
- [cis-eQTL](cis-eqtl.md): one genomic window at a time with
  [`iter_regions(...)`](../api/reading.md#genoio.Dataset.iter_regions).
