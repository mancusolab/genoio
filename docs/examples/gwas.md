# GWAS

GWAS scans usually use one phenotype vector and one covariate matrix across
many fixed-width genotype blocks.

This example keeps the analysis setup simple: sample IDs are assumed to match
between genotype, phenotype, and covariate files. In a production scan, you
would also check missing phenotypes, duplicate IDs, and covariate encoding.

The important `genoio` pattern is the block iterator.
[`iter_blocks(...)`](../api/reading.md#genoio.Dataset.iter_blocks) reads at
most `size` retained variants per iteration. The variant filter is not a Python
callback. It's a small expression object that Rust evaluates while reading the
source.

```python
import polars as pl
import genoio

ds = genoio.pfile("data/chr22_hg38")

# `samples` defines the row order of every matrix returned by this dataset.
samples = ds.samples()

phenotypes = pl.read_csv("phenotypes.tsv", separator="\t")
covariates = pl.read_csv("covariates.tsv", separator="\t")

design = (
    # Start from genotype sample order so y and C match matrix rows.
    samples.select("iid")
    .join(phenotypes.select("iid", "trait"), on="iid", how="left")
    .join(covariates.select("iid", "age", "sex", "PC1", "PC2"), on="iid", how="left")
)

y = design["trait"].to_numpy()
C = design.select("age", "sex", "PC1", "PC2").to_numpy()

# Filters are combined with Python operators, then evaluated by the Rust reader.
# This keeps the returned matrices limited to common biallelic SNPs.
variant_filter = genoio.maf(max=0.05) & genoio.snp() & genoio.biallelic()

for X, variants in ds.iter_blocks(10_000, variants=variant_filter, return_variants=True):
    # X has shape (samples, variants_in_this_block).
    # variants describes the columns in X for this block only.
    association_scan(X, y, C, variants=variants)
```

Use [`iter_blocks(...)`](../api/reading.md#genoio.Dataset.iter_blocks) when your
analysis can process any fixed number of retained variants at a time. That is a
good fit for genome-wide scans because the block boundary is only a memory and
throughput choice.
