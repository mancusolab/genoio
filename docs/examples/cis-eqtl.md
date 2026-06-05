# cis-eQTL

cis-eQTL scans use genomic windows as the iteration unit.
[`iter_regions(...)`](../api/reading.md#genoio.Dataset.iter_regions) keeps each
region result separate while still using indexed region reads when the source
format supports them.

This example assumes the expression matrix has one row per sample and one
column per gene. The `genes` table supplies the gene ID and its cis window:

| gene_id | cis_region |
|---|---|
| ENSG000001 | 22:20000000-21000000 |
| ENSG000002 | 22:22000000-23000000 |

The unusual part is the region list. A
[`region(...)`](../api/filters.md#genoio.region) is also a filter object. To
scan only common biallelic SNPs in each cis window, build one combined filter
per gene: `region & maf & snp & biallelic`. That lets `genoio` use region
pushdown first, then apply the remaining predicates while reading candidate
variants.

```python
import polars as pl
import genoio

ds = genoio.bgen("data/chr22_hg38.bgen")

# `samples` defines the row order of X_region for every returned region.
samples = ds.samples()

expression = pl.read_csv("expression.tsv", separator="\t")
covariates = pl.read_csv("covariates.tsv", separator="\t")
genes = pl.read_csv("cis_windows.tsv", separator="\t")

design = (
    # Start from genotype sample order so expression, covariates, and genotypes align.
    samples.select("iid")
    .join(expression, on="iid", how="left")
    .join(covariates.select("iid", "age", "sex", "PC1", "PC2"), on="iid", how="left")
)
C = design.select("age", "sex", "PC1", "PC2").to_numpy()

# This shared filter is applied inside every cis window.
variant_filter = genoio.maf(max=0.05) & genoio.snp() & genoio.biallelic()

# Each list element is a complete read filter for one gene's cis window.
# The region predicate can use indexed VCF/BGEN reads when an index exists.
regions = [genoio.region(region) & variant_filter for region in genes["cis_region"]]

for region_index, (region, (X_region, variants)) in enumerate(
    ds.iter_regions(regions, dosage="dosage", return_variants=True)
):
    # Use the region index to pull the matching gene metadata and phenotype.
    gene = genes.row(region_index, named=True)
    y_region = design[gene["gene_id"]].to_numpy()

    # X_region contains only variants retained for this gene's cis window.
    cis_scan(X_region, y_region, C, gene=gene["gene_id"], region=region, variants=variants)
```

Use [`iter_regions(...)`](../api/reading.md#genoio.Dataset.iter_regions) when
each genomic interval has its own phenotype or analysis context. Unlike
[`iter_blocks(...)`](../api/reading.md#genoio.Dataset.iter_blocks), the
iterator boundary is part of the scientific question, not just a matrix chunk
size.
