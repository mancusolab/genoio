# genoio

`genoio` is a genotype matrix IO layer for statistical genetics software. It
provides one Python API for VCF/BCF, PLINK1, and PLINK2 inputs, with Rust
readers for parsing, filtering, and matrix construction.

`genoio` is designed for tool developers and researchers who need stable
behavior across genotype formats:

- matrices are returned with samples on rows and variants on columns
- sample metadata rows match matrix rows
- variant metadata rows match matrix columns
- block reads stream retained variant chunks with matching metadata

## Installation

From a local checkout:

```bash
pip install .
```

For development, use the Makefile to create `.venv`, install development
dependencies, and build the Rust extension:

```bash
make build
```

Useful development targets:

```bash
make help
make build-release
make build-wheel
make verify
```

## Documentation

Available at [https://mancusolab.github.io/genoio](https://mancusolab.github.io/genoio).

## Quick Examples

Read sample metadata, align phenotypes and covariates by `iid`, then stream
genotype blocks into an association or QTL scanner:

```python
import polars as pl
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()

phenotypes = pl.read_csv("phenotypes.tsv", separator="\t")
covariates = pl.read_csv("covariates.tsv", separator="\t")

design = (
    samples.select("iid")
    .join(phenotypes.select("iid", "expression"), on="iid", how="left")
    .join(covariates.select("iid", "age", "sex", "PC1", "PC2"), on="iid", how="left")
)

y = design["expression"].to_numpy()
C = design.select("age", "sex", "PC1", "PC2").to_numpy()

variant_filter = genoio.region("22:20000000-21000000") & genoio.missing_rate(max=0.1)

for X, variants in ds.blocks(5_000, variants=variant_filter, return_variants=True):
    run_association_scan(X, y, covariates=C, samples=samples, variants=variants)
```

Open supported genotype sources with the matching constructor:

```python
vcf_ds = genoio.vcf("cohort.vcf.gz")
bed_ds = genoio.bfile("cohort")       # .bed/.bim/.fam
pgen_ds = genoio.pfile("cohort.pgen") # .pgen/.pvar/.psam
```

Use serializable filters when reading whole matrices or blocks:

```python
rare = (
    genoio.chrom("22")
    & genoio.maf(max=0.01)
    & genoio.missing_rate(max=0.1)
)

X, variants = pgen_ds.read(variants=rare, return_variants=True)
```
