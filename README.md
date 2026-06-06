# genoio

`genoio` is a genotype matrix IO layer for statistical genetics software. It
provides one Python API for VCF/BCF, BGEN, PLINK1, and PLINK2 inputs, with Rust
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

For development, install `uv`, then use the Makefile to sync `.venv` from
`uv.lock` and build the Rust extension:

```bash
make build
```

Useful development targets:

```bash
make help
make build-release
make build-wheel
make lock
make verify
```

## Documentation

Available at [https://mancusolab.github.io/genoio](https://mancusolab.github.io/genoio).

## Quick Examples

### GWAS

```python
import polars as pl
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()

phenotypes = pl.read_csv("phenotypes.tsv", separator="\t")
covariates = pl.read_csv("covariates.tsv", separator="\t")

design = (
    samples.select("iid")
    .join(phenotypes.select("iid", "trait"), on="iid", how="left")
    .join(covariates.select("iid", "age", "sex", "PC1", "PC2"), on="iid", how="left")
)

y = design["trait"].to_numpy()
C = design.select("age", "sex", "PC1", "PC2").to_numpy()

variant_filter = genoio.maf(max=0.05) & genoio.snp() & genoio.biallelic()

for X, variants in ds.iter_blocks(10_000, variants=variant_filter, return_variants=True):
    association_scan(X, y, C, variants=variants)
```

### cis-eQTL

```python
import polars as pl
import genoio

ds = genoio.bgen("data/chr22_hg38.bgen")
samples = ds.samples()

expression = pl.read_csv("expression.tsv", separator="\t")
covariates = pl.read_csv("covariates.tsv", separator="\t")
genes = pl.read_csv("cis_windows.tsv", separator="\t")

design = (
    samples.select("iid")
    .join(expression, on="iid", how="left")
    .join(covariates.select("iid", "age", "sex", "PC1", "PC2"), on="iid", how="left")
)
C = design.select("age", "sex", "PC1", "PC2").to_numpy()

variant_filter = genoio.maf(max=0.05) & genoio.snp() & genoio.biallelic()
regions = [genoio.region(region) & variant_filter for region in genes["cis_region"]]

for region_index, (region, (X_region, variants)) in enumerate(
    ds.iter_regions(regions, dosage="dosage", return_variants=True)
):
    gene = genes.row(region_index, named=True)
    y_region = design[gene["gene_id"]].to_numpy()
    cis_scan(X_region, y_region, C, gene=gene["gene_id"], region=region, variants=variants)
```

Open supported genotype sources with the matching constructor:

```python
vcf_ds = genoio.vcf("cohort.vcf.gz")
bed_ds = genoio.bfile("cohort")       # .bed/.bim/.fam
pgen_ds = genoio.pfile("cohort.pgen") # .pgen/.pvar[.zst]/.psam
bgen_ds = genoio.bgen("cohort.bgen")  # .bgen plus optional .sample
```

Use `dosage="dosage"` for stored dosage values. BGEN v1.2+ Layout 2 biallelic
diploid dosage records are returned as expected A1 allele counts. Genotype reads
of phased BGEN records sum source haplotype probabilities to expected diploid A1
dosage; `kind="haplo", dosage="dosage"` returns expected A1 dosage per
haplotype row:

```python
X = bgen_ds.read(dosage="dosage")
H = bgen_ds.read(kind="haplo", dosage="dosage")
```

PLINK2 haplotype reads support source-encoded explicit phased hardcalls and
explicit phased full dosages. Explicit phased hardcall haplotypes can also be
read sparsely when retained calls are non-missing:

```python
H_hardcall = pgen_ds.read(kind="haplo", dosage="hardcall")
H_hardcall_sparse = pgen_ds.read(kind="haplo", dosage="hardcall", sparse=True)
H_dosage = pgen_ds.read(kind="haplo", dosage="dosage")
```

For BGEN region reads, place a bgenix SQLite index beside the source as
`cohort.bgen.bgi`. Concrete region filters use that index when present and
fall back to a sequential BGEN scan when it is absent:

```python
X, variants = bgen_ds.read(
    dosage="dosage",
    variants=genoio.region("22:20000000-21000000"),
    return_variants=True,
)
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
