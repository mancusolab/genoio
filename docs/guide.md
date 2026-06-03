# Getting started

Most downstream tools want one stable sample axis and many chunks of variant
columns. `genoio` makes that pattern explicit: resolve the source once, read
sample metadata once, then process matrix blocks with their matching variant
metadata.

```python
import genoio

ds = genoio.open("data/chr22_hg38", format="plink2")
samples = ds.samples()

for X, variants in ds.blocks(10_000, return_variants=True):
    run_association_scan(X, samples=samples, variants=variants)
```

`samples` describes the rows of every block. `variants` describes the columns
of the current block only. That avoids loading all variant metadata globally and
matching it back to columns after the fact.

---

## Reading one matrix

Use `read` when the analysis needs one matrix in memory.

```python
X = genoio.read("data/chr22_hg38", format="plink2")
```

The returned array has samples on rows and variants on columns:

```python
n_samples, n_variants = X.shape
```

Ask for metadata when you need to preserve row and column labels:

```python
X, samples, variants = genoio.read(
    "data/chr22_hg38",
    format="plink2",
    return_samples=True,
    return_variants=True,
)
```

!!! note "Matrix orientation"
    `genoio` always returns matrices with samples as rows and variants as
    columns. This is true for VCF, PLINK1, PLINK2, dense reads, sparse reads,
    and block reads.

---

## Opening a dataset

Use `open` when code will reuse the same source for metadata, whole reads, or
blocks.

```python
ds = genoio.open("data/chr22_hg38.pgen")

samples = ds.samples()
variants = ds.variants()
X = ds.read()
```

PLINK sources can be passed as a shared prefix or as one member file. For
example, `data/chr22_hg38`, `data/chr22_hg38.pgen`, and
`data/chr22_hg38.psam` all resolve to the same PLINK2 dataset when the
companion files are present.

Pass `format=...` when a prefix could resolve to more than one format:

```python
ds = genoio.open("data/chr22_hg38", format="plink1")
```

---

## Streaming blocks

Block reads yield consecutive retained variants after filtering. A block size
of `5_000` means up to 5,000 variants that passed the filter, not 5,000 raw
source records.

```python
rare = genoio.maf(max=0.05) & genoio.missing_rate(max=0.1)

for X, variants in genoio.blocks(
    "data/chr22_hg38.vcf.gz",
    5_000,
    variants=rare,
    return_variants=True,
):
    run_association_scan(X, variants=variants)
```

The top-level `genoio.blocks(...)` function is a convenience wrapper around
`genoio.open(...).blocks(...)`.

---

## Missing data

Dense reads support three missing-data policies:

```python
genoio.read(path, missing="nan")     # default for dense reads
genoio.read(path, missing="raise")   # fail if retained calls are missing
genoio.read(path, missing="impute")  # per-variant mean imputation
```

Sparse reads currently require `missing="raise"` because this release does not
store sparse missing-value masks.

---

## Sparse matrices

Use `sparse=True` for SciPy CSC output, or `sparse="csr"` for CSR output.

```python
X_csc = genoio.read("data/chr22_hg38", format="plink1", sparse=True)
X_csr = genoio.read("data/chr22_hg38", format="plink1", sparse="csr")
```

Sparse genotype columns are oriented to the minor allele by default to reduce
stored nonzeros.

---

## Haplotype rows

Phased VCF genotypes can be read as haplotype rows:

```python
H = genoio.read("phased.vcf.gz", kind="haplo")
```

Each retained sample contributes two output rows. Haplotype reads require
phased diploid genotypes in retained variants. PLINK1 and PLINK2 haplotype
reads are not implemented in this release.
