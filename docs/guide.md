# Getting started

Most downstream tools want one stable sample axis and many chunks of variant
columns. `genoio` makes that pattern explicit: resolve the source once, read
sample metadata once, then process matrix blocks with their matching variant
metadata.

```python
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()
y = load_phenotype_vector(samples["iid"])

for X, variants in ds.blocks(10_000, return_variants=True):
    run_association_scan(X, y, samples=samples, variants=variants)
```

`samples` describes the rows of every block. `variants` describes the columns
of the current block only. That avoids loading all variant metadata globally and
matching it back to columns after the fact.

---

## Metadata frames

Sample and variant metadata are returned as Polars DataFrames. `samples()` is
source ordered and uses `iid` as the sample identifier for phenotype or
covariate joins. Variant frames returned from `read` or `blocks` are ordered to
match matrix columns after filtering. See [Reading](api/reading.md#metadata-frames)
for the metadata schemas.

---

## Reading one matrix

Use a dataset's `read` method when the analysis needs one matrix in memory.

```python
X = genoio.pfile("data/chr22_hg38").read()
```

The returned array has samples on rows and variants on columns:

```python
n_samples, n_variants = X.shape
```

Ask for metadata when you need to preserve row and column labels:

```python
X, samples, variants = genoio.pfile("data/chr22_hg38").read(
    return_samples=True,
    return_variants=True,
)
```

!!! note "Matrix orientation"
    `genoio` always returns matrices with samples as rows and variants as
    columns. This is true for VCF, PLINK1, PLINK2, dense reads, sparse reads,
    and block reads.

---

## Constructing a dataset

Use the constructor that matches the source format:

```python
vcf_ds = genoio.vcf("data/chr22_hg38.vcf.gz")
bfile_ds = genoio.bfile("data/chr22_hg38")
pfile_ds = genoio.pfile("data/chr22_hg38.pgen")

samples = pfile_ds.samples()
variants = pfile_ds.variants()
X = pfile_ds.read()
```

PLINK sources can be passed as a shared prefix or as one member file. For
example, `data/chr22_hg38`, `data/chr22_hg38.pgen`, and
`data/chr22_hg38.psam` all resolve to the same PLINK2 dataset when passed to
`pfile(...)` and the companion files are present.

---

## Streaming blocks

Block reads yield consecutive retained variants after filtering. A block size
of `5_000` means up to 5,000 variants that passed the filter, not 5,000 raw
source records.

```python
rare = genoio.maf(max=0.05) & genoio.missing_rate(max=0.1)

ds = genoio.vcf("data/chr22_hg38.vcf.gz")
samples = ds.samples()
y = load_phenotype_vector(samples["iid"])

for X, variants in ds.blocks(
    5_000,
    variants=rare,
    return_variants=True,
):
    run_association_scan(X, y, samples=samples, variants=variants)
```

---

## Missing data

Dense reads support three missing-data policies:

```python
ds.read(missing="nan")     # default for dense reads
ds.read(missing="raise")   # fail if retained calls are missing
ds.read(missing="impute")  # per-variant mean imputation
```

Sparse reads currently require `missing="raise"` because this release does not
store sparse missing-value masks.

---

## Sparse matrices

Use `sparse=True` for SciPy CSC output, or `sparse="csr"` for CSR output.

```python
X_csc = genoio.bfile("data/chr22_hg38").read(sparse=True)
X_csr = genoio.bfile("data/chr22_hg38").read(sparse="csr")
```

Sparse genotype columns are oriented to the minor allele by default to reduce
stored nonzeros.

---

## Haplotype rows

Phased VCF genotypes can be read as haplotype rows:

```python
H = genoio.vcf("phased.vcf.gz").read(kind="haplo")
```

Each retained sample contributes two output rows. Haplotype reads require
phased diploid genotypes in retained variants. PLINK1 and PLINK2 haplotype
reads are not implemented in this release.
