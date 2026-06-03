# genoio

`genoio` reads VCF, PLINK1, and PLINK2 genotype files into Python matrices.
The Python API resolves sources and assembles results. Rust readers parse
records, apply filters, and build matrices.

```python
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()

for X, variants in ds.blocks(10_000, return_variants=True):
    # X has shape (samples, variants_in_this_block).
    run_association_scan(X, samples=samples, variants=variants)
```

Three constructors resolve supported sources:

- [`vcf`](api/reading.md#genoio.vcf) for VCF/BCF files.
- [`bfile`](api/reading.md#genoio.bfile) for PLINK1 `.bed/.bim/.fam` file sets.
- [`pfile`](api/reading.md#genoio.pfile) for PLINK2 `.pgen/.pvar/.psam` file sets.

Each constructor returns a reusable [`Dataset`](api/reading.md#genoio.Dataset)
with `read`, `blocks`, `samples`, and `variants` methods.

Dense reads return NumPy arrays with shape `(samples, variants)`. Sparse reads
return SciPy sparse matrices. Metadata is returned as Polars DataFrames.

See [Getting started](guide.md) for the basic workflow, or [Filtering](filtering.md)
for the expression system used to select variants while reading.
