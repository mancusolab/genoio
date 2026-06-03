# genoio

`genoio` reads VCF, PLINK1, and PLINK2 genotype files into Python matrices.
The Python API resolves sources and assembles results. Rust readers parse
records, apply filters, and build matrices.

```python
import genoio

ds = genoio.open("data/chr22_hg38", format="plink2")
samples = ds.samples()

for X, variants in ds.blocks(10_000, return_variants=True):
    # X has shape (samples, variants_in_this_block).
    run_association_scan(X, samples=samples, variants=variants)
```

Three entry points cover the common workflows:

- [`open`](api/reading.md#genoio.open) resolves a source once and returns a
  reusable dataset.
- [`read`](api/reading.md#genoio.read) reads one matrix now.
- [`blocks`](api/reading.md#genoio.blocks) streams variant blocks with matrix
  columns and variant metadata kept in the same order.

Dense reads return NumPy arrays with shape `(samples, variants)`. Sparse reads
return SciPy sparse matrices. Metadata is returned as Polars DataFrames.

See [Getting started](guide.md) for the basic workflow, or [Filtering](filtering.md)
for the expression system used to select variants while reading.
