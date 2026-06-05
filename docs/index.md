# Getting started

`genoio` reads VCF, PLINK1, PLINK2, and BGEN genotype files into Python
matrices. The Python API resolves sources and assembles results. Rust readers
parse records, apply filters, and build matrices.

## Installation

Install the development version from this repository:

```bash
pip install git+https://github.com/mancusolab/genoio.git
```

For local development, use the project environment and build the Rust extension
in place:

```bash
make build-dev
```

## Quick example

```python
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()
y = load_phenotype_vector(samples["iid"])
C = load_covariates(samples["iid"])

for X, variants in ds.iter_blocks(10_000, return_variants=True):
    # X has shape (samples, variants_in_this_block).
    # `y` and `C` must be aligned to the rows described by `samples`.
    association_scan(X, y, C, variants=variants)
```

Four constructors resolve supported sources:

- [`vcf`](api/reading.md#genoio.vcf) for VCF/BCF files.
- [`bfile`](api/reading.md#genoio.bfile) for PLINK1 `.bed/.bim/.fam` file sets.
- [`pfile`](api/reading.md#genoio.pfile) for PLINK2 `.pgen/.pvar[.zst]/.psam` file sets.
- [`bgen`](api/reading.md#genoio.bgen) for BGEN `.bgen` files.

Each constructor returns a reusable [`Dataset`](api/reading.md#genoio.Dataset)
with [`read`](api/reading.md#genoio.Dataset.read),
[`iter_blocks`](api/reading.md#genoio.Dataset.iter_blocks),
[`iter_regions`](api/reading.md#genoio.Dataset.iter_regions),
[`samples`](api/reading.md#genoio.Dataset.samples), and
[`variants`](api/reading.md#genoio.Dataset.variants) methods.

Dense reads return NumPy arrays with shape `(samples, variants)`. Sparse reads
return SciPy sparse matrices. Metadata is returned as Polars DataFrames.

## Next steps

Use [Examples](examples/index.md) for GWAS and cis-eQTL scan sketches. Read
[Filtering](filtering.md) for variant selection and region pushdown,
[Format support](formats.md) for source-specific behavior, or
[API Reading](api/reading.md) for matrix options, missing-data handling, sparse
output, and iterator contracts.
