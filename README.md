# genoio

Genotype matrix IO for Python. One API for VCF, PLINK1, and PLINK2; Rust
readers for parsing, filtering, and matrix construction.

## Example

```python
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()
y = load_phenotype_vector(samples["iid"])

for X, variants in ds.blocks(10_000, return_variants=True):
    # X has shape (samples, variants_in_this_block).
    # `variants` describes the columns of X in the same order.
    # `y` must be aligned to the rows described by `samples`.
    run_association_scan(X, y, samples=samples, variants=variants)
```

This is the main use case: samples are fixed, variants stream by block, and
each block carries the metadata needed to interpret its columns. Downstream
tools do not need to load all variant metadata globally and match it back to
matrix columns after the fact.

`vcf(...)`, `bfile(...)`, and `pfile(...)` construct reusable datasets:

- `vcf(...)`: VCF/BCF source.
- `bfile(...)`: PLINK1 `.bed/.bim/.fam` source.
- `pfile(...)`: PLINK2 `.pgen/.pvar/.psam` source.

Each dataset has `read(...)`, `blocks(...)`, `samples()`, and `variants()`.
Block reads stream variant chunks while keeping each block's matrix and variant
metadata aligned.

Dense reads return NumPy arrays with shape `(samples, variants)`. Sparse reads
return SciPy sparse matrices. Metadata is returned as Polars DataFrames.

## What genoio Does

`genoio` turns genotype files into matrices:

- rows are samples
- columns are variants
- values are diploid allele counts `0`, `1`, or `2`
- missing calls are handled by an explicit policy
- sample and variant metadata can be returned with the matrix

The package is meant for analysis pipelines that need a consistent matrix API
across file formats. It is not a full variant annotation framework.

```python
X_vcf = genoio.vcf("data/chr22_hg38.vcf.gz").read()
X_bed = genoio.bfile("data/chr22_hg38").read()
X_pgen = genoio.pfile("data/chr22_hg38").read()
```

PLINK sources can be passed as a prefix or as one member file such as
`data/chr22_hg38.pgen`; `genoio` resolves companion files from the shared
prefix. The constructor chooses the format, so unrelated files with the same
prefix do not affect resolution.

## Filtering Expressions

Filters are serializable expression objects. Python builds the expression; Rust
evaluates it while reading records.

```python
rare_high_quality = (
    genoio.region("22:20000000-21000000")
    & genoio.qual(min=20)
    & genoio.maf(max=0.05)
)

ds = genoio.vcf("data/chr22_hg38.vcf.gz")
samples = ds.samples()
y = load_phenotype_vector(samples["iid"])

for X, variants in ds.blocks(
    5_000,
    variants=rare_high_quality,
    return_variants=True,
):
    run_association_scan(X, y, samples=samples, variants=variants)
```

Expressions compose with Python operators:

```python
genoio.chrom("22") & genoio.snp()
genoio.maf(max=0.01) | genoio.id_in(["rs123", "rs456"])
~genoio.missing_rate(max=0.1)
```

There are two kinds of predicates:

- **Metadata predicates** use fields already present in the source record:
  chromosome, position, ID, REF/ALT structure, and `QUAL`.
- **Genotype predicates** require decoding retained genotypes first: MAF, MAC,
  missing rate, and polymorphism.

That distinction matters for speed. Metadata predicates can drop records before
matrix decoding. A concrete VCF/BCF `region(...)` can also use a `.tbi` or
`.csi` index when the compressed source has one. Genotype predicates must decode
the candidate variant before deciding whether to keep it.

Available filters:

- `chrom("22")`
- `region("22:20000000-21000000")`
- `snp()`
- `biallelic()`
- `qual(min=..., max=...)`
- `maf(min=..., max=...)`
- `mac(min=..., max=...)`
- `missing_rate(max=...)`
- `polymorphic()`
- `id_in([...])`

For compressed VCF/BCF sources, unindexed region reads are rejected. This avoids
silently doing a full compressed-file scan when the user asked for a region.

## Performance

The Rust backend avoids Python loops over variants and samples. For block reads,
matrix construction happens in Rust and crosses the Python boundary once per
block.

On a local Apple Silicon M1 Mac (`arm64`, Python 3.11), release builds on the
first 1,000 variants show:

| Read | genoio | comparison |
|---|---:|---:|
| VCF | 0.106 s | `cyvcf2` matrix construction 0.238 s |
| PLINK1 | 0.397 s | `pandas_plink` matrix construction 2.438 s |
| PLINK2 | 0.026 s | `pgenlib` matrix construction 0.008 s |

These numbers are workload- and machine-dependent. They are included to set
expectations, not as a universal benchmark. The main point is that `genoio`
builds matrices without Python-level variant iteration.

The benchmark fixture is a local `data/chr22_hg38` dataset with 3,202 samples.
It isn't distributed with the repository. It comes from the PLINK 2
[1000 Genomes phase 3 hg38 resources](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg):
the chromosome 22 PLINK 2 files were used as the source, then converted with
`plink2` to VCF and PLINK1 `.bed/.bim/.fam` files so each script reads the same
underlying genotypes.

Run benchmarks on your machine:

```bash
python -m maturin develop --release
python scripts/benchmark_vcf.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink2.py --max-variants 1000 --repeats 3
```

Optional comparison packages are used when installed: `cyvcf2` for VCF,
`pandas_plink` for PLINK1, and `pgenlib` for PLINK2.

## Metadata

```python
ds = genoio.pfile("data/chr22_hg38")

samples = ds.samples()
variants = ds.variants()
```

`samples` and `variants` are Polars DataFrames. Variant metadata identifies
matrix columns and dosage orientation.

See [Reading](docs/api/reading.md#metadata-frames) for the metadata schemas.

Return metadata alongside a whole-matrix read:

```python
X, sample_metadata, variant_metadata = genoio.pfile("data/chr22_hg38").read(
    return_samples=True,
    return_variants=True,
)
```

For block reads, prefer reading samples once and returning variants per block:

```python
ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()
y = load_phenotype_vector(samples["iid"])

for X, variants in ds.blocks(10_000, return_variants=True):
    run_association_scan(X, y, samples=samples, variants=variants)
```

## Samples and Missing Data

Pass sample IDs with `samples=...` to keep a subset of rows:

```python
X, sample_metadata = genoio.bfile("data/chr22_hg38").read(
    samples=["HG00096", "HG00097"],
    return_samples=True,
)
```

Rows remain in source order, not in the order of the requested list. Duplicate
sample IDs in the request are rejected.

Dense reads support three missing-data policies:

```python
ds.read(missing="nan")     # default for dense reads
ds.read(missing="raise")   # fail if retained calls are missing
ds.read(missing="impute")  # per-variant mean imputation
```

Sparse reads currently require `missing="raise"` because this release does not
store sparse missing-value masks.

## Sparse and Haplotype Reads

```python
X = genoio.bfile("data/chr22_hg38").read(sparse=True)
X_csr = genoio.bfile("data/chr22_hg38").read(sparse="csr")
```

Sparse matrices are built with variants as columns and samples as rows. By
default, sparse genotype columns are oriented to the minor allele to reduce
stored nonzeros.

Phased VCF genotypes can be read as haplotype rows:

```python
H = genoio.vcf("phased.vcf.gz").read(kind="haplo")
```

Each retained sample contributes two output rows. Haplotype reads require phased
diploid genotypes in retained variants. PLINK1 and PLINK2 haplotype reads are
not implemented in this release.

## Installation

From a local checkout:

```bash
pip install -e ".[dev]"
python -m maturin develop --release
```

For a quick editable development build, omit `--release`:

```bash
python -m maturin develop
```

## Format Support

| Format | Inputs | Genotype reads | Haplotype reads | Notes |
|---|---|---:|---:|---|
| VCF/BCF | `.vcf`, `.vcf.gz`, `.bcf` | yes | phased VCF only | Indexed region filters use `.tbi` or `.csi` when available. |
| PLINK1 | `.bed` + `.bim` + `.fam` | yes | no | Variant-major BED files are supported. |
| PLINK2 | `.pgen` + `.pvar` + `.psam` | yes | no | Biallelic hard-call PGEN records are supported. |

## Current Limitations

- PLINK2 support is limited to biallelic hard-call records. Dosage tracks are
  not implemented yet.
- Sparse reads do not preserve missing-value masks.
- Haplotype reads are currently VCF-only.
- Region pushdown is implemented for concrete indexed VCF/BCF region filters,
  not for arbitrary filter expressions.

## Development

```bash
pip install -e ".[dev]"
pre-commit install
python -m maturin develop
pytest -q
cargo test --manifest-path rust/Cargo.toml --workspace
```

Run the full hygiene checks:

```bash
pre-commit run --all-files
```
