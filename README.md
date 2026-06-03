# genoio

`genoio` reads genotype files into Python matrices with one API across VCF,
PLINK1, and PLINK2. The Python layer handles source resolution and result
assembly; Rust readers do the format parsing, filtering, and matrix
construction.

The common pattern is to resolve a source once, keep sample metadata fixed, and
stream variant blocks with the metadata for each block:

```python
import genoio

ds = genoio.open("data/chr22_hg38", format="plink2")
samples = ds.samples()

for X, variants in ds.blocks(10_000, return_variants=True):
    # X has shape (samples, variants_in_this_block).
    # `variants` describes the columns of X in the same order.
    run_association_scan(X, samples=samples, variants=variants)
```

Dense reads return NumPy arrays. Sparse reads return SciPy sparse matrices.
Metadata comes back as Polars DataFrames when requested.

`genoio` is useful when you want to write analysis code once and point it at a
VCF, PLINK1, or PLINK2 source without changing the rest of the pipeline.

## What It Does

`genoio` turns genotype sources into matrices:

- rows are samples
- columns are variants
- values are diploid allele counts `0`, `1`, or `2`
- missing calls are handled by an explicit policy
- sample and variant metadata can be returned with the matrix

It also supports:

- variant filters such as `maf(max=0.05)`, `qual(min=20)`, and genomic regions
- sample keep lists
- dense or sparse output
- block-wise reads for larger scans
- phased VCF haplotype reads

The package does not try to be a full variant annotation framework. It focuses
on the part many downstream tools need first: getting a correctly oriented
genotype matrix and enough metadata to keep rows and columns interpretable.

## Install

From a local checkout:

```bash
pip install -e ".[dev]"
python -m maturin develop --release
```

For a quick editable development build, omit `--release`:

```bash
python -m maturin develop
```

Use release builds for benchmarks and real performance comparisons.

## Reading Genotypes

Open a dataset once when downstream code will reuse the same source:

```python
import genoio

ds = genoio.open("data/chr22_hg38.vcf.gz")
samples = ds.samples()

for X, variants in ds.blocks(5_000, return_variants=True):
    analyze_block(X, samples=samples, variants=variants)
```

Use `read(...)` when you only need one matrix:

```python
X_vcf = genoio.read("data/chr22_hg38.vcf.gz")
X_bed = genoio.read("data/chr22_hg38", format="plink1")
X_pgen = genoio.read("data/chr22_hg38", format="plink2")
```

You can pass a PLINK prefix or one member file such as
`data/chr22_hg38.pgen`; `genoio` resolves companion files from the shared
prefix.

## Metadata

```python
samples = genoio.samples("data/chr22_hg38", format="plink2")
variants = genoio.variants("data/chr22_hg38", format="plink2")
```

`samples` and `variants` are Polars DataFrames. Variant metadata includes source
alleles, normalized `a0`/`a1` alleles, optional `qual`, and genotype-derived
statistics when a read computes them for filtering.

You can return metadata alongside a matrix:

```python
X, sample_metadata, variant_metadata = genoio.read(
    "data/chr22_hg38",
    format="plink2",
    return_samples=True,
    return_variants=True,
)
```

## Filtering Expressions

Filters are serializable expression objects. Python builds the expression;
Rust evaluates it while reading records.

```python
rare_high_quality = (
    genoio.region("22:20000000-21000000")
    & genoio.qual(min=20)
    & genoio.maf(max=0.05)
)

X, variants = genoio.read(
    "data/chr22_hg38.vcf.gz",
    variants=rare_high_quality,
    return_variants=True,
)
```

Expressions compose with Python operators:

```python
genoio.chrom("22") & genoio.snp()
genoio.maf(max=0.01) | genoio.id_in(["rs123", "rs456"])
~genoio.missing_rate(max=0.1)
```

There are two kinds of predicates:

- metadata predicates use fields already present in the source record, such as
  chromosome, position, ID, REF/ALT structure, and `QUAL`
- genotype predicates require decoding retained genotypes first, such as MAF,
  MAC, missing rate, and polymorphism

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

The Rust backend is designed to avoid Python loops over variants and samples.
For common block reads, matrix construction happens in Rust and crosses the
Python boundary once per block.

On the local `data/chr22_hg38` fixture with 3,202 samples, release builds show:

| Read | genoio | comparison |
|---|---:|---:|
| VCF, first 1,000 variants | ~0.10 s | `cyvcf2` matrix construction ~0.23 s |
| PLINK2, first 100 variants | ~0.003 s | `pgenlib` matrix construction ~0.004 s |
| PLINK2, first 1,000 variants | ~0.02 s | `pgenlib` matrix construction ~0.007 s |
| PLINK2, first 10,000 variants | ~0.30 s | `pgenlib` matrix construction ~0.05 s |

These numbers are workload- and machine-dependent. They are included to set
expectations, not as a universal benchmark. The main point is that `genoio`
does not require Python-level variant iteration to build the matrix.

Run the benchmark scripts on your machine:

```bash
python -m maturin develop --release
python scripts/benchmark_vcf.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink2.py --max-variants 1000 --repeats 3
```

Optional comparison packages are used when installed:

- `cyvcf2` for VCF
- `pandas_plink` for PLINK1
- `pgenlib` for PLINK2

## Block Reads

Use `blocks(...)` when you want to process variants incrementally:

```python
for X in genoio.blocks("data/chr22_hg38", 10_000, format="plink2"):
    print(X.shape)
```

Block reads accept the same read options as `read(...)`:

```python
for X, variants in genoio.blocks(
    "data/chr22_hg38.vcf.gz",
    5_000,
    variants=genoio.maf(max=0.05),
    return_variants=True,
):
    ...
```

Blocks are defined in retained-variant order after filtering. If a filter keeps
only every tenth variant, a block of size `5_000` contains up to `5_000`
retained variants, not `5_000` raw source records.

## Samples

Pass sample IDs with `samples=...` to keep a subset of rows:

```python
X, sample_metadata = genoio.read(
    "data/chr22_hg38",
    format="plink1",
    samples=["HG00096", "HG00097"],
    return_samples=True,
)
```

Rows remain in source order, not in the order of the requested list. Duplicate
sample IDs in the request are rejected.

## Missing Data

Dense reads support three missing-data policies:

```python
genoio.read(path, missing="nan")     # default for dense reads
genoio.read(path, missing="raise")   # fail if retained calls are missing
genoio.read(path, missing="impute")  # per-variant mean imputation
```

Sparse reads currently require `missing="raise"` because this release does not
store sparse missing-value masks.

## Sparse Reads

```python
X = genoio.read("data/chr22_hg38", format="plink1", sparse=True)
X_csr = genoio.read("data/chr22_hg38", format="plink1", sparse="csr")
```

Sparse matrices are built with variants as columns and samples as rows. By
default, sparse genotype columns are oriented to the minor allele to reduce
stored nonzeros.

## Haplotype Reads

Phased VCF genotypes can be read as haplotype rows:

```python
H = genoio.read("phased.vcf.gz", kind="haplo")
```

Each retained sample contributes two output rows. Haplotype reads require phased
diploid genotypes in retained variants. PLINK1 and PLINK2 haplotype reads are
not implemented in this release.

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
