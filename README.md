# genoio

`genoio` is a Python package for reading genotype matrices from common genetics
file formats. It exposes a small Python API and uses Rust readers for format
parsing.

The package is aimed at workflows that need the same matrix API across VCF,
PLINK1, and PLINK2 inputs.

## Supported Formats

| Format | Inputs | Genotype reads | Haplotype reads | Notes |
|---|---|---:|---:|---|
| VCF/BCF | `.vcf`, `.vcf.gz`, `.bcf` | yes | phased VCF only | Indexed region filters use `.tbi` or `.csi` when available. |
| PLINK1 | `.bed` + `.bim` + `.fam` | yes | no | Variant-major BED files are supported. |
| PLINK2 | `.pgen` + `.pvar` + `.psam` | yes | no | Biallelic hard-call PGEN records are supported. |

Dense genotype reads return a NumPy array with shape `(samples, variants)`.
Sparse reads return a SciPy sparse matrix. Metadata is returned as Polars
DataFrames.

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

## Quick Start

```python
import genoio

ds = genoio.open("data/chr22_hg38.vcf.gz")

X = ds.read()
print(X.shape)
```

Read a PLINK prefix the same way:

```python
X = genoio.read("data/chr22_hg38", format="plink2")
```

You can also pass a member file such as `data/chr22_hg38.pgen`; `genoio`
resolves the companion files from the shared prefix.

## Metadata

```python
samples = genoio.samples("data/chr22_hg38", format="plink2")
variants = genoio.variants("data/chr22_hg38", format="plink2")
```

`samples` and `variants` are Polars DataFrames. Variant metadata includes source
alleles, normalized `a0`/`a1` alleles, optional `qual`, and genotype-derived
statistics when a read computes them for filtering.

## Filters

Variant filters are serializable expressions evaluated by the Rust backend.
They can be combined with `&`, `|`, and `~`.

```python
import genoio

variant_filter = (
    genoio.region("22:20000000-21000000")
    & genoio.qual(min=20)
    & genoio.maf(max=0.05)
)

X, variants = genoio.read(
    "data/chr22_hg38.vcf.gz",
    variants=variant_filter,
    return_variants=True,
)
```

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

For compressed VCF/BCF sources, concrete region filters use a `.tbi` or `.csi`
index when present. Unindexed compressed region reads are rejected because a
full scan would be misleading for region-specific workflows.

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

Blocks are defined in retained-variant order after filtering.

## Haplotype Reads

Phased VCF genotypes can be read as haplotype rows:

```python
H = genoio.read("phased.vcf.gz", kind="haplo")
```

Each retained sample contributes two output rows. Haplotype reads require phased
diploid genotypes in retained variants. PLINK1 and PLINK2 haplotype reads are
not implemented in this release.

## Benchmarks

Benchmark scripts are available in `scripts/`:

```bash
python scripts/benchmark_vcf.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink2.py --max-variants 1000 --repeats 3
```

Optional comparison packages are used when installed:

- `cyvcf2` for VCF
- `pandas_plink` for PLINK1
- `pgenlib` for PLINK2

Run benchmarks against a release build:

```bash
python -m maturin develop --release
```

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
