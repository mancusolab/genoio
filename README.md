# genoio

`genoio` is a genotype matrix IO layer for statistical genetics tools. It gives
pipeline authors one Python API for VCF/BCF, PLINK1, and PLINK2 inputs, with
Rust readers handling parsing, filtering, and matrix construction.

The package is aimed at libraries and workflows that need stable matrix and
metadata contracts, not at one-off genotype analysis. It is useful when an
association, QTL, fine-mapping, or simulation tool needs to ingest several
genotype formats without owning file-format parsing itself.

## Association Workflow

The common pattern is: read sample metadata once, align phenotypes and
covariates to those samples, then stream genotype blocks with matching variant
metadata.

```python
import polars as pl
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()

phenotypes = pl.read_csv("phenotypes.tsv", separator="\t")
covariates = pl.read_csv("covariates.tsv", separator="\t")

design = (
    # Start from genoio sample order so y and C match genotype matrix rows.
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

Extra individuals in phenotype or covariate tables are ignored because the join
starts from genotype samples. If you pass `samples=...` to `read` or `blocks`,
build the design from the returned filtered sample frame:

```python
for X, block_samples, variants in ds.blocks(
    5_000,
    samples=["sample_1", "sample_2", "sample_3"],
    variants=variant_filter,
    return_samples=True,
    return_variants=True,
):
    design = (
        block_samples.select("iid")
        .join(phenotypes.select("iid", "expression"), on="iid", how="left")
        .join(covariates.select("iid", "age", "sex", "PC1", "PC2"), on="iid", how="left")
    )
    y = design["expression"].to_numpy()
    C = design.select("age", "sex", "PC1", "PC2").to_numpy()
    run_association_scan(X, y, covariates=C, samples=block_samples, variants=variants)
```

## What Genoio Provides

`genoio` exposes a small set of contracts that downstream tools can build on:

- genotype matrices have samples on rows and variants on columns
- dense reads return NumPy arrays
- sparse reads return SciPy sparse matrices
- sample and variant metadata are Polars DataFrames
- returned sample rows match matrix rows
- returned variant rows match matrix columns
- filters preserve source order among retained variants
- sample keep lists preserve source order among retained samples
- filters that retain no variants return shape `(n_samples, 0)`, not an error

The core variant metadata schema is:

| Column | Meaning |
|---|---|
| `chrom` | Chromosome or contig label. |
| `pos` | 1-based variant position. |
| `id` | Variant ID from the source. |
| `a0` | Allele counted as dosage `0` in returned matrices. |
| `a1` | Allele counted as dosage `1` or `2` in returned matrices. |

`a1` is the counted allele. It is not guaranteed to be the VCF ALT allele after
format normalization or sparse minor-allele orientation.

## Sources

Use the constructor that matches the file format:

```python
vcf_ds = genoio.vcf("cohort.vcf.gz")
bed_ds = genoio.bfile("cohort")       # .bed/.bim/.fam
pgen_ds = genoio.pfile("cohort.pgen") # .pgen/.pvar/.psam
```

Each constructor returns a reusable `Dataset` with:

- `read(...)` for one dense or sparse matrix
- `blocks(...)` for streaming retained variant blocks
- `samples()` for sample metadata
- `variants()` for variant metadata

PLINK sources can be passed as a shared prefix or as one member file. The
constructor determines the file-set type, so same-stem files from other formats
do not affect resolution.

## Filtering

Filters are serializable expression objects. Python builds the expression, and
Rust evaluates it while reading records.

```python
rare = (
    genoio.chrom("22")
    & genoio.maf(max=0.01)
    & genoio.missing_rate(max=0.1)
)

for X, variants in ds.blocks(10_000, variants=rare, return_variants=True):
    run_association_scan(X, y, covariates=C, samples=samples, variants=variants)
```

Expressions compose with Python operators:

```python
genoio.chrom("22") & genoio.snp()
genoio.maf(max=0.01) | genoio.id_in(["rs123", "rs456"])
~genoio.missing_rate(max=0.1)
```

Metadata predicates such as `chrom`, `region`, `id_in`, `snp`, `biallelic`, and
`qual` can often be applied before genotype decoding. Genotype-stat predicates
such as `maf`, `mac`, `missing_rate`, and `polymorphic` require decoded
genotypes for candidate variants.

For compressed VCF/BCF sources, concrete `region(...)` filters use `.tbi` or
`.csi` indexes when present. Unindexed compressed region reads are rejected
instead of silently scanning the whole file.

## Missing Data, Sparse Reads, and Haplotypes

Dense reads support explicit missing-data policies:

```python
ds.read(missing="nan")     # default for dense reads
ds.read(missing="raise")   # fail if retained calls are missing
ds.read(missing="impute")  # per-variant mean imputation
```

Sparse reads currently require `missing="raise"` because this release does not
store sparse missing-value masks.

```python
X_csc = ds.read(sparse=True)
X_csr = ds.read(sparse="csr")
```

Phased VCF genotypes can be read as haplotype rows:

```python
H = genoio.vcf("phased.vcf.gz").read(kind="haplo")
```

Each retained diploid sample contributes two haplotype rows. PLINK1 and PLINK2
haplotype reads are not implemented in this release.

## Format Support

| Format | Inputs | Genotype reads | Haplotype reads | Notes |
|---|---|---:|---:|---|
| VCF/BCF | `.vcf`, `.vcf.gz`, `.bcf` | yes | phased VCF only | Indexed region filters use `.tbi` or `.csi` when available. |
| PLINK1 | `.bed` + `.bim` + `.fam` | yes | no | Variant-major BED files are supported. |
| PLINK2 | `.pgen` + `.pvar` + `.psam` | yes | no | Biallelic hard-call PGEN records are supported. |

Current limitations:

- PLINK2 support is limited to biallelic hard-call records. Dosage tracks are
  not implemented yet.
- Sparse reads do not preserve missing-value masks.
- Haplotype reads are currently VCF-only.
- Region pushdown is implemented for concrete indexed VCF/BCF region filters,
  not for arbitrary filter expressions.

## Performance

The Rust backend avoids Python loops over variants and samples. For block reads,
matrix construction happens in Rust and crosses the Python boundary once per
block.

Benchmark scripts are included for local comparisons against `cyvcf2`,
`pandas_plink`, and `pgenlib` when those packages are installed:

```bash
python -m maturin develop --release
python scripts/benchmark_vcf.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink2.py --max-variants 1000 --repeats 3
```

See [docs/performance.md](docs/performance.md) for the current local benchmark
context and commands.

## Installation

From a local checkout:

```bash
pip install -e ".[dev,docs]"
python -m maturin develop --release
```

For a quick editable development build, omit `--release`:

```bash
python -m maturin develop
```

## Development

Run the full local verification suite before publishing or opening a PR:

```bash
make verify
```

`make verify` builds the extension, runs Rust formatting, Clippy, Rust tests,
Pyright, Python tests, MkDocs strict mode, and Rust documentation with warnings
as errors.

More documentation lives in [docs/](docs/index.md).
