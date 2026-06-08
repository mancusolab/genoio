# genoio

`genoio` reads common statistical-genetics file formats into Python matrices.
It exposes one API for VCF/BCF, PLINK1, PLINK2, and BGEN inputs, with Rust
readers underneath for parsing, filtering, and matrix construction.

Use it when downstream code needs predictable matrix contracts across formats:

- samples are rows and variants are columns
- sample metadata rows match matrix rows
- variant metadata rows match matrix columns
- block and region reads return metadata aligned to each matrix chunk

The full documentation is at
[mancusolab.github.io/genoio](https://mancusolab.github.io/genoio).

## Installation

Install the development version from GitHub:

```bash
pip install git+https://github.com/mancusolab/genoio.git
```

From a local checkout:

```bash
pip install .
```

For development, build the Rust extension in the project environment:

```bash
make build-dev
```

## Quick Example

```python
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()

y = load_phenotype_vector(samples["iid"])
C = load_covariates(samples["iid"])

rare = genoio.maf(max=0.01) & genoio.missing_rate(max=0.1)

for X, variants in ds.iter_blocks(10_000, variants=rare, return_variants=True):
    # X has shape (samples, variants_in_this_block).
    association_scan(X, y, C, variants=variants)
```

`read(...)` loads one matrix. `iter_blocks(...)` streams retained variant
chunks. `iter_regions(...)` reads one result per genomic interval.

```python
X, samples, variants = ds.read(return_samples=True, return_variants=True)
```

## Supported Inputs

| Format | Constructor | Inputs | Notes |
|---|---|---|---|
| VCF/BCF | `genoio.vcf(...)` | `.vcf`, `.vcf.gz`, `.bcf` | Hardcall reads by default; dense `FORMAT/DS` dosage reads are supported. |
| PLINK1 | `genoio.bfile(...)` | `.bed/.bim/.fam` | Variant-major BED hardcall reads. |
| PLINK2 | `genoio.pfile(...)` | `.pgen/.pvar[.zst]/.psam` | Hardcalls, dense genotype dosage, and explicit phased haplotype records. |
| BGEN | `genoio.bgen(...)` | `.bgen` plus optional `.sample` | Dense dosage-backed reads for supported BGEN v1.2+ Layout 2 records. |

Dense reads return NumPy arrays. Sparse reads return SciPy sparse matrices.
Metadata is returned as Polars DataFrames.

For source-specific behavior and current limitations, read
[Format support](docs/formats.md).

## Filtering And Read Options

Filters are serializable expressions:

```python
variants = (
    genoio.chrom("22")
    & genoio.region("22:20000000-21000000")
    & genoio.maf(max=0.05)
)
X, variants = ds.read(variants=variants, return_variants=True)
```

By default, genotype reads return hardcall A1 allele counts. Use
`dosage="dosage"` for source dosage or probability values when the format
supports them, `sparse=True` for SciPy CSC output, and `kind="haplo"` for
haplotype rows.

```python
X = genoio.bgen("cohort.bgen").read(dosage="dosage")
H = genoio.pfile("phased").read(kind="haplo", dosage="hardcall")
```

See [Filtering](docs/filtering.md) for filter expressions and pushdown rules,
[Reading](docs/api/reading.md) for matrix options, and
[Format support](docs/formats.md) for source-specific limitations.

## Development

Useful targets:

```bash
make help
make build-dev
make build-wheel
make verify
```

The project uses Python for the public API and Rust for file parsing and matrix
construction. Build configuration lives in `pyproject.toml`; Rust crates live
under `rust/`.

## License

`genoio` is distributed under the MIT license. See [LICENSE](LICENSE).
