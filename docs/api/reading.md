# Reading

Use these entry points to open genotype sources and read matrices. `vcf`,
`bfile`, `pfile`, and `bgen` identify the input format and return a `Dataset`; the
dataset methods then read dense NumPy arrays, sparse SciPy matrices, metadata
frames, or streaming blocks. Dense and sparse genotype matrices use samples as
rows and variants as columns, so block reads can be concatenated by variant when
sample ordering is unchanged.

Filters may retain no variants. In that case, dense reads return an array with
shape `(n_samples, 0)`, sparse reads return a SciPy sparse matrix with the same
shape, and `return_variants=True` returns an empty variant frame with the normal
schema. Block reads yield no blocks when no variant passes the filter.

By default, `kind="geno"` returns diploid hardcall allele counts: `0`, `1`, or
`2` copies of the A1 allele, with missing calls handled by the selected
missing-data policy. The `dosage` option controls this genotype value source:
`dosage="hardcall"` uses source hard calls, while `dosage="dosage"` requires
dosage-backed values. Dense genotype dosage reads are supported for VCF
`FORMAT/DS`, PLINK2 unphased biallelic dosage records, and BGEN v1.2+ Layout 2
biallelic diploid dosage records. Genotype reads of phased BGEN records
collapse the two source haplotype probabilities to expected diploid A1 dosage.

Use `kind="haplo"` for one row per source haplotype. VCF/BCF haplotype reads
use phased hardcall `FORMAT/GT`. PLINK2 haplotype reads support source-encoded
explicit phased hardcalls with `dosage="hardcall"` and explicit phased full
dosages with `dosage="dosage"`. BGEN haplotype reads support source-encoded
phased BGEN v1.2+ Layout 2 biallelic diploid probabilities with
`dosage="dosage"` and return expected A1 dosage per haplotype row. Sparse
PLINK2 explicit phased hardcall haplotypes are supported. Sparse dosage, sparse
PLINK2 dosage haplotypes, sparse BGEN haplotypes, and PLINK1 dosage are not
implemented yet and raise `genoio.UnsupportedRepresentation`.

Genotype-stat filters such as `maf`, `mac`, and `missing_rate` use the selected
value source, so dosage reads compute those statistics from expected allele
dosages rather than hardcall counts. For haplotype reads, those filters are
evaluated on the collapsed diploid expected A1 dosage for each sample, while the
returned matrix remains haplotype-level.

## Reading one matrix

Use [`Dataset.read(...)`](#genoio.Dataset.read) when the requested matrix fits
in memory.

```python
X = genoio.pfile("cohort").read()
```

Ask for metadata when downstream code needs row and column labels:

```python
X, samples, variants = genoio.pfile("cohort").read(
    return_samples=True,
    return_variants=True,
)
```

The returned matrix always has samples on rows and variants on columns:

```python
n_samples, n_variants = X.shape
```

## Metadata frames

[`Dataset.samples()`](#genoio.Dataset.samples) returns a Polars DataFrame in
source sample order. Genotype reads and blocks return the same sample columns when
`return_samples=True`. BGEN reads require real sample IDs, either embedded in
the `.bgen` file or provided by the same-prefix `.sample` file.

Haplotype reads return one row per haplotype. Their returned sample frames add
`source_sample_index` and `haplotype_index` so each haplotype row can be mapped
back to its original diploid sample.

??? info "Sample metadata schema"
    | Column | Meaning |
    |---|---|
    | `fid` | Family ID when present in the source; otherwise null. |
    | `iid` | Sample ID. Use this column to align phenotypes and covariates. |
    | `father` | Paternal ID from PLINK-style sample metadata, or null. |
    | `mother` | Maternal ID from PLINK-style sample metadata, or null. |
    | `sex` | Source sex code as text when present, or null. |
    | `phenotype` | Source phenotype field as text when present, or null. |
    | `source_sample_index` | Haplotype-read-only source sample row index. |
    | `haplotype_index` | Haplotype-read-only haplotype number within the source sample. |

[`Dataset.variants()`](#genoio.Dataset.variants) returns a Polars DataFrame in
source variant order. [`read(..., return_variants=True)`](#genoio.Dataset.read),
[`iter_blocks(..., return_variants=True)`](#genoio.Dataset.iter_blocks), and
[`iter_regions(..., return_variants=True)`](#genoio.Dataset.iter_regions) return
the same five columns for the retained variants, in matrix-column order. The
schema is format-neutral: `a1` is the allele counted by returned genotype
values, not necessarily the VCF ALT allele.

??? info "Variant metadata schema"
    | Column | Meaning |
    |---|---|
    | `chrom` | Chromosome or contig label. |
    | `pos` | 1-based variant position. |
    | `id` | Variant ID from the source. |
    | `a0` | Allele counted as `0` in returned genotype matrices. |
    | `a1` | Allele counted as `1` or `2` in returned genotype matrices. |

## Block iteration

Use [`Dataset.iter_blocks(...)`](#genoio.Dataset.iter_blocks) for scans that
apply the same operation across many variant columns. A block size of `5_000`
means up to 5,000 retained variants, not 5,000 raw source records.

```python
rare = genoio.maf(max=0.05) & genoio.missing_rate(max=0.1)

for X, variants in ds.iter_blocks(
    5_000,
    variants=rare,
    return_variants=True,
):
    association_scan(X, y, variants=variants)
```

Each block keeps the dataset sample order. Variant metadata returned with a
block describes that block's columns only.

## BGEN region reads

BGEN dosage reads use a same-path bgenix SQLite index for concrete region
filters when one is present. For `cohort.bgen`, `genoio` looks for
`cohort.bgen.bgi`. Without that index, BGEN region filters fall back to the
normal sequential scan.

```python
ds = genoio.bgen("cohort.bgen")
X, variants = ds.read(
    dosage="dosage",
    variants=genoio.region("22:20000000-21000000"),
    return_variants=True,
)
```

If a read only needs the matrix, leave `return_samples=False` and
`return_variants=False`. BGEN matrix-only reads avoid returning metadata frames.

## Region iteration

Use [`Dataset.iter_regions(...)`](#genoio.Dataset.iter_regions) when each
genomic interval is an analysis unit, as in cis-eQTL scans. Each yielded item is
`(region, result)`, where `region` is the original filter object and `result`
follows the normal [`read(...)`](#genoio.Dataset.read) return contract.

```python
regions = [
    genoio.region("22:20000000-21000000"),
    genoio.region("22:22000000-23000000"),
]

for region, (X, variants) in ds.iter_regions(regions, return_variants=True):
    scan_region(region, X, variants)
```

[`iter_regions(...)`](#genoio.Dataset.iter_regions) supplies `variants` from the
`regions` argument, so it rejects a separate `variants=` read option.

## Missing data

Dense reads support three missing-data policies:

```python
ds.read(missing="nan")     # default for dense reads
ds.read(missing="raise")   # fail if retained calls are missing
ds.read(missing="impute")  # per-variant mean imputation
```

Sparse reads currently require `missing="raise"` because this release does not
store sparse missing-value masks.

## Sparse matrices

Use `sparse=True` for SciPy CSC output, or `sparse="csr"` for CSR output.

```python
X_csc = genoio.bfile("cohort").read(sparse=True)
X_csr = genoio.bfile("cohort").read(sparse="csr")
```

Sparse genotype columns are oriented to the minor allele by default to reduce
stored nonzeros.

## Haplotype rows

Phased VCF genotypes can be read as haplotype rows:

```python
H = genoio.vcf("phased.vcf.gz").read(kind="haplo")
```

PLINK2 explicit phased hardcall records can be read with the default haplotype
value source:

```python
H = genoio.pfile("phased_hardcall").read(kind="haplo", dosage="hardcall")
```

PLINK2 explicit phased full-dosage records are dosage-backed haplotype rows:

```python
H = genoio.pfile("phased_dosage").read(kind="haplo", dosage="dosage")
```

BGEN haplotype reads require BGEN v1.2+ Layout 2 phased biallelic diploid
probability records and return source-encoded phased dosage as expected A1
dosage per haplotype row:

```python
H = genoio.bgen("phased.bgen").read(kind="haplo", dosage="dosage")
```

Each retained sample contributes two output rows. Haplotype reads require
phased diploid records in retained variants. PLINK2 hardcall haplotypes require
explicit phased hardcall records; PLINK2 dosage haplotypes require explicit
phased full dosage records. BGEN haplotypes require phased probability records.
PLINK2 hardcall haplotypes can be read sparsely when retained calls are
non-missing. PLINK1 haplotype reads, sparse PLINK2 dosage haplotypes, sparse
BGEN haplotypes, and hardcall-from-dosage conversion are not implemented in
this release.

Unsupported retained records fail the read. Metadata-only filters, such as
explicit variant ID lists and concrete regions, can skip unsupported records
before decoding their genotype or haplotype payloads.

::: genoio.Dataset

::: genoio.vcf

::: genoio.bfile

::: genoio.pfile

::: genoio.bgen
