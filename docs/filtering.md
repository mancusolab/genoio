# Filtering

Filters are serializable expression objects. Python builds the expression, and
Rust evaluates it while reading records.

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

---

## Metadata and genotype predicates

There are two kinds of predicates.

Metadata predicates use fields already present in the source record:
chromosome, position, ID, REF/ALT structure, and `QUAL`. These predicates can
drop records before genotype decoding. Examples include `chrom(...)`,
`region(...)`, `id_in(...)`, `snp()`, `biallelic()`, and `qual(...)`.

Genotype predicates require retained genotypes to be decoded first. MAF, MAC,
missing rate, and polymorphism are genotype predicates. Examples include
`maf(...)`, `mac(...)`, `missing_rate(...)`, and `polymorphic()`.

This distinction matters for speed. A filter like `qual(min=20) & snp()` can
discard records before matrix construction. A filter like `maf(max=0.05)` must
decode candidate genotypes before deciding whether to keep the variant.

Mixed expressions use both paths. In
`region("22:20000000-21000000") & qual(min=20) & maf(max=0.05)`, the reader can
use the region and quality predicates to reduce the candidate set before
computing MAF on the remaining variants.

---

## Region filters

Use a concrete region string with 1-based inclusive coordinates:

```python
region = genoio.region("22:20000000-21000000")
```

For compressed VCF/BCF sources, region reads require an index. `genoio` rejects
unindexed compressed region reads instead of silently scanning the full file.
When a `.tbi` or `.csi` index is present, the reader uses it to seek to the
requested region.

For BGEN dosage reads, concrete region filters use a same-path bgenix SQLite
index when present. For `cohort.bgen`, `genoio` looks for `cohort.bgen.bgi`.
When that index is absent, BGEN reads fall back to the normal sequential scan.

---

## Available predicates

```python
genoio.chrom("22")
genoio.region("22:20000000-21000000")
genoio.snp()
genoio.biallelic()
genoio.qual(min=20)
genoio.maf(max=0.05)
genoio.mac(min=10)
genoio.missing_rate(max=0.1)
genoio.polymorphic()
genoio.id_in(["rs123", "rs456"])
```

Threshold predicates are inclusive. For example, `maf(max=0.05)` keeps variants
with minor allele frequency less than or equal to `0.05`.

Filters are allowed to retain zero variants. Whole reads return an empty
variant axis with shape `(n_samples, 0)`, returned variant metadata keeps the
standard columns with zero rows, and block reads yield no blocks.

---

## Advanced: filter IR

Filter expressions are converted to a small intermediate representation before
they cross into Rust:

```python
expr = genoio.chrom("22") & genoio.maf(max=0.05)
expr.to_ir()
```

The IR records predicate names, validated parameters, and boolean structure. It
does not contain Python callbacks. That keeps filtering portable across whole
reads, block reads, and Rust reader implementations.

Rust uses the IR to separate cheap metadata decisions from data-dependent
genotype decisions. Metadata predicates can be evaluated before matrix
construction, and concrete VCF/BCF or BGEN region predicates can use an index
when one is available. Before reading, Rust also normalizes simple boolean expressions:
overlapping conjoined regions are reduced to their intersection, repeated
threshold predicates are tightened, conjoined `id_in` predicates are
intersected, and contradictory predicates become an empty result without
scanning variant records. Genotype predicates are delayed until the candidate
variant's genotypes are decoded, then the same retained-variant decision
controls both matrix columns and returned variant metadata.

Treat `to_ir()` as an inspection aid rather than a stable wire format. Build
filters with the Python constructors so validation stays consistent.

See [Filter API](api/filters.md) for signatures and validation rules.
