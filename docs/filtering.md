# Filtering

Use filters to keep the variants you want before they become a matrix. Python
builds the expression; Rust evaluates it while reading records.

```python
rare_high_quality = (
    genoio.region("22:20000000-21000000")
    & genoio.qual(min=20)
    & genoio.maf(max=0.05)
)

ds = genoio.vcf("data/chr22_hg38.vcf.gz")
samples = ds.samples()
y = load_phenotype_vector(samples["iid"])

with ds.iter_blocks(
    5_000,
    variants=rare_high_quality,
    return_variants=True,
) as blocks:
    for X, variants in blocks:
        association_scan(X, y, variants=variants)
```

Filters compose with Python operators:

```python
genoio.chrom("22") & genoio.snp()
genoio.maf(max=0.01) | genoio.id_in(["rs123", "rs456"])
~genoio.missing_rate(max=0.1)
```

---

## Metadata filters and genotype filters

Some predicates use metadata already stored on the variant record: chromosome,
position, ID, REF/ALT structure, and `QUAL`. These can remove records before
genotypes are decoded.

Other predicates need the genotype values. `maf(...)`, `mac(...)`,
`missing_rate(...)`, and `polymorphic()` are evaluated only for records that
survive the metadata part of the expression. For common genotype-stat filters,
the Rust reader can often evaluate a backend-specific filter plan from native
counts or packed genotypes instead of building full variant statistics first.

!!! tip "Combine metadata and genotype filters"
    `region("22:20000000-21000000") & qual(min=20) & maf(max=0.05)` lets the
    reader narrow the candidate set before evaluating MAF. You do not need to
    hand-order simple expressions for performance: Rust simplifies the filter IR
    and partially evaluates metadata predicates before genotype decoding.

---

## Region filters

Use a concrete region string with 1-based inclusive coordinates:

```python
region = genoio.region("22:20000000-21000000")
```

!!! note "Indexed regions"
    Compressed text VCF region reads require a `.tbi` or `.csi` index. `genoio`
    rejects unindexed compressed text VCF region reads instead of silently
    scanning the whole file.

    BGEN dosage reads use a same-path bgenix SQLite index when one exists. For
    `cohort.bgen`, the expected path is `cohort.bgen.bgi`. Without that index,
    BGEN falls back to a sequential scan.

```python
bgen_ds = genoio.bgen("cohort.bgen")
region = genoio.region("22:20000000-21000000")

with bgen_ds.iter_blocks(
    1_000,
    dosage="dosage",
    variants=region,
    return_variants=True,
) as blocks:
    for X, variants in blocks:
        analyze_region(X, variants)
```

---

## Common predicates

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

Threshold predicates are inclusive. `maf(max=0.05)` keeps variants with minor
allele frequency less than or equal to `0.05`.

??? warning "When a filter matches no variants"
    Empty filters are valid. Whole reads return shape `(n_samples, 0)`, returned
    variant metadata has the standard columns with zero rows, and block reads
    yield no blocks.

??? info "Advanced: filter IR"
    Filter expressions are converted to a small intermediate representation
    before they cross into Rust:

    ```python
    expr = genoio.chrom("22") & genoio.maf(max=0.05)
    expr.to_ir()
    ```

    The IR records predicate names, validated parameters, and boolean structure.
    It does not contain Python callbacks. That keeps filtering portable across
    whole reads, block reads, and Rust reader implementations.

    Rust uses the IR to separate metadata decisions from genotype-dependent
    decisions. Metadata predicates can be evaluated before matrix construction,
    and concrete compressed text VCF or BGEN region predicates can use an index
    when one is available. Before reading, Rust also normalizes simple boolean
    expressions: overlapping conjoined regions are reduced to their
    intersection, repeated threshold predicates are tightened, conjoined `id_in`
    predicates are intersected, and contradictory predicates become an empty
    result without scanning variant records.

    Treat `to_ir()` as an inspection aid rather than a stable wire format. Build
    filters with the Python constructors so validation stays consistent.

See [Filter API](api/filters.md) for signatures and validation rules.
