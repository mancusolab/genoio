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
drop records before genotype decoding.

Genotype predicates require retained genotypes to be decoded first. MAF, MAC,
missing rate, and polymorphism are genotype predicates.

This distinction matters for speed. A filter like `qual(min=20) & snp()` can
discard records before matrix construction. A filter like `maf(max=0.05)` must
decode candidate genotypes before deciding whether to keep the variant.

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

See [Filter API](api/filters.md) for signatures and validation rules.
