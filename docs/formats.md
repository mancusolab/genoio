# Format support

`genoio` exposes one matrix API across VCF, PLINK1, PLINK2, and BGEN sources.
Genotype reads return hardcall allele counts by default. Dense VCF, PLINK2, and
BGEN reads can instead use stored dosage values with `dosage="dosage"`.

| Format | Inputs | Genotype reads | Haplotype reads | Notes |
|---|---|---:|---:|---|
| VCF/BCF | `.vcf`, `.vcf.gz`, `.bcf` | yes; dense `FORMAT/DS` dosage supported | phased hardcall `FORMAT/GT` records | Indexed region filters use `.tbi` or `.csi` when available. |
| PLINK1 | `.bed` + `.bim` + `.fam` | yes | no | Variant-major BED files are supported. |
| PLINK2 | `.pgen` + `.pvar` or `.pvar.zst` + `.psam` | yes; dense unphased biallelic dosage supported | dense explicit phased hardcalls with `dosage="hardcall"` and explicit phased dosages with `dosage="dosage"` | Biallelic hard-call PGEN records are supported. Sparse PLINK2 haplotypes are not implemented. |
| BGEN | `.bgen` plus optional same-prefix `.sample` | dense `kind="geno", dosage="dosage"` only | dense Layout 2 phased biallelic diploid probabilities with `kind="haplo", dosage="dosage"` | Dosage-backed BGEN reads use expected A1 dosage values. Concrete region filters use a same-path `.bgen.bgi` index when present. |

---

## Source resolution

VCF and BCF inputs are single files:

```python
X = genoio.vcf("cohort.vcf.gz").read()
```

PLINK inputs are file sets. Pass either the shared prefix or one member file:

```python
X = genoio.pfile("cohort").read()
X = genoio.pfile("cohort.pgen").read()
```

Use `bfile(...)` for PLINK1 prefixes and `pfile(...)` for PLINK2 prefixes.
The constructor chooses the file-set type, so same-stem files from other
formats are ignored.

For PLINK2, `genoio` accepts either an uncompressed `.pvar` or a zstd-compressed
`.pvar.zst`. If both exist for the same prefix, `.pvar` is used.

BGEN inputs are `.bgen` files with an optional same-prefix `.sample` file:

```python
X = genoio.bgen("cohort.bgen").read(dosage="dosage")
```

BGEN reads require real sample IDs, either embedded in the `.bgen` file or
provided by the same-prefix `.sample` file. Layout 2 biallelic diploid dosage
records are returned as expected A1 allele dosages. Genotype reads of phased
BGEN records collapse the two source haplotype probabilities to expected
diploid A1 dosage; haplotype reads return expected A1 dosage per haplotype row.
Matrix-only BGEN reads avoid returning sample and variant metadata unless
`return_samples=True` or `return_variants=True` is requested.

For concrete region filters such as
[`genoio.region("22:20000000-21000000")`](api/filters.md#genoio.region), BGEN
dosage reads use a same-path bgenix SQLite index when present. For
`cohort.bgen`, the expected index path is `cohort.bgen.bgi`. If the index is
absent, reads fall back to the normal sequential scan. The index is used only
for concrete region pushdown; other predicates still run through the normal
metadata or genotype filter path after candidate records are read.

---

## Current limitations

- `dosage="dosage"` currently supports dense VCF `FORMAT/DS` reads, dense
  PLINK2 unphased biallelic genotype dosage reads, dense PLINK2 explicit
  phased full-dosage haplotype reads, and dense BGEN Layout 2 biallelic diploid
  dosage reads. Phased BGEN records can return haplotype rows or collapse to
  expected diploid A1 dosage. Sparse dosage reads are not implemented.
- PLINK1 has no dosage representation in BED files.
- PLINK2 support is limited to biallelic hard-call, unphased genotype dosage,
  explicit phased hardcall haplotype, and explicit phased full-dosage haplotype
  records.
- BGEN support is limited to dense dosage-backed genotype and haplotype reads.
  Hardcall conversion, sparse reads, multiallelic BGEN, variable ploidy, and
  unsupported compression and layout values are not supported. `.bgi` pushdown
  is limited to concrete region filters.
- Sparse reads do not preserve missing-value masks.
- Haplotype reads are implemented for phased VCF hardcalls, explicit phased
  PLINK2 hardcall/full-dosage records, and phased BGEN dosage records.
- Hardcall-from-dosage conversion is not performed. Requests for hardcall
  haplotypes must use hardcalls encoded by the source format, not probabilities.
- Unsupported retained records fail the read. Records excluded by metadata-only
  filters, such as explicit ID or region filters, can be skipped before their
  genotype or haplotype payload is decoded.
- Region pushdown is implemented for concrete indexed VCF/BCF and BGEN region
  filters, not for arbitrary filter expressions.
