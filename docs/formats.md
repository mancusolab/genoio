# Format support

`genoio` exposes one matrix API across VCF, PLINK1, PLINK2, and BGEN sources.
Genotype reads return hardcall allele counts by default. Where a format stores
dense dosages, pass `dosage="dosage"` to read those values instead.

| Format | Inputs | Genotype reads | Haplotype reads | Notes |
|---|---|---:|---:|---|
| VCF/BCF | `.vcf`, `.vcf.gz`, `.vcf.bgz`, `.bcf` | yes; dense `FORMAT/DS` dosage supported | phased hardcall `FORMAT/GT` records | Indexed region filters use `.tbi` or `.csi` when available. |
| PLINK1 | `.bed` + `.bim` + `.fam` | yes | no | Variant-major BED files are supported. |
| PLINK2 | `.pgen` + `.pvar` or `.pvar.zst` + `.psam` | yes; dense unphased biallelic dosage supported | explicit phased hardcalls with `dosage="hardcall"`; dense explicit phased dosages with `dosage="dosage"` | Biallelic hard-call PGEN records are supported. Sparse PLINK2 hardcall haplotypes are supported; sparse PLINK2 dosage haplotypes are intentionally unsupported. |
| BGEN | `.bgen` plus optional same-prefix `.sample` | dense `kind="geno", dosage="dosage"` only | dense BGEN v1.2+ Layout 2 phased biallelic diploid probabilities with `kind="haplo", dosage="dosage"` | Dosage-backed BGEN reads use expected A1 dosage values. Concrete region filters use a same-path `.bgen.bgi` index when present. |

!!! tip "Choose the source by what you need to read"
    For hardcall genotype matrices, VCF/BCF, PLINK1, and PLINK2 are the usual
    choices. For dosage matrices, use VCF `FORMAT/DS`, PLINK2 dosage records, or
    BGEN dosage records. For haplotypes, use phased VCF hardcalls, phased PLINK2
    records, or phased BGEN dosage records.

---

## Source resolution

VCF and BCF inputs are single files.

```python
X = genoio.vcf("cohort.vcf.gz").read()
```

PLINK inputs are file sets. Pass either the shared prefix or one member file.

```python
X = genoio.pfile("cohort").read()
X = genoio.pfile("cohort.pgen").read()
```

Use `bfile(...)` for PLINK1 prefixes and `pfile(...)` for PLINK2 prefixes.
The constructor chooses the file-set type, so same-stem files from other
formats are ignored.

For PLINK2, `genoio` accepts either an uncompressed `.pvar` or a zstd-compressed
`.pvar.zst`. If both exist for the same prefix, `.pvar` is used.

BGEN inputs are `.bgen` files with an optional same-prefix `.sample` file.

```python
X = genoio.bgen("cohort.bgen").read(dosage="dosage")
```

!!! note "BGEN sample IDs"
    BGEN reads require real sample IDs, either embedded in the `.bgen` file or
    provided by the same-prefix `.sample` file.

??? info "BGEN dosage details"
    BGEN v1.2+ Layout 2 biallelic diploid dosage records are returned as
    expected A1 allele dosages. Genotype reads of phased BGEN records collapse
    the two source haplotype probabilities to expected diploid A1 dosage.
    Haplotype reads return expected A1 dosage per haplotype row.

    Matrix-only BGEN reads avoid returning sample and variant metadata unless
    `return_samples=True` or `return_variants=True` is requested.

    For concrete region filters such as
    [`genoio.region("22:20000000-21000000")`](api/filters.md#genoio.region),
    BGEN dosage reads use `cohort.bgen.bgi` when that index exists. The index is
    used only for concrete region pushdown; other predicates still run through
    the normal filter path after candidate records are read.

---

## Boundaries to know

The table above is the contract most users need. These details matter when a
file mixes record encodings or when you request sparse, dosage, or haplotype
output.

- **Dosages.** `dosage="dosage"` reads dense VCF `FORMAT/DS`, dense PLINK2
  biallelic dosages, and dense BGEN v1.2+ Layout 2 biallelic diploid dosages.
  PLINK1 dosages, sparse dosages, and hardcall-from-dosage conversion are not
  supported.
- **Haplotypes.** Haplotype reads use phased VCF hardcalls, explicit phased
  PLINK2 hardcall/full-dosage records, or phased BGEN dosage records. Hardcall
  haplotypes must come from source hardcalls, not probabilities or dosages.
- **Record encodings.** PLINK2 support is limited to biallelic hard-call,
  unphased genotype dosage, explicit phased hardcall haplotype, and explicit
  phased full-dosage haplotype records. BGEN support is limited to dense
  dosage-backed genotype and haplotype reads.
- **Sparse and indexed reads.** Sparse reads don't preserve missing-value masks.
  Region pushdown is implemented for concrete indexed VCF/BCF and BGEN region
  filters, not arbitrary filter expressions.

!!! warning "Unsupported retained records fail the read"
    `genoio` skips records removed by metadata-only filters, such as explicit ID
    or region filters, before decoding their genotype payload. If an unsupported
    record survives filtering, the read fails instead of silently changing the
    data.
