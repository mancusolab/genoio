# Format support

`genoio` exposes one matrix API across VCF, PLINK1, and PLINK2 sources.
Genotype reads return hardcall allele counts by default. Dense VCF and PLINK2
reads can instead use stored dosage values with `dosage="dosage"`.

| Format | Inputs | Genotype reads | Haplotype reads | Notes |
|---|---|---:|---:|---|
| VCF/BCF | `.vcf`, `.vcf.gz`, `.bcf` | yes; dense `FORMAT/DS` dosage supported | phased VCF only | Indexed region filters use `.tbi` or `.csi` when available. |
| PLINK1 | `.bed` + `.bim` + `.fam` | yes | no | Variant-major BED files are supported. |
| PLINK2 | `.pgen` + `.pvar` or `.pvar.zst` + `.psam` | yes; dense unphased biallelic dosage supported | no | Biallelic hard-call PGEN records are supported. |

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

---

## Current limitations

- `dosage="dosage"` currently supports dense VCF `FORMAT/DS` reads and dense
  PLINK2 unphased biallelic dosage reads. Sparse dosage reads are not
  implemented.
- PLINK1 has no dosage representation in BED files.
- PLINK2 support is limited to biallelic hard-call and unphased dosage records.
- Sparse reads do not preserve missing-value masks.
- Haplotype reads are currently VCF-only.
- Region pushdown is implemented for concrete indexed VCF/BCF region filters,
  not for arbitrary filter expressions.
