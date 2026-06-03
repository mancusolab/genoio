# Format support

`genoio` exposes one matrix API across VCF, PLINK1, and PLINK2 sources.

| Format | Inputs | Genotype reads | Haplotype reads | Notes |
|---|---|---:|---:|---|
| VCF/BCF | `.vcf`, `.vcf.gz`, `.bcf` | yes | phased VCF only | Indexed region filters use `.tbi` or `.csi` when available. |
| PLINK1 | `.bed` + `.bim` + `.fam` | yes | no | Variant-major BED files are supported. |
| PLINK2 | `.pgen` + `.pvar` + `.psam` | yes | no | Biallelic hard-call PGEN records are supported. |

---

## Source resolution

VCF and BCF inputs are single files:

```python
X = genoio.read("cohort.vcf.gz")
```

PLINK inputs are file sets. Pass either the shared prefix or one member file:

```python
X = genoio.read("cohort", format="plink2")
X = genoio.read("cohort.pgen")
```

Use `format=...` when a shared prefix could refer to more than one supported
format.

---

## Current limitations

- PLINK2 support is limited to biallelic hard-call records. Dosage tracks are
  not implemented yet.
- Sparse reads do not preserve missing-value masks.
- Haplotype reads are currently VCF-only.
- Region pushdown is implemented for concrete indexed VCF/BCF region filters,
  not for arbitrary filter expressions.
