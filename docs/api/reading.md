# Reading

Use these entry points to open genotype sources and read matrices. `vcf`,
`bfile`, and `pfile` identify the input format and return a `Dataset`; the
dataset methods then read dense NumPy arrays, sparse SciPy matrices, metadata
frames, or streaming blocks. Dense and sparse genotype matrices use samples as
rows and variants as columns, so block reads can be concatenated by variant when
sample ordering is unchanged.

## Metadata frames

`Dataset.samples()` returns a Polars DataFrame in source sample order. The same
schema is returned when `read(..., return_samples=True)` or
`blocks(..., return_samples=True)` is used.

Haplotype reads return one row per haplotype. Those returned sample frames also
include `source_sample_index` and `haplotype_index` columns so each haplotype row
can be mapped back to its original diploid sample.

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

`Dataset.variants()` returns a Polars DataFrame in source variant order.
`read(..., return_variants=True)` and `blocks(..., return_variants=True)` return
the same columns for the retained variants, in matrix-column order.

??? info "Variant metadata schema"
    | Column | Meaning |
    |---|---|
    | `chrom` | Chromosome or contig label. |
    | `pos` | 1-based variant position. |
    | `id` | Variant ID from the source. |
    | `a0` | Allele counted as dosage `0` in returned matrices. |
    | `a1` | Allele counted as dosage `1` or `2` in returned matrices. |
    | `ref` | Reference allele when the source records one; otherwise null. |
    | `alt` | Alternate allele string when the source records one; otherwise null. |
    | `source_a0` | Original source allele corresponding to `a0`. |
    | `source_a1` | Original source allele corresponding to `a1`. |
    | `flipped` | Whether `genoio` flipped source allele orientation for normalized output. |
    | `qual` | VCF `QUAL` value when available; otherwise null. |
    | `af` | Alternate-allele frequency when genotype statistics were computed; otherwise null. |
    | `maf` | Minor-allele frequency when genotype statistics were computed; otherwise null. |
    | `mac` | Minor-allele count when genotype statistics were computed; otherwise null. |
    | `missing_rate` | Missing-call fraction when genotype statistics were computed; otherwise null. |
    | `n_called` | Number of called genotypes when genotype statistics were computed; otherwise null. |

::: genoio.Dataset

::: genoio.vcf

::: genoio.bfile

::: genoio.pfile
