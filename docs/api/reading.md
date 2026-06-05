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
`2` copies of `a1`, with missing calls handled by the selected missing-data
policy. The `dosage` option controls this genotype value source:
`dosage="hardcall"` uses source hard calls, while `dosage="dosage"` requires
dosage-backed values. This release supports dense VCF `FORMAT/DS`, PLINK2
unphased biallelic dosage reads, and BGEN Layout 2 biallelic diploid dosage
reads. Phased BGEN records are collapsed to expected diploid A1 dosage. BGEN
dosage values are expected copies of `a1`. Sparse dosage, haplotype dosage, and
PLINK1 dosage are not implemented yet and raise
`genoio.UnsupportedRepresentation`.
Genotype-stat filters such as `maf`, `mac`, and `missing_rate` use the selected
value source, so dosage reads compute those statistics from expected allele
dosages rather than hardcall counts.

## Metadata frames

`Dataset.samples()` returns a Polars DataFrame in source sample order. Genotype
reads and blocks return the same sample columns when
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

`Dataset.variants()` returns a Polars DataFrame in source variant order.
`read(..., return_variants=True)` and `blocks(..., return_variants=True)` return
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

::: genoio.Dataset

::: genoio.vcf

::: genoio.bfile

::: genoio.pfile

::: genoio.bgen
