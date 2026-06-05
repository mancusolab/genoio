# Model Symbol Table

## Context
- Plan slug: `haplotype-io`
- Generated date: `2026-06-05`

| Symbol | Meaning | Domain/Support | Shape/Type | Defined In Source | Notes |
| --- | --- | --- | --- | --- | --- |
| `X_haplo` | Dense haplotype matrix returned by `Dataset.read(kind="haplo", ...)`. | Supported source-encoded haplotype representations only. | NumPy array with shape `(2 * retained_samples, retained_variants)` for diploid reads. | `src/genoio/_api.py`; `rust/genoio-py/src/lib.rs` | Rows are ordered by source sample, then haplotype index. |
| `a1` | Counted allele in returned matrices. | Biallelic retained variants. | Variant metadata string column. | `rust/genoio-core/src/metadata.rs`; `src/genoio/_assembly.py` | Values count expected copies of `a1`. |
| `h` | Haplotype row index within a source sample. | Diploid retained samples. | Integer `0` or `1`. | Existing haplotype metadata contract. | Returned as `haplotype_index`. |
| `source_sample_index` | Source-order sample index backing a haplotype row. | Retained samples after sample filtering. | Integer metadata column. | Existing haplotype metadata contract. | Used to map haplotype rows back to diploid samples. |
| `variant_window` | Retained-variant window used by `iter_blocks(...)`. | Non-negative start/length over retained variants. | Rust `VariantWindow`. | `rust/genoio-core/src/metadata.rs`; `src/genoio/_api.py` | Windowing occurs after filters. |
| `region` | Concrete genomic interval filter. | 1-based inclusive `chrom:start-end`. | `FilterExpr` / Rust `VariantFilter`. | `src/genoio/_filters.py`; `rust/genoio-core/src/filter.rs` | May use indexed VCF/BGEN pushdown. |

## Checks
- [x] No undefined symbols.
- [x] No conflicting symbol reuse.
- [x] Support/domain constraints are explicit.
