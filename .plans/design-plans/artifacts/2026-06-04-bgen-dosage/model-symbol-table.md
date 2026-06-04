# Model Symbol Table

## Context
- Plan slug: `bgen-dosage`
- Generated date: `2026-06-04`

| Symbol | Meaning | Domain/Support | Shape/Type | Defined In Source | Notes |
| --- | --- | --- | --- | --- | --- |
| `N` | Number of samples in the BGEN file | Positive integer, must match sample identifier source | `u32` in BGEN header, `usize` internally | BGEN v1.3 header and sample block | Used for sample count validation. |
| `M` | Number of variant data blocks | Non-negative integer | `u32` in BGEN header, `usize` internally | BGEN v1.3 header | Used for metadata/read loop validation. |
| `K` | Number of alleles in one variant | First implementation requires `K = 2` | `u16` in variant block | BGEN v1.3 variant block | Multiallelic variants are rejected initially. |
| `P(AA), P(AB), P(BB)` | Biallelic unphased diploid genotype probabilities | Each probability in `[0, 1]`; last probability may be inferred by Layout 2 | `f32`/`f64` while decoding | BGEN v1.3 Layout 2 probability block | Used to compute `a1` dosage. |
| `dosage_a1` | Expected count of allele `a1` | `[0, 2]` or missing | `f32` matrix value plus missing mask | genoio dense dosage contract | `dosage_a1 = P(AB) + 2 * P(BB)`. |

## Checks
- [x] No undefined symbols.
- [x] No conflicting symbol reuse.
- [x] Support/domain constraints are explicit.
