# Equation To Code Map

## Context
- Plan slug: `bgen-dosage`
- Generated date: `2026-06-04`

| Equation ID | Equation (LaTeX or text) | Intended Computation | Target Module/Function | Test ID | Status |
| --- | --- | --- | --- | --- | --- |
| EQ-1 | `dosage_a1 = P(AB) + 2 * P(BB)` | Convert biallelic unphased diploid genotype probabilities into expected copies of `a1`. | `genoio-io` BGEN dense dosage decoder | `bgen-dosage.AC3.*` | planned |

## Checks
- [x] Objective sign and optimization direction are not applicable; this is deterministic file decoding.
- [x] Update rules are not applicable.
- [x] Every mapped equation has a corresponding test target.
