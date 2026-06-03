# genoio Implementation Plan

This directory is the implementation-plan source for the initial `genoio` package build.

## Scope

- Package skeleton, public API, and source resolution.
- Metadata reads and source capabilities.
- Dense VCF and PLINK1 genotype reads.
- Rust-evaluated variant filters.
- Sparse reads and variant block iteration.
- PLINK2 parser strategy deferral and phased VCF haplotype reads.

## Plan Files

- `phase_01.md` through `phase_06.md`: phase tasks, acceptance criteria, and verification commands.
- `test-requirements.md`: final acceptance-criteria coverage requirements and required verification commands.
- `phase_*_review_remediation.md`: review findings and remediation evidence captured during execution.

## Status

Implementation and phase-level reviews completed. Final verification uses the commands listed in `test-requirements.md`.
