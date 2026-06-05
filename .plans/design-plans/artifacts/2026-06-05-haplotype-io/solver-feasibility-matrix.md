# Solver Feasibility Matrix

## Context
- Plan slug: `haplotype-io`
- Generated date: `2026-06-05`

| Candidate | Problem Form Fit (root/least-squares/minimize) | AD Compatibility | Constraint Handling | Status/Error Mapping | Feasible (yes/no) | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Optimistix | n/a | n/a | n/a | n/a | no | No solver is introduced by this I/O design. |
| Lineax | n/a | n/a | n/a | n/a | no | No linear solve is introduced by this I/O design. |
| Custom Solver | n/a | n/a | n/a | Source decoding errors map to existing `genoio` exceptions. | no | Binary decoding belongs in existing Rust readers, not a numerical solver. |

## Decision
- Preferred solver path: n/a
- Reason: This design decodes source-encoded haplotype representations; it does
  not define an optimization or inference problem.
- Benchmark or validation requirement before implementation: Benchmarks are
  added after correctness tests for the new I/O paths.
