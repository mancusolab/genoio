# Solver Feasibility Matrix

## Context
- Plan slug: `bgen-dosage`
- Generated date: `2026-06-04`

| Candidate | Problem Form Fit (root/least-squares/minimize) | AD Compatibility | Constraint Handling | Status/Error Mapping | Feasible (yes/no) | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Optimistix | Not applicable | Not applicable | Not applicable | Not applicable | no | BGEN support is file decoding, not numerical optimization. |
| Lineax | Not applicable | Not applicable | Not applicable | Not applicable | no | BGEN support is file decoding, not linear algebra solving. |
| Custom Solver | Not applicable | Not applicable | Not applicable | Not applicable | no | No solver is required. |

## Decision
- Preferred solver path: none
- Reason: deterministic binary parsing and dosage computation only.
- Benchmark or validation requirement before implementation: decoding tests and read-path benchmarks, not solver validation.
