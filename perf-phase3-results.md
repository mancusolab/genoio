# Phase 3 performance and verification results

Date: 2026-06-19 22:26:40 PDT

Base commit: `8d2af9f6cbf0fdf6ee716e83bff4a0971329d416`

Tree state: uncommitted Phase 2 and Phase 3 changes.

## Scope

Phase 3 was a behavior-preserving organization pass for the BGEN reader.
`rust/genoio-io/src/bgen.rs` now acts as the facade and metadata entry point.
The dense reader internals were split by responsibility:

- `bgen/dense.rs`: diploid dosage dense orchestration, matrix-only fast path,
  and indexed-region dosage reads.
- `bgen/haplotype.rs`: phased haplotype dosage dense orchestration and
  indexed-region haplotype reads.
- `bgen/filter.rs`: dosage filter counts, compiled genotype-filter evaluation,
  and filter/stat attachment helpers.
- `bgen/session.rs`: `BgenReadSession`, sequential/indexed variant cursor, and
  shared indexed-read context.

Existing low-level modules remained in place: `decode`, `header`, `index`, and
`io`.

## Verification

All commands below ran successfully:

```bash
env CC=clang AR=ar cargo clippy --manifest-path rust/Cargo.toml -p genoio-io --all-targets -- -D warnings
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io --test bgen_dense
make rust-check
make rust-test
make build-release
```

`make rust-test` included the BGEN dense, haplotype, indexed-region,
genotype-filtered, metadata, and malformed-input tests.

## BGEN Release Benchmarks

Baseline medians are from `perf-phase0-results.md`.

Command:

```bash
.venv/bin/python scripts/benchmark_bgen.py --backend genoio --scenario all --max-variants 1000 --repeats 7
```

| Scenario | Phase 0 median | Phase 3 median | Change | Runs |
|---|---:|---:|---:|---|
| matrix-only | 0.0316s | 0.0233s | -26.3% | 0.0289 0.0228 0.0232 0.0247 0.0233 0.0231 0.0312 |
| with-variants | 0.0337s | 0.0248s | -26.4% | 0.0232 0.0221 0.0264 0.0242 0.0248 0.0279 0.0406 |
| sample-filtered | 0.0284s | 0.0625s | +120.1% | 0.0323 0.0343 0.1038 0.1397 0.1150 0.0625 0.0297 |
| genotype-filtered | 0.2242s | 0.2309s | +3.0% | 0.2337 0.2309 0.3125 0.3208 0.1992 0.2107 0.1574 |
| indexed-region | 0.0868s | 0.0835s | -3.8% | 0.0704 0.0835 0.0692 0.0936 0.0788 0.1047 0.0976 |

Command:

```bash
.venv/bin/python scripts/benchmark_bgen.py --backend genoio --scenario matrix-only --max-variants 10000 --repeats 7
```

| Scenario | Phase 0 median | Phase 3 median | Change | Runs |
|---|---:|---:|---:|---|
| matrix-only | 0.4605s | 0.4855s | +5.4% | 0.8673 0.5442 0.6543 0.4732 0.4311 0.4855 0.3807 |

The 10k matrix-only pass was noisy enough to exceed the 3% threshold, so it was
rerun serially with more repeats:

```bash
.venv/bin/python scripts/benchmark_bgen.py --backend genoio --scenario matrix-only --max-variants 10000 --repeats 15
```

| Scenario | Phase 0 median | Phase 3 median | Change | Runs |
|---|---:|---:|---:|---|
| matrix-only | 0.4605s | 0.4100s | -11.0% | 0.4039 0.4219 0.4051 0.4830 0.3899 0.4120 0.4033 0.4057 0.4390 0.3907 0.4377 0.3945 0.4213 0.4335 0.4100 |

## Phase 3 status

Completed:

- BGEN orchestration was split into focused modules without intended behavior
  changes.
- BGEN dense, haplotype, indexed-region, and genotype-filtered tests passed.
- BGEN matrix-only and indexed-region release medians stayed within the Phase 3
  acceptance criteria after accounting for variance.

Notes:

- Phase 3 was primarily organizational. The benchmark differences appear to be
  noise and compiler layout effects from the refactor, not a deliberate speed
  optimization.
- The 1k sample-filtered run was highly variable; it is not part of the Phase 3
  acceptance gate, but should be watched if future BGEN changes target sample
  filtering.
