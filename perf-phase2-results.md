# Phase 2 performance and verification results

Date: 2026-06-19 22:18:37 PDT

Base commit: `8d2af9f6cbf0fdf6ee716e83bff4a0971329d416`

Tree state: uncommitted Phase 2 PLINK2 genotype-filter optimization.

## Profile

Phase 2 profiling used `samply` against the release benchmark workload:

```bash
env CARGO_PROFILE_RELEASE_DEBUG=true make build-release
samply record --save-only --unstable-presymbolicate \
  -o /private/tmp/genoio-plink2-genotype-filtered-presym.json \
  -- .venv/bin/python scripts/benchmark_plink2.py \
    --backend genoio \
    --scenario genotype-filtered \
    --max-variants 10000 \
    --repeats 8 \
    --no-compare
```

Profile output:

- Benchmark median during profiling: `0.7307s` for 10k genotype-filtered.
- Saved profile: `/private/tmp/genoio-plink2-genotype-filtered-presym.json`.
- Saved symbols: `/private/tmp/genoio-plink2-genotype-filtered-presym.syms.json`.

Top relevant inclusive Rust samples:

| Symbol | Inclusive samples | Share |
|---|---:|---:|
| `genoio_io::plink::plink2::pgen::header::read_variable_width_header_body` | 4773 | 64.8% |
| `genoio_io::plink::plink2::pgen::io::read_plink2_variant_packed` | 898 | 12.2% |
| `genoio_io::hardcall::HardcallBatch::expand_into_sample_major` | 514 | 7.0% |
| `genoio_io::plink::plink2::pgen::main_track::decode_one_bit_record_with_cursor` | 170 | 2.3% |
| `genoio_io::plink::plink2::pgen::main_track::decode_difflist` | 94 | 1.3% |
| `genoio_io::hardcall::evaluate_packed_hardcall_filter` | 90 | 1.2% |

The profile showed full variable-width PGEN header parsing as the dominant cost
for matrix-only genotype-filtered block reads.

## Change

Matrix-only genotype-stat-filtered PLINK2 block reads now initialize with a
prefix PGEN header instead of a full variable-width header. If the retained
window is not satisfied within that prefix, the read retries with a larger
prefix until the window is satisfied or the full PGEN is decoded.

This keeps full-header behavior for ordinary reads and metadata-returning
reads. The optimized path still checks that the companion PVAR exists, but it
does not parse PVAR rows.

## Tests

Added focused coverage in `rust/genoio-io/tests/block_windows.rs`:

- `plink2_matrix_only_genotype_filter_extends_variable_width_prefix`
- `plink2_matrix_only_genotype_filter_prefix_ignores_later_unsupported_records`

The existing PLINK2 matrix tests cover fixed-width, variable-width,
LD-compressed, one-bit, and difflist records.

Verification commands:

```bash
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io --test block_windows plink2_matrix_only_genotype_filter
make rust-check
make rust-test
make build-release
```

All verification commands succeeded.

## PLINK2 Release Benchmarks

Baseline medians are from `perf-phase0-results.md` and
`perf-phase1-results.md`.

Command:

```bash
.venv/bin/python scripts/benchmark_plink2.py --backend genoio --scenario genotype-filtered --max-variants 1000 --repeats 21 --no-compare
```

| Scenario | Phase 0 median | Phase 1 median | Phase 2 median | Change vs Phase 1 | Runs |
|---|---:|---:|---:|---:|---|
| genotype-filtered | 0.5871s | 0.5943s | 0.0602s | -89.9% | 0.0635 0.0636 0.0651 0.0632 0.0580 0.0560 0.0617 0.0588 0.0562 0.0576 0.0569 0.0560 0.0566 0.0642 0.0635 0.0618 0.0599 0.0602 0.0633 0.0634 0.0585 |

Non-target checks:

| Scenario | Phase 0 median | Phase 1 median | Phase 2 median | Notes |
|---|---:|---:|---:|---|
| matrix-only, 1k | 0.0093s | 0.0091s | 0.0091s | 21-repeat serial rerun within tolerance |
| matrix-only, 10k | 0.0960s | 0.0987s | 0.0985s | 21-repeat serial rerun within tolerance |
| sample-filtered, 1k | 0.0090s | 0.0084s | 0.0085s | 21-repeat serial rerun within tolerance |
| with-variants, 1k | 0.0127s | 0.0129s | 0.0134s | short-read median remained above 3%; min stayed close to baseline and this path was not changed |

The initial 7-repeat benchmark pass showed noisy non-target medians; serial
higher-repeat reruns brought matrix-only and sample-filtered checks back within
tolerance. The with-variants median remained above the 3% threshold, but the
optimized code path is guarded by `matrix_only` and genotype-stat-only filters,
so this appears to be short-read variance rather than a touched-path regression.

Haplotype benchmark commands were attempted and skipped for the same local
fixture limitation recorded in Phase 0:

- `--kind haplo-hardcall` failed with `unphased pgen hardcall record retained in haplotype read`.
- `--kind haplo-dosage` failed with `pgen record does not contain explicit phased dosage values`.

## Phase 2 status

Completed:

- Sampled profile identified the changed code path.
- Matrix-only genotype-filtered PLINK2 block reads avoid full variable-width
  header parsing when a prefix is enough.
- Focused prefix-extension and unused-later-record tests were added.
- Rust check/test gates passed.
- PLINK2 genotype-filtered median improved by roughly 90% on the 1k benchmark.

Open:

- A valid large explicit phased PLINK2 fixture is still needed before haplotype
  benchmark regression gates can be enforced.
- The 1k with-variants benchmark remains noisy enough to sit above the 3%
  median threshold despite no touched-path change.
