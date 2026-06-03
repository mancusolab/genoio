# Phase 6 TDD and Regression Evidence

## Provenance

- Date: 2026-06-03 HST.
- Phase range under review: `9a47fde1b3252bde0ebbdb8d184de15803cff1da..75477309634e75c605ecf48230798d0bdd5584bc`.
- Review-fix scope: update stale performance-document assertions and add missing packed-batch integration regression confidence.
- Process note: the original Phase 6 commits did not include a Phase 6 TDD evidence artifact, and the implementation commits did not change the plan-named integration test files. This artifact records post-review red/green and mutation/regression evidence; it is not a fabricated record of the original implementation process.
- Toolchain note: Rust commands use `env CC=clang AR=ar` because this repository's local build environment expects those tool overrides for Cargo verification.

## Performance Documentation Red/Green

Initial red command:

```bash
pytest tests/test_performance_docs.py -q
```

Outcome: failed as expected, `1 failed`. The stale test still asserted Phase 4 direct-fill provenance and timing text:

```text
AssertionError: assert '0.0145-0.0147 s' in ...
```

Fix applied: `tests/test_performance_docs.py` now asserts the Phase 6 packed-batch decision text, provenance commit `3bb767085c43c8a39687fa93e4b238c305d3c5bc`, 1,000-variant medians `0.0100 s` and `0.0069 s`, 10,000-variant medians `0.0941 s` and `0.0497 s`, and the decision to keep packed batches for unfiltered dense source windows.

Green command:

```bash
pytest tests/test_performance_docs.py -q
```

Outcome: passed, `1 passed`.

## Packed-Batch Regression Tests

Existing baseline green commands before mutation:

```bash
env CC=clang AR=ar cargo test -p genoio-io packed_variant_batch_expands_like_variant_at_a_time
```

Outcome: passed, `1 passed; 0 failed; 4 filtered out` for the library unit tests.

```bash
env CC=clang AR=ar cargo test -p genoio-io plink2_dense_fixed_width_source_window_matches_full_read_slice --test plink2_dense
```

Outcome: passed, `1 passed; 0 failed; 17 filtered out`.

```bash
env CC=clang AR=ar cargo test -p genoio-io plink2_dense_variable_width_source_window_matches_full_read_slice --test plink2_dense
```

Outcome: passed, `1 passed; 0 failed; 17 filtered out`.

Added regression test:

```text
rust/genoio-io/tests/plink2_dense.rs::plink2_dense_fixed_width_source_window_crosses_packed_batch_boundary
```

The test builds a 70-variant fixed-width PLINK2 fixture, reads a 66-variant source window that crosses the private packed-batch size of 64, and checks both metadata-bearing and matrix-only source-window output against the full-read sample-major slice. It also checks the first and last returned variant IDs (`rs2` and `rs67`) to keep metadata alignment tied to the same window.

Green command after adding the regression:

```bash
env CC=clang AR=ar cargo test -p genoio-io plink2_dense_fixed_width_source_window_crosses_packed_batch_boundary --test plink2_dense
```

Outcome: passed, `1 passed; 0 failed; 18 filtered out`.

## Temporary Mutation

Applied a temporary mutation in `rust/genoio-io/src/plink2.rs`:

```rust
let variant_index = batch_variant_index;
```

This intentionally ignored `variant_start` in `PackedVariantBatch::expand_into_sample_major`, so the final partial batch overwrote columns at the start of the output instead of writing after the first full batch. The mutation was reverted after the red checks below.

Red command:

```bash
env CC=clang AR=ar cargo test -p genoio-io packed_variant_batch_expands_like_variant_at_a_time
```

Outcome: failed as expected, `0 passed; 1 failed; 4 filtered out`.

Key failure:

```text
test plink2::tests::packed_variant_batch_expands_like_variant_at_a_time ... FAILED
assertion `left == right` failed
```

Red command:

```bash
env CC=clang AR=ar cargo test -p genoio-io plink2_dense_fixed_width_source_window_crosses_packed_batch_boundary --test plink2_dense
```

Outcome: failed as expected, `0 passed; 1 failed; 18 filtered out`.

Key failure:

```text
test plink2_dense_fixed_width_source_window_crosses_packed_batch_boundary ... FAILED
assertion `left == right` failed
```

## Restored Green

Restored the packed-batch output index calculation:

```rust
let variant_index = variant_start + batch_variant_index;
```

Then reran the focused checks.

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io packed_variant_batch_expands_like_variant_at_a_time
```

Outcome: passed, `1 passed; 0 failed; 4 filtered out`.

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io plink2_dense_fixed_width_source_window_crosses_packed_batch_boundary --test plink2_dense
```

Outcome: passed, `1 passed; 0 failed; 18 filtered out`.

Command:

```bash
git diff -- rust/genoio-io/src/plink2.rs
```

Outcome: no output; the temporary mutation was fully reverted.
