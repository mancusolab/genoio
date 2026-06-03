# Phase 5 TDD and Mutation Evidence

## Provenance

- Date: 2026-06-03 11:14:36 HST
- Base/head reviewed: `7c443b988b31f939b58e37cb7663fd0c327b4b0e..7cd69b66edf937585f4c8caec0f0fde19bf7900f`
- Evidence type: mutation/regression evidence for the behavior-changing packed PLINK2 decoder rewrite.
- Toolchain note: commands use `env CC=clang AR=ar` because the active shell environment points Cargo's default `AR` at a missing `arm64-apple-darwin20.0.0-ar`.

## Baseline Green

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io plink2_dense --test plink2_dense
```

Outcome: passed, `17 passed; 0 failed; 1 filtered out`.

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io variable_record_helpers_write_packed_genotypes
```

Outcome: passed, `1 passed; 0 failed; 3 filtered out`.

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io packed_genotypes_copy_and_invert_0_2
```

Outcome: passed, `1 passed; 0 failed; 3 filtered out`.

## Temporary Mutation

Applied a temporary mutation in `rust/genoio-io/src/plink2.rs`:

```rust
fn invert_0_2(&mut self) {
    for sample_index in 0..self.sample_ct {
        match self.get(sample_index) {
            0 => self.set(sample_index, 0),
            2 => self.set(sample_index, 2),
            _ => {}
        }
    }
}
```

This intentionally made packed LD inversion a no-op. The mutation was reverted after the red checks below.

## Red Evidence

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io packed_genotypes_copy_and_invert_0_2
```

Outcome: failed as expected, `0 passed; 1 failed; 3 filtered out`.

Key failure:

```text
test plink2::tests::packed_genotypes_copy_and_invert_0_2 ... FAILED
left: [0, 1, 2, 3, 3]
right: [2, 1, 0, 3, 3]
```

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io variable_record_helpers_write_packed_genotypes
```

Outcome: failed as expected, `0 passed; 1 failed; 3 filtered out`.

Key failure:

```text
test plink2::tests::variable_record_helpers_write_packed_genotypes ... FAILED
left: [0, 1, 0, 3]
right: [2, 1, 2, 3]
```

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io plink2_dense_decodes_variable_width_hardcall_records --test plink2_dense
```

Outcome: failed as expected, `0 passed; 1 failed; 17 filtered out`.

Key failure:

```text
test plink2_dense_decodes_variable_width_hardcall_records ... FAILED
left: [0.0, 0.0, 2.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
right: [0.0, 0.0, 2.0, 2.0, 0.0, 1.0, 1.0, 0.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0]
```

## Restored Green

Restored `PackedGenotypes::invert_0_2` to swap categories `0` and `2`, then reran the same focused checks.

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io packed_genotypes_copy_and_invert_0_2
```

Outcome: passed, `1 passed; 0 failed; 3 filtered out`.

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io variable_record_helpers_write_packed_genotypes
```

Outcome: passed, `1 passed; 0 failed; 3 filtered out`.

Command:

```bash
env CC=clang AR=ar cargo test -p genoio-io plink2_dense_decodes_variable_width_hardcall_records --test plink2_dense
```

Outcome: passed, `1 passed; 0 failed; 17 filtered out`.

Command:

```bash
git diff -- rust/genoio-io/src/plink2.rs
```

Outcome: no output; the temporary mutation was fully reverted.
