# Phase 4 targeted decoder tests

Date: 2026-06-19
Base commit before Phase 4 changes: `c5be20b`

## Coverage added

- Added PGEN bitpack helper tests for base-128 varint decoding, truncated
  continuation errors, out-of-range varints, fixed-width sample IDs, byte/bit
  bounds checks, and difflist sample-id ordering/bounds validation.
- Added PGEN main-track tests for four packed difflist values, duplicate and
  out-of-range difflist sample IDs, and one-bit records with and without an
  extra difflist overlay.
- Added BGEN dosage filter count tests for monomorphic variants with missing
  calls, all-missing variants, MAC/MAF/missing-rate decisions, and invalid
  dosage inputs.

## Existing parity coverage

- `tests/block_windows.rs` covers PLINK2 matrix-only genotype-filter windows,
  including variable-width prefixes and skipping PVAR rows.
- `tests/plink2_dense.rs` covers fixed-width, variable-width, LD-compressed,
  one-bit, difflist, sample-filtered, and matrix-only source-window behavior.
- `tests/bgen_dense.rs` covers BGEN matrix-only reads, sample-filtered reads,
  genotype-stat filters, indexed windows, and haplotype dosage paths.

## Verification

Focused checks:

```bash
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io plink::plink2::pgen::bitpack
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io plink::plink2::pgen::main_track
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io bgen::filter
```

Result:

- PGEN bitpack: 6 passed.
- PGEN main-track: 4 passed.
- BGEN filter: 8 passed.

Full gates:

```bash
make rust-fmt
make rust-check
make rust-test
make build-release
```

Result: all passed.
