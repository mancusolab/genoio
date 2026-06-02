# Phase 4 Baseline Review Remediation Evidence

## Issues Fixed

- `variants=` generator inputs were consumed during validation and then passed to Rust as an empty `id_in([])` filter.
- Indexed VCF region pushdown was missing; region reads used a full scan even when a tabix-indexed VCF was available.
- Phase 4 behavior changes now have red/green regression evidence.

## Red Evidence

### Generator Variant IDs

Command:

```bash
pytest -q tests/test_filters.py::test_generator_variant_id_selection_is_consumed_once
```

Result before fix:

```text
FAILED tests/test_filters.py::test_generator_variant_id_selection_is_consumed_once
AssertionError: assert [] == ['rs1', 'rs2']
```

### Indexed VCF Region Pushdown

Command:

```bash
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io --test filter_metadata indexed_vcf_region_filter_fetches_exact_start_and_end_positions
```

Result before fix:

```text
test indexed_vcf_region_filter_fetches_exact_start_and_end_positions ... FAILED
assertion `left == right` failed
  left: 4
 right: 2
```

The test asserts `candidate_variants == 2` for region `1:20-30` on a four-record indexed VCF. The pre-fix full scan counted all four records.

## Green Evidence

### Targeted Regression Tests

Commands:

```bash
pytest -q tests/test_filters.py::test_generator_variant_id_selection_is_consumed_once
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io --test filter_metadata
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-core --test filter_ir
pytest -q tests/test_filters.py
```

Results after fix:

```text
1 passed in 0.10s
4 passed; 0 failed
4 passed; 0 failed
15 passed in 0.12s
```

