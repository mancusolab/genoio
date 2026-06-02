# Phase 4 Baseline Review Remediation Evidence

## Issues Fixed

- `variants=` generator inputs were consumed during validation and then passed to Rust as an empty `id_in([])` filter.
- Indexed VCF region pushdown was missing; region reads used a full scan even when a tabix-indexed VCF was available.
- Phase 4 behavior changes now have red/green regression evidence.
- `NaN` rate thresholds constructed non-JSON-compliant filter IR instead of raising `InvalidOptionError`.
- Non-pushdown-safe region expressions on compressed unindexed VCFs were rejected instead of falling back to full-scan Rust evaluation.

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

## CLI/API Specialist Remediation

### NaN Rate Threshold Validation

Command:

```bash
pytest -q tests/test_filters.py
```

Result before fix:

```text
FAILED tests/test_filters.py::test_filter_constructors_reject_invalid_values[<lambda>6]
FAILED tests/test_filters.py::test_filter_constructors_reject_invalid_values[<lambda>7]
FAILED tests/test_filters.py::test_filter_constructors_reject_invalid_values[<lambda>9]
Failed: DID NOT RAISE <class 'genoio._errors.InvalidOptionError'>
```

### Compressed Unindexed Non-Pushdown Region Fallback

Command:

```bash
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io compressed_vcf_non_pushdown_region_filter_falls_back_to_full_scan
```

Result before fix:

```text
test compressed_vcf_non_pushdown_region_filter_falls_back_to_full_scan ... FAILED
Parse { message: "region filter on compressed VCF requires an index" }
```

### Focused Green Evidence

Commands:

```bash
pytest -q tests/test_filters.py
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io compressed_vcf_non_pushdown_region_filter_falls_back_to_full_scan
```

Results after fix:

```text
18 passed in 0.41s
test compressed_vcf_non_pushdown_region_filter_falls_back_to_full_scan ... ok
```
