# Phase 7 TDD and Mutation Evidence

Date: 2026-06-03 Pacific/Honolulu

## Scope

Phase 7 changed PLINK2 parser/filter behavior by computing genotype-filter
statistics from packed hard-call categories before expanding retained variants.

The original Phase 7 commit history did not preserve a clean red-before-green
TDD sequence: the first implementation commit combined production changes and
tests. This file records honest post-review evidence instead of claiming commit
history proves the original red step.

## Regression Added Post-Review

Added `filter_genotype_stats_plink2_variable_width_selected_samples_attach_stats`
in `rust/genoio-io/tests/filter_genotype_stats.rs`.

The test covers:
- variable-width PGEN records;
- selected samples requested out of source order and returned in source order;
- a missing call in the retained filtered output;
- MAF and missing-rate filter decisions from packed stats;
- attached `af`, `maf`, `mac`, `missing_rate`, and `n_called` metadata.

Initial green check against the current implementation:

```text
Command: env CC=clang AR=ar cargo test -p genoio-io filter_genotype_stats --test filter_genotype_stats
Exit: 0
Outcome: 4 tests passed, including filter_genotype_stats_plink2_variable_width_selected_samples_attach_stats.
```

## Mutation Check: Packed Missing Calls Broken

Temporary mutation applied to `rust/genoio-io/src/plink2.rs`:

```diff
-                3 => missing_count += 1,
+                3 => hom_ref_count += 1,
```

This deliberately made packed category `3` count as homozygous reference instead
of missing.

Red check:

```text
Command: env CC=clang AR=ar cargo test -p genoio-io filter_genotype_stats --test filter_genotype_stats
Exit: 101
Outcome: 1 passed, 3 failed.

Failed tests:
- filter_genotype_stats_plink2_match_expanded_stats_and_attach_metadata
  - attached missing_rate changed from expected Some(0.5) to Some(0.33333334)
- filter_genotype_stats_plink2_variable_width_selected_samples_attach_stats
  - attached missing_rate changed from expected Some(0.5) to Some(0.25)
- filter_genotype_stats_plink2_sparse_keeps_dense_filter_semantics
  - sparse path attempted to retain a variant with missing values and failed
    with "sparse missing values are not stored in this release"
```

The new variable-width selected-sample test fails specifically when packed
missing-call accounting is broken.

## Restoration Green Check

Restored `rust/genoio-io/src/plink2.rs`:

```diff
-                3 => hom_ref_count += 1,
+                3 => missing_count += 1,
```

Green check:

```text
Command: env CC=clang AR=ar cargo test -p genoio-io filter_genotype_stats --test filter_genotype_stats
Exit: 0
Outcome: 4 tests passed.

Passing tests:
- filter_genotype_stats_plink2_match_expanded_stats_and_attach_metadata
- filter_genotype_stats_plink2_variable_width_selected_samples_attach_stats
- filter_genotype_stats_plink2_sparse_keeps_dense_filter_semantics
- filter_genotype_stats_use_called_genotypes_before_missing_imputation
```

## Conclusion

The post-review mutation check demonstrates that the packed-stat/filter tests
detect incorrect packed missing-call handling and return green after restoration.
This is durable evidence for Phase 7 review closure, not a retroactive claim
about the original commit sequence.
