# PLINK2 Speed Optimization Implementation Plan

**Goal:** Use packed genotype counts for genotype-stat filters so dropped variants do not require full float expansion.

**Architecture:** Add packed hard-call counting in `rust/genoio-io/src/plink2.rs` and reuse the existing `VariantStats` and `VariantFilter` contracts from `genoio-core`. Keep filter diagnostics and attached stats identical to the expanded path.

**Tech Stack:** Rust `VariantStats`, `VariantFilter`, packed PLINK2 decoder, existing filter tests, Python public filter expressions.

**Scope:** 7 phases from original design.

**Codebase verified:** 2026-06-03 Pacific/Honolulu.

---

## Acceptance Criteria Coverage

This phase implements and tests:

### plink2-speed.AC4: Packed genotype architecture is correct
- **plink2-speed.AC4.4 Success:** Packed genotype counts produce the same
  MAF/MAC/missingness/polymorphism filter decisions and attached stats as the
  current expanded path.

### plink2-speed.AC5: Public behavior remains stable
- **plink2-speed.AC5.2 Success:** Existing Python constructors, dataset methods,
  matrix orientation, dtype/missing behavior, and metadata alignment remain
  unchanged.
- **plink2-speed.AC5.3 Failure:** Unsupported dosage, multiallelic, or malformed
  PGEN records continue to fail explicitly rather than silently taking a fast
  path.

### Phase Evidence
- Post-review TDD/mutation evidence is recorded in
  `phase_07_tdd_evidence.md`.

---

<!-- START_TASK_1 -->
### Task 1: Compute VariantStats From Packed Genotypes

**Verifies:** plink2-speed.AC4.4, plink2-speed.AC5.2

**Files:**
- Modify: `rust/genoio-io/src/plink2.rs:1024-1161`
- Test: `rust/genoio-io/tests/plink2_dense.rs:80-197`

**Implementation:**
Add private method:

```rust
impl PackedGenotypes {
    fn stats_for_selected(&self, source_indices: &[usize]) -> Result<genoio_core::VariantStats>
}
```

Count hard-call categories directly:
- category `0`: homozygous reference,
- category `1`: heterozygous,
- category `2`: homozygous alternate,
- category `3`: missing.

Match `compute_variant_stats` semantics exactly:
- missing calls excluded from allele-frequency and MAC denominators,
- `missing_rate` based on all selected samples,
- `polymorphic` based on called genotypes,
- all-missing variants keep optional frequency/MAC fields as `None` where current behavior does.

If `compute_variant_stats` has subtle behavior, add a small shared core helper rather than duplicating divergent formulas. Keep any new helper in `rust/genoio-core/src/filter.rs` only if both VCF and PLINK2 can use it.

**Testing:**
Tests must verify each AC listed above:
- plink2-speed.AC4.4: packed stats match `compute_variant_stats` on representative fixed-width and variable-width fixtures.
- plink2-speed.AC5.2: attached variant stats are unchanged.

**Verification:**
Run: `cargo test -p genoio-io plink2_dense --test plink2_dense`
Expected: PLINK2 stats parity tests pass.

**Commit:** `perf: compute plink2 stats from packed genotypes`
<!-- END_TASK_1 -->

<!-- START_TASK_2 -->
### Task 2: Use Packed Stats In PLINK2 Filter Decisions

**Verifies:** plink2-speed.AC4.4, plink2-speed.AC5.2

**Files:**
- Modify: `rust/genoio-io/src/plink2.rs:97-176`
- Test: `rust/genoio-io/tests/filter_genotype_stats.rs:1-56`
- Test: `tests/test_filters.py`

**Implementation:**
In `read_plink2_dense_windowed`, when `PartialFilterDecision::NeedGenotypes` is returned:
- decode the current variant to packed genotypes,
- compute `VariantStats` from packed genotypes,
- evaluate the filter using those stats,
- attach stats to retained variants exactly as before,
- expand to float values only after the variant is retained and belongs in the requested window.

Keep sparse PLINK2 filtering semantically correct. If sparse still needs expanded values to build CSC output, it may expand retained variants after packed filter acceptance.

**Testing:**
Tests must verify each AC listed above:
- plink2-speed.AC4.4: MAF, MAC, missingness, and polymorphic filters retain the same PLINK2 variants and attached stats as the current expanded path.
- plink2-speed.AC5.2: diagnostics distinguish metadata and genotype drops exactly as before.

**Verification:**
Run: `cargo test -p genoio-io filter_genotype_stats --test filter_genotype_stats`
Expected: genotype-stat filter tests pass.

Run: `pytest tests/test_filters.py -q`
Expected: Python public filter tests pass.

**Commit:** `perf: filter plink2 variants with packed stats`
<!-- END_TASK_2 -->

<!-- START_TASK_3 -->
### Task 3: Preserve Explicit Failure Paths

**Verifies:** plink2-speed.AC5.3

**Files:**
- Modify: `rust/genoio-io/tests/plink2_dense.rs:133-144`
- Modify: `rust/genoio-io/tests/plink2_dense.rs:283-347`

**Implementation:**
Add or update tests that malformed and unsupported records still fail explicitly after packed counting:
- unsupported PGEN mode,
- unsupported compression type,
- LD-compressed record before non-LD state,
- non-increasing difflist sample IDs,
- variable-width block offset mismatch,
- truncated record.

Do not route malformed records through the count fast path if the normal packed decoder would reject them.

**Testing:**
Tests must verify each AC listed above:
- plink2-speed.AC5.3: unsupported or malformed records raise errors with source context instead of silently returning matrices.

**Verification:**
Run: `cargo test -p genoio-io plink2_dense --test plink2_dense`
Expected: malformed-record tests pass.

**Commit:** `test: preserve plink2 packed failure semantics`
<!-- END_TASK_3 -->

<!-- START_TASK_4 -->
### Task 4: Benchmark Genotype-Filtered Reads

**Verifies:** plink2-speed.AC4.4, plink2-speed.AC5.2, plink2-speed.AC5.3

**Files:**
- Modify: `.plans/implementation-plans/2026-06-03-plink2-speed/benchmark-baseline.md`
- Modify: `docs/performance.md:15-30`

**Implementation:**
Run the genotype-filtered benchmark scenario and compare it to the Phase 1 baseline. Record retained variant counts, diagnostics when available, and median time.

If genotype-filtered reads do not improve, document the dominant cost. Keep the packed-count implementation only if parity is strong and it does not regress other scenarios.

**Testing:**
Tests must verify each AC listed above:
- plink2-speed.AC4.4: genotype-filtered benchmark exercises packed counts.
- plink2-speed.AC5.2: retained variants and matrices match the public contract.
- plink2-speed.AC5.3: malformed fixture tests still pass after the fast path.

**Verification:**
Run: `cargo test`
Expected: all Rust tests pass.

Run: `pytest -q`
Expected: all Python tests pass.

Run: `python scripts/benchmark_plink2.py --scenario genotype-filtered --prefix data/chr22_hg38 --max-variants 1000 --repeats 5 --no-compare`
Expected: genotype-filtered benchmark completes and output is recorded.

**Commit:** `bench: record plink2 packed filter performance`
<!-- END_TASK_4 -->
