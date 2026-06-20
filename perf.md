# genoio-io performance plan

This plan targets Rust reader performance in `genoio-io` without changing the
public Python API. The main goal is to make future speed work measurable and
low-risk. Organization cleanup is included only where it lowers the cost of
profiling or reduces duplicate hot-path logic.

## Ground rules

- Measure with release builds before and after each change.
- Keep PyO3 out of `genoio-io`; Python boundary work stays in `genoio-py`.
- Do not add `unsafe` for speed unless profiling proves the need and the safety
  contract is documented.
- Keep existing matrix, metadata, filter, missing-value, and ordering contracts.
- Prefer small refactors that keep decode loops readable.

## Scope guardrails

VCF and BCF are out of scope for this plan unless release benchmarks or profiles
identify a specific bottleneck. Their current structure is mostly sound: BCF
already shares dense GT/DS traversal, and text VCF keeps hot noodles record
loops explicit. If VCF/BCF becomes a target, write a separate focused plan with
its own baseline and acceptance criteria.

## Phase 0: establish the baseline

Use the existing benchmark scripts and local `data/chr22_hg38` fixture.

1. Build release mode:

   ```bash
   make build-release
   ```

2. Record the commit, CPU, OS, Python version, Rust version, and fixture path.

3. Run the current baseline:

   ```bash
   python scripts/benchmark_plink2.py --scenario all --max-variants 1000 --repeats 7
   python scripts/benchmark_plink2.py --scenario matrix-only --max-variants 10000 --repeats 7 --no-compare
   python scripts/benchmark_plink2.py --kind haplo-hardcall --scenario matrix-only --backend genoio --max-variants 1000 --repeats 7
   python scripts/benchmark_plink2.py --kind haplo-dosage --scenario matrix-only --backend genoio --max-variants 1000 --repeats 7
   python scripts/benchmark_bgen.py --scenario all --max-variants 1000 --repeats 7
   python scripts/benchmark_bgen.py --scenario matrix-only --backend both --max-variants 10000 --repeats 7
   python scripts/benchmark_vcf.py --max-variants 1000 --repeats 7
   ```

4. Add a lightweight result format if needed. The current scripts print human
   output; a `--json` option would make branch comparisons less error-prone.

Acceptance criteria:

- Baseline results are saved with commit hashes.
- Benchmarks run against release-mode Rust.
- The comparison method can be repeated on `main` and the feature branch.

## Phase 1: reduce PGEN refactor risk

`rust/genoio-io/src/plink/plink2/pgen.rs` is the largest maintenance and
performance risk. Split it by responsibility before doing deeper hot-path work.

Proposed modules:

- `pgen/header.rs`: header parsing, mode validation, record offsets.
- `pgen/io.rs`: file seeks, fixed-width record reads, variable-width record reads.
- `pgen/main_track.rs`: hardcall main-track decode.
- `pgen/dosage_track.rs`: dosage overlays and fixed-width dosage decode.
- `pgen/haplotype_track.rs`: phase and haplotype dosage decode.
- `pgen/bitpack.rs`: bit helpers, varints, sample-id width, bounds checks.

Cleanups to include:

- Share full-header and prefix-header parsing instead of keeping parallel match
  arms.
- Keep one implementation of one-bit record decode.
- Remove duplicate bit helpers.
- Add short module docs that state each module's binary-format responsibility.

Acceptance criteria:

- No behavior changes.
- `make rust-check` and `make rust-test` pass.
- PLINK2 release benchmarks stay within 3% of baseline median unless variance
  explains the difference.

Result: completed in `perf-phase1-results.md`.

## Phase 2: profile and improve PGEN genotype-filtered reads

Current benchmark data shows PLINK2 genotype-filtered reads as the largest
PLINK2 cost surface. Profile that case before changing the decoder.

Focus areas:

- Replace `decode_difflist` allocation with a visitor-style path where it shows
  up in profiles.
- Reuse decoded first-id buffers in variable-width records where lifetimes are
  clear.
- Check whether packed hardcall expansion does avoidable work for sample-filtered
  reads.
- Keep LD-compressed record state explicit; correctness matters more than a small
  allocation win.

Acceptance criteria:

- Flamegraph or sampled profile points to the changed code.
- Matrix parity tests cover fixed-width, variable-width, LD-compressed, one-bit,
  and difflist records.
- PLINK2 genotype-filtered median improves without regressing matrix-only,
  sample-filtered, haplotype hardcall, or haplotype dosage cases by more than 3%.

Result: completed in `perf-phase2-results.md`.

## Phase 3: give BGEN the same structure as PLINK2

`rust/genoio-io/src/bgen.rs` repeats the retained-variant loop across
sequential/indexed and dosage/haplotype paths. Split organization first, then
measure.

Proposed modules:

- `bgen/session.rs`: `BgenReadSession`, cursor types, index position checks.
- `bgen/dense.rs`: diploid dosage dense orchestration.
- `bgen/haplotype.rs`: phased haplotype dense orchestration.
- `bgen/filter.rs`: dosage filter counts and filter evaluation helpers.

Potential performance work after the split:

- Share retention and index-validation flow without hiding decode work.
- Keep matrix-only unfiltered reads on their current fast path.
- Reuse probability and haplotype decode buffers across variants.
- Measure whether indexed-region reads pay extra validation or seek costs.

Acceptance criteria:

- BGEN dense, haplotype, indexed-region, and genotype-filtered tests pass.
- BGEN matrix-only and indexed-region release medians stay within 3% of baseline.
- Any BGEN speed claim includes before/after medians and fixture details.

Result: completed in `perf-phase3-results.md`.

## Phase 4: add targeted decoder tests

Add tests where refactors would otherwise rely on benchmark fixtures alone.

Test targets:

- PGEN varint bounds and malformed input errors.
- PGEN difflist ordering, duplicate sample IDs, out-of-range sample IDs, and
  packed value decoding.
- PGEN one-bit records with and without extra difflist entries.
- BGEN dosage-filter count invariants for missing, monomorphic, polymorphic,
  MAC, MAF, and missing-rate filters.
- Matrix parity for sample-filtered and matrix-only paths.

Property tests are useful for pure helpers if they stay fast and deterministic.

## Phase 5: update public performance docs

After each measured improvement, update `docs/performance.md` with:

- Commit or branch used.
- Fixture path and workload.
- Release build command.
- Median, min, and repeat count.
- Comparison backend version when used.

Do not update headline numbers from debug builds or one-off timing runs.

## Verification checklist for each performance PR

- `make rust-fmt`
- `make rust-check`
- `make rust-test`
- `make build-release`
- Relevant release benchmark commands from Phase 0
- Before/after summary with medians and percentage change
- Notes on any benchmark variance or skipped comparison backend
