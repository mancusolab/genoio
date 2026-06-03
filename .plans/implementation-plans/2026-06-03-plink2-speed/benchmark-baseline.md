# PLINK2 Baseline Benchmark

## Provenance

- Date: 2026-06-03 09:12:22 HST
- Git commit: `257554e8711c80ede05ce8fe10eadb92d36e3b7e`
- Build status: Rust extension rebuilt in release mode with
  `env CC=clang AR=ar python -m maturin develop --release`
- Machine note: benchmark script did not emit a machine note. Local machine
  recorded separately as `macOS-26.3.1-arm64-arm-64bit`, `arm64`;
  `uname -a` reported `Darwin nmancuso861.local 25.3.0 ... RELEASE_ARM64_T8103 arm64`.
- Data prefix: `data/chr22_hg38`
- Optional comparison dependency: `pgenlib` importable from
  `/Users/nicholas/micromamba/lib/python3.11/site-packages/pgenlib.cpython-311-darwin.so`
- `.pvar` handling: local fixture has `chr22_hg38.pvar.zst`; `zstd` was available
  at `/Users/nicholas/micromamba/bin/zstd`.

## TDD Evidence

Task 1 introduced behavior-changing benchmark CLI coverage.

Red command:

```bash
pytest -q tests/test_benchmark_plink2_cli.py
```

Red result before implementation: failed because `--scenario` was unrecognized
and `scripts/benchmark_plink2.py` did not expose `read_genoio_matrix_only`.

Green command:

```bash
pytest -q tests/test_benchmark_plink2_cli.py
```

Green result after implementation: `3 passed`.

Task 2 added PLINK2 public contract tests only. The new Python and Rust tests
passed against the existing implementation, so no production parser change was
required. The Rust contract-test run initially needed an explicit toolchain
environment because `ar` was not otherwise available:

```bash
env CC=clang AR=ar cargo test
```

## Matrix-Only With pgenlib Comparison

Command:

```bash
python scripts/benchmark_plink2.py --scenario matrix-only --prefix data/chr22_hg38 --max-variants 1000 --repeats 5
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0220s min=0.0208s runs=0.0208 0.0220 0.0239 0.0215 0.0239
pgenlib_pgenreader
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0074s min=0.0070s runs=0.0094 0.0074 0.0070 0.0073 0.0075
comparison
  genoio_plink2_matrix_only.shape=(3202, 1000) pgenlib_pgenreader.shape=(3202, 1000)
  allclose=True max_abs_diff=0
```

Same-variant-count matrix-only summary:

| Reader | Median | Shape | Sum | Missing |
| --- | ---: | --- | ---: | ---: |
| `genoio_plink2_matrix_only` | 0.0220 s | `(3202, 1000)` | 72577 | 0 |
| `pgenlib_pgenreader` | 0.0074 s | `(3202, 1000)` | 72577 | 0 |

`pgenlib` was 2.97x faster for this matrix-only 1,000-variant read.

## All Scenarios

Command:

```bash
python scripts/benchmark_plink2.py --scenario all --prefix data/chr22_hg38 --max-variants 1000 --repeats 5 --no-compare
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0220s min=0.0216s runs=0.0219 0.0216 0.0224 0.0222 0.0220
pgenlib_pgenreader
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0070s min=0.0068s runs=0.0078 0.0068 0.0070 0.0074 0.0069
genoio_plink2_with_variants
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0232s min=0.0229s runs=0.0266 0.0262 0.0232 0.0230 0.0229
  variant_metadata length=1000
skipped pgenlib comparison for with-variants: pgenlib does not provide the same metadata/filter contract
genoio_plink2_sample_filtered
  matrix shape=(1601, 1000) dtype=float32 sum=34913 missing=0
  time median=0.0158s min=0.0152s runs=0.0159 0.0154 0.0158 0.0152 0.0162
skipped pgenlib comparison for sample-filtered: pgenlib does not provide the same metadata/filter contract
genoio_plink2_genotype_filtered
  matrix shape=(3202, 1000) dtype=float32 sum=447988 missing=0
  time median=16.1535s min=16.0589s runs=16.1535 16.2978 16.1233 16.0589 18.3733
skipped pgenlib comparison for genotype-filtered: pgenlib does not provide the same metadata/filter contract
```

Scenario summary:

| Scenario | Reader | Median | Shape | Notes |
| --- | --- | ---: | --- | --- |
| Matrix-only | `genoio_plink2_matrix_only` | 0.0220 s | `(3202, 1000)` | Dense matrix block only |
| Matrix-only | `pgenlib_pgenreader` | 0.0070 s | `(3202, 1000)` | Same variant count; comparison disabled only for value diff |
| Metadata-returning | `genoio_plink2_with_variants` | 0.0232 s | `(3202, 1000)` | `variant_metadata length=1000` |
| Sample-filtered | `genoio_plink2_sample_filtered` | 0.0158 s | `(1601, 1000)` | First half of source samples retained |
| Genotype-filtered | `genoio_plink2_genotype_filtered` | 16.1535 s | `(3202, 1000)` | `maf(min=0.01)` filter |

For non-matrix-only scenarios, pgenlib comparison was intentionally skipped by
the benchmark script because pgenlib does not provide the same metadata/filter
contract.

## Post Metadata Short-Circuit Benchmark

Provenance:

- Date: 2026-06-03 09:40:34 HST
- Git commit before documentation commit: `bd1a4b2a679a0ea9c5fbc9b81dc8eb1d7505c09a`
- Build status: Rust extension rebuilt in release mode with
  `env CC=clang AR=ar python -m maturin develop --release`
- Machine note: benchmark script did not emit a machine note. Local machine
  recorded separately as `macOS-26.3.1-arm64-arm-64bit`, `arm64`;
  `uname -a` reported `Darwin nmancuso861.local 25.3.0 ... RELEASE_ARM64_T8103 arm64`.
- Data prefix: `data/chr22_hg38`
- Optional comparison dependency: `pgenlib` importable from
  `/Users/nicholas/micromamba/lib/python3.11/site-packages/pgenlib.cpython-311-darwin.so`
- `.pvar` handling: local fixture has `chr22_hg38.pvar.zst`; `zstd` was available
  at `/Users/nicholas/micromamba/bin/zstd`.

Command:

```bash
python scripts/benchmark_plink2.py --scenario all --prefix data/chr22_hg38 --max-variants 1000 --repeats 5 --no-compare
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0178s min=0.0173s runs=0.0178 0.0174 0.0178 0.0173 0.0178
pgenlib_pgenreader
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0068s min=0.0065s runs=0.0080 0.0067 0.0065 0.0069 0.0068
genoio_plink2_with_variants
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0229s min=0.0226s runs=0.0253 0.0239 0.0226 0.0228 0.0229
  variant_metadata length=1000
skipped pgenlib comparison for with-variants: pgenlib does not provide the same metadata/filter contract
genoio_plink2_sample_filtered
  matrix shape=(1601, 1000) dtype=float32 sum=34913 missing=0
  time median=0.0157s min=0.0149s runs=0.0157 0.0153 0.0149 0.0160 0.0161
skipped pgenlib comparison for sample-filtered: pgenlib does not provide the same metadata/filter contract
genoio_plink2_genotype_filtered
  matrix shape=(3202, 1000) dtype=float32 sum=447988 missing=0
  time median=16.0566s min=15.8323s runs=16.0243 16.1316 16.0566 15.8323 16.0804
skipped pgenlib comparison for genotype-filtered: pgenlib does not provide the same metadata/filter contract
```

Before/after matrix-only median comparison:

| Scenario | Phase 1 baseline median | Post-short-circuit median | Change |
| --- | ---: | ---: | ---: |
| Matrix-only | 0.0220 s | 0.0178 s | 1.24x faster |
| Metadata-returning | 0.0232 s | 0.0229 s | 1.01x faster |
| Sample-filtered | 0.0158 s | 0.0157 s | 1.01x faster |
| Genotype-filtered | 16.1535 s | 16.0566 s | 1.01x faster |

The matrix-only path is the benchmark scenario that exercises the metadata
short-circuit. The metadata-returning scenario still returned
`variant_metadata length=1000`, confirming that metadata paths continue to
work. Sample-filtered and genotype-filtered medians remained effectively
unchanged, which is consistent with those scenarios not using the matrix-only
fast path; gating is covered by the contract tests.

## Phase 3 Scratch-Reuse and Sequential-Read Benchmark

Provenance:

- Date: 2026-06-03 10:01:09 HST
- Git commit before documentation commit: `a6e06c311a160c434b637658bf3776e2c9b59805`
- Build status: Rust extension rebuilt in release mode with
  `env CC=clang AR=ar python -m maturin develop --release`
- Machine note: benchmark script did not emit a machine note. Local machine
  recorded separately as `macOS-26.3.1-arm64-arm-64bit`, `arm64`;
  `uname -a` reported `Darwin nmancuso861.local 25.3.0 ... RELEASE_ARM64_T8103 arm64`.
- Data prefix: `data/chr22_hg38`
- Optional comparison dependency: `pgenlib` importable from
  `/Users/nicholas/micromamba/lib/python3.11/site-packages/pgenlib.cpython-311-darwin.so`
- `.pvar` handling: local fixture has `chr22_hg38.pvar.zst`; `zstd` was available
  at `/Users/nicholas/micromamba/bin/zstd`.

Command:

```bash
python scripts/benchmark_plink2.py --scenario all --prefix data/chr22_hg38 --max-variants 1000 --repeats 5 --no-compare
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0148s min=0.0144s runs=0.0144 0.0155 0.0148 0.0146 0.0152
pgenlib_pgenreader
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0065s min=0.0064s runs=0.0080 0.0065 0.0066 0.0064 0.0064
genoio_plink2_with_variants
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0185s min=0.0182s runs=0.0226 0.0197 0.0182 0.0185 0.0185
  variant_metadata length=1000
skipped pgenlib comparison for with-variants: pgenlib does not provide the same metadata/filter contract
genoio_plink2_sample_filtered
  matrix shape=(1601, 1000) dtype=float32 sum=34913 missing=0
  time median=0.0129s min=0.0129s runs=0.0129 0.0129 0.0129 0.0130 0.0130
skipped pgenlib comparison for sample-filtered: pgenlib does not provide the same metadata/filter contract
genoio_plink2_genotype_filtered
  matrix shape=(3202, 1000) dtype=float32 sum=447988 missing=0
  time median=13.2681s min=13.0629s runs=13.3502 13.1876 13.0629 13.2944 13.2681
skipped pgenlib comparison for genotype-filtered: pgenlib does not provide the same metadata/filter contract
```

Phase 1 baseline to Phase 3 median comparison:

| Scenario | Phase 1 baseline median | Phase 3 median | Change |
| --- | ---: | ---: | ---: |
| Matrix-only | 0.0220 s | 0.0148 s | 1.49x faster |
| Metadata-returning | 0.0232 s | 0.0185 s | 1.25x faster |
| Sample-filtered | 0.0158 s | 0.0129 s | 1.22x faster |
| Genotype-filtered | 16.1535 s | 13.2681 s | 1.22x faster |

The scratch-reuse and sequential fixed-width read changes are visible on the
small fixture across the non-genotype-filter hot paths: matrix-only, metadata-
returning, and sample-filtered scenarios all reduced their medians relative to
the Phase 1 baseline. The genotype-filtered scenario also improved, but its
absolute runtime is still dominated by filter evaluation over genotypes: it
remains a 13-second path while the direct matrix and metadata scenarios are
below 20 milliseconds. The remaining gap to `pgenlib` in matrix-only reads is
still visible (`pgenlib` median 0.0065 s versus `genoio` median 0.0148 s).

## Phase 4 Direct Sample-Major Fill Decision

Provenance:

- Date: 2026-06-03 10:41:33 HST
- Git commit before documentation commit: `135f3d2a165882263ed3520872998473bfd9615b`
- Build status: Rust extension rebuilt in release mode with
  `env CC=clang AR=ar python -m maturin develop --release`
- Direct-fill status: current build uses direct sample-major fill for
  unfiltered source windows in the matrix-only dense path.
- Data prefix: `data/chr22_hg38`
- Benchmark script label decision: no script label was added. The benchmark
  scenario label remains `genoio_plink2_matrix_only`; this artifact records the
  implementation state and decision without changing benchmark output format.

Required 1,000-variant command:

```bash
python scripts/benchmark_plink2.py --scenario matrix-only --prefix data/chr22_hg38 --max-variants 1000 --repeats 5 --no-compare
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0147s min=0.0147s runs=0.0147 0.0147 0.0147 0.0150 0.0149
pgenlib_pgenreader
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0066s min=0.0065s runs=0.0079 0.0067 0.0065 0.0066 0.0066
```

Repeat 1,000-variant command used to classify the small absolute delta against
the Phase 3 baseline:

```bash
python scripts/benchmark_plink2.py --scenario matrix-only --prefix data/chr22_hg38 --max-variants 1000 --repeats 5 --no-compare
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0145s min=0.0144s runs=0.0151 0.0148 0.0144 0.0145 0.0145
pgenlib_pgenreader
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0066s min=0.0064s runs=0.0077 0.0064 0.0064 0.0066 0.0066
```

Required 10,000-variant command:

```bash
python scripts/benchmark_plink2.py --scenario matrix-only --prefix data/chr22_hg38 --max-variants 10000 --repeats 5 --no-compare
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 10000) dtype=float32 sum=802432 missing=0
  time median=0.2746s min=0.2708s runs=0.2753 0.2708 0.2718 0.2766 0.2746
pgenlib_pgenreader
  matrix shape=(3202, 10000) dtype=float32 sum=802432 missing=0
  time median=0.0467s min=0.0465s runs=0.0476 0.0467 0.0467 0.0465 0.0522
```

Direct-fill comparison:

| Window | Phase 3 `genoio` median | Phase 4 direct-fill `genoio` median | `pgenlib` median | Decision signal |
| --- | ---: | ---: | ---: | --- |
| 1,000 variants | 0.0148 s | 0.0145-0.0147 s | 0.0066 s | Neutral to slightly improved across repeat runs |
| 10,000 variants | Not recorded | 0.2746 s | 0.0467 s | New scale evidence; pgenlib remains materially faster |

Decision: keep the direct sample-major fill path. The 1,000-variant direct-fill
measurements are neutral to slightly improved relative to the Phase 3
scratch-reuse baseline while removing the variant-major accumulation and
transpose copy for unfiltered matrix-only source windows. Proceed to packed
batch transpose work because `pgenlib` remains materially faster: 2.16x faster
on the best direct-fill 1,000-variant run and 5.88x faster at 10,000 variants.

## Phase 6 Packed Batch Transpose Decision

Provenance:

- Date: 2026-06-03 11:33:16 HST
- Git commit before documentation commit: `3bb767085c43c8a39687fa93e4b238c305d3c5bc`
- Build status: Rust extension rebuilt in release mode with
  `env CC=clang AR=ar python -m maturin develop --release`
- Packed-batch status: current build uses private packed variant batches for
  unfiltered dense source-window construction in both matrix-only and
  metadata-bearing paths.
- Data prefix: `data/chr22_hg38`
- Batch size: 64 packed variants.
- Decision: keep the packed-batch production source-window path.

Required 1,000-variant command:

```bash
python scripts/benchmark_plink2.py --scenario matrix-only --prefix data/chr22_hg38 --max-variants 1000 --repeats 5
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0100s min=0.0095s runs=0.0095 0.0109 0.0100 0.0101 0.0095
pgenlib_pgenreader
  matrix shape=(3202, 1000) dtype=float32 sum=72577 missing=0
  time median=0.0069s min=0.0068s runs=0.0102 0.0073 0.0069 0.0068 0.0068
comparison
  genoio_plink2_matrix_only.shape=(3202, 1000) pgenlib_pgenreader.shape=(3202, 1000)
  allclose=True max_abs_diff=0
```

Required 10,000-variant command:

```bash
python scripts/benchmark_plink2.py --scenario matrix-only --prefix data/chr22_hg38 --max-variants 10000 --repeats 5 --no-compare
```

Output:

```text
genoio_plink2_matrix_only
  matrix shape=(3202, 10000) dtype=float32 sum=802432 missing=0
  time median=0.0941s min=0.0926s runs=0.1121 0.0935 0.0941 0.0926 0.0950
pgenlib_pgenreader
  matrix shape=(3202, 10000) dtype=float32 sum=802432 missing=0
  time median=0.0497s min=0.0472s runs=0.0472 0.0491 0.0508 0.0504 0.0497
```

Packed-batch comparison:

| Window | Phase 4 direct-fill `genoio` median | Phase 6 packed-batch `genoio` median | `pgenlib` median | Decision signal |
| --- | ---: | ---: | ---: | --- |
| 1,000 variants | 0.0145-0.0147 s | 0.0100 s | 0.0069 s | Improved; remaining pgenlib gap is 1.45x |
| 10,000 variants | 0.2746 s | 0.0941 s | 0.0497 s | Improved; remaining pgenlib gap is 1.89x |

Decision rationale: packed batches materially improve the unfiltered matrix-only
source-window path at both required window sizes while preserving exact matrix
agreement with `pgenlib` on the 1,000-variant comparison. Keep the production
path and retain the parity tests.

## Phase 7 Packed Genotype-Stat Filtering Decision

Provenance:

- Date: 2026-06-03 12:04:59 HST
- Git commit before documentation commit: `c673ec535a566f8f24b53147d7022ac654c73431`
- Build status: Rust extension rebuilt in release mode with
  `env CC=clang AR=ar python -m maturin develop --release`
- Packed genotype-stat status: current build computes PLINK2 genotype-filter
  statistics from packed hard-call categories before expanding retained
  variants to float matrix values.
- Data prefix: `data/chr22_hg38`
- Decision: keep the packed-count production path.

Required genotype-filtered command:

```bash
python scripts/benchmark_plink2.py --scenario genotype-filtered --prefix data/chr22_hg38 --max-variants 1000 --repeats 5 --no-compare
```

Output:

```text
genoio_plink2_genotype_filtered
  matrix shape=(3202, 1000) dtype=float32 sum=447988 missing=0
  time median=7.8559s min=7.6775s runs=7.8559 8.1471 7.9423 7.7117 7.6775
skipped pgenlib comparison for genotype-filtered: pgenlib does not provide the same metadata/filter contract
```

Genotype-filtered comparison:

| Scenario | Phase 1 baseline median | Phase 6 median | Phase 7 packed-count median | Matrix summary |
| --- | ---: | ---: | ---: | --- |
| Genotype-filtered, 1,000 retained variants | 16.1535 s | Not remeasured | 7.8559 s | shape `(3202, 1000)`, sum `447988`, missing `0` |

Retained-count signal: the benchmark matrix shape remained `(3202, 1000)`,
which is the same retained variant count and sample count previously recorded
for the genotype-filtered scenario. The matrix sum and missing count also match
the Phase 1 and Phase 3 genotype-filtered benchmark summaries, so the packed
stats path preserved the public retained-matrix contract on this fixture.

Decision rationale: keep the packed-count implementation. Targeted Rust and
Python filter tests covered MAF, MAC, missingness, polymorphic decisions,
attached stats, sparse filtering, and explicit malformed-record failures. The
benchmark median improved from the Phase 1 genotype-filtered baseline of
16.1535 s to 7.8559 s while retaining the same matrix shape, sum, and missing
count. The dominant remaining cost is scanning enough source variants and
building the retained matrix to produce 1,000 variants after the genotype
filter; the packed-count path removes the previous full float expansion cost
for variants dropped by genotype statistics.
