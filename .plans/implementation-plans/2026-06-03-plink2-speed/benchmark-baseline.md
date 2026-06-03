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
unchanged, which confirms those scenarios did not accidentally use the
matrix-only fast path.
