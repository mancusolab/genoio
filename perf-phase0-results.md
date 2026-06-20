# Phase 0 performance baseline

Date: 2026-06-19 21:39:03 PDT

## Environment

- Branch: `main`
- Commit: `f8553abc2fb6aea2635686bd2d034da07a96eb2b`
- OS: macOS 26.3.1 build 25D2128
- Kernel: `Darwin nmancuso861.local 25.3.0 Darwin Kernel Version 25.3.0: Wed Jan 28 20:53:31 PST 2026; root:xnu-12377.91.3~2/RELEASE_ARM64_T8103 arm64`
- CPU detail: unavailable from `sysctl` in this managed shell (`Operation not permitted`)
- Python: `Python 3.11.15`
- Rust: `rustc 1.95.0 (59807616e 2026-04-14)`
- Cargo: `cargo 1.95.0 (f2d3ce0bd 2026-03-21)`
- `genoio`: 0.1.0
- NumPy: 2.4.6
- SciPy: 1.17.1
- Polars: 1.41.2

## Fixture

Prefix: `data/chr22_hg38`

Files present:

- `data/chr22_hg38.pgen` - 126M
- `data/chr22_hg38.psam` - 76K
- `data/chr22_hg38.pvar.zst` - 41M
- `data/chr22_hg38.bgen` - 290M
- `data/chr22_hg38.bgen.bgi`
- `data/chr22_hg38.sample`
- `data/chr22_hg38.vcf.gz` - 430M
- `data/chr22_hg38.vcf.gz.tbi`
- `data/chr22_hg38.bed` - 815M
- `data/chr22_hg38.bim` - 31M
- `data/chr22_hg38.fam` - 67K

## Release build

Command:

```bash
make build-release
```

Result: succeeded. `maturin develop --release` built and installed
`genoio-0.1.0-cp311-cp311-macosx_11_0_arm64.whl`.

## Benchmark results

All successful runs used the release-mode extension installed by
`make build-release`.

### PLINK2

Command:

```bash
.venv/bin/python scripts/benchmark_plink2.py --backend genoio --scenario all --max-variants 1000 --repeats 7
```

| Scenario | Shape | Median | Min | Runs |
|---|---:|---:|---:|---|
| matrix-only | `(3202, 1000)` | 0.0093s | 0.0090s | 0.0090 0.0095 0.0098 0.0093 0.0101 0.0091 0.0091 |
| with-variants | `(3202, 1000)` | 0.0127s | 0.0124s | 0.0144 0.0127 0.0127 0.0132 0.0124 0.0126 0.0142 |
| sample-filtered | `(1601, 1000)` | 0.0090s | 0.0086s | 0.0090 0.0089 0.0086 0.0109 0.0087 0.0097 0.0116 |
| genotype-filtered | `(3202, 1000)` | 0.5871s | 0.5566s | 0.5971 0.5948 0.5871 0.5566 0.6924 0.5688 0.5681 |

Variant metadata length for `with-variants`: 1000.

Command:

```bash
.venv/bin/python scripts/benchmark_plink2.py --backend genoio --scenario matrix-only --max-variants 10000 --repeats 7 --no-compare
```

| Scenario | Shape | Median | Min | Runs |
|---|---:|---:|---:|---|
| matrix-only | `(3202, 10000)` | 0.0960s | 0.0918s | 0.1103 0.1186 0.1055 0.0933 0.0922 0.0918 0.0960 |

Skipped comparison backend:

- `pgenlib` is not importable in `.venv`.
- `/Users/nicholas/Projects/plink-ng/2.0/Python` contains source files, but no importable built extension was found.

Skipped PLINK2 haplotype commands:

- `--kind haplo-hardcall` failed on `data/chr22_hg38` with `unphased pgen hardcall record retained in haplotype read`.
- `--kind haplo-dosage` failed on `data/chr22_hg38` with `pgen record does not contain explicit phased dosage values`.
- No large local phased PLINK2 benchmark prefix was found under this repository.

### BGEN

Command:

```bash
.venv/bin/python scripts/benchmark_bgen.py --backend genoio --scenario all --max-variants 1000 --repeats 7
```

| Scenario | Shape | Median | Min | Runs |
|---|---:|---:|---:|---|
| matrix-only | `(3202, 1000)` | 0.0316s | 0.0284s | 0.0316 0.0357 0.0339 0.0284 0.0395 0.0295 0.0294 |
| with-variants | `(3202, 1000)` | 0.0337s | 0.0303s | 0.0303 0.0342 0.0327 0.0359 0.0337 0.0366 0.0319 |
| sample-filtered | `(1601, 1000)` | 0.0284s | 0.0197s | 0.0284 0.0232 0.0365 0.0197 0.0328 0.0300 0.0254 |
| genotype-filtered | `(3202, 1000)` | 0.2242s | 0.1458s | 0.2242 0.2273 0.2208 0.1699 0.3037 0.2801 0.1458 |
| indexed-region | `(3202, 1000)` | 0.0868s | 0.0713s | 0.0713 0.0738 0.0831 0.0868 0.1242 0.0995 0.0886 |

Variant metadata length for `with-variants`: 1000.

Command:

```bash
.venv/bin/python scripts/benchmark_bgen.py --backend genoio --scenario matrix-only --max-variants 10000 --repeats 7
```

| Scenario | Shape | Median | Min | Runs |
|---|---:|---:|---:|---|
| matrix-only | `(3202, 10000)` | 0.4605s | 0.3879s | 0.4247 0.5275 0.4605 0.5318 0.4167 0.4826 0.3879 |

Skipped comparison backend:

- `bgen_reader` is not installed in `.venv`.

### VCF

Command:

```bash
.venv/bin/python scripts/benchmark_vcf.py --backend genoio --max-variants 1000 --repeats 7
```

| Scenario | Shape | Median | Min | Runs |
|---|---:|---:|---:|---|
| matrix-only | `(3202, 1000)` | 0.0422s | 0.0415s | 0.0418 0.0417 0.0415 0.0422 0.0429 0.0432 0.0442 |

Skipped comparison backend:

- `cyvcf2` is not installed in `.venv`.

## Phase 0 status

Completed:

- Release-mode build succeeded.
- Benchmark environment and fixture details were recorded.
- Genoio release-mode baselines were recorded for PLINK2, BGEN, and VCF.
- Optional comparison backend gaps were recorded.

Open:

- Add a machine-readable benchmark output mode, such as `--json`, before repeated
  branch comparisons.
- Install or build optional comparison backends if cross-library comparisons are
  required: `pgenlib`, `bgen_reader`, and `cyvcf2`.
- Provide a large explicit phased PLINK2 benchmark prefix before treating
  haplotype PLINK2 timings as part of the baseline.
