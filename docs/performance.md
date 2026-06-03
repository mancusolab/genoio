# Performance

`genoio` builds matrices in Rust and crosses the Python boundary once per read
or block. This avoids Python loops over variants and samples.

The current local benchmark was run on an Apple Silicon M1 Mac (`arm64`) with
Python 3.11. The Rust extension was built in release mode. Each run reads the
first 1,000 variants and constructs a dense `float32` matrix with shape
`(3202, 1000)`.

| Source | genoio median | comparison median | Result |
|---|---:|---:|---|
| VCF vs `cyvcf2` | 0.1064 s | 0.2383 s | genoio 2.24x faster |
| PLINK1 vs `pandas_plink` | 0.3968 s | 2.4382 s | genoio 6.15x faster |
| PLINK2 matrix-only vs `pgenlib` | 0.0148 s | 0.0065 s | `pgenlib` 2.28x faster |

These numbers are workload- and machine-dependent. Treat them as local
benchmarks, not universal claims.

---

## Benchmark data

The benchmark scripts default to `data/chr22_hg38`, but that directory is a
local fixture and isn't distributed with the repository.

The fixture comes from the PLINK 2
[1000 Genomes phase 3 hg38 resources](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg).
The chromosome 22 PLINK 2 files were used as the source, then converted with
`plink2` to VCF and PLINK1 `.bed/.bim/.fam` files. That keeps the comparisons
focused on reader behavior rather than differences in samples or variants. The
VCF header records `##source=PLINKv2.0`.

For the PLINK2 timing above, the raw Phase 1 baseline, post-short-circuit
benchmark, and Phase 3 scratch-reuse/sequential-read benchmark are recorded in
`.plans/implementation-plans/2026-06-03-plink2-speed/benchmark-baseline.md`.
The Phase 3 PLINK2 run was collected on 2026-06-03 at git commit
`a6e06c311a160c434b637658bf3776e2c9b59805`, after rebuilding the Rust extension
in release mode with `env CC=clang AR=ar python -m maturin develop --release`,
on macOS arm64. The same run measured PLINK2 metadata-returning,
sample-filtered, and genotype-filtered scenarios with 1,000 variants. The
matrix-only median improved from the Phase 1 baseline of 0.0220 s to 0.0148 s;
the metadata-returning median was 0.0185 s, the sample-filtered median was
0.0129 s, and the genotype-filtered median was 13.2681 s.

For the timings above, the VCF, PLINK1, and PLINK2 matrix-only comparisons
produced the same matrix summary:
shape `(3202, 1000)`, `float32` dtype, genotype sum `72577`, zero missing
values, and exact agreement with the comparison reader.

---

## Run local benchmarks

Build the Rust extension in release mode before benchmarking:

```bash
python -m maturin develop --release
```

Then run the benchmark scripts:

```bash
python scripts/benchmark_vcf.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink2.py --scenario all --max-variants 1000 --repeats 5
python scripts/benchmark_plink2.py --scenario matrix-only --max-variants 10000 --repeats 5 --no-compare
```

Optional comparison packages are used when installed:

- `cyvcf2` for VCF
- `pandas_plink` for PLINK1
- `pgenlib` for PLINK2

---

## What affects speed

Metadata filters are cheaper than genotype filters because they can run before
matrix decoding. Region filters on indexed compressed VCF/BCF sources can also
avoid scanning unrelated records.

PLINK2 performance depends on how much metadata the caller requests. Block
reads that only need matrix data can avoid parsing full variant metadata.
Returning `variants` per block costs more, but it keeps matrix columns
interpretable for downstream tools.
