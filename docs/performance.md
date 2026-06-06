# Performance

`genoio` builds matrices in Rust and crosses the Python boundary once per read
or block. This avoids Python loops over variants and samples.

The current local benchmark was run on an Apple Silicon M1 Mac (`arm64`) with
Python 3.11. The Rust extension was built in release mode. The headline
comparisons read the first 1,000 variants into a dense `float32` matrix.

| Source | genoio median | comparison median | Result |
|---|---:|---:|---|
| VCF vs `cyvcf2` | 0.1064 s | 0.2383 s | genoio 2.24x faster |
| PLINK1 vs `pandas_plink` | 0.0163 s | 2.5449 s | genoio 156.13x faster |
| PLINK2 matrix-only vs `pgenlib` | 0.0094 s | 0.0069 s | `pgenlib` 1.36x faster |

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

The headline PLINK2 timing above is the matrix-only case. Additional PLINK2
scenarios measure the cost of returning metadata and applying filters:

| Scenario | What it measures | Median |
|---|---|---:|
| matrix-only | Read only the genotype matrix. | 0.0094 s |
| with variants | Return the matrix plus variant metadata. | 1.6101 s |
| sample-filtered | Read half the samples, preserving source sample order. | 1.5683 s |
| genotype-filtered | Apply genotype-stat filters before returning retained variants. | 0.5727 s |

The matrix-only comparison produced exact agreement with `pgenlib` for shape
`(3202, 1000)`, genotype sum `72577`, and zero missing values. At 10,000
variants, the matrix-only medians were 0.0879 s for `genoio` and 0.0500 s for
`pgenlib`.

The genotype-filtered scenario returned shape `(3202, 1000)`, genotype sum
`447988`, and zero missing values. It computes statistics for candidate
variants, but bounded block reads stop once the requested number of retained
variants has been returned.

The local BGEN fixture stores phased BGEN v1.2+ Layout 2 biallelic diploid
dosage records. `genoio` collapses those haplotype probabilities to expected
diploid A1 dosage. On the same machine and release build, BGEN scenario medians
were:

| Scenario | What it measures | Median |
|---|---|---:|
| matrix-only | Read only the dosage matrix. | 0.0532 s |
| with variants | Return the matrix plus variant metadata. | 0.0699 s |
| sample-filtered | Read half the samples, preserving source sample order. | 0.0471 s |
| genotype-filtered | Apply genotype-stat filters before returning retained variants. | 0.3102 s |
| indexed-region | Read a bounded region through a same-path `.bgen.bgi` index. | 0.0627 s |

At 10,000 variants, BGEN matrix-only median time was 0.6285 s.

For a direct BGEN matrix-only comparison, `scripts/benchmark_bgen.py` can also
compute expected dosages through `bgen_reader`/`cbgen`. The high-level
`bgen_reader.read(slice(...))` path is not used for this mixed-width local
fixture because it raises worker-thread broadcast errors; the benchmark instead
uses `open_bgen` metadata with `cbgen.read_probability()` per variant. The
comparison package was unavailable during the latest run, so the table below
records the previous direct comparison.

| Variants | genoio median | bgen_reader/cbgen median | Value check |
|---:|---:|---:|---|
| 1,000 | 0.1175 s | 0.1603 s | exact `allclose`, max diff 0 |
| 10,000 | 1.1968 s | 1.1133 s | exact `allclose`, max diff 0 |

---

## Run local benchmarks

Build the Rust extension in release mode before benchmarking:

```bash
make build-release
```

Then run the benchmark scripts:

```bash
python scripts/benchmark_vcf.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink2.py --scenario all --max-variants 1000 --repeats 5
python scripts/benchmark_plink2.py --scenario matrix-only --max-variants 10000 --repeats 5 --no-compare
```

For BGEN, use a BGEN v1.2+ Layout 2 biallelic diploid dosage fixture:

```bash
python scripts/benchmark_bgen.py --scenario all --max-variants 1000 --repeats 5
python scripts/benchmark_bgen.py --scenario matrix-only --backend both --max-variants 1000 --repeats 5
python scripts/benchmark_bgen.py --scenario indexed-region --region 22:20000000-21000000 --max-variants 1000 --repeats 5
```

Optional comparison packages are used when installed:

- `cyvcf2` for VCF
- `pandas_plink` for PLINK1
- `pgenlib` for PLINK2
- `bgen_reader` for BGEN

---

## What affects speed

Metadata filters are cheaper than genotype filters because they can run before
matrix decoding. Region filters on indexed compressed VCF/BCF sources and BGEN
sources with a same-path `.bgen.bgi` index can also avoid scanning unrelated
records.

PLINK2 performance depends on how much metadata the caller requests. Block
reads that only need matrix data can avoid parsing full variant metadata.
Returning `variants` per block costs more, but it keeps matrix columns
interpretable for downstream tools.
