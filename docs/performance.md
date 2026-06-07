# Performance

`genoio` builds matrices in Rust and crosses the Python boundary once per read
or block. This avoids Python loops over variants and samples.

The current local benchmark was run on an Apple Silicon M1 Mac (`arm64`) with
Python 3.11. The Rust extension was built in release mode. The headline
comparisons read the first 1,000 variants into a dense `float32` matrix. PLINK1
and PLINK2 comparison benchmarks used `numpy 1.26.4` because `pandas-plink`
currently declares `numpy<2.0`.

| Source | genoio median | comparison median | Result |
|---|---:|---:|---|
| VCF vs `cyvcf2` | 0.1064 s | 0.2383 s | genoio 2.24x faster |
| PLINK1 vs `pandas_plink` | 0.0110 s | 2.5486 s | genoio 231.69x faster |
| PLINK2 matrix-only vs `pgenlib` | 0.0109 s | 0.0058 s | `pgenlib` 1.88x faster |

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

The latest 10,000-variant PLINK scenario sweep used the same release build with
`numpy 2.4.6`. It measures block reads with and without returned metadata,
sample filtering, and genotype-stat filtering:

| Source | Scenario | What it measures | Median |
|---|---|---|---:|
| PLINK1 | matrix-only | Read only the genotype matrix. | 0.1030 s |
| PLINK1 | with variants | Return the matrix plus variant metadata. | 0.1160 s |
| PLINK1 | sample-filtered | Read half the samples, preserving source sample order. | 0.0487 s |
| PLINK1 | genotype-filtered | Apply genotype-stat filters before returning retained variants. | 0.1701 s |
| PLINK2 | matrix-only | Read only the genotype matrix. | 0.1060 s |
| PLINK2 | with variants | Return the matrix plus variant metadata. | 0.1543 s |
| PLINK2 | sample-filtered | Read half the samples, preserving source sample order. | 0.0820 s |
| PLINK2 | genotype-filtered | Apply genotype-stat filters before returning retained variants. | 0.8845 s |

The matrix-only comparison produced exact agreement with `pgenlib` for shape
`(3202, 1000)`, genotype sum `72577`, and zero missing values. The PLINK1
comparison produced the same shape and sum, with exact agreement against
`pandas_plink`.

The PLINK1 matrix-only path and genotype-stat filter path work directly from
packed BED hard calls when variant metadata is not requested. They still require
the companion `.bim` file to exist, but they do not parse BIM rows on those
matrix-only block reads. Genotype-filtered block reads compute statistics for
candidate variants, then stop once the requested number of retained variants has
been returned.

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
