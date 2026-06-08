# Performance

`genoio` reads genotype data in Rust and returns NumPy, SciPy, and Polars
objects to Python. The goal is simple: keep variant and sample loops out of
Python, and cross the Python boundary once per full read or block.

On local 1000 Genomes chromosome 22 benchmarks, that makes `genoio` especially
fast for VCF and PLINK1 reads, competitive for BGEN dosage reads, and close to
`pgenlib` for PLINK2 matrix-only reads. Absolute timings matter more than
relative ratios here: a 2x gap at a few milliseconds is rarely the bottleneck in
an analysis pipeline.

## Headline Results

These benchmarks were run on an Apple Silicon M1 Mac (`arm64`) with Python 3.11
and a release-mode Rust extension. The headline comparison reads the first
1,000 variants into a dense `float32` matrix.

| Source | genoio median | Comparison | Result |
|---|---:|---:|---|
| VCF | 0.1064 s | `cyvcf2`: 0.2383 s | genoio 2.24x faster |
| PLINK1 | 0.0110 s | `pandas_plink`: 2.5486 s | genoio 231.69x faster |
| PLINK2 matrix-only | 0.0109 s | `pgenlib`: 0.0058 s | `pgenlib` 1.88x faster |

Treat these as local benchmarks, not universal speed guarantees. File layout,
storage, filters, metadata requests, and Python environment all affect runtime.

## Larger Block Reads

The 10,000-variant sweep measures common block-read patterns: matrix-only reads,
returning variant metadata, sample filtering, and genotype-stat filtering.

| Source | Scenario | Median | Notes |
|---|---|---:|---|
| PLINK1 | matrix-only | 0.1030 s | Fast packed BED hard-call path. |
| PLINK1 | with variants | 0.1160 s | Adds variant metadata. |
| PLINK1 | sample-filtered | 0.0487 s | Reads half the samples. |
| PLINK1 | genotype-filtered | 0.1701 s | Computes stats before returning retained variants. |
| PLINK2 | matrix-only | 0.1060 s | Fast path when metadata is not needed. |
| PLINK2 | with variants | 0.1543 s | Metadata parsing is visible but modest. |
| PLINK2 | sample-filtered | 0.0820 s | Reads half the samples. |
| PLINK2 | genotype-filtered | 0.8845 s | Current largest PLINK2 cost surface. |

The main practical lesson: metadata and sample filters are cheap enough for
routine use. Genotype-stat filters do more work because they must inspect
genotypes before deciding which variants to keep.

## BGEN Dosage Reads

The local BGEN fixture stores phased BGEN v1.2+ Layout 2 biallelic diploid
dosage records. `genoio` collapses haplotype probabilities to expected diploid
A1 dosage.

| Scenario | Median | Notes |
|---|---:|---|
| matrix-only | 0.0532 s | Reads only the dosage matrix. |
| with variants | 0.0699 s | Adds variant metadata. |
| sample-filtered | 0.0471 s | Reads half the samples. |
| genotype-filtered | 0.3102 s | Computes dosage-based stats before returning variants. |
| indexed-region | 0.0627 s | Uses a same-path `.bgen.bgi` index. |

At 10,000 variants, BGEN matrix-only median time was 0.6285 s.

A previous direct matrix-only comparison against `bgen_reader`/`cbgen` produced:

| Variants | genoio median | `bgen_reader`/`cbgen` median |
|---:|---:|---:|
| 1,000 | 0.1175 s | 0.1603 s |
| 10,000 | 1.1968 s | 1.1133 s |

## Benchmark Data

The benchmark scripts default to `data/chr22_hg38`. That directory is a local
1000 Genomes-derived fixture and is not distributed with the repository.

The fixture starts from the PLINK 2
[1000 Genomes phase 3 hg38 resources](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg).
Chromosome 22 PLINK2 files are used as the source, then converted with `plink2`
to VCF and PLINK1. This keeps format comparisons focused on reader behavior
rather than differences in samples or variants.

## Run Local Benchmarks

Build the Rust extension in release mode first:

```bash
make build-release
```

Then run the relevant benchmark:

```bash
python scripts/benchmark_vcf.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink2.py --scenario all --max-variants 1000 --repeats 5
python scripts/benchmark_plink2.py --scenario matrix-only --max-variants 10000 --repeats 5 --no-compare
```

For BGEN:

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

## What Affects Speed

Matrix-only reads are fastest because `genoio` can skip metadata work that the
caller did not request. Returning `variants` costs more, but it keeps matrix
columns interpretable.

Metadata filters are cheaper than genotype filters because they can run before
matrix decoding. Region filters on indexed compressed VCF/BCF sources and BGEN
sources with a same-path `.bgen.bgi` index can also skip unrelated records.
