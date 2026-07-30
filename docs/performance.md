# Performance

`genoio` reads genotype data in Rust and returns NumPy, SciPy, and Polars
objects to Python. The goal is simple: keep variant and sample loops out of
Python, and cross the Python boundary once per full read or block.

On local 1000 Genomes chromosome 22 benchmarks, that makes `genoio` fast for
VCF, PLINK1, and BGEN reads, and close to `pgenlib` for PLINK2 matrix-only
reads. Absolute timings matter more than relative ratios here: a 2x gap at a
few milliseconds is rarely the bottleneck in an analysis pipeline.

## Headline Results

These benchmarks were run on an Apple Silicon M1 Mac (`arm64`) with Python 3.11
and a release-mode Rust extension. The headline comparison reads dense
`float32` matrices.

| Source | Variants | genoio median | Comparison | Result |
|---|---:|---:|---:|---|
| VCF | 1,000 | 0.0445 s | `cyvcf2`: 0.2141 s | genoio 4.8x faster |
| PLINK1 | 1,000 | 0.0103 s | `pandas_plink`: 1.4575 s | genoio 141.5x faster |
| PLINK2 matrix-only | 1,000 | 0.0100 s | `pgenlib`: 0.0074 s | `pgenlib` 1.4x faster |
| BGEN dosage | 1,000 | 0.0235 s | `bgen_reader`/`cbgen`: 0.1076 s | genoio 4.6x faster |
| BGEN dosage | 1,000 | 0.0235 s | `bgen`: 0.0240 s | about parity |
| VCF | 10,000 | 0.5355 s | `cyvcf2`: 2.2634 s | genoio 4.2x faster |
| PLINK1 | 10,000 | 0.0796 s | `pandas_plink`: 1.5335 s | genoio 19.3x faster |
| PLINK2 matrix-only | 10,000 | 0.0886 s | `pgenlib`: 0.0490 s | `pgenlib` 1.8x faster |
| BGEN dosage | 10,000 | 0.3869 s | `bgen_reader`/`cbgen`: 1.0550 s | genoio 2.7x faster |
| BGEN dosage | 10,000 | 0.3869 s | `bgen`: 0.3028 s | `bgen` 1.3x faster |

Treat these as local benchmarks, not universal speed guarantees. File layout,
storage, filters, metadata requests, and Python environment all affect runtime.

## Single-block read scenarios

The sweep below measures common first-block patterns: matrix-only reads,
returning variant metadata, sample filtering, and genotype-stat filtering.
These scenarios call `next(dataset.iter_blocks(...))`; they don't measure
sustained iteration. VCF rows use 1,000 variants because the compressed VCF
scenarios are more expensive than the packed binary formats. PLINK rows use
10,000 variants.

| Source | Scenario | Variants | Median | Notes |
|---|---|---:|---:|---|
| VCF | metadata scan | all | 8.5002 s | Reads sample and variant metadata. |
| VCF | matrix-only | 1,000 | 0.0445 s | Dense genotype hardcalls. |
| VCF | with variants | 1,000 | 0.0452 s | Adds variant metadata. |
| VCF | sample-filtered | 1,000 | 0.0273 s | Reads half the samples. |
| VCF | genotype-filtered | 1,000 | 0.2210 s | Computes stats before returning retained variants. |
| VCF | indexed-region | 1,000 | 0.0592 s | Uses the `.tbi` index. |
| VCF | indexed-region sample-filtered | 1,000 | 0.0440 s | Combines region and sample filters. |
| VCF | haplotype matrix-only | 1,000 | 0.0619 s | Returns phased hardcall haplotype rows. |
| VCF | sparse matrix-only | 1,000 | 0.0458 s | Returns sparse CSC hardcalls. |
| PLINK1 | matrix-only | 10,000 | 0.0804 s | Fast packed BED hard-call path. |
| PLINK1 | with variants | 10,000 | 0.0971 s | Adds variant metadata. |
| PLINK1 | sample-filtered | 10,000 | 0.0479 s | Reads half the samples. |
| PLINK1 | genotype-filtered | 10,000 | 0.1569 s | Computes stats before returning retained variants. |
| PLINK2 | matrix-only | 10,000 | 0.0872 s | Fast path when metadata is not needed. |
| PLINK2 | with variants | 10,000 | 0.1184 s | Adds variant metadata. |
| PLINK2 | sample-filtered | 10,000 | 0.0569 s | Reads half the samples. |
| PLINK2 | genotype-filtered | 10,000 | 0.4925 s | Current largest PLINK2 cost surface. |

The main practical lesson: metadata and sample filters are cheap enough for
routine use. Genotype-stat filters do more work because they must inspect
genotypes before deciding which variants to keep.

## Sustained `iter_blocks` throughput

Use `scripts/benchmark_iter_blocks.py` to measure more than the first yielded
block. It consumes the same exact retained-variant prefix for every requested
block size, so the sweep changes boundary frequency without changing decoded
work. The timed path only checks shapes and counts. It doesn't concatenate the
matrices.

```bash
python scripts/benchmark_iter_blocks.py \
  --source-format vcf \
  --path data/chr22_hg38.vcf.gz \
  --label candidate-05bc6c1 \
  --block-sizes 128,512,2048 \
  --max-variants 16384 \
  --scenario all \
  --warmups 1 \
  --repeats 5 \
  --output-json /tmp/genoio-iter-blocks.json
```

The total must be divisible by every block size. The script reports actual
blocks, variants, rows, per-run times, median time, and variants per second.
It reuses a single `Dataset` and closes every iterator explicitly, including
bounded scans that stop before end of file.

For an A/B comparison, use distinct labels such as `baseline-69d02b3` and
`candidate-05bc6c1` while keeping the source, command, release build, and
machine fixed. Each JSON report records the label, resolved source path,
`genoio` version, Python version, platform, and machine so results can be
audited after the run.

## BGEN Dosage Reads

The local BGEN fixture stores BGEN v1.2+ Layout 2 biallelic diploid dosage
records. `genoio` returns expected diploid A1 dosage values for this fixture.

| Scenario | Variants | Median | Notes |
|---|---:|---:|---|
| matrix-only | 1,000 | 0.0235 s | Reads only the dosage matrix. |
| with variants | 1,000 | 0.0247 s | Adds variant metadata. |
| sample-filtered | 1,000 | 0.0182 s | Reads half the samples. |
| genotype-filtered | 1,000 | 0.1511 s | Computes dosage-based stats before returning variants. |
| indexed-region | 1,000 | 0.0702 s | Uses a same-path `.bgen.bgi` index. |
| matrix-only | 10,000 | 0.3869 s | Larger block read. |

Direct matrix-only comparisons against optional BGEN readers produced:

| Variants | genoio median | `bgen_reader`/`cbgen` median | `bgen` median |
|---:|---:|---:|---:|
| 1,000 | 0.0235 s | 0.1076 s | 0.0240 s |
| 10,000 | 0.3869 s | 1.0550 s | 0.3028 s |

## Benchmark Data

The benchmark scripts default to `data/chr22_hg38`. That directory is a local
1000 Genomes-derived fixture and is not distributed with the repository.

The fixture starts from the PLINK 2
[1000 Genomes phase 3 hg38 resources](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg).
Chromosome 22 PLINK2 files are used as the source, then converted with `plink2`
to VCF and PLINK1. This keeps format comparisons focused on reader behavior
rather than differences in samples or variants.

The converted VCF fixture does not contain `FORMAT/DS`, so VCF dosage reads are
not included. The default PLINK2 and BGEN fixtures also do not contain the
phased records needed for haplotype dosage benchmarks.

## Run Local Benchmarks

Build the Rust extension in release mode first:

```bash
make build-release
```

Then run the relevant benchmark:

```bash
python scripts/benchmark_iter_blocks.py --source-format vcf --path data/chr22_hg38.vcf.gz --label candidate-05bc6c1 --block-sizes 128,512,2048 --max-variants 16384 --scenario all --repeats 5
python scripts/benchmark_vcf.py --scenario all --max-variants 1000 --repeats 5
python scripts/benchmark_vcf.py --scenario matrix-only --max-variants 10000 --repeats 5
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 5
python scripts/benchmark_plink1.py --max-variants 10000 --repeats 5
python scripts/benchmark_plink2.py --scenario all --max-variants 1000 --repeats 5
python scripts/benchmark_plink2.py --scenario matrix-only --max-variants 10000 --repeats 5 --no-compare
python scripts/benchmark_plink2.py --scenario all --max-variants 10000 --repeats 5 --backend genoio --no-compare
```

For BGEN:

```bash
python scripts/benchmark_bgen.py --scenario all --backend all --max-variants 1000 --repeats 5
python scripts/benchmark_bgen.py --scenario matrix-only --backend all --max-variants 10000 --repeats 5
python scripts/benchmark_bgen.py --scenario indexed-region --region 22:20000000-21000000 --max-variants 1000 --repeats 5
```

For the filter benchmark:

```bash
python scripts/benchmark_filter_perf.py --source-format bfile --path data/chr22_hg38 --max-variants 10000 --repeats 5 --window-mode retained --filter-shape maf --maf-min 0.01
```

Optional comparison packages are used when installed:

- `cyvcf2` for VCF
- `pandas_plink` for PLINK1
- `pgenlib` for PLINK2
- `bgen_reader`, `cbgen`, and `bgen` for BGEN

## What Affects Speed

Matrix-only reads are fastest because `genoio` can skip metadata work that the
caller did not request. Returning `variants` costs more, but it keeps matrix
columns interpretable.

Metadata filters are cheaper than genotype filters because they can run before
matrix decoding. Region filters on indexed compressed text VCF sources and BGEN
sources with a same-path `.bgen.bgi` index can also skip unrelated records.

For genotype-stat filters, pushing the filter into the Rust reader avoids a
Python-side full-window read followed by NumPy post-filtering. This retained
window benchmark reads up to 10,000 variants passing `maf(min=0.01)`.

| Source | Rust-side filter | NumPy post-filter | Result |
|---|---:|---:|---|
| VCF | 2.4315 s | 5.8654 s | Rust 2.4x faster |
| PLINK1 | 0.1675 s | 1.3927 s | Rust 8.3x faster |
| PLINK2 | 0.4740 s | 1.8106 s | Rust 3.8x faster |
| BGEN | 1.3510 s | 3.5479 s | Rust 2.6x faster |
