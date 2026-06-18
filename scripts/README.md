# Benchmark Scripts

These scripts compare bounded dense matrix construction across `genoio` and common format-specific readers.

Defaults use `data/chr22_hg38` and read the first 1000 variants. Use `--max-variants` to adjust the workload.

`data/chr22_hg38` is a local benchmark fixture, not a repository fixture. It
comes from the PLINK 2
[1000 Genomes phase 3 hg38 resources](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg).
The chromosome 22 PLINK 2 files were used as the source, then converted with
`plink2` to VCF and PLINK1 `.bed/.bim/.fam` files so each script reads the same
underlying genotypes.

```bash
python scripts/benchmark_vcf.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink1.py --max-variants 1000 --repeats 3
python scripts/benchmark_plink2.py --scenario all --max-variants 1000 --repeats 5
python scripts/benchmark_plink2.py --scenario matrix-only --max-variants 10000 --repeats 5 --no-compare
python scripts/benchmark_plink2.py --kind haplo-hardcall --scenario matrix-only --backend genoio --max-variants 1000
python scripts/benchmark_plink2.py --kind haplo-dosage --scenario matrix-only --backend genoio --max-variants 1000
```

Backend-specific runs are useful when an optional comparison package is unavailable:

```bash
python scripts/benchmark_vcf.py --backend cyvcf2 --max-variants 1000
python scripts/benchmark_plink1.py --backend pandas_plink --max-variants 1000
python scripts/benchmark_plink2.py --backend pgenlib --max-variants 1000
python scripts/benchmark_bgen.py --backend bgen_reader --max-variants 1000
```

Comparison backends are optional. Install `cyvcf2`, `pandas-plink`, `pgenlib`,
or `bgen-reader` only when you need those specific comparisons.

The PLINK2 `genoio` benchmark accepts either `.pvar` or `.pvar.zst`; compressed
PVAR metadata is decompressed by the Rust reader. The default `--kind geno`
continues to time genotype hardcall reads. Use `--kind haplo-hardcall` for
explicit phased hardcall PGEN records and `--kind haplo-dosage` for explicit
phased full-dosage PGEN records. The `pgenlib` comparison is skipped for
haplotype modes because the benchmark comparison path only reads genotype
hardcalls.

The BGEN benchmark reads dosage values from `<prefix>.bgen` and uses
`<prefix>.sample` for the sample-filtered scenario when the BGEN file does not
embed sample identifiers. The default `--kind geno` reads expected diploid A1
dosage. Use `--kind haplo` for source-encoded phased Layout 2 biallelic diploid
probabilities returned as expected A1 dosage per haplotype row. The
`indexed-region` scenario uses `--region` and a same-path `<prefix>.bgen.bgi`
index when present.

```bash
python scripts/benchmark_bgen.py --scenario all --max-variants 1000 --repeats 5
python scripts/benchmark_bgen.py --scenario matrix-only --backend both --max-variants 1000 --repeats 5
python scripts/benchmark_bgen.py --scenario matrix-only --backend all --max-variants 1000 --repeats 5
python scripts/benchmark_bgen.py --scenario indexed-region --region 22:20000000-21000000 --max-variants 1000 --repeats 5
python scripts/benchmark_bgen.py --kind haplo --scenario matrix-only --backend genoio --max-variants 1000
```

The BGEN `--backend both` comparison computes expected dosage through
`bgen_reader`/`cbgen` and checks matrix parity with `genoio`. The comparison
path reads probabilities variant-by-variant because the high-level
`bgen_reader.read(slice(...))` path can fail on mixed-width BGEN probability
records. Use `--backend bgen` or `--backend all` to compare against the
optional Cython/C++ `bgen` package, which reads per-variant `alt_dosage`.
Haplotype scenarios skip comparison backends because this script does not
reshape phased probabilities into haplotype rows for those backends.

`pgenlib` must be importable for the PLINK2 comparison backend. If it is built in the symlinked PLINK repository but not installed, pass:

```bash
python scripts/benchmark_plink2.py --pgenlib-path plink-ng/2.0/Python
```
