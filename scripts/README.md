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
```

Backend-specific runs are useful when an optional comparison package is unavailable:

```bash
python scripts/benchmark_vcf.py --backend cyvcf2 --max-variants 1000
python scripts/benchmark_plink1.py --backend pandas_plink --max-variants 1000
python scripts/benchmark_plink2.py --backend pgenlib --max-variants 1000
```

The PLINK2 `genoio` benchmark needs an uncompressed `.pvar`. If only `.pvar.zst` is present, the script decompresses it into a temporary directory with `zstd`.

`pgenlib` must be importable for the PLINK2 comparison backend. If it is built in the symlinked PLINK repository but not installed, pass:

```bash
python scripts/benchmark_plink2.py --pgenlib-path plink-ng/2.0/Python
```
