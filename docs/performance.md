# Performance

`genoio` builds matrices in Rust and crosses the Python boundary once per read
or block. This avoids Python loops over variants and samples.

On the local `data/chr22_hg38` fixture with 3,202 samples, release builds show:

| Read | genoio | comparison |
|---|---:|---:|
| VCF, first 1,000 variants | ~0.10 s | `cyvcf2` matrix construction ~0.23 s |
| PLINK2, first 100 variants | ~0.003 s | `pgenlib` matrix construction ~0.004 s |
| PLINK2, first 1,000 variants | ~0.02 s | `pgenlib` matrix construction ~0.007 s |
| PLINK2, first 10,000 variants | ~0.30 s | `pgenlib` matrix construction ~0.05 s |

These numbers are workload- and machine-dependent. Treat them as local
benchmarks, not universal claims.

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
python scripts/benchmark_plink2.py --max-variants 1000 --repeats 3
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
