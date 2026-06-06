# FAQ

## Can I use genoio output with JAX?

Yes. Dense reads return NumPy arrays, so convert each block with
`jax.numpy.asarray` or `jax.device_put` before calling JAX code.

```python
import jax.numpy as jnp

for X, variants in ds.iter_blocks(
    block_size=4096,
    variants=genoio.maf(max=0.05) & genoio.snp() & genoio.biallelic(),
    return_variants=True,
):
    # Convert inside the loop so the analysis does not need one whole-genome
    # genotype matrix in device memory.
    X_jax = jnp.asarray(X)
    scan_block(X_jax, variants)
```

Keep the host-to-device transfer cost in mind. `genoio` reads data on the CPU;
JAX kernels run fastest when each block is large enough to amortize that transfer.

## Can I use genoio output with PyTorch?

Yes. Dense reads return NumPy arrays, and `torch.from_numpy` converts compatible
CPU arrays without copying.

```python
import torch

for X, variants in ds.iter_blocks(
    block_size=4096,
    variants=genoio.maf(max=0.05) & genoio.snp() & genoio.biallelic(),
    return_variants=True,
):
    # from_numpy shares CPU memory with X. Clone if downstream code mutates it.
    X_tensor = torch.from_numpy(X).to(device)
    scan_block(X_tensor, variants)
```

If you move tensors to a GPU, PyTorch will copy the data to that device. For
large scans, convert and transfer one block at a time.

## What about sparse matrices?

Sparse reads return SciPy sparse matrices. That is useful for methods that
already accept SciPy inputs. JAX and PyTorch have their own sparse APIs, so
convert deliberately if your model expects those representations.

For many association scans, dense blocks are simpler. They also keep the same
shape across formats and make it easier to call NumPy, JAX, PyTorch, or compiled
association kernels.

## Why does genoio use filter objects instead of Python callbacks?

Filter objects describe the selection before records are read. Python builds an
expression such as
[`genoio.region("22:20000000-21000000")`](api/filters.md#genoio.region) and
[`genoio.maf(max=0.05)`](api/filters.md#genoio.maf); Rust evaluates that
expression while reading the file.

That matters for region reads. A concrete region filter can be pushed into
indexed VCF/BCF and BGEN reads, so `genoio` can jump to the requested genomic
interval instead of scanning the full file.

## When should I use read, iter_blocks, or iter_regions?

Use [`read(...)`](api/reading.md#genoio.Dataset.read) when the requested matrix
fits in memory and you want one array.

Use [`iter_blocks(...)`](api/reading.md#genoio.Dataset.iter_blocks) for
GWAS-style scans. It returns fixed-width chunks of retained variants, which is a
good fit for algorithms that apply the same test to many independent variant
columns.

Use [`iter_regions(...)`](api/reading.md#genoio.Dataset.iter_regions) for cis
scans or other local analyses. It takes a sequence of region filters and returns
one matrix per requested interval.

## Are samples reordered?

No. Retained samples stay in source order. Use
[`ds.samples()`](api/reading.md#genoio.Dataset.samples) or
[`read(return_samples=True)`](api/reading.md#genoio.Dataset.read) to align
phenotype and covariate tables before a scan.

## Are variants reordered?

No. Variants are returned in source order after filtering. For region iteration,
each result contains the variants retained inside that region.

## Which haplotype representations are supported?

VCF/BCF haplotype reads use phased hardcall `GT` records. PLINK2 haplotype
reads support explicit phased hardcalls with `dosage="hardcall"` and explicit
phased full dosages with `dosage="dosage"`. BGEN haplotype reads support Layout
2 phased biallelic diploid probabilities with `kind="haplo",
dosage="dosage"`, returned as expected A1 dosage per haplotype row.

`genoio` does not convert probabilities into hardcalls. Sparse PLINK2 dosage
haplotypes and sparse BGEN haplotypes are not supported, while PLINK2 explicit
phased hardcall haplotypes can be read sparsely when retained calls are
non-missing. If an unsupported record is retained, the read fails;
metadata-only filters can skip unsupported records before decode.

## Universal Standard

[Is this a new universal standard](https://xkcd.com/927/)? No. Shh.
