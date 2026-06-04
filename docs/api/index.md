# API overview

The public API is intentionally small. Start here if you want to understand how
the pieces fit together before reading the generated reference pages: source
constructors create `Dataset` objects, datasets read dense or sparse genotype
matrices, and filters describe which variants should be retained. Matrices are
returned with samples on rows and variants on columns, with optional Polars
metadata frames for sample and variant alignment.

- [`vcf`](reading.md#genoio.vcf), [`bfile`](reading.md#genoio.bfile), and
  [`pfile`](reading.md#genoio.pfile) return reusable [`Dataset`](reading.md#genoio.Dataset)
  objects.
- [`Dataset.read`](reading.md#genoio.Dataset.read) reads one matrix.
- [`Dataset.blocks`](reading.md#genoio.Dataset.blocks) streams matrix blocks.
- [`Dataset.samples`](reading.md#genoio.Dataset.samples) and
  [`Dataset.variants`](reading.md#genoio.Dataset.variants) return metadata
  frames.
- Filter constructors in [Filters](filters.md) build serializable variant
  predicates.

Public errors are listed in [Errors](errors.md).
