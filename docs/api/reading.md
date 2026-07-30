# Reading

For worked examples, see [GWAS](../examples/gwas.md),
[cis-eQTL](../examples/cis-eqtl.md), and [Filtering](../filtering.md).

::: genoio.Dataset

::: genoio.BlockIterator

## `iter_blocks()` return-type compatibility

`Dataset.iter_blocks()` now returns `BlockIterator` instead of a Python
generator. Normal iteration, `next()`, `list()`, `close()`, and
`contextlib.closing()` remain supported. Generator-only methods such as
`send()` and `throw()`, along with generator introspection attributes, were
intentionally removed as a minor-release compatibility change.

Use ordinary loop control for iteration and `with dataset.iter_blocks(...)` for
deterministic cleanup. There is no direct `throw()` equivalent.

::: genoio.vcf

::: genoio.bfile

::: genoio.pfile

::: genoio.bgen
