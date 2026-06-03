# API overview

The public API is intentionally small.

- [`open`](reading.md#genoio.open) returns a reusable [`Dataset`](reading.md#genoio.Dataset).
- [`read`](reading.md#genoio.read) reads one matrix.
- [`blocks`](reading.md#genoio.blocks) streams matrix blocks.
- [`samples`](reading.md#genoio.samples) and [`variants`](reading.md#genoio.variants)
  return metadata frames.
- Filter constructors in [Filters](filters.md) build serializable variant
  predicates.

Public errors are listed in [Errors](errors.md).
