# genoio

`genoio` reads VCF/BCF, PLINK1, PLINK2, and BGEN genotype sources into Python
matrices. The public API is Python; the file readers are Rust. This pairing
keeps the interface simple while moving parsing and matrix construction into
compiled, efficient code.

It is built for downstream genetics tools that need stable matrix contracts:
samples on rows, variants on columns, and metadata aligned to returned matrices,
blocks, or regions.

## Install

```bash
pip install git+https://github.com/mancusolab/genoio.git
```

From a local checkout:

```bash
pip install .
```

## Documentation

For complete documentation and examples, see the
[documentation site](https://mancusolab.github.io/genoio).

## Quick Example

```python
import genoio

ds = genoio.pfile("data/chr22_hg38")
samples = ds.samples()
y = load_phenotype_vector(samples["iid"])
C = load_covariates(samples["iid"])

for X, variants in ds.iter_blocks(10_000, return_variants=True):
    # X has shape (samples, variants_in_this_block).
    # `y` and `C` must be aligned to the rows described by `samples`.
    association_scan(X, y, C, variants=variants)
```

Use `read(...)` for one matrix, `iter_blocks(...)` for streaming scans, and
`iter_regions(...)` for interval-based workflows.

## License

MIT. See [LICENSE](LICENSE).

---

## Notes

`genoio` was developed by members of the Mancuso Lab with assistance from Codex,
following the practices described in the
[scientific-software-playbook](https://github.com/mancusolab/scientific-software-playbook)
and [coding-skills](https://github.com/mancusolab/coding-skills) repositories.
