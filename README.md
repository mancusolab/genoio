# genoio

`genoio` reads VCF/BCF, PLINK1, PLINK2, and BGEN genotype sources into Python
matrices. The public API is Python; the file readers are Rust. This pairing
keeps the interface simple while moving parsing and matrix construction into
compiled, efficient code.

It is built for downstream genetics tools that need stable matrix contracts:
samples on rows, variants on columns, and metadata aligned to returned matrices,
blocks, or regions.

Documentation: [mancusolab.github.io/genoio](https://mancusolab.github.io/genoio)

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

This sketch/mockup shows how to perform a blockwise scan for GWAS using `genoio`.

```python
import genoio

# load PLINK2 genotype data
ds = genoio.pfile("data/chr22_hg38")
y = load_phenotypes()

# set up filters
common = genoio.maf(min=0.01) & genoio.missing_rate(max=0.1)

# iterate blockwise
for X, variants in ds.iter_blocks(10_000, variants=common, return_variants=True):
    association_scan(X, y, variants=variants)
```

Phenotypes and covariates should be aligned to `ds.samples()`.

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
