# Phase 5 Review Remediation Evidence

## Issues Fixed

- Sparse reads now build CSC directly in VCF and PLINK1 readers instead of converting through a dense matrix.
- PLINK1 block reads now seek/read per-variant BED payloads instead of loading the full `.bed` file.
- VCF and PLINK1 block reads apply variant-window gating before genotype decode when filters do not require genotype statistics.
- Omitted sparse `missing` now defaults to `raise`, while explicit `missing="nan"` and `missing="impute"` remain rejected.
- Invalid sparse missing options now raise `InvalidOptionError` instead of leaking raw `TypeError`.

## Red Evidence

### Sparse Defaults And Structured Errors

Command:

```bash
pytest -q tests/test_sparse_read.py tests/test_blocks.py
```

Result before fix:

```text
FAILED tests/test_sparse_read.py::test_sparse_true_returns_csc_and_preserves_tuple_metadata
FAILED tests/test_sparse_read.py::test_sparse_invalid_missing_policy_raises_structured_error
FAILED tests/test_blocks.py::test_sparse_blocks_work_with_default_missing_policy
```

### Structural Architecture Gaps

Command:

```bash
rg "sparse_from_dense_minor_flipped|fs::read\(bed\)|let dense = read_.*dense" rust/genoio-io/src rust/genoio-core/src
```

Result before fix showed sparse I/O paths calling dense readers and PLINK1 using `fs::read(bed)`.

## Green Evidence

Commands:

```bash
env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml --workspace
env CC=clang AR=ar python -m maturin develop
pytest -q tests/test_sparse_read.py tests/test_blocks.py tests/test_dense_read.py tests/test_filters.py
pytest -q
ruff check src tests
python - <<'PY'
import genoio
print(genoio.__name__)
PY
```

Results after fix:

```text
Rust workspace tests passed.
Editable extension rebuilt.
43 passed in targeted Phase 5 Python suite.
73 passed in full Python suite.
ruff reported no diagnostics.
Smoke import printed genoio.
```
