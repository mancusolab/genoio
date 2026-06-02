# Phase 3 Code Review Remediation Evidence

## TDD Evidence Finding

Reviewer finding: Phase 3 behavior changes were committed in a single implementation
commit, so the review could not audit failing-test-first evidence from repository
history alone.

Remediation:
- Preserve the implementor-reported transcript evidence for the original Phase 3 work.
- Add fresh failing-test-first evidence for the review-fix regressions in commit
  `00ab16b`.

Original Phase 3 transcript evidence reported by implementor:
- RED Rust command: `env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml --workspace vcf_dense`
- RED Rust result: exit 101. Key failures: missing `read_vcf_dense`,
  `read_plink1_dense`, `DenseGenotypeMatrix`, and `DenseDiagnostics`.
- RED Python command: `pytest -q tests/test_dense_read.py`
- RED Python result: exit 1. Key failures: `Dataset.read()` raised
  `NotImplementedError`; `genoio.SampleFilterError` did not exist.
- GREEN targeted result: VCF dense tests passed, PLINK1 dense tests passed, and
  `tests/test_dense_read.py` passed.

Fresh review-fix red/green evidence:
- RED Rust command: `env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml --workspace vcf_dense`
- RED Rust result: exit 101. The new multi-ALT regression showed dense VCF accepted
  a record with `ALT=G,T`, producing `DenseGenotypeMatrix { values: [1.0],
  ... alt_allele: Some("G,T") ... }`.
- GREEN Rust command: `env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml --workspace`
- GREEN Rust result: exit 0. Rust workspace tests passed, including
  `vcf_dense_rejects_multi_alt_records_even_when_gt_uses_first_alt`.

- RED Python command: `pytest -q tests/test_dense_read.py`
- RED Python result: exit 1. Key failures: multi-ALT dense VCF did not raise a
  structured source error, and validation reported duplicate samples before
  dtype/missing incompatibility.
- GREEN Python command: `pytest -q tests/test_dense_read.py`
- GREEN Python result: exit 0. Nine dense-read tests passed.

- RED all-missing impute check: `pytest -q tests/test_dense_read.py` after
  temporarily disabling the all-missing guard.
- RED all-missing impute result: exit 1. The new regression failed, proving it
  catches missing all-missing imputation behavior.
- GREEN all-missing impute result: the restored implementation raised
  `MissingDataError` and the focused dense-read test suite passed.

Final verification after Phase 3 review fixes:
- `env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml --workspace`
  -> exit 0.
- `env CC=clang AR=ar python -m maturin develop` -> exit 0.
- `pytest -q tests/test_dense_read.py tests/test_metadata.py` -> exit 0,
  16 passed.
- `pytest -q tests/test_dense_read.py` -> exit 0, 9 passed.
- `pytest -q` -> exit 0, 39 passed.
- `ruff check src tests` -> exit 0.
