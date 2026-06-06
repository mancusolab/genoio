# Haplotype I/O Human Test Plan

## Preconditions

- Worktree at `b468e09730efb33d29eeca3e3a95bb598fc238d9`.
- Local Rust/Python development environment available.
- Fixture-writing tests can create temporary VCF, BGEN, and PLINK2 files.

## Phase Checks

| Step | Action | Expected Result |
| --- | --- | --- |
| 1 | Run `cargo test -p genoio-io --test bgen_dense` | BGEN haplotype dosage, missingness, region pushdown, unsupported record, and filter tests pass. |
| 2 | Run `cargo test -p genoio-io --test plink2_dense` | PLINK2 hardcall/dosage haplotype decode, ordering, unsupported retained records, and regressions pass. |
| 3 | Run `cargo test -p genoio-io --test vcf_haplotype` | Existing VCF hardcall haplotype behavior remains supported. |
| 4 | Run `cargo test -p genoio-io --test filter_genotype_stats` | Genotype-stat filtering regressions pass. |
| 5 | Run `pytest tests/test_public_api.py tests/test_haplotype.py tests/test_dense_read.py tests/test_blocks.py tests/test_filters.py -q` | Python API, dense reads, blocks, regions, sparse rejection, and metadata tests pass. |
| 6 | Run `pytest tests/test_benchmark_bgen_cli.py tests/test_benchmark_plink2_cli.py -q` | Benchmark CLI accepts and dispatches dense haplotype modes. |
| 7 | Run `make verify` | Full project verification passes. |

## End-to-End Scenarios

| Scenario | Steps | Expected Result |
| --- | --- | --- |
| VCF hardcall haplotypes | Open phased VCF, call `read(kind="haplo")`, request samples/variants. | Matrix has one row per haplotype; metadata includes `source_sample_index` and `haplotype_index`. |
| BGEN dosage haplotypes | Open phased Layout 2 BGEN, call `read(kind="haplo", dosage="dosage")`. | Returns expected A1 dosage per haplotype row; missing samples become `nan` by default. |
| PLINK2 hardcall haplotypes | Open explicit phased hardcall PGEN, call `read(kind="haplo", dosage="hardcall")`. | Returns `0`/`1` haplotype rows in source sample then haplotype order. |
| PLINK2 dosage haplotypes | Open explicit phased full-dosage PGEN, call `read(kind="haplo", dosage="dosage")`. | Returns per-haplotype A1 dosages and preserves variant column order. |
| Unsupported sparse modes | Request sparse PLINK2/BGEN haplotypes. | Raises `UnsupportedRepresentation` with dense-mode guidance. |
| Unsupported retained records | Retain unphased/multiallelic/invalid-ploidy records during haplotype read. | Read fails before returning invalid haplotype values. |
| Filtered unsupported records | Apply metadata filter excluding unsupported records. | Unsupported payload is skipped before decode and supported retained records read successfully. |
| Region/block reads | Use `iter_blocks(...)` and `iter_regions(...)` for BGEN and PLINK2 haplotypes. | Chunks concatenate to full reads; region results preserve requested region association and source variant order. |
| Genotype-stat filters | Apply `maf`, `mac`, and `missing_rate` to haplotype reads. | Filters use collapsed diploid expected A1 dosage; returned matrix remains haplotype-level. |

## Traceability

| AC ID | Automated Evidence | Human Step |
| --- | --- | --- |
| AC1.1-AC1.4 | `test_public_api.py`, `test_haplotype.py`, `test_dense_read.py`, `vcf_haplotype`, PLINK2 Rust tests | Phase checks 2, 3, 5 |
| AC2.1-AC2.5 | `bgen_dense`, `test_dense_read.py`, `test_haplotype.py` | Phase checks 1, 5 |
| AC3.1-AC3.5 | `plink2_dense`, `filter_genotype_stats`, `test_haplotype.py`, `test_dense_read.py` | Phase checks 2, 4, 5 |
| AC4.1-AC4.5 | Rust BGEN/PLINK2 tests plus Python haplotype/block/filter tests | Phase checks 1, 2, 4, 5 |
| AC5.1-AC5.4 | Targeted suites, docs grep surface, benchmark CLI tests, `make verify` | Phase checks 6, 7 |
