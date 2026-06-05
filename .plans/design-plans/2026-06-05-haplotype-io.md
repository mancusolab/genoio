# Source-faithful Haplotype I/O Design

## Status
Approved

## Handoff Decision
- Current decision: ready
- Ready for implementation: yes
- Blocking items: none

## Metadata
- Date: 2026-06-05
- Slug: haplotype-io
- Artifact Directory: `.plans/design-plans/artifacts/2026-06-05-haplotype-io`

## Summary
Extend `genoio` with dense, source-faithful haplotype reads for PLINK2 and
BGEN. The design preserves the existing public `Dataset.read(...)` contract:
`kind="haplo"` requests haplotype rows, `dosage` selects the encoded source
representation, and unsupported retained records raise instead of being
converted or silently skipped.

The implementation follows the current Python/Rust split. Python validates
options and assembles NumPy/metadata outputs. PyO3 dispatches to format-specific
Rust readers. Rust performs source-native decoding, filtering, retained-variant
windowing, sample selection, and metadata row construction.

This plan explicitly avoids hardcall-from-dosage conversion and sparse
haplotype reads. Correctness, error behavior, docs, and parity tests come before
benchmarks.

## Problem Statement
`genoio` currently supports VCF haplotype reads, PLINK2 genotype hardcalls and
unphased dosages, and BGEN genotype dosages. Researchers building new genetics
methods also need haplotype-level matrices from PLINK2 and BGEN when those
haplotype representations are explicitly encoded in the source data. The I/O
contract should stay source-faithful: `genoio` should read encoded hardcalls,
dosages, and phased haplotypes, but it should not invent hardcalls from dosage
probabilities or apply analysis-policy thresholds during input.

## Definition of Done
Design a source-faithful haplotype I/O expansion for `genoio`. The plan should
cover dense haplotype reads for PLINK2 and BGEN without deriving hardcalls from
dosages: PLINK2 supports explicit hardcall phase via `kind="haplo",
dosage="hardcall"` and explicit phased dosage via `kind="haplo",
dosage="dosage"`; BGEN supports phased probability records as expected A1
dosage per haplotype row via `kind="haplo", dosage="dosage"`. Unsupported or
unphased retained records fail the read with `UnsupportedRepresentation` or
actionable source/representation errors rather than being skipped or converted.
Sparse haplotype reads, hardcall-from-dosage conversion, and researcher
workflow helpers are out of scope; benchmarks are included after correctness
and API/docs/tests are designed.

## Goals and Non-Goals
### Goals
- Add a precise design for dense PLINK2 haplotype reads from explicit phase
  information.
- Add a precise design for dense PLINK2 phased dosage haplotype reads when the
  PGEN record explicitly contains phased dosage data.
- Add a precise design for dense BGEN phased probability reads as expected A1
  dosage per haplotype row.
- Preserve the existing matrix orientation and metadata conventions: rows are
  samples or haplotypes, columns are variants, and haplotype sample metadata
  maps rows back to source samples.
- Specify tests, docs, errors, and benchmarks needed to finalize the I/O
  contract.

### Non-Goals
- Sparse haplotype matrices for PLINK2 or BGEN.
- Hardcall genotype or haplotype calls derived from dosage probabilities.
- Downstream researcher workflow helpers such as sample alignment, LD,
  residualization, association scans, or fine-mapping utilities.
- Broad BGEN hardcall support, multiallelic BGEN support, or variable ploidy
  support unless required to reject unsupported retained records clearly.

## Existing Patterns
Confirmed codebase patterns:

- Public API lives in `src/genoio/_api.py`; constructors resolve source files
  and `Dataset.read(...)` validates options before calling the Rust extension.
- Python matrix assembly lives in `src/genoio/_assembly.py`; dense reads return
  NumPy arrays and optional Polars metadata frames.
- PyO3 dispatch lives in `rust/genoio-py/src/lib.rs`; it maps `format`,
  `kind`-specific entrypoints, and `dosage` options to `genoio_io` readers.
- Format-specific binary readers live in `rust/genoio-io/src/`. PLINK2 and
  BGEN parsing should remain in `plink2.rs` and `bgen.rs` unless the
  implementation grows enough to justify a stable submodule boundary.
- Shared matrix, metadata, sparse, filter, and capability contracts live in
  `rust/genoio-core/src/`.
- Existing VCF haplotype behavior is the semantic template: retained unphased
  records fail; filtered-out unphased records can be dropped before haplotype
  decoding.
- Existing block iteration uses retained-variant windows. New haplotype readers
  must preserve `read(...)`, `iter_blocks(...)`, and `iter_regions(...)`
  behavior where the format supports the requested representation.

## Model Acquisition Path
- Path: `existing-codebase-port`
- Why this path: This is an extension of the existing `genoio` Python/Rust I/O
  architecture, not a new statistical model. The design should preserve current
  source resolution, filter, matrix assembly, metadata, and error patterns.
- User selection confirmation: The user requested the design in the current
  repository and current branch, with scope centered on finalizing I/O.

## Required Workflow States
- model_path_decided: yes
- codebase_investigation_complete_if_port: yes
- simulation_contract_complete_if_in_scope: n/a

## Model Specification Sources
| Source ID | Path/Link | Type | Notes | Confidence (high/med/low) |
| --- | --- | --- | --- | --- |
| SRC-1 | `src/genoio/_api.py` | existing-code | Public read validation and Dataset methods. | high |
| SRC-2 | `rust/genoio-py/src/lib.rs` | existing-code | PyO3 dispatch for dense, sparse, and VCF haplotype reads. | high |
| SRC-3 | `rust/genoio-io/src/plink2.rs` | existing-code | PGEN hardcall/dosage decoding and current phased-track rejection points. | high |
| SRC-4 | `rust/genoio-io/src/bgen.rs` | existing-code | BGEN Layout 2 dosage decoding and indexed region pushdown. | high |
| SRC-5 | `rust/genoio-io/tests/vcf_haplotype.rs` | existing-code | Existing haplotype behavior for phased VCF and retained unphased records. | high |
| SRC-6 | https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html | format-spec | BGEN Layout 2 phased probability storage and probability scaling. | high |
| SRC-7 | https://www.cog-genomics.org/plink/2.0/formats | official-doc | PLINK2 file format overview and `.pgen`/`.psam`/`.pvar` roles. | high |
| SRC-8 | https://deepwiki.com/chrchang/plink-ng/4.1-pgen-format | secondary-reference | PGEN mode and record-type summary including phased hardcalls and phased dosages. Use as orientation, not sole source for bit-level implementation. | med |

## Model Option Analysis (Required When `suggested-model`)
| Candidate ID | Model Family | When It Fits | Key Assumptions | Failure Modes | Supporting Citation(s) | Selection Status |
| --- | --- | --- | --- | --- | --- | --- |
| n/a | n/a | n/a | n/a | n/a | n/a | n/a |

## Existing Codebase Port Contract (Required When `existing-codebase-port`)
- Porting objective: Extend the existing `genoio` I/O architecture with dense
  PLINK2 and BGEN haplotype reads while preserving current matrix orientation,
  metadata, filtering, windowing, and error patterns.
- Source selection confirmation: The user requested work in the current
  repository and branch. The implementation should follow the existing native
  Rust reader pattern instead of introducing external reader dependencies.

### Source Pin
| Source ID | Source Type (`local-directory` or `github-url`) | Path/URL | Commit/Tag | Notes |
| --- | --- | --- | --- | --- |
| PORT-SRC-1 | local-directory | `/Users/nicholas/Projects/genoio` | current `main` branch at design time | Existing Python/Rust package. |

### Behavior Inventory And Parity Targets
| Behavior ID | Surface (`cli`/`api`/`numerics`/`io`) | Current Behavior | Target Behavior | Evidence Plan (tests/golden outputs) |
| --- | --- | --- | --- | --- |
| PORT-BHV-1 | api | `Dataset.read(kind="haplo")` currently supports phased VCF only. | Add dense PLINK2 and BGEN haplotype support for source-encoded representations. | Python and Rust tests for supported and unsupported combinations. |
| PORT-BHV-2 | io | PLINK2 genotype reads support hardcalls and unphased dosages; phased tracks are currently rejected. | Decode explicit PLINK2 hardcall phase and phased dosage tracks only on haplotype paths. | PGEN fixtures covering phased hardcalls, phased dosages, unphased retained records, filters, sample selection, and windows. |
| PORT-BHV-3 | io | BGEN genotype dosage reads collapse phased records to diploid expected A1 dosage. | Add BGEN haplotype dosage reads returning one row per haplotype as expected A1 dosage. | BGEN fixtures covering phased records, missing samples, filters, sample selection, windows, and indexed region reads. |
| PORT-BHV-4 | io | Unsupported retained haplotype records fail for VCF. | Preserve fail-on-retained-unsupported behavior for PLINK2 and BGEN. | Negative tests for unphased retained records and unsupported multiallelic/ploidy cases. |
| PORT-BHV-5 | cli | Benchmark scripts cover format read performance. | Add benchmark coverage after correctness for new haplotype paths. | Benchmark scripts run on representative phased PLINK2/BGEN data when fixtures are provided. |

## Codebase Investigation Findings (Required When `existing-codebase-port`)
- Investigation mode: `local-directory`
- Investigation completion: yes
- Investigator: direct Codex investigation using `rg`, source reads, and test reads.

| Finding ID | Source Scope | Summary | Evidence (file:line or commit:path:line) | Status (`confirmed`/`discrepancy`/`addition`/`missing`) |
| --- | --- | --- | --- | --- |
| PORT-INV-1 | Python API | `Dataset.read(...)` already supports `kind="haplo"` as a public option but Python validation currently disallows dosage-backed haplotype reads. | `src/genoio/_api.py` | confirmed |
| PORT-INV-2 | PyO3 dispatch | `read_haplotypes_dense` currently dispatches VCF/BCF only and rejects BGEN; PLINK2 is not handled. | `rust/genoio-py/src/lib.rs` | confirmed |
| PORT-INV-3 | PLINK2 reader | Current PGEN validation rejects hardcall-phase-with-dosage and phased-dosage tracks. These are the key implementation seams for the new haplotype paths. | `rust/genoio-io/src/plink2.rs` | confirmed |
| PORT-INV-4 | BGEN reader | Current BGEN dosage reader decodes phased records by collapsing haplotype probabilities to diploid expected A1 dosage. | `rust/genoio-io/src/bgen.rs`; `rust/genoio-io/tests/bgen_dense.rs` | confirmed |
| PORT-INV-5 | Haplotype metadata | Existing haplotype sample frames use `source_sample_index` and `haplotype_index` to map rows back to source samples. | `src/genoio/_api.py`; `rust/genoio-io/tests/vcf_haplotype.rs` | confirmed |

## External Research Findings (When Triggered)
| Claim ID | Claim | Source URL | Source Type | Access Date | Confidence (high/med/low) |
| --- | --- | --- | --- | --- | --- |
| EXT-1 | BGEN Layout 2 supports phased data; when `Phased=1`, probabilities are stored by haplotype and allele, and stored integer values scale to probabilities. | https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html | official-spec | 2026-06-05 | high |
| EXT-2 | PLINK2 `.pgen` is the binary genotype table paired with `.pvar` and `.psam`; the format supports richer representations than PLINK1 BED. | https://www.cog-genomics.org/plink/2.0/formats | official-doc | 2026-06-05 | high |
| EXT-3 | PGEN storage modes/record types include phased hardcalls and phased dosages. | https://deepwiki.com/chrchang/plink-ng/4.1-pgen-format | secondary-reference | 2026-06-05 | med |

## Mathematical Sanity Checks
- Summary: This is an I/O representation design, not a new statistical model.
  Numeric sanity checks are limited to source-encoded probability scaling,
  missing-value propagation, row/column shape invariants, and allele-count
  orientation.
- Blocking issues: PLINK2 phased dosage bit-level details must be verified
  against the PGEN specification/reference behavior before implementation.
- Accepted risks: The design intentionally remains biallelic diploid for this
  milestone. Multiallelic and variable-ploidy retained records will fail.

Detailed artifacts:
- `.plans/design-plans/artifacts/2026-06-05-haplotype-io/model-symbol-table.md`
- `.plans/design-plans/artifacts/2026-06-05-haplotype-io/equation-to-code-map.md`

## Solver Strategy Decision
- User preference: Native, source-faithful I/O readers; avoid analysis-policy
  conversions.
- Chosen strategy: Incremental native readers in the existing Rust modules.
- Why this strategy: It matches current architecture, keeps dependencies stable,
  and leaves room to factor shared haplotype decoding later if repeated logic
  becomes substantial.

## Solver Translation Feasibility
- Summary: No numerical solver is introduced. The translation problem is binary
  source decoding into existing dense matrix contracts.
- Blocking constraints: PLINK2 phased dosage record decoding must be specified
  precisely enough for tests before coding.
- Custom-solver rationale (if chosen): n/a

Detailed artifact:
- `.plans/design-plans/artifacts/2026-06-05-haplotype-io/solver-feasibility-matrix.md`

## Layer Contracts
### Ingress
- Contract: Existing constructors resolve VCF/BCF, PLINK1, PLINK2, and BGEN
  sources. `Dataset.read(kind="haplo", dosage=...)` selects dense haplotype
  decoding when supported by the source representation.
- Rejection rules: Reject unsupported retained records, unsupported source
  combinations, sparse haplotype requests for PLINK2/BGEN, hardcall-from-dosage
  requests, multiallelic retained records, variable ploidy, and missing
  companion metadata required for sample IDs.

### Pipeline
- Contract: Python validates public options; PyO3 dispatches by source format
  and dosage source; Rust decodes source-native records, applies filters,
  retained-variant windows, and sample selection.
- Validation-first checks: Validate `kind`, `dosage`, `sparse`, filter IR,
  sample filters, source support, PGEN layout/record type, BGEN Layout 2 block
  structure, sample count, ploidy, phasedness, and allele count before
  producing output.

### Numerics
- Contract: Haplotype matrices keep samples as row groups and variants as
  columns. Each retained diploid sample contributes two rows, ordered
  haplotype 0 then haplotype 1 in source sample order. Values count expected
  A1 dosage per haplotype: `0/1` for hardcall phase, fractional values for
  source-encoded phased dosages/probabilities.
- Result/status semantics: Missing haplotypes propagate through the existing
  dense missing-mask handling. Unsupported retained representations raise
  typed public errors rather than returning partial results.

### Egress
- Contract: Python returns the existing dense read shapes: matrix alone or
  matrix plus requested `samples`/`variants` Polars frames. Haplotype sample
  frames include `source_sample_index` and `haplotype_index`.
- Output/exit-code mapping: n/a; no new CLI command is introduced. Benchmark
  scripts may add haplotype modes for timing.

## Data Conversion and Copy Strategy
For each source format, record copy mode (`zero-copy`, `mmap`, `single-copy fallback`) and rationale.

- PLINK2: single-copy fallback. Rust decodes packed PGEN records into
  variant-major buffers, transposes to sample/haplotype-major output, and PyO3
  transfers values to Python for NumPy assembly.
- BGEN: single-copy fallback. Rust decodes packed probability blocks into
  selected haplotype values and transposes/assembles dense output as existing
  BGEN genotype dosage reads do.
- Python output: NumPy arrays remain the public dense matrix boundary; no
  zero-copy contract is introduced in this design.

## Multi-Input Reconciliation Contract (Required When Multiple Tabular Sources Feed Numerics)
- Sources: n/a
- Entity key(s) (for example subject/sample ID): n/a
- Join type and rationale: n/a
- Duplicate-key policy: n/a
- Missing-key policy: n/a
- Row-order freeze policy: n/a
- Reconciliation accounting (matched/dropped/retained counts): n/a
- Conversion boundary (where reconciled tabular data becomes arrays/PyTrees): n/a

## Validation Strategy
- Boundary checks: Public API validation rejects invalid `kind`, `dosage`, and
  sparse combinations before dispatch. Rust validates source headers, record
  dimensions, sample counts, phasedness, ploidy, biallelic constraints, and
  payload lengths.
- Shape/range/domain checks: Output shape is `(2 * retained_samples,
  retained_variants)` for diploid haplotype reads. BGEN probabilities must be
  in the encoded range and convert through the BGEN bit-depth scale. PGEN
  dosage missing sentinels must set the missing mask.
- Multi-input alignment checks (key uniqueness, overlap expectations, deterministic row ordering):
  Sample filters preserve source order. PLINK2 `.psam` and `.pvar` counts must
  match `.pgen`. BGEN `.sample` or embedded sample IDs must match header sample
  count.
- Failure semantics: Unsupported retained records fail the whole read. Records
  rejected by metadata filters may be skipped before decoding. Records rejected
  by genotype-stat filters may fail only when decoding is required to evaluate
  the retained candidate.

## Testing and Verification Strategy
- TDD scope: Add failing tests first for Python API support/rejection,
  PyO3 dispatch, PLINK2 hardcall haplotypes, PLINK2 phased dosage haplotypes,
  BGEN haplotype dosages, filters, sample selection, region/window iteration,
  and metadata row mapping.
- Regression strategy: Preserve all existing VCF haplotype, BGEN genotype
  dosage, PLINK2 genotype hardcall/dosage, cross-format parity, docs, and
  benchmark CLI tests. Add negative tests for hardcall-from-dosage requests and
  unsupported retained records.
- Verification commands: `make verify`; targeted Rust tests under
  `rust/genoio-io/tests/`; targeted Python tests under `tests/`; benchmark
  scripts after correctness.

## Implementation Phases
<!-- START_PHASE_1 -->
### Phase 1: Public Contract and Dispatch
**Goal:** Make the public API and PyO3 dispatch express the intended dense
haplotype combinations without implementing all decoders yet.

**Components:**
- `src/genoio/_api.py` option validation and source-support rules for
  `kind="haplo"` with `dosage="hardcall"` or `dosage="dosage"`.
- `rust/genoio-py/src/lib.rs` haplotype dense dispatch for VCF/BCF, PLINK2,
  and BGEN, while keeping sparse haplotype support limited to current behavior.
- Public docs and API reference text that state source-faithful representation
  rules.

**Dependencies:** Confirmed Definition of Done and existing VCF haplotype
behavior.

**Done when:** Python API tests show supported requests route to the correct
backend surface, unsupported sparse/hardcall-from-dosage combinations fail with
actionable errors, and existing VCF haplotype behavior is unchanged.
<!-- END_PHASE_1 -->

<!-- START_PHASE_2 -->
### Phase 2: BGEN Dense Haplotype Dosage Reads
**Goal:** Decode BGEN Layout 2 phased probability records into dense haplotype
dosage rows.

**Components:**
- `rust/genoio-io/src/bgen.rs` BGEN haplotype dosage reader following existing
  BGEN dosage filtering, sample selection, windowing, and indexed region
  pushdown patterns.
- BGEN tests for phased biallelic diploid probability records, missing samples,
  sample filters, metadata filters, genotype-stat filters where applicable,
  retained-variant windows, and `.bgi` region reads.
- Python tests proving `genoio.bgen(...).read(kind="haplo",
  dosage="dosage")`, `iter_blocks(...)`, and `iter_regions(...)` return the
  expected shapes and metadata.

**Dependencies:** Phase 1 dispatch contract.

**Done when:** BGEN dense haplotype dosage reads return `(2 * samples,
variants)` matrices with source-ordered haplotype rows, phased probabilities
are converted to expected A1 haplotype dosage, unsupported retained BGEN
records fail, and all targeted Rust/Python tests pass.
<!-- END_PHASE_2 -->

<!-- START_PHASE_3 -->
### Phase 3: PLINK2 Dense Haplotype Reads
**Goal:** Decode explicit PLINK2 hardcall phase and phased dosage records into
dense haplotype rows.

**Components:**
- `rust/genoio-io/src/plink2.rs` PGEN haplotype hardcall and phased dosage
  readers integrated with existing PVAR/PSAM parsing, record-window logic,
  sample filters, and retained-variant filtering.
- PLINK2 tests for fixed-width and variable-width phased representations where
  supported by the PGEN record layout, including unphased retained-record
  failures and metadata-filtered skip behavior.
- Python tests proving `genoio.pfile(...).read(kind="haplo",
  dosage="hardcall")` and `genoio.pfile(...).read(kind="haplo",
  dosage="dosage")` preserve shape, metadata, missing handling, and iteration
  contracts.

**Dependencies:** Phase 1 dispatch contract and PGEN bit-level design validated
against the PGEN specification/reference behavior.

**Done when:** PLINK2 dense haplotype hardcall and phased dosage reads pass
targeted Rust/Python tests, unsupported retained records fail, and existing
genotype hardcall/dosage behavior is unchanged.
<!-- END_PHASE_3 -->

<!-- START_PHASE_4 -->
### Phase 4: Documentation, Benchmarks, and Release Checks
**Goal:** Make the new I/O contract visible, reproducible, and ready for
release consideration.

**Components:**
- Docs updates in `docs/formats.md`, `docs/api/reading.md`, `docs/faq.md`, and
  relevant examples or performance notes.
- Benchmark script updates for dense haplotype paths after correctness is
  established.
- Release-readiness checks covering CI, wheel builds, docs, and benchmark
  instructions.

**Dependencies:** Phases 2 and 3.

**Done when:** Docs build strictly, benchmarks can time the new haplotype read
paths when representative inputs are available, `make verify` passes, and
release notes clearly state supported and unsupported haplotype representations.
<!-- END_PHASE_4 -->

## Simulation And Inference-Consistency Validation
- In scope: no
- Simulate entrypoint/signature: n/a
- Inputs: n/a
- Outputs: n/a
- Seed/RNG policy: n/a

### Assumption Alignment
| Inference Assumption | Simulation Rule | Mismatch Risk | Mitigation |
| --- | --- | --- | --- |
| n/a | n/a | n/a | n/a |

### Planned Validation Experiments
| Experiment ID | Type (recovery/SBC/PPC) | Success Criterion | Notes |
| --- | --- | --- | --- |
| n/a | n/a | n/a | Simulation is out of scope for this I/O design. |

## Risks and Open Questions
| ID | Risk or Question | Severity | Mitigation or Next Step | Owner |
| --- | --- | --- | --- | --- |
| R1 | PLINK2 phased dosage record details may be more complex than current genotype dosage overlays. | high | Validate bit-level interpretation against PGEN spec/reference fixtures before implementation. | implementation owner |
| R2 | Supporting both BGEN indexed region reads and haplotype dosage may duplicate genotype dosage code. | med | Start incremental in `bgen.rs`; factor shared helpers only after tests expose meaningful duplication. | implementation owner |
| R3 | Genotype-stat filters on haplotype dosage reads may be ambiguous if existing stats assume diploid genotype columns. | resolved | Compute `maf`, `mac`, and `missing_rate` from collapsed diploid expected A1 dosage for filter evaluation, then return source-faithful haplotype rows for retained variants. | user decision |
| R4 | Public error mapping may blur malformed source vs unsupported representation. | med | Add tests for exception classes and messages at Python API boundaries. | implementation owner |

## Additional Considerations
- Source-faithful I/O means the library may expose fractional haplotype dosage
  values. Users who want hard haplotypes from probabilities must threshold
  downstream in their analysis code.
- The design leaves capability introspection optional. It may become useful
  later, but it is not required to implement the new read paths.
- Genotype-stat filters on haplotype reads use diploid per-sample expected A1
  dosage for filter evaluation. This keeps MAF/MAC/missingness semantics aligned
  with normal variant QC while preserving haplotype values in the returned
  matrix.

## Acceptance Criteria
### haplotype-io.AC1: Public API Contract
- `haplotype-io.AC1.1`: `Dataset.read(kind="haplo", dosage="hardcall")`
  supports VCF and PLINK2 when retained records encode phased hardcalls.
- `haplotype-io.AC1.2`: `Dataset.read(kind="haplo", dosage="dosage")`
  supports PLINK2 phased dosage and BGEN phased probability records.
- `haplotype-io.AC1.3`: Sparse PLINK2/BGEN haplotype requests fail with
  `UnsupportedRepresentation` and an actionable message.
- `haplotype-io.AC1.4`: Hardcall-from-dosage conversion is never performed or
  implied by defaults.

### haplotype-io.AC2: BGEN Haplotype Dosage
- `haplotype-io.AC2.1`: BGEN phased biallelic diploid Layout 2 records return
  expected A1 dosage per haplotype row.
- `haplotype-io.AC2.2`: BGEN haplotype rows are ordered by source sample, then
  haplotype index.
- `haplotype-io.AC2.3`: Missing BGEN samples propagate through the existing
  dense missing-data policies.
- `haplotype-io.AC2.4`: BGEN region filters use `.bgen.bgi` pushdown when
  available and preserve source variant order.
- `haplotype-io.AC2.5`: Unsupported retained BGEN records fail the read.

### haplotype-io.AC3: PLINK2 Haplotype Reads
- `haplotype-io.AC3.1`: Explicit PLINK2 hardcall phase records return `0`/`1`
  A1 haplotype rows.
- `haplotype-io.AC3.2`: Explicit PLINK2 phased dosage records return expected
  A1 dosage per haplotype row.
- `haplotype-io.AC3.3`: Unphased retained PLINK2 records fail haplotype reads.
- `haplotype-io.AC3.4`: Metadata-filtered unsupported PLINK2 records may be
  skipped before decoding.
- `haplotype-io.AC3.5`: Existing PLINK2 genotype hardcall and unphased dosage
  reads remain unchanged.

### haplotype-io.AC4: Metadata, Iteration, and Filtering
- `haplotype-io.AC4.1`: Haplotype sample metadata includes
  `source_sample_index` and `haplotype_index`.
- `haplotype-io.AC4.2`: `iter_blocks(...)` yields retained-variant windows for
  supported haplotype reads.
- `haplotype-io.AC4.3`: `iter_regions(...)` yields one haplotype read per
  region filter for supported sources.
- `haplotype-io.AC4.4`: Returned variant metadata remains in matrix-column
  order after filtering.
- `haplotype-io.AC4.5`: On haplotype reads, `maf`, `mac`, and `missing_rate`
  filters are evaluated from collapsed diploid expected A1 dosage, while the
  returned matrix remains haplotype-level.

### haplotype-io.AC5: Verification, Docs, and Benchmarks
- `haplotype-io.AC5.1`: Rust and Python tests cover success and failure cases
  for PLINK2 and BGEN haplotype reads.
- `haplotype-io.AC5.2`: Docs state supported haplotype representations and
  non-goals clearly.
- `haplotype-io.AC5.3`: Benchmark scripts can time the new dense haplotype
  paths after correctness is complete.
- `haplotype-io.AC5.4`: `make verify` passes before the feature is considered
  complete.

## Glossary
- **A1**: The allele counted by returned `genoio` genotype or haplotype values.
- **Dense read**: A read returning a NumPy matrix rather than a SciPy sparse
  matrix.
- **Haplotype read**: A read with one row per haplotype instead of one row per
  diploid sample.
- **Hardcall**: Discrete genotype or haplotype allele calls encoded in the
  source file.
- **Phased dosage**: Source-encoded dosage/probability data that distinguishes
  haplotypes.
- **Retained record**: A source variant record that passes metadata/genotype
  filters and is eligible for matrix output.
- **Source-faithful I/O**: Reading representations encoded in the source data
  without applying downstream analysis conversions such as dosage thresholding.

## Status Transition Log
| Date | From | To | Why | By |
| --- | --- | --- | --- | --- |
| 2026-06-05 | N/A | Draft | Plan created | |
| 2026-06-05 | Draft | Approved | User approved incremental native-reader design and resolved genotype-stat filter behavior. | Codex |
