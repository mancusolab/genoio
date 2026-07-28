# Persistent Block Reader Sessions Design

## Status
Approved for Implementation

## Handoff Decision
- Current decision: approved
- Ready for implementation: yes
- Blocking items: None.

## Metadata
- Date: 2026-07-28
- Slug: persistent-block-reader
- Artifact Directory: `.plans/design-plans/artifacts/2026-07-28-persistent-block-reader`

## Summary
The design replaces `iter_blocks()`'s repeated retained-window reads with a
stateful Rust `BlockReader`. A backend-neutral facade dispatches to persistent
VCF, BCF, BGEN, PLINK1, or PLINK2 sessions. Each session owns its file handles,
parsed header, sample and filter state, source cursor, decoder state,
diagnostics, and reusable buffers, allowing successive blocks to continue from
the prior position.

A private PyO3 class exposes the reader to `Dataset._block_iterator()`.
Construction remains inside the Python generator to preserve lazy source
opening; Rust performs opening and decoding without holding the GIL, then
existing adapters create the public Python results. Rust ownership also ties
file cleanup to generator exhaustion, closure, failure, or garbage collection.

## Problem Statement
`Dataset.iter_blocks()` currently produces bounded retained-variant chunks by
starting a new Rust read for every block. Each call rebuilds options, reopens
source and companion files, reparses headers and sample selection, and—on
sequential or stateful formats—rescans or redecodes an increasingly long
prefix. Memory remains bounded, but total work can grow quadratically with the
number of variants. The implementation needs persistent per-iterator native
state without changing the public Python read contract.

## Definition of Done
- All existing `iter_blocks()` formats and modes use one persistent native
  session per iterator, eliminating repeated source reopening and prefix
  decoding.
- Outputs, filtering, metadata, lazy error timing, independent-iterator
  behavior, indexed pushdown, and bounded memory remain compatible.
- Tests verify read equivalence, resource cleanup, open-once behavior, and
  at-most-once candidate decoding.
- `read()` and `iter_regions()` remain unchanged. The latter's scaling and
  inaccurate BCF-index documentation become separate follow-up work.

## Goals and Non-Goals
### Goals
- Give every supported `Dataset.iter_blocks()` path one persistent native
  reader session per iterator.
- Reduce total source traversal from repeated-prefix work to one pass over
  candidate records.
- Preserve the current Python API, block contents, filtering, metadata,
  diagnostics, lazy failures, and bounded-memory behavior.
- Prove persistence with deterministic open and decode counters.

### Non-Goals
- Changing `Dataset.read()` or its stateless Rust entry points.
- Improving `Dataset.iter_regions()` scaling or adding BCF index pushdown.
- Adding a public session API, background prefetch, or concurrent block
  decoding.
- Changing genotype decoding, filter definitions, missing-value semantics, or
  output types.

## Existing Patterns
- `src/genoio/_api.py` owns public option validation, lazy Python generators,
  and conversion of native buffers to NumPy, SciPy, and Polars objects.
- `rust/genoio-py/src/reads.rs` parses Python arguments into Rust-owned values,
  releases the GIL around I/O, and delegates output conversion to
  `rust/genoio-py/src/output.rs`.
- `rust/genoio-py/src/errors.rs` contains the panic boundary and maps structured
  `GenoioError` variants to stable Python exception classes.
- `rust/genoio-io` owns source-format dispatch, parsing, filtering, decoding,
  and matrix construction. The Python adapter does not contain format logic.
- `rust/genoio-io/src/bgen/session.rs` already models an owned BGEN file/header
  session, but its lifetime currently ends after one read call.
- `rust/genoio-io/src/retention.rs` defines retained-variant window semantics.
  The persistent reader keeps the same ordering and filter decisions while
  replacing increasing retained-window offsets with a forward cursor.

This design adds `rust/genoio-io/src/blocks.rs` because the cross-format
stateful reader is a stable shared contract and lifecycle boundary. It adds
`rust/genoio-py/src/blocks.rs` because a mutable native iterator is distinct
from the stateless functions in `reads.rs`. No new Python module is justified;
`Dataset.iter_blocks()` remains in `src/genoio/_api.py`. Format sessions stay
inside the existing VCF, BGEN, and PLINK module trees; their exact file split is
decided by cohesion and file size during implementation.

## Architecture

`genoio-io` exposes a backend-neutral contract:

```text
BlockReader::open(source, options, block_size) -> Result<BlockReader>
BlockReader::next_block(&mut self) -> Result<Option<BlockOutput>>

BlockOutput = Dense(DenseGenotypeMatrix) | Sparse(SparseGenotypeMatrix)
```

`BlockReader` contains a private enum selecting a text VCF, BCF, BGEN, PLINK1,
or PLINK2 session. The facade performs dispatch only. Each format session owns
its files, parsed header and sample selection, variant filter, source position,
cumulative diagnostics, decoder state, reusable record buffers, and scratch
space. Existing stateless `read_*_windowed` functions remain available to
`read()` and `iter_regions()`.

The private PyO3 `BlockReader` class owns the `genoio-io` reader behind
synchronized mutable access. Construction and `next_block()` release the GIL.
Conversion into Python-owned arrays and Arrow-backed metadata uses the existing
output adapters under the GIL.

`Dataset._block_iterator()` constructs the native reader inside the Python
generator body. Source opening therefore remains lazy until the first
`next()`. The generator repeatedly requests native blocks and applies the
existing public result conversion. EOF, `generator.close()`, exceptions, and
garbage collection drop the native reader and close its files through Rust
ownership. Separate Python iterators create independent sessions.

## Model Acquisition Path
- Path: `provided-model`
- Why this path: This infrastructure-only refactor has no scientific model.
  The readiness schema's `provided-model` category represents the provided
  behavioral specification in the existing API, Rust readers, and regression
  tests. No statistical or inferential model is added or selected.
- User selection confirmation: The user approved the complete persistent-reader
  scope and the shared Rust facade.

## Required Workflow States
- model_path_decided: yes
- codebase_investigation_complete_if_port: n/a
- simulation_contract_complete_if_in_scope: n/a

## Model Specification Sources
| Source ID | Path/Link | Type | Notes | Confidence (high/med/low) |
| --- | --- | --- | --- | --- |
| SRC-1 | `src/genoio/_api.py`, `rust/genoio-io`, `rust/genoio-py`, `tests/test_blocks.py` | Local behavioral specification | Existing public behavior and reader contracts; no scientific model. | high |

## Model Option Analysis (Required When `suggested-model`)
Not applicable.

## Existing Codebase Port Contract (Required When `existing-codebase-port`)
Not applicable; no external implementation is being ported.

## Codebase Investigation Findings (Required When `existing-codebase-port`)
- Investigation mode: `local-directory`
- Investigation completion: yes
- Investigator: `scientific-codebase-investigation-pass`

| Finding ID | Source Scope | Summary | Evidence (file:line or commit:path:line) | Status (`confirmed`/`discrepancy`/`addition`/`missing`) |
| --- | --- | --- | --- | --- |
| INV-1 | Python blocks | Each block invokes a fresh retained-window read with an increasing start offset. | `src/genoio/_api.py:593` | confirmed |
| INV-2 | Python regions | `iter_regions()` independently calls `read()` per region and remains outside this design. | `src/genoio/_api.py:584` | confirmed |
| INV-3 | PyO3 boundary | Stateless reads already parse Python values before releasing the GIL and convert outputs afterward. | `rust/genoio-py/src/reads.rs:68`, `rust/genoio-py/src/errors.rs:29` | confirmed |
| INV-4 | BGEN | `BgenReadSession` already owns the reader and header, but only for one read call. | `rust/genoio-io/src/bgen/session.rs:39` | confirmed |
| INV-5 | Retention | Current window state counts variants after filters rather than raw source records. | `rust/genoio-io/src/retention.rs:26` | confirmed |
| INV-6 | Indexed text VCF | The CSI query borrows its BGZF reader, so the current query cannot be stored directly in an owned session. | `rust/genoio-io/src/vcf/text/source.rs:201` | confirmed |
| INV-7 | PLINK2 | Variable-width PGEN decoding advances `PgenDecoderState` sequentially before metadata rejection. | `rust/genoio-io/src/plink/plink2/dense.rs:150` | confirmed |
| INV-8 | BCF documentation | The Python docstring claims BCF index pushdown, but BCF bypasses the indexed text-VCF route. | `src/genoio/_api.py:565`, `rust/genoio-io/src/vcf.rs:71` | discrepancy |

## External Research Findings (When Triggered)
Not triggered. The design depends on local ownership and behavior contracts, not
an external API or new dependency.

## Mathematical Sanity Checks
- Summary: Let `C` be candidate source records, `R` retained variants, and `S`
  block size. Repeated retained-window reads can require
  `O(C * ceil(R / S))` work and become `O(R² / S)` when all candidates are
  retained. A persistent cursor reduces traversal to `O(C)`. Genotype output
  memory remains `O(samples * S)` plus reusable decode scratch and an existing
  index-record list.
- Blocking issues: None.
- Accepted risks: Native indexes may retain `O(index records)` positions, as
  they do in the current indexed BGEN and VCF paths.

Detailed artifacts:
- `.plans/design-plans/artifacts/2026-07-28-persistent-block-reader/model-symbol-table.md`
- `.plans/design-plans/artifacts/2026-07-28-persistent-block-reader/equation-to-code-map.md`

## Solver Strategy Decision
- User preference: Not applicable.
- Chosen strategy: Not applicable.
- Why this strategy: No numerical solver is involved.

## Solver Translation Feasibility
- Summary: Not applicable.
- Blocking constraints: None.
- Custom-solver rationale (if chosen): Not applicable.

Detailed artifact:
- `.plans/design-plans/artifacts/2026-07-28-persistent-block-reader/solver-feasibility-matrix.md`

## Layer Contracts
### Ingress
- Contract: `Dataset.iter_blocks(size, **read_options)` keeps its current public
  signatures and return types. Python validates block size, read options, and
  representation support before returning the generator.
- Rejection rules: Invalid options fail eagerly. Source opening, header errors,
  and record errors remain lazy and occur on the `next()` call that reaches
  them.

### Pipeline
- Contract: One `genoio-io::BlockReader` owns one forward-moving format session.
  `next_block()` returns at most `block_size` variants that survive all filters,
  in source order.
- Validation-first checks: Each backend validates its header and companion
  relationships once. Metadata predicates run before genotype decoding when
  the format allows it. Genotype-stat predicates run before a variant consumes
  an output slot.

### Numerics
- Contract: Existing hardcall, dosage, haplotype, missing-value, and sparse
  encoding routines remain authoritative. The session changes their lifetime,
  not their numerical behavior.
- Result/status semantics: A retained record is fully consumed before a block
  returns. Filter rejects do not consume a retained-variant slot. Cumulative
  diagnostics retain their current prefix meaning.

### Egress
- Contract: `BlockOutput` contains the existing dense or sparse core matrix
  structure. PyO3 transfers owned buffers to NumPy and Arrow adapters, and
  Python constructs the same NumPy/SciPy/Polars result shape as `read()`.
- Output/exit-code mapping: `Some(block)` yields one public result; `None`
  produces `StopIteration`. Existing structured Rust errors map to the existing
  Python exception classes. Panics map to `RustInternalError`.

## Data Conversion and Copy Strategy
All formats decode one block into Rust-owned dense or sparse buffers. PyO3 moves
those buffers into NumPy without allocating a second full block buffer.
Metadata uses the existing Arrow C stream adapter before Polars consumes it.
No format allocates a full-dataset genotype matrix.

| Source format | Source decode strategy | Python transfer |
| --- | --- | --- |
| Text VCF / BCF | Reuse record and genotype scratch; write accepted variants into the current block. | Move dense or sparse vectors to NumPy; transfer requested metadata through Arrow. |
| BGEN | Reuse probability/decompression buffers in `BgenReadSession`; write accepted variants into the current block. | Same owned-vector and Arrow transfer. |
| PLINK1 | Reuse packed hardcall and batch buffers while BED/BIM cursors advance. | Same owned-vector and Arrow transfer. |
| PLINK2 | Reuse `PgenDecoderState`, packed/dosage buffers, and PVAR record state. | Same owned-vector and Arrow transfer. |

## Multi-Input Reconciliation Contract (Required When Multiple Tabular Sources Feed Numerics)
- Sources: Existing PLINK BED/BIM/FAM, PGEN/PVAR/PSAM, and optional BGEN sample
  companions.
- Entity key(s) (for example subject/sample ID): Existing source sample IDs and
  source variant order.
- Join type and rationale: No new join is introduced. Existing parsers and
  `select_samples_source_order` remain authoritative.
- Duplicate-key policy: Preserve current parser and selector behavior.
- Missing-key policy: Preserve current sample-filter errors and diagnostics.
- Row-order freeze policy: Selected samples and retained variants remain in
  source order.
- Reconciliation accounting (matched/dropped/retained counts): Preserve
  cumulative `DenseDiagnostics`.
- Conversion boundary (where reconciled tabular data becomes arrays/PyTrees):
  Format sessions construct core matrix and metadata buffers before the PyO3
  boundary.

## Validation Strategy
- Boundary checks: Preserve eager Python option validation and Rust source
  capability validation. The native constructor accepts Rust-owned paths,
  filters, and options only.
- Shape/range/domain checks: Preserve block-size validation, source header
  validation, payload-length checks, allele/ploidy constraints, and sparse
  missing-value rejection.
- Multi-input alignment checks (key uniqueness, overlap expectations, deterministic row ordering):
  Preserve current companion counts, requested-sample selection, and
  source-order alignment checks.
- Failure semantics: Do not prefetch or validate past the current block. A
  native or conversion error terminates the Python generator and drops the
  session. EOF is sticky and performs no further I/O.

## Testing and Verification Strategy
- TDD scope: Add failing Rust contract tests before each backend session and
  failing Python integration tests before the PyO3 cutover.
- Regression strategy: Concatenate session blocks and compare them with a
  whole read for every supported representation. Include empty, exact-width,
  partial-final, filtered, indexed, metadata, missing-value, and delayed-error
  cases. Use test-only counters for source opens and candidate decodes instead
  of timing thresholds.
- Verification commands:
  - `cargo test --workspace`
  - `.venv/bin/pytest -p no:capture`
  - `pre-commit run --all-files`

## Implementation Phases
<!-- START_PHASE_1 -->
### Phase 1: Shared block-reader contract
**Goal:** Establish the backend-neutral lifecycle and output contract without
changing Python behavior.

**Components:**
- `rust/genoio-io/src/blocks.rs` — owns `BlockReader`, `BlockOutput`, common
  options, sticky EOF, and cross-format dispatch.
- Existing retention and matrix components — provide retained-variant and
  output semantics without duplicating passive types in new modules.
- Rust contract tests — cover block capacity, retained boundaries, cumulative
  diagnostics, sticky EOF, and bounded allocation.

**Dependencies:** None.

**Done when:** The shared contract compiles, its lifecycle tests pass, and no
public Python path uses it yet. Covers
`persistent-block-reader.AC2.2`, `persistent-block-reader.AC2.6`,
`persistent-block-reader.AC3.5`, and `persistent-block-reader.AC4.1`.
<!-- END_PHASE_1 -->

<!-- START_PHASE_2 -->
### Phase 2: BGEN and PLINK1 sessions
**Goal:** Make BGEN and PLINK1 block traversal persistent across every supported
representation.

**Components:**
- Existing BGEN modules — extend `BgenReadSession` to own sample selection,
  sequential or `.bgi` cursor state, diagnostics, and reusable decode buffers.
- Existing PLINK1 modules — retain BED/BIM/FAM readers, selected samples,
  packed decoder state, variant position, and batch scratch.
- Backend contract tests — verify whole-read parity, filters, metadata,
  open-once behavior, at-most-once decoding, and early drop.

**Dependencies:** Phase 1.

**Done when:** All supported BGEN and PLINK1 session modes pass their parity and
work-counter tests. Covers `persistent-block-reader.AC1.2`,
`persistent-block-reader.AC1.3`, `persistent-block-reader.AC2.1`,
`persistent-block-reader.AC2.3`, and `persistent-block-reader.AC4.2` for these
formats.
<!-- END_PHASE_2 -->

<!-- START_PHASE_3 -->
### Phase 3: PLINK2 sessions
**Goal:** Preserve PGEN/PVAR decoder state across blocks for every supported
PLINK2 representation.

**Components:**
- Existing PLINK2 modules — retain PGEN/PVAR/PSAM readers, the parsed PGEN
  header, sample selection, current variant position, `PgenDecoderState`, and
  reusable hardcall/dosage/haplotype buffers.
- PLINK2 contract tests — cover fixed-width, variable-width, LD-compressed,
  filtered, dense, sparse, dosage, haplotype, metadata, and delayed PVAR error
  paths.

**Dependencies:** Phase 1.

**Done when:** PLINK2 records that depend on prior decoder state cross block
boundaries correctly, and all parity and work-counter tests pass. Covers
`persistent-block-reader.AC1.2`, `persistent-block-reader.AC1.3`,
`persistent-block-reader.AC2.1`, `persistent-block-reader.AC2.4`,
`persistent-block-reader.AC3.2`, and `persistent-block-reader.AC4.2` for
PLINK2.
<!-- END_PHASE_3 -->

<!-- START_PHASE_4 -->
### Phase 4: VCF and BCF sessions
**Goal:** Make text VCF and BCF traversal persistent while preserving current
indexed and sequential routing.

**Components:**
- Existing text VCF source modules — own sequential reader state and add an
  owned BGZF chunk cursor for indexed region filters without a self-referential
  query object.
- Existing BCF modules — retain the BGZF reader, parsed header, reusable record,
  selected samples, diagnostics, and decode buffers.
- VCF/BCF contract tests — cover uncompressed, BGZF sequential, tabix/CSI
  indexed, BCF, filters, representations, delayed record errors, and early
  drop.

**Dependencies:** Phase 1.

**Done when:** VCF and BCF session modes pass parity and work-counter tests, and
indexed text VCF reads traverse only their existing query chunks. Covers
`persistent-block-reader.AC1.2`, `persistent-block-reader.AC1.3`,
`persistent-block-reader.AC2.1`, `persistent-block-reader.AC2.3`,
`persistent-block-reader.AC3.2`, and `persistent-block-reader.AC4.2` for VCF
and BCF.
<!-- END_PHASE_4 -->

<!-- START_PHASE_5 -->
### Phase 5: PyO3 adapter and Python cutover
**Goal:** Route `Dataset.iter_blocks()` through the complete native session
matrix without changing its public contract.

**Components:**
- `rust/genoio-py/src/blocks.rs` — owns the synchronized private native class,
  GIL release, panic/error boundary, and existing output conversion calls.
- PyO3 module registration and `src/genoio/_rust.pyi` — expose the private class
  consistently to runtime and typing checks.
- `src/genoio/_api.py` — lazily constructs the native reader inside
  `_block_iterator()` and removes repeated retained-window calls.
- Python integration tests — cover every supported format/mode, output parity,
  eager option validation, lazy source and record failures, independent
  iterators, close/error cleanup, and no stateless fallback.

**Dependencies:** Phases 2, 3, and 4.

**Done when:** The complete Python suite passes with `pytest -p no:capture`,
Rust and stub-parity checks pass, deterministic counters prove the traversal
contract, and `read()` plus `iter_regions()` remain on their existing paths.
Covers all acceptance criteria.
<!-- END_PHASE_5 -->

## Simulation And Inference-Consistency Validation
- In scope: no
- Reason: This refactor does not change an inferential model or numerical
  estimator. Fixture-based I/O parity tests are the relevant validation.

## Risks and Open Questions
| ID | Risk or Question | Severity | Mitigation or Next Step | Owner |
| --- | --- | --- | --- | --- |
| R1 | Extracting sessions duplicates existing decode loops and allows behavior to diverge. | High | Refactor each existing loop around shared state and require whole-read parity tests in the same phase. | Implementer |
| R2 | The current indexed text-VCF query borrows its BGZF reader. | High | Use an owned chunk cursor with the same seek and chunk-boundary semantics; do not add a self-referential dependency. | Implementer |
| R3 | PLINK2 LD-compressed records lose prior decoder state at a block boundary. | High | Store one `PgenDecoderState` for the session and add a fixture whose dependent record starts a later block. | Implementer |
| R4 | A backend silently retains the stateless fallback. | High | Require the complete format/mode matrix and deterministic open/decode counters before the Python cutover. | Implementer |
| R5 | A session allocates metadata or genotype storage proportional to the full dataset. | Medium | Assert block capacities and monitor allocation shape in Rust contract tests. Existing index-record lists are the only allowed variant-count-sized state. | Implementer |
| R6 | Persistent diagnostics change private native results. | Medium | Keep cumulative counters in the session and compare each block's diagnostics with the equivalent retained-window read. | Implementer |
| R7 | BCF indexed-pushdown documentation remains inaccurate. | Low | Track as a separate follow-up; do not expand this refactor. | Maintainer |

## Additional Considerations
**Error timing:** Option validation stays eager. Source opening occurs on the
first `next()`. Record validation stops at the current block boundary, so a
malformed later record fails only when its block is requested.

**Threading:** The native object serializes mutable access and releases the GIL
during open and decode. It does not start a worker thread or prefetch blocks.

**Follow-up scope:** `iter_regions()` also repeats read setup and can rescan
unindexed sources once per region. Efficient arbitrary-region iteration needs
separate indexing and memory-policy decisions, so it is not part of this plan.

## Acceptance Criteria
### persistent-block-reader.AC1: Persistent traversal replaces stateless pagination
- **persistent-block-reader.AC1.1 Success:** Every currently supported format
  and representation uses `BlockReader`; no `iter_blocks()` path falls back to
  repeated retained-window reads.
- **persistent-block-reader.AC1.2 Success:** Each required source and companion
  file is opened once per iterator.
- **persistent-block-reader.AC1.3 Success:** Each candidate record is visited
  once and each required genotype payload is decoded at most once per iterator.
- **persistent-block-reader.AC1.4 Edge:** Two iterators over the same `Dataset`
  own independent sessions and can advance independently.

### persistent-block-reader.AC2: Public results remain compatible
- **persistent-block-reader.AC2.1 Success:** Concatenated blocks equal
  `Dataset.read()` for every supported dense/sparse, genotype/haplotype,
  hardcall/dosage, filtering, and metadata combination.
- **persistent-block-reader.AC2.2 Edge:** Blocks contain at most `size` retained
  variants in source order; empty, exact-width, partial-final, and all-filtered
  inputs behave correctly.
- **persistent-block-reader.AC2.3 Success:** Sample/variant filters and existing
  text-VCF/BGEN index pushdown produce the same results as current reads.
- **persistent-block-reader.AC2.4 Edge:** Variable-width and LD-compressed
  PLINK2 records decode correctly when dependencies cross block boundaries.
- **persistent-block-reader.AC2.5 Success:** Requested sample metadata is
  present on every block, and variant metadata remains aligned with matrix
  columns.
- **persistent-block-reader.AC2.6 Success:** Native cumulative diagnostics
  match the equivalent retained-window reads.
- **persistent-block-reader.AC2.7 Failure:** Missing-value and
  unsupported-representation errors retain their existing behavior and
  exception classes.

### persistent-block-reader.AC3: Laziness and lifecycle remain compatible
- **persistent-block-reader.AC3.1 Failure:** Invalid block sizes and read
  options fail eagerly without opening source files.
- **persistent-block-reader.AC3.2 Failure:** Source/header errors occur on the
  first `next()`, while malformed later records fail only when the affected
  block is requested.
- **persistent-block-reader.AC3.3 Failure:** Structured Rust errors retain
  their Python exception mapping; panics become `RustInternalError`.
- **persistent-block-reader.AC3.4 Edge:** Exhaustion, `generator.close()`, read
  errors, and garbage collection release session resources without affecting
  independent iterators.
- **persistent-block-reader.AC3.5 Edge:** Once EOF is observed, subsequent
  native calls return `None` without further I/O.

### persistent-block-reader.AC4: Memory and scaling are bounded
- **persistent-block-reader.AC4.1 Success:** Genotype output allocation is
  proportional to rows × block size, plus reusable decode scratch and existing
  index-position storage; no full genotype matrix is allocated.
- **persistent-block-reader.AC4.2 Success:** Deterministic open/decode counters
  demonstrate one setup and linear candidate traversal for every backend/mode.
  Timing benchmarks are informative, not pass/fail gates.

### persistent-block-reader.AC5: Excluded APIs remain unchanged
- **persistent-block-reader.AC5.1 Success:** `Dataset.read()` continues using
  its existing stateless paths and retains its output behavior.
- **persistent-block-reader.AC5.2 Success:** `Dataset.iter_regions()` behavior
  and routing remain unchanged; its scaling and BCF-index documentation
  discrepancy remain follow-up work.

## Glossary
- **Arrow C stream adapter:** An interoperability mechanism that transfers
  columnar metadata through the Arrow C interface without serializing it into
  an intermediate format.
- **BCF:** Binary Call Format, the binary counterpart to VCF for variant and
  genotype records.
- **BED/BIM/FAM:** The PLINK1 file trio containing binary genotypes, variant
  metadata, and sample or family metadata, respectively.
- **BGEN:** A binary genotype format designed to store variant data and genotype
  probabilities or dosages.
- **BGZF:** Blocked gzip compression that supports virtual offsets and random
  access within compressed genomic files.
- **`BlockReader`:** The proposed backend-neutral, stateful Rust facade that
  owns one format session and advances it block by block.
- **Candidate record:** A source record visited and evaluated by the reader; it
  may later be rejected by filters.
- **CSI:** Coordinate-Sorted Index, an index format that maps genomic regions to
  chunks in a coordinate-sorted source.
- **Dosage:** The expected alternate-allele count for a sample, often
  represented as a fractional value.
- **GIL:** Python's Global Interpreter Lock, which the native adapter releases
  while opening files and decoding blocks.
- **Hardcall:** A discrete genotype assignment rather than a probability or
  expected allele count.
- **Haplotype:** Allele states represented separately for individual chromosome
  copies rather than as a combined genotype.
- **Index pushdown:** Applying a region filter through an index so unrelated
  source chunks are not parsed or decoded.
- **LD-compressed:** A PGEN encoding that represents a genotype record relative
  to an earlier record using linkage disequilibrium, requiring decoder state to
  persist across blocks.
- **PGEN/PVAR/PSAM:** The PLINK2 file trio containing genotypes, variant
  metadata, and sample metadata, respectively.
- **PLINK1 / PLINK2:** Two generations of the PLINK genotype dataset formats,
  using BED/BIM/FAM and PGEN/PVAR/PSAM file sets.
- **Polars:** A columnar DataFrame library used here to consume Arrow-backed
  metadata.
- **PyO3:** A Rust framework for implementing Python extension modules and
  converting values across the Python-Rust boundary.
- **Retained variant:** A candidate variant that passes all applicable filters
  and occupies a block output slot.
- **Retained-window read:** The current stateless operation that selects a
  positional window of retained variants, potentially rescanning earlier
  records for every block.
- **Sticky EOF:** End-of-file state that remains recorded so later calls return
  immediately without further I/O.
- **tabix:** A genomic indexing and query convention for coordinate-sorted,
  BGZF-compressed text files.
- **VCF:** Variant Call Format, a text format for variant records and optional
  per-sample genotype fields.

## Status Transition Log
| Date | From | To | Why | By |
| --- | --- | --- | --- | --- |
| 2026-07-28 | N/A | Draft | Plan created | |
| 2026-07-28 | Draft | In Review | Architecture and acceptance criteria approved by the user. | Codex |
| 2026-07-28 | In Review | Approved for Implementation | Automated readiness validation passed. | Codex |
