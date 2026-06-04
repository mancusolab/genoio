# BGEN Dosage Support Design

## Status
Approved

## Handoff Decision
- Current decision: ready-for-implementation-planning
- Ready for implementation: yes
- Blocking items: none

## Metadata
- Date: 2026-06-04
- Slug: bgen-dosage
- Artifact Directory: `.plans/design-plans/artifacts/2026-06-04-bgen-dosage`

## Summary
BGEN support will add a fourth source format to `genoio` using the same dense
dosage contract already used for VCF `FORMAT/DS` and PLINK2 dosage tracks. The
first implementation deliberately targets a narrow, beta-suitable subset of the
BGEN v1.3 format: Layout 2, biallelic, diploid, unphased probability data with
real sample IDs.

The reader will return expected copies of `a1`, support sample and variant
metadata alignment, and reuse the existing filter, block-read, missing-data,
and dense matrix contracts. Broader BGEN features such as hardcall conversion,
phased data, multiallelic data, variable ploidy, sparse dosage, and index
pushdown are deferred.

## Problem Statement
Statistical genetics tools often receive imputed genotypes in BGEN format.
`genoio` currently offers a unified Python API for VCF, PLINK1, and PLINK2, but
BGEN users must rely on separate libraries and then manually reconcile sample
IDs, variant metadata, dosage orientation, filtering semantics, and matrix
shape. That creates avoidable integration work for downstream tools such as
eQTL mappers, where genotype matrices must be deterministically aligned with
phenotypes and covariates.

The design adds BGEN as a first-class `genoio` source while preserving the
existing library contract: samples on rows, variants on columns, allele-count
values for `a1`, explicit missing masks, deterministic metadata frames, and
validation errors for unsupported representations.

## Definition of Done
Design a first BGEN implementation for `genoio` that supports dense dosage
reads matching existing VCF/PLINK2 dosage behavior. The design requires real
sample IDs from an embedded BGEN sample block or a companion `.sample` file,
rejects anonymous/generated sample IDs, returns dosage as copies of `a1`,
supports metadata and genotype-stat filters, and preserves block-read behavior.
The first scope excludes hardcall conversion, sparse reads, haplotype reads,
phased BGEN, multiallelic BGEN, variable ploidy, and index/region pushdown.

## Goals and Non-Goals
### Goals
- Add BGEN source resolution with required sample IDs from either embedded
  sample identifiers or a same-prefix companion `.sample` file.
- Add dense BGEN dosage reads through the existing `kind="geno"` API when
  `dosage="dosage"` is requested.
- Compute dosage as expected copies of `a1` for biallelic unphased diploid
  probabilities.
- Reuse existing sample filtering, metadata filtering, genotype-stat filtering,
  block reads, missing-data handling, and dense matrix return contracts.
- Reject unsupported BGEN features with clear `UnsupportedRepresentation` or
  source/metadata errors.
- Keep the first implementation Rust-native inside `genoio-io`, with PyO3 only
  dispatching to the Rust reader.

### Non-Goals
- No BGEN hardcall conversion.
- No sparse BGEN reads.
- No haplotype reads from BGEN.
- No phased BGEN support.
- No multiallelic BGEN support.
- No variable-ploidy BGEN support.
- No `.bgi`/bgenix index support or region pushdown in the first slice.
- No generated anonymous sample IDs.
- No dependency on the symlinked Cython `bgen` package at runtime.

## Existing Patterns
- Source resolution lives in `src/genoio/_source.py` and returns a
  `ResolvedSource` with logical members. BGEN should follow this by adding
  `SourceFormat.BGEN` and member keys such as `"bgen"` and optionally
  `"sample"`.
- The Python `Dataset` API delegates dense reads to the private Rust extension
  through `rust/genoio-py/src/lib.rs`. BGEN should add a dense dispatch branch
  rather than introducing Python-side parsing.
- Format-specific Rust readers live in `rust/genoio-io/src/{vcf,plink1,plink2}.rs`
  and return `genoio-core` contracts. BGEN warrants a new module
  `rust/genoio-io/src/bgen.rs` because BGEN is an independent binary format
  boundary with substantial parsing logic.
- Dense reads return `DenseGenotypeMatrix`, with sample-major output after
  format decoding. BGEN should use the same container and diagnostics fields.
- Sample keep-lists are used for membership only; matrix row order remains
  source order through `select_samples_source_order`.
- Variant filters use `VariantFilter::partial_decision` to skip
  metadata-rejected records before genotype decoding. BGEN should follow this
  for metadata-only filters and compute dosage stats only when genotype-derived
  predicates require them.
- Missing values are represented by a dense value plus a missing mask in Rust,
  then mapped by Python according to the requested missing-data policy.

## Model Acquisition Path
- Path: `provided-model`
- Why this path: BGEN decoding is specified by the official BGEN v1.3 format
  contract, not a statistical model selection exercise.
- User selection confirmation: User requested BGEN support matching existing
  `genoio` behavior, required sample IDs, and provided the official BGEN v1.3
  specification URL plus a local Cython package reference.

## Required Workflow States
- model_path_decided: yes
- codebase_investigation_complete_if_port: n/a
- simulation_contract_complete_if_in_scope: n/a

## Model Specification Sources
| Source ID | Path/Link | Type | Notes | Confidence (high/med/low) |
| --- | --- | --- | --- | --- |
| SRC-1 | https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html | official specification | BGEN v1.3 header, sample block, variant block, Layout 2 probability encoding, compression, and flags. | high |
| SRC-2 | `bgen/` symlink to local Cython package | reference implementation | Useful behavioral reference for header/sample parsing and dosage APIs; not a runtime dependency. | medium |

## Model Option Analysis (Required When `suggested-model`)
| Candidate ID | Model Family | When It Fits | Key Assumptions | Failure Modes | Supporting Citation(s) | Selection Status |
| --- | --- | --- | --- | --- | --- | --- |
| N/A | N/A | N/A | N/A | N/A | N/A | N/A |

## Existing Codebase Port Contract (Required When `existing-codebase-port`)
- Porting objective: N/A
- Source selection confirmation: N/A

### Source Pin
| Source ID | Source Type (`local-directory` or `github-url`) | Path/URL | Commit/Tag | Notes |
| --- | --- | --- | --- | --- |
| N/A | N/A | N/A | N/A | N/A |

### Behavior Inventory And Parity Targets
| Behavior ID | Surface (`cli`/`api`/`numerics`/`io`) | Current Behavior | Target Behavior | Evidence Plan (tests/golden outputs) |
| --- | --- | --- | --- | --- |
| N/A | N/A | N/A | N/A | N/A |

## Codebase Investigation Findings (Required When `existing-codebase-port`)
- Investigation mode: N/A
- Investigation completion: N/A
- Investigator: direct local scan during design

| Finding ID | Source Scope | Summary | Evidence (file:line or commit:path:line) | Status (`confirmed`/`discrepancy`/`addition`/`missing`) |
| --- | --- | --- | --- | --- |
| INV-1 | Source resolution | Existing formats are represented by `SourceFormat` and `ResolvedSource` member maps. | `src/genoio/_source.py` | confirmed |
| INV-2 | PyO3 dispatch | Dense reads dispatch by source format and dosage source. | `rust/genoio-py/src/lib.rs` | confirmed |
| INV-3 | IO module shape | Each existing binary/text source format has a cohesive Rust reader module returning `genoio-core` containers. | `rust/genoio-io/src/lib.rs` | confirmed |
| INV-4 | Local BGEN reference | The symlinked Cython package parses headers, sample blocks, variants, and probability data, but generates integer IDs when no sample IDs exist. `genoio` must intentionally diverge by rejecting anonymous IDs. | `bgen/src/samples.cpp` | confirmed |

## External Research Findings (When Triggered)
| Claim ID | Claim | Source URL | Source Type | Access Date | Confidence (high/med/low) |
| --- | --- | --- | --- | --- | --- |
| EXT-1 | A BGEN file contains a header, optional sample identifier block, then consecutive variant data blocks. | https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html | official-doc | 2026-06-04 | high |
| EXT-2 | Header flags encode compression, layout, and sample identifier presence; Layout 2 is recommended for new files. | https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html | official-doc | 2026-06-04 | high |
| EXT-3 | Layout 2 probability blocks include sample count, allele count, min/max ploidy, per-sample ploidy/missingness bytes, phased flag, bit depth, and packed probabilities. | https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html | official-doc | 2026-06-04 | high |
| EXT-4 | Layout 2 unphased probabilities omit the final genotype probability, which is inferred as one minus the sum of stored probabilities. | https://www.chg.ox.ac.uk/~gav/bgen_format/spec/latest.html | official-doc | 2026-06-04 | high |

## Mathematical Sanity Checks
- Summary: The only computation is deterministic dosage conversion for
  biallelic unphased diploid probabilities:
  `dosage_a1 = P(AB) + 2 * P(BB)`.
- Blocking issues: None for the first supported subset.
- Accepted risks: Probability bit-depth quantization means decoded values may
  differ from source floating-point examples by the BGEN encoding resolution.

Detailed artifacts:
- `.plans/design-plans/artifacts/2026-06-04-bgen-dosage/model-symbol-table.md`
- `.plans/design-plans/artifacts/2026-06-04-bgen-dosage/equation-to-code-map.md`

## Solver Strategy Decision
- User preference: N/A
- Chosen strategy: No solver.
- Why this strategy: BGEN support is deterministic file parsing and data
  conversion, not statistical inference or optimization.

## Solver Translation Feasibility
- Summary: Not applicable.
- Blocking constraints: None.
- Custom-solver rationale (if chosen): Not applicable.

Detailed artifact:
- `.plans/design-plans/artifacts/2026-06-04-bgen-dosage/solver-feasibility-matrix.md`

## Layer Contracts
### Ingress
- Contract: Resolve `.bgen` sources into a `ResolvedSource` with a required
  BGEN member and an optional `.sample` companion member. The reader accepts
  either embedded sample identifiers or a companion sample file, but must
  produce real sample IDs before any matrix read or metadata return.
- Rejection rules: Reject missing BGEN files, unsupported extensions, absent
  sample IDs, mismatched sample counts, invalid BGEN magic/header values,
  unsupported layout values, unsupported compression values, malformed sample
  blocks, and malformed companion `.sample` files.

### Pipeline
- Contract: Python validates API options, then PyO3 dispatches BGEN dense
  dosage reads into `genoio-io`. The Rust BGEN module parses metadata,
  evaluates filters, decodes retained dosage values, and returns
  `DenseGenotypeMatrix`.
- Validation-first checks: `dosage="dosage"` is required for BGEN dense reads.
  Sparse and haplotype BGEN reads fail before decode. Hardcall BGEN reads fail
  before decode until a hardcall conversion policy is explicitly designed.

### Numerics
- Contract: Decode biallelic unphased diploid probabilities into expected
  copies of `a1`. Missing samples remain missing in the dense missing mask.
  Genotype-stat filters use the decoded dosage values, matching VCF/PLINK2
  dosage semantics.
- Result/status semantics: Unsupported BGEN representation errors are
  user-actionable and should map to public `UnsupportedRepresentation` where
  possible. Malformed files or inconsistent counts should map to source or
  metadata errors consistent with existing readers.

### Egress
- Contract: Return the same Python shapes as other dense genotype reads:
  `np.ndarray` with samples on rows and variants on columns, optional sample
  and variant Polars DataFrames, and block iterators that concatenate to the
  full read.
- Output/exit-code mapping: Library API only; no CLI exit-code contract.

## Data Conversion and Copy Strategy
For each source format, record copy mode (`zero-copy`, `mmap`, `single-copy fallback`) and rationale.

| Source | Copy Mode | Rationale |
| --- | --- | --- |
| BGEN header/sample metadata | Single-copy into Rust-owned structs | Header and sample blocks are small and validated once. |
| BGEN variant metadata | Streaming parse into retained `VariantRecord`s | Allows metadata filters to drop variants before genotype decode. |
| BGEN compressed probability block | Single-copy compressed read, decompressed scratch buffer per variant | Compression requires materialization before bit unpacking. Scratch buffers should be reused across variants. |
| Dense output matrix | Rust-owned vectors converted through existing PyO3 dense conversion | Matches existing VCF/PLINK dense behavior. |

## Multi-Input Reconciliation Contract (Required When Multiple Tabular Sources Feed Numerics)
- Sources: BGEN file plus optional companion `.sample` file.
- Entity key(s) (for example subject/sample ID): sample ID.
- Join type and rationale: no join; sample IDs are an ordered source of truth
  for genotype rows.
- Duplicate-key policy: reject duplicate sample IDs.
- Missing-key policy: reject absent sample IDs.
- Row-order freeze policy: source sample order is preserved. Keep-lists filter
  membership only and do not reorder rows.
- Reconciliation accounting (matched/dropped/retained counts): reuse
  `DenseDiagnostics` sample counts and filtering diagnostics.
- Conversion boundary (where reconciled tabular data becomes arrays/PyTrees):
  Rust reader returns `DenseGenotypeMatrix`; Python converts to NumPy/Polars.

## Validation Strategy
- Boundary checks: validate header length, offset, magic bytes, flag bits,
  sample identifier availability, sample counts, variant count progression,
  block lengths, decompression lengths, and EOF/truncation behavior.
- Shape/range/domain checks: first slice requires Layout 2, `K = 2`,
  unphased flag, min/max ploidy both 2, per-sample ploidy 2 or missing, and
  probability bit depth in the spec-supported range. Decoded probabilities must
  produce dosage in `[0, 2]`.
- Multi-input alignment checks (key uniqueness, overlap expectations, deterministic row ordering):
  embedded sample IDs or companion sample IDs must have exactly `N` unique IDs.
  Requested sample filters use existing source-order selection semantics.
- Failure semantics: unsupported but valid BGEN features raise
  unsupported-representation errors; malformed BGEN data raises source/metadata
  parse errors with file context.

## Testing and Verification Strategy
- TDD scope: source resolution, sample ID requirements, header parsing,
  sample block parsing, companion `.sample` parsing, biallelic Layout 2 dosage
  decoding, missingness, metadata filters, genotype-stat filters, blocks, and
  public API error mapping.
- Regression strategy: build small deterministic BGEN fixtures in Rust tests,
  and use selected fixtures from the local `bgen` package only as references
  when their licensing and stability are appropriate. Tests should assert
  behavior rather than matching private parser internals.
- Verification commands:
  - `env CC=clang AR=ar cargo test --manifest-path rust/Cargo.toml -p genoio-io --test bgen_dense`
  - `.venv/bin/pytest tests/test_source_resolution.py tests/test_dense_read.py -q`
  - `make verify`

## Implementation Phases
<!-- START_PHASE_1 -->
### Phase 1: Source Resolution and API Contract
**Goal:** Introduce BGEN as a recognized source format without decoding genotype
probabilities yet.

**Components:** Extend `src/genoio/_source.py` with `SourceFormat.BGEN` and
`resolve_bgen`; expose the constructor in the public API; add PyO3 dispatch
branches that return clear unsupported errors for unimplemented BGEN read
paths.

**Dependencies:** Existing source resolution and public API patterns.

**Done when:** `.bgen` files and optional same-prefix `.sample` files resolve
correctly, missing/anonymous sample-ID scenarios are represented in tests, and
unsupported read paths fail with stable public errors.
<!-- END_PHASE_1 -->

<!-- START_PHASE_2 -->
### Phase 2: BGEN Header, Sample, and Metadata Reader
**Goal:** Parse enough BGEN metadata to return source capabilities, samples,
and variant rows without decoding dosage probabilities.

**Components:** Add `rust/genoio-io/src/bgen.rs` as the BGEN reader boundary;
parse BGEN v1.3 header fields and flags; parse embedded sample identifier
blocks; parse companion `.sample` IDs when embedded IDs are absent; stream
Layout 2 variant identifying data into `VariantRecord`s.

**Dependencies:** Phase 1 source resolution and existing `MetadataOutput`
contract.

**Done when:** BGEN metadata reads return sample and variant frames with the
core schema, reject unsupported layouts/compression/features clearly, and
validate sample count consistency.
<!-- END_PHASE_2 -->

<!-- START_PHASE_3 -->
### Phase 3: Dense Layout 2 Dosage Decode
**Goal:** Decode first-scope BGEN probability blocks into dense `a1` dosage
matrices.

**Components:** Implement compressed/uncompressed probability block handling,
zlib decompression and zstd if feasible with existing dependencies, Layout 2
header validation, packed-bit probability extraction, missing mask handling,
and `dosage_a1 = P(AB) + 2 * P(BB)`.

**Dependencies:** Phase 2 metadata parsing.

**Done when:** Dense reads produce sample-by-variant `float32` dosage matrices
with correct missing masks for biallelic unphased diploid Layout 2 BGEN
fixtures.
<!-- END_PHASE_3 -->

<!-- START_PHASE_4 -->
### Phase 4: Filters, Blocks, and Public Behavior
**Goal:** Make BGEN dosage reads behave like existing dense dosage sources
under filtering and block iteration.

**Components:** Reuse `VariantFilter` partial decisions; decode probability
blocks only when needed for retained variants or genotype-stat predicates;
compute dosage stats with `compute_dosage_variant_stats`; preserve
`VariantWindow` block semantics; map errors through Python public exception
types.

**Dependencies:** Phase 3 dense decode.

**Done when:** Metadata filters, genotype-stat filters, sample filters, empty
variant results, and block concatenation match full-read semantics.
<!-- END_PHASE_4 -->

<!-- START_PHASE_5 -->
### Phase 5: Documentation, Performance Baseline, and Hardening
**Goal:** Document the BGEN beta contract and establish baseline reliability.

**Components:** Update format support docs, API reading docs, source resolution
tests, public API tests, malformed-file tests, and benchmark scripts if useful.
Document unsupported BGEN features and sample-ID requirements explicitly.

**Dependencies:** Phase 4 public behavior.

**Done when:** Documentation states the exact supported BGEN subset, full
`make verify` passes, and any BGEN benchmark numbers are clearly scoped to the
first implementation.
<!-- END_PHASE_5 -->

## Simulation And Inference-Consistency Validation
- In scope: no
- Simulate entrypoint/signature: N/A
- Inputs: N/A
- Outputs: N/A
- Seed/RNG policy: N/A

### Assumption Alignment
| Inference Assumption | Simulation Rule | Mismatch Risk | Mitigation |
| --- | --- | --- | --- |
| N/A | N/A | N/A | N/A |

### Planned Validation Experiments
| Experiment ID | Type (recovery/SBC/PPC) | Success Criterion | Notes |
| --- | --- | --- | --- |
| N/A | N/A | N/A | N/A |

## Risks and Open Questions
| ID | Risk or Question | Severity | Mitigation or Next Step | Owner |
| --- | --- | --- | --- | --- |
| R1 | Packed-bit Layout 2 decoding is easy to get subtly wrong. | high | Use small hand-built fixtures across bit depths and compare against probability/dosage expectations. | implementer |
| R2 | BGEN files without embedded IDs may rely on external `.sample` conventions. | medium | Require explicit same-prefix `.sample` resolution and reject missing or count-mismatched IDs. | implementer |
| R3 | Zstandard support may require an additional Rust dependency or reuse of existing zstd support. | medium | Check current dependency tree during implementation; zlib support is required first, zstd is included only if straightforward. | implementer |
| R4 | Without index pushdown, large BGEN metadata filters may require full scans. | medium | Document beta limitation; preserve block reads and add `.bgi`/region pushdown later. | maintainer |
| R5 | Public `mac` metadata remains integer-only while dosage MAC can be fractional. | low | Filtering is internally exact; defer public fractional stats-column schema decision to a separate design if needed. | maintainer |

## Additional Considerations
- The local symlinked Cython `bgen` package should remain a reference, not a
  vendored dependency. Its behavior of generating integer sample IDs when no
  IDs exist is intentionally not adopted.
- BGEN `.sample` files may contain additional columns. The first reader needs
  only the first ID column and count validation, matching the minimal alignment
  requirement.
- The first implementation should prefer correctness and clear errors over
  broad format support. Adding phased, multiallelic, or indexed reads later is
  safer after the narrow Layout 2 path is well tested.

## Acceptance Criteria
### Source Resolution and Sample IDs
- `bgen-dosage.AC1.1`: `.bgen` paths resolve as `SourceFormat.BGEN` with a
  logical `"bgen"` member.
- `bgen-dosage.AC1.2`: same-prefix `.sample` files are discovered as optional
  sample ID companions.
- `bgen-dosage.AC1.3`: BGEN files without embedded sample IDs and without a
  companion `.sample` file are rejected.
- `bgen-dosage.AC1.4`: embedded or companion sample ID counts must match the
  BGEN header sample count.
- `bgen-dosage.AC1.5`: duplicate sample IDs are rejected.

### Supported and Unsupported Representations
- `bgen-dosage.AC2.1`: dense `kind="geno", dosage="dosage"` reads are supported
  for Layout 2 biallelic unphased diploid BGEN records.
- `bgen-dosage.AC2.2`: BGEN hardcall reads reject with
  `UnsupportedRepresentation`.
- `bgen-dosage.AC2.3`: sparse and haplotype BGEN reads reject with
  `UnsupportedRepresentation`.
- `bgen-dosage.AC2.4`: phased, multiallelic, variable-ploidy, unsupported
  compression, unsupported layout, and malformed probability blocks reject with
  clear errors.

### Dosage Semantics
- `bgen-dosage.AC3.1`: returned BGEN matrices use samples on rows and variants
  on columns.
- `bgen-dosage.AC3.2`: dosage values are expected copies of `a1`.
- `bgen-dosage.AC3.3`: biallelic unphased diploid probabilities decode as
  `P(AB) + 2 * P(BB)`.
- `bgen-dosage.AC3.4`: missing BGEN sample calls are preserved in the dense
  missing mask and obey Python missing-data policy.

### Filtering and Blocks
- `bgen-dosage.AC4.1`: sample filters retain source sample order.
- `bgen-dosage.AC4.2`: metadata-only variant filters can drop variants before
  probability decode.
- `bgen-dosage.AC4.3`: `maf`, `mac`, `missing_rate`, and `polymorphic` filters
  evaluate from decoded dosage values.
- `bgen-dosage.AC4.4`: filters retaining zero variants return normal empty
  matrix and metadata shapes.
- `bgen-dosage.AC4.5`: BGEN block reads concatenate to the full dense read for
  supported options.

### Documentation and Verification
- `bgen-dosage.AC5.1`: docs describe the supported BGEN subset and required
  sample ID behavior.
- `bgen-dosage.AC5.2`: docs list unsupported BGEN features explicitly.
- `bgen-dosage.AC5.3`: `make verify` passes after implementation.

## Glossary
- **BGEN**: Binary genotype format commonly used for imputed genotype
  probabilities.
- **Layout 2**: Modern BGEN variant block layout with explicit ploidy,
  missingness, phased/unphased flag, bit depth, and packed probabilities.
- **Sample identifier block**: Optional BGEN block after the header containing
  one sample ID per sample.
- **Companion `.sample` file**: External sample metadata file used when BGEN
  sample identifiers are not embedded.
- **`a1` dosage**: Expected count of the `a1` allele in a genotype, matching
  the existing `genoio` dosage convention.
- **Genotype-stat filter**: Variant filter requiring decoded genotype values,
  such as `maf`, `mac`, `missing_rate`, or `polymorphic`.

## Status Transition Log
| Date | From | To | Why | By |
| --- | --- | --- | --- | --- |
| 2026-06-04 | N/A | Draft | Plan created | |
| 2026-06-04 | Draft | Approved | User approved acceptance criteria and design scope. | Nicholas |
