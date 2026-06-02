# PLINK2 Parser Strategy

genoio detects complete PLINK2 `.pgen/.pvar/.psam` file sets, but full decode is deferred for the current release. Public read calls raise `UnsupportedFormatError` with a message pointing back to this strategy.

## Options

### From-scratch Rust `.pgen` parser

- Licensing: The parser could live under genoio's own license if written without copying PLINK2 code, but the implementation would need careful provenance review.
- Packaging: A Rust-native parser fits the existing PyO3 wheel pipeline and avoids shipping a second native dependency.
- Platform wheels: This option is wheel-friendly once implemented because it compiles with the existing Rust workspace.
- Phased and dosage support: `.pgen` supports hardcalls, phase, dosage, and compact encodings. Covering those correctly would require a substantial conformance suite.
- Multiallelic support: PLINK2 supports richer variant representations than the current biallelic genotype path. A partial parser risks silently narrowing data semantics.
- Build reproducibility: Reproducibility is strong because the parser would be compiled from the repo, but correctness would depend on maintaining bit-level compatibility with PLINK2.

This path has the cleanest packaging story and the highest long-term maintenance cost.

### Bind to pgenlib

- Licensing: PLINK2/pgenlib licensing and redistribution terms must be reviewed before bundling or linking it into genoio wheels.
- Packaging: pgenlib introduces a native dependency outside the current Rust workspace. The binding path could be FFI from Rust or a Python-side adapter, but either adds build and ABI surface.
- Platform wheels: macOS, Linux, and Windows wheels would need reliable pgenlib builds. CI must prove stable linking for all supported Python versions and architectures.
- Phased and dosage support: pgenlib is the most likely path to correct coverage of PLINK2 phase and dosage semantics because it is maintained with the format.
- Multiallelic support: pgenlib is also the safest option for current PLINK2 multiallelic behavior, subject to binding API coverage.
- Build reproducibility: Reproducible builds depend on vendoring or pinning pgenlib source and documenting the toolchain used for wheels.

This path has the best correctness outlook and the most packaging risk.

### Metadata-only/deferred support

- Licensing: Detection and structured rejection do not require bundling PLINK2 code.
- Packaging: No new dependency is introduced, and existing wheels remain unchanged.
- Platform wheels: No additional wheel matrix work is required.
- Phased and dosage support: Reads are rejected, so genoio does not misrepresent unsupported phase or dosage encodings.
- Multiallelic support: Reads are rejected, avoiding accidental loss of multiallelic information.
- Build reproducibility: Existing reproducibility is preserved because no new native dependency is added.

This path makes the support boundary explicit while leaving room for a later parser spike.

## Decision

Full PLINK2 decode is deferred. genoio will detect complete `.pgen/.pvar/.psam` sources and raise a structured deferred-decode error for public reads until a follow-up spike selects either a pinned pgenlib binding or a tested Rust parser. The next implementation decision should prefer pgenlib if licensing and wheel reproducibility are acceptable; otherwise, a Rust parser needs a dedicated conformance plan before any partial decode ships.
