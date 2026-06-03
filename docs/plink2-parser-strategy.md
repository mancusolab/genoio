# PLINK2 Parser Strategy

genoio detects complete PLINK2 `.pgen/.pvar/.psam` file sets and supports Rust-native decode paths for biallelic hard-call PGEN files.

Supported today:

- `.psam` sample metadata with a header line containing `IID`, plus optional `FID`, `PAT`/`MAT`, `SEX`, and phenotype columns.
- `.pvar` variant metadata with explicit `#CHROM ... POS ID REF ALT`-style headers, plus PLINK1-compatible no-header five- or six-column fallback.
- `.pgen` storage mode `0x02`, the fixed-width unphased hard-call format where all records are PGEN category-code type `0`.
- `.pgen` storage mode `0x10`, the variable-width format, when records are biallelic hardcalls using main-track compression types `0`, `1`, `2`, `3`, `4`, `6`, or `7`.
- Hardcall phase auxiliary tracks are allowed for genotype reads and ignored by `kind="geno"`.
- Public `kind="geno"` dense reads, sparse reads with the existing missing-value restriction, metadata returns, filters, and blocks.

Unsupported PLINK2 encodings are rejected instead of approximated:

- Variable-width modes with ignorable extensions or external indexes: `0x11`, `0x20`, and `0x21`.
- Fixed-width dosage modes `0x03` and `0x04`.
- Multiallelic patch sets, dosage tracks, phased-dosage tracks, and external PGEN indexes.
- Compressed `.pgen.zst`/`.pvar.zst` source members are not parsed directly; use decompressed `.pgen`/`.pvar` members.

## Options

### From-scratch Rust `.pgen` parser

- Licensing: The parser could live under genoio's own license if written without copying PLINK2 code, but the implementation would need careful provenance review.
- Packaging: A Rust-native parser fits the existing PyO3 wheel pipeline and avoids shipping a second native dependency.
- Platform wheels: This option is wheel-friendly once implemented because it compiles with the existing Rust workspace.
- Phased and dosage support: `.pgen` supports hardcalls, phase, dosage, and compact encodings. Covering those correctly requires a substantial conformance suite.
- Multiallelic support: PLINK2 supports richer variant representations than the current biallelic genotype path. Partial parsing must reject unsupported semantics.
- Build reproducibility: Reproducibility is strong because the parser would be compiled from the repo, but correctness would depend on maintaining bit-level compatibility with PLINK2.

This path has the cleanest packaging story and the highest long-term maintenance cost. genoio now uses this path for fixed-width and selected variable-width hard-call subsets.

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
- Phased and dosage support: Unsupported richer reads are rejected, so genoio does not misrepresent unsupported phase or dosage encodings.
- Multiallelic support: Unsupported richer reads are rejected, avoiding accidental loss of multiallelic information.
- Build reproducibility: Existing reproducibility is preserved because no new native dependency is added.

This remains the fallback for encodings outside the supported hard-call subset.

## Decision

Implement a Rust-native hard-call parser for `.pgen` mode `0x02` and the biallelic hard-call subset of mode `0x10`, guided by the official PLINK2 PGEN specification and pgenlib Python API behavior. Keep dosage, multiallelic patches, external indexes, and compressed source members explicitly unsupported until follow-up spikes extend the Rust parser with conformance tests or introduce a pinned pgenlib binding with reviewed licensing and reproducible wheel builds.
