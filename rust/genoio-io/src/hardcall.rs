// pattern: Functional Core

//! Shared hard-call storage, filtering, and dense batching.
//!
//! VCF and PLINK readers import this facade; implementation details stay split
//! by responsibility beneath it.

mod batch;
mod counts;
mod packed;

pub(crate) use batch::{flush_hardcall_batch_into_sample_major, HardcallBatch};
pub(crate) use counts::{evaluate_hardcall_counts_filter, HardcallCounts};
pub(crate) use packed::{evaluate_packed_hardcall_filter, PackedHardcalls};

pub(crate) const HARDCALL_BATCH_SIZE: usize = 64;
