// pattern: Mixed
//
// Imperative payload I/O delegates to pure generic and fast dosage decoders.

//! BGEN Layout 2 probability payload reading and dosage decoding.

mod fast;
mod generic;
mod payload;

pub(super) use fast::SampleMajorSlotMut;
pub(super) use generic::{
    decode_buffered_dosage_values, decode_buffered_haplotype_values,
    try_decode_buffered_dosage_values_into_sample_major_slot,
    try_decode_buffered_dosage_values_with_counts, DosageDecodeBuffers, HaplotypeDecodeBuffers,
};
pub(super) use payload::{
    read_layout2_probability_payload_into, skip_layout2_probability_payload_raw,
    ProbabilityPayloadBuffers,
};
