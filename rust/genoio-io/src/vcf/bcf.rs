// pattern: Mixed
//
// The facade exposes BCF reader operations while implementation modules
// separate record decoding from stateful source iteration.

//! BCF reader facade.

mod decode;
mod haplotype;
mod record;
mod source;

pub(super) use source::{
    read_dense_windowed, read_dosage_dense_windowed, read_haplotypes_dense_windowed,
    read_haplotypes_sparse_windowed, read_metadata, read_sparse_windowed,
};
