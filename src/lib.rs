//! Core library for the Virustic educational sequence assembler.
//!
//! The implementation deliberately uses only Rust's standard library. That
//! keeps the graph and parsing algorithms visible, auditable, and easy to
//! build in constrained environments.

pub mod assembly;
pub mod dna;
pub mod graph;
pub mod io;

pub use assembly::{assemble_unitigs, write_contigs, Contig};
pub use graph::{DeBruijnGraph, GraphError, GraphStats};
pub use io::{read_sequences, SequenceRecord};
