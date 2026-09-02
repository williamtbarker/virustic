use std::io::Cursor;
use virustic::{assemble_unitigs, read_sequences, DeBruijnGraph};

#[test]
fn reads_build_a_reproducible_assembly() {
    let input = b">r1\nAACGTT\n>r2\nAACGTT\n";
    let reads = read_sequences(Cursor::new(input)).unwrap();
    let mut graph = DeBruijnGraph::new(4).unwrap();
    for read in reads {
        graph.add_sequence(&read.sequence).unwrap();
    }
    graph.prune(2);
    let contigs = assemble_unitigs(&graph, 0);
    assert_eq!(contigs.len(), 1);
    assert_eq!(contigs[0].sequence, "AACGTT");
    assert_eq!(contigs[0].minimum_coverage, 2);
}
