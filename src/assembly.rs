use crate::graph::DeBruijnGraph;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};

#[derive(Debug, Clone, PartialEq)]
pub struct Contig {
    pub sequence: String,
    pub edge_count: usize,
    pub minimum_coverage: u32,
    pub mean_coverage: f64,
    pub circular: bool,
}

type Edge = (String, String);

/// Assemble maximal non-branching paths (unitigs) from a de Bruijn graph.
///
/// Every retained edge is emitted exactly once. Results are sorted by length
/// and sequence to keep output reproducible across machines and runs.
pub fn assemble_unitigs(graph: &DeBruijnGraph, minimum_length: usize) -> Vec<Contig> {
    let adjacency = graph.adjacency();
    let mut indegree: BTreeMap<&str, usize> =
        adjacency.keys().map(|node| (node.as_str(), 0)).collect();
    for outgoing in adjacency.values() {
        for target in outgoing.keys() {
            *indegree.entry(target.as_str()).or_default() += 1;
        }
    }

    let outdegree = |node: &str| adjacency.get(node).map_or(0, BTreeMap::len);
    let mut used: BTreeSet<Edge> = BTreeSet::new();
    let mut contigs = Vec::new();

    // Paths beginning at a source, sink-adjacent node, or branch.
    for (source, outgoing) in adjacency {
        let source_in = *indegree.get(source.as_str()).unwrap_or(&0);
        let source_out = outgoing.len();
        if source_out == 0 || (source_in == 1 && source_out == 1) {
            continue;
        }
        for (target, coverage) in outgoing {
            if used.contains(&(source.clone(), target.clone())) {
                continue;
            }
            let (nodes, coverages) = extend_path(
                source, target, *coverage, adjacency, &indegree, &outdegree, &mut used, None,
            );
            push_contig(&mut contigs, nodes, coverages, false, minimum_length);
        }
    }

    // Any unused edges belong to isolated one-in/one-out cycles.
    for (source, outgoing) in adjacency {
        for (target, coverage) in outgoing {
            if used.contains(&(source.clone(), target.clone())) {
                continue;
            }
            let (nodes, coverages) = extend_path(
                source,
                target,
                *coverage,
                adjacency,
                &indegree,
                &outdegree,
                &mut used,
                Some(source),
            );
            push_contig(&mut contigs, nodes, coverages, true, minimum_length);
        }
    }

    contigs.sort_by(|left, right| {
        right
            .sequence
            .len()
            .cmp(&left.sequence.len())
            .then_with(|| left.sequence.cmp(&right.sequence))
    });
    contigs
}

#[allow(clippy::too_many_arguments)]
fn extend_path<'a, F>(
    source: &'a str,
    target: &'a str,
    first_coverage: u32,
    adjacency: &'a BTreeMap<String, BTreeMap<String, u32>>,
    indegree: &BTreeMap<&'a str, usize>,
    outdegree: &F,
    used: &mut BTreeSet<Edge>,
    cycle_origin: Option<&str>,
) -> (Vec<String>, Vec<u32>)
where
    F: Fn(&str) -> usize,
{
    let mut nodes = vec![source.to_owned(), target.to_owned()];
    let mut coverages = vec![first_coverage];
    used.insert((source.to_owned(), target.to_owned()));
    let mut current = target;

    loop {
        if cycle_origin == Some(current) {
            break;
        }
        if indegree.get(current).copied().unwrap_or(0) != 1 || outdegree(current) != 1 {
            break;
        }
        let Some((next, coverage)) = adjacency.get(current).and_then(|edges| edges.iter().next())
        else {
            break;
        };
        let edge = (current.to_owned(), next.clone());
        if used.contains(&edge) {
            break;
        }
        used.insert(edge);
        nodes.push(next.clone());
        coverages.push(*coverage);
        current = next;
    }
    (nodes, coverages)
}

fn push_contig(
    contigs: &mut Vec<Contig>,
    nodes: Vec<String>,
    coverages: Vec<u32>,
    circular: bool,
    minimum_length: usize,
) {
    let Some(first) = nodes.first() else {
        return;
    };
    let mut sequence = first.clone();
    for node in nodes.iter().skip(1) {
        if let Some(base) = node.as_bytes().last() {
            sequence.push(char::from(*base));
        }
    }
    if sequence.len() < minimum_length || coverages.is_empty() {
        return;
    }
    let minimum_coverage = coverages.iter().copied().min().unwrap_or(0);
    let total: u64 = coverages.iter().map(|value| u64::from(*value)).sum();
    let mean_coverage = total as f64 / coverages.len() as f64;
    contigs.push(Contig {
        sequence,
        edge_count: coverages.len(),
        minimum_coverage,
        mean_coverage,
        circular,
    });
}

pub fn write_contigs<W: Write>(mut writer: W, contigs: &[Contig]) -> io::Result<()> {
    for (index, contig) in contigs.iter().enumerate() {
        writeln!(
            writer,
            ">contig_{:04} length={} edges={} min_cov={} mean_cov={:.2} circular={}",
            index + 1,
            contig.sequence.len(),
            contig.edge_count,
            contig.minimum_coverage,
            contig.mean_coverage,
            contig.circular
        )?;
        for chunk in contig.sequence.as_bytes().chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assembles_a_linear_path() {
        let mut graph = DeBruijnGraph::new(3).unwrap();
        graph.add_sequence("AACGT").unwrap();
        let contigs = assemble_unitigs(&graph, 0);
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].sequence, "AACGT");
        assert_eq!(contigs[0].edge_count, 3);
        assert!(!contigs[0].circular);
    }

    #[test]
    fn splits_paths_at_a_branch() {
        let mut graph = DeBruijnGraph::new(3).unwrap();
        graph.add_sequence("AACG").unwrap();
        graph.add_sequence("AACT").unwrap();
        let sequences: BTreeSet<String> = assemble_unitigs(&graph, 0)
            .into_iter()
            .map(|contig| contig.sequence)
            .collect();
        assert!(sequences.contains("AAC"));
        assert!(sequences.contains("ACG"));
        assert!(sequences.contains("ACT"));
    }

    #[test]
    fn emits_isolated_cycles() {
        let mut graph = DeBruijnGraph::new(3).unwrap();
        graph.add_sequence("ATAT").unwrap();
        let contigs = assemble_unitigs(&graph, 0);
        assert_eq!(contigs.len(), 1);
        assert!(contigs[0].circular);
        assert_eq!(contigs[0].edge_count, 2);
    }

    #[test]
    fn filters_short_contigs() {
        let mut graph = DeBruijnGraph::new(3).unwrap();
        graph.add_sequence("AACGT").unwrap();
        assert!(assemble_unitigs(&graph, 6).is_empty());
    }
}
