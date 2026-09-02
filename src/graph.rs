use crate::dna::{normalize_dna, DnaError};
use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt::{Display, Formatter};

type EdgeMap = BTreeMap<String, u32>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphError {
    InvalidK(usize),
    InvalidDna(DnaError),
}

impl Display for GraphError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidK(k) => write!(formatter, "k must be at least 2, received {k}"),
            Self::InvalidDna(error) => Display::fmt(error, formatter),
        }
    }
}

impl Error for GraphError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidDna(error) => Some(error),
            Self::InvalidK(_) => None,
        }
    }
}

impl From<DnaError> for GraphError {
    fn from(error: DnaError) -> Self {
        Self::InvalidDna(error)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GraphStats {
    pub k: usize,
    pub nodes: usize,
    pub edges: usize,
    pub observations: u64,
}

/// A deterministic directed de Bruijn graph.
///
/// Nodes are `(k - 1)`-mers and edges are observed k-mers. Edge values retain
/// observation counts so low-support edges can be removed before assembly.
#[derive(Debug, Clone)]
pub struct DeBruijnGraph {
    k: usize,
    adjacency: BTreeMap<String, EdgeMap>,
}

impl DeBruijnGraph {
    pub fn new(k: usize) -> Result<Self, GraphError> {
        if k < 2 {
            return Err(GraphError::InvalidK(k));
        }
        Ok(Self {
            k,
            adjacency: BTreeMap::new(),
        })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    /// Add all unambiguous k-mers in a read. Windows containing `N` are
    /// skipped, preventing an ambiguous symbol from creating graph edges.
    pub fn add_sequence(&mut self, sequence: &str) -> Result<usize, GraphError> {
        let sequence = normalize_dna(sequence)?;
        if sequence.len() < self.k {
            return Ok(0);
        }

        let mut added = 0;
        for window in sequence.as_bytes().windows(self.k) {
            if window.contains(&b'N') {
                continue;
            }
            let prefix = std::str::from_utf8(&window[..self.k - 1])
                .expect("validated DNA is ASCII")
                .to_owned();
            let suffix = std::str::from_utf8(&window[1..])
                .expect("validated DNA is ASCII")
                .to_owned();
            let count = self
                .adjacency
                .entry(prefix)
                .or_default()
                .entry(suffix.clone())
                .or_default();
            *count = count.saturating_add(1);
            self.adjacency.entry(suffix).or_default();
            added += 1;
        }
        Ok(added)
    }

    /// Remove edges observed fewer than `minimum_coverage` times, followed by
    /// nodes that no longer participate in any edge.
    pub fn prune(&mut self, minimum_coverage: u32) {
        let threshold = minimum_coverage.max(1);
        for outgoing in self.adjacency.values_mut() {
            outgoing.retain(|_, coverage| *coverage >= threshold);
        }

        let referenced: BTreeSet<String> = self
            .adjacency
            .values()
            .flat_map(|outgoing| outgoing.keys().cloned())
            .collect();
        self.adjacency
            .retain(|node, outgoing| !outgoing.is_empty() || referenced.contains(node));
    }

    pub fn stats(&self) -> GraphStats {
        GraphStats {
            k: self.k,
            nodes: self.adjacency.len(),
            edges: self.adjacency.values().map(BTreeMap::len).sum(),
            observations: self
                .adjacency
                .values()
                .flat_map(BTreeMap::values)
                .map(|count| u64::from(*count))
                .sum(),
        }
    }

    pub(crate) fn adjacency(&self) -> &BTreeMap<String, EdgeMap> {
        &self.adjacency
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_repeated_edges() {
        let mut graph = DeBruijnGraph::new(3).unwrap();
        graph.add_sequence("ACGACG").unwrap();
        let stats = graph.stats();
        assert_eq!(stats.observations, 4);
        assert_eq!(graph.adjacency["AC"]["CG"], 2);
    }

    #[test]
    fn skips_ambiguous_windows() {
        let mut graph = DeBruijnGraph::new(3).unwrap();
        assert_eq!(graph.add_sequence("ACNTA").unwrap(), 0);
        assert_eq!(graph.stats().edges, 0);
    }

    #[test]
    fn prunes_low_support_edges() {
        let mut graph = DeBruijnGraph::new(3).unwrap();
        graph.add_sequence("ACGT").unwrap();
        graph.add_sequence("ACG").unwrap();
        graph.prune(2);
        assert_eq!(graph.stats().edges, 1);
        assert_eq!(graph.stats().nodes, 2);
    }
}
