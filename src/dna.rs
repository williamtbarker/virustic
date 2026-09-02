use std::error::Error;
use std::fmt::{Display, Formatter};

/// Error returned when a sequence contains a symbol outside the supported
/// DNA alphabet (`A`, `C`, `G`, `T`, and `N`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaError {
    pub position: usize,
    pub symbol: char,
}

impl Display for DnaError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "unsupported DNA symbol {:?} at zero-based position {}",
            self.symbol, self.position
        )
    }
}

impl Error for DnaError {}

/// Normalize a DNA sequence to uppercase and validate its alphabet.
pub fn normalize_dna(sequence: &str) -> Result<String, DnaError> {
    let mut normalized = String::with_capacity(sequence.len());
    for (position, symbol) in sequence.chars().enumerate() {
        let upper = symbol.to_ascii_uppercase();
        if matches!(upper, 'A' | 'C' | 'G' | 'T' | 'N') {
            normalized.push(upper);
        } else {
            return Err(DnaError { position, symbol });
        }
    }
    Ok(normalized)
}

/// Return the reverse complement of a DNA sequence.
pub fn reverse_complement(sequence: &str) -> Result<String, DnaError> {
    let normalized = normalize_dna(sequence)?;
    Ok(normalized
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'N' => 'N',
            _ => unreachable!("normalize_dna validates the alphabet"),
        })
        .collect())
}

/// Select a strand-independent representation of a DNA k-mer.
pub fn canonical_kmer(kmer: &str) -> Result<String, DnaError> {
    let forward = normalize_dna(kmer)?;
    let reverse = reverse_complement(&forward)?;
    Ok(if forward <= reverse { forward } else { reverse })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_lowercase_dna() {
        assert_eq!(normalize_dna("acgtn").unwrap(), "ACGTN");
    }

    #[test]
    fn reports_invalid_symbol_position() {
        let error = normalize_dna("ACXG").unwrap_err();
        assert_eq!(error.position, 2);
        assert_eq!(error.symbol, 'X');
    }

    #[test]
    fn computes_reverse_complement() {
        assert_eq!(reverse_complement("AACGTN").unwrap(), "NACGTT");
    }

    #[test]
    fn canonical_form_is_strand_independent() {
        assert_eq!(
            canonical_kmer("AACG").unwrap(),
            canonical_kmer("CGTT").unwrap()
        );
    }
}
