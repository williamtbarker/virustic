use crate::dna::{normalize_dna, DnaError};
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::io::{self, BufRead, Write};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    pub id: String,
    pub sequence: String,
}

#[derive(Debug)]
pub enum ParseError {
    Io(io::Error),
    Format { line: usize, message: String },
    InvalidDna { id: String, source: DnaError },
}

impl Display for ParseError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "I/O error: {error}"),
            Self::Format { line, message } => write!(formatter, "line {line}: {message}"),
            Self::InvalidDna { id, source } => {
                write!(formatter, "invalid sequence in record {id:?}: {source}")
            }
        }
    }
}

impl Error for ParseError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::InvalidDna { source, .. } => Some(source),
            Self::Format { .. } => None,
        }
    }
}

impl From<io::Error> for ParseError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

/// Read FASTA or four-line FASTQ records from a buffered reader.
///
/// The format is detected from the first non-empty line. FASTA sequences may
/// span multiple lines. FASTQ sequence and quality data must each occupy one
/// line, which is the overwhelmingly common interchange form for short reads.
pub fn read_sequences<R: BufRead>(reader: R) -> Result<Vec<SequenceRecord>, ParseError> {
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;
    let first = lines
        .iter()
        .position(|line| !line.trim().is_empty())
        .ok_or_else(|| ParseError::Format {
            line: 1,
            message: "input is empty".to_owned(),
        })?;

    match lines[first].as_bytes().first() {
        Some(b'>') => parse_fasta(&lines, first),
        Some(b'@') => parse_fastq(&lines, first),
        _ => Err(ParseError::Format {
            line: first + 1,
            message: "expected a FASTA '>' or FASTQ '@' header".to_owned(),
        }),
    }
}

fn checked_record(id: String, sequence: String) -> Result<SequenceRecord, ParseError> {
    if sequence.is_empty() {
        return Err(ParseError::Format {
            line: 1,
            message: format!("record {id:?} has an empty sequence"),
        });
    }
    let sequence = normalize_dna(&sequence).map_err(|source| ParseError::InvalidDna {
        id: id.clone(),
        source,
    })?;
    Ok(SequenceRecord { id, sequence })
}

fn parse_fasta(lines: &[String], start: usize) -> Result<Vec<SequenceRecord>, ParseError> {
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut sequence = String::new();

    for (index, raw_line) in lines.iter().enumerate().skip(start) {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                records.push(checked_record(id, std::mem::take(&mut sequence))?);
            }
            let id = header.trim();
            if id.is_empty() {
                return Err(ParseError::Format {
                    line: index + 1,
                    message: "FASTA header has no identifier".to_owned(),
                });
            }
            current_id = Some(id.to_owned());
        } else if current_id.is_some() {
            sequence.push_str(line);
        } else {
            return Err(ParseError::Format {
                line: index + 1,
                message: "sequence appears before the first FASTA header".to_owned(),
            });
        }
    }

    if let Some(id) = current_id {
        records.push(checked_record(id, sequence)?);
    }
    Ok(records)
}

fn parse_fastq(lines: &[String], start: usize) -> Result<Vec<SequenceRecord>, ParseError> {
    let mut records = Vec::new();
    let mut index = start;

    while index < lines.len() {
        if lines[index].trim().is_empty() {
            index += 1;
            continue;
        }
        if index + 3 >= lines.len() {
            return Err(ParseError::Format {
                line: index + 1,
                message: "incomplete four-line FASTQ record".to_owned(),
            });
        }

        let header = lines[index].trim();
        let sequence = lines[index + 1].trim();
        let separator = lines[index + 2].trim();
        let quality = lines[index + 3].trim();
        let id = header.strip_prefix('@').ok_or_else(|| ParseError::Format {
            line: index + 1,
            message: "FASTQ header must begin with '@'".to_owned(),
        })?;
        if id.is_empty() {
            return Err(ParseError::Format {
                line: index + 1,
                message: "FASTQ header has no identifier".to_owned(),
            });
        }
        if !separator.starts_with('+') {
            return Err(ParseError::Format {
                line: index + 3,
                message: "FASTQ separator must begin with '+'".to_owned(),
            });
        }
        if sequence.len() != quality.len() {
            return Err(ParseError::Format {
                line: index + 4,
                message: format!(
                    "quality length {} differs from sequence length {}",
                    quality.len(),
                    sequence.len()
                ),
            });
        }
        records.push(checked_record(id.to_owned(), sequence.to_owned())?);
        index += 4;
    }
    Ok(records)
}

/// Write sequence records as wrapped FASTA.
pub fn write_fasta<W: Write>(mut writer: W, records: &[SequenceRecord]) -> io::Result<()> {
    for record in records {
        writeln!(writer, ">{id}", id = record.id)?;
        for chunk in record.sequence.as_bytes().chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn reads_multiline_fasta() {
        let input = b">read-1\nacgt\nac\n>read-2 description\nTTNN\n";
        let records = read_sequences(Cursor::new(input)).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].sequence, "ACGTAC");
        assert_eq!(records[1].id, "read-2 description");
    }

    #[test]
    fn reads_fastq_and_validates_quality_length() {
        let input = b"@read-1\nACGT\n+\n!!!!\n";
        let records = read_sequences(Cursor::new(input)).unwrap();
        assert_eq!(records[0].sequence, "ACGT");

        let bad = b"@read-1\nACGT\n+\n!!!\n";
        assert!(read_sequences(Cursor::new(bad)).is_err());
    }

    #[test]
    fn writes_wrapped_fasta() {
        let records = vec![SequenceRecord {
            id: "example".to_owned(),
            sequence: "A".repeat(81),
        }];
        let mut output = Vec::new();
        write_fasta(&mut output, &records).unwrap();
        let text = String::from_utf8(output).unwrap();
        assert!(text.contains(&format!("{}\nA\n", "A".repeat(80))));
    }
}
