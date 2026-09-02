use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use virustic::{assemble_unitigs, read_sequences, write_contigs, DeBruijnGraph};

const HELP: &str = "Virustic 0.1.0
A small, deterministic de Bruijn graph assembler.

USAGE:
    virustic --input READS --output CONTIGS [OPTIONS]

OPTIONS:
    -i, --input PATH          Input FASTA or four-line FASTQ
    -o, --output PATH         Output FASTA
    -k, --kmer-size INTEGER   K-mer size, at least 2 [default: 31]
        --min-coverage N      Remove edges observed fewer than N times [default: 2]
        --min-length N        Do not emit contigs shorter than N bases [default: 0]
    -h, --help                Print this help
";

#[derive(Debug, PartialEq, Eq)]
struct Config {
    input: PathBuf,
    output: PathBuf,
    kmer_size: usize,
    minimum_coverage: u32,
    minimum_length: usize,
}

impl Config {
    fn parse<I>(arguments: I) -> Result<Self, String>
    where
        I: IntoIterator<Item = String>,
    {
        let mut input = None;
        let mut output = None;
        let mut kmer_size = 31;
        let mut minimum_coverage = 2;
        let mut minimum_length = 0;
        let mut arguments = arguments.into_iter();

        while let Some(argument) = arguments.next() {
            let mut value = || {
                arguments
                    .next()
                    .ok_or_else(|| format!("missing value after {argument}"))
            };
            match argument.as_str() {
                "-i" | "--input" => input = Some(PathBuf::from(value()?)),
                "-o" | "--output" => output = Some(PathBuf::from(value()?)),
                "-k" | "--kmer-size" => {
                    kmer_size = value()?
                        .parse()
                        .map_err(|_| "k-mer size must be an integer".to_owned())?;
                }
                "--min-coverage" => {
                    minimum_coverage = value()?
                        .parse()
                        .map_err(|_| "minimum coverage must be an integer".to_owned())?;
                }
                "--min-length" => {
                    minimum_length = value()?
                        .parse()
                        .map_err(|_| "minimum length must be an integer".to_owned())?;
                }
                unknown => return Err(format!("unknown argument {unknown:?}")),
            }
        }

        if kmer_size < 2 {
            return Err("k-mer size must be at least 2".to_owned());
        }
        if minimum_coverage == 0 {
            return Err("minimum coverage must be at least 1".to_owned());
        }
        Ok(Self {
            input: input.ok_or_else(|| "--input is required".to_owned())?,
            output: output.ok_or_else(|| "--output is required".to_owned())?,
            kmer_size,
            minimum_coverage,
            minimum_length,
        })
    }
}

fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let input = BufReader::new(File::open(&config.input)?);
    let records = read_sequences(input)?;
    let mut graph = DeBruijnGraph::new(config.kmer_size)?;
    let mut accepted_windows = 0usize;
    for record in &records {
        accepted_windows += graph.add_sequence(&record.sequence)?;
    }
    if accepted_windows == 0 {
        return Err(format!(
            "no unambiguous {}-mers were found; choose a smaller k or inspect the reads",
            config.kmer_size
        )
        .into());
    }

    graph.prune(config.minimum_coverage);
    let stats = graph.stats();
    let contigs = assemble_unitigs(&graph, config.minimum_length);
    let output = BufWriter::new(File::create(&config.output)?);
    write_contigs(output, &contigs)?;
    eprintln!(
        "processed {} reads; retained {} nodes and {} edges; wrote {} contigs",
        records.len(),
        stats.nodes,
        stats.edges,
        contigs.len()
    );
    Ok(())
}

fn main() {
    let arguments: Vec<String> = env::args().skip(1).collect();
    if arguments
        .iter()
        .any(|argument| argument == "-h" || argument == "--help")
    {
        print!("{HELP}");
        return;
    }
    let config = match Config::parse(arguments) {
        Ok(config) => config,
        Err(message) => {
            eprintln!("error: {message}\n\n{HELP}");
            std::process::exit(2);
        }
    };
    if let Err(error) = run(config) {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_cli_configuration() {
        let config = Config::parse(
            [
                "--input",
                "reads.fasta",
                "--output",
                "contigs.fasta",
                "--kmer-size",
                "21",
                "--min-coverage",
                "3",
            ]
            .into_iter()
            .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(config.kmer_size, 21);
        assert_eq!(config.minimum_coverage, 3);
    }

    #[test]
    fn rejects_too_small_k() {
        let error = Config::parse(
            ["-i", "reads.fasta", "-o", "out.fasta", "-k", "1"]
                .into_iter()
                .map(str::to_owned),
        )
        .unwrap_err();
        assert!(error.contains("at least 2"));
    }
}
