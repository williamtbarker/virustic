# Virustic

Virustic is a compact, deterministic de Bruijn-graph assembler written in
Rust. It turns FASTA or conventional four-line FASTQ reads into unitigs while
retaining edge coverage statistics.

The project is intentionally dependency-free. Its purpose is to make the core
data structures and graph traversal easy to inspect, test, and extend—not to
compete with production assemblers.

## What it demonstrates

- validated FASTA and FASTQ parsing
- directed de Bruijn graph construction
- coverage-aware edge pruning
- deterministic maximal non-branching path extraction
- explicit handling of isolated cycles
- library and command-line interfaces
- unit, integration, formatting, and lint checks in CI

## Quick start

```bash
cargo run --release -- \
  --input examples/reads.fasta \
  --output contigs.fasta \
  --kmer-size 7 \
  --min-coverage 2
```

Inspect the result:

```bash
cat contigs.fasta
```

Run the quality gates:

```bash
./scripts/verify.sh
```

Mac-specific installation and publishing steps are in
[`MACBOOK_SETUP.md`](MACBOOK_SETUP.md).

## Algorithm

For each unambiguous k-mer, Virustic adds an edge from its `(k - 1)`-base
prefix to its `(k - 1)`-base suffix and increments the edge's observation
count. Edges below the requested coverage are discarded. The assembler then
walks each maximal non-branching path once; disconnected one-in/one-out cycles
are emitted in a second pass.

Ordered maps and stable output sorting make results reproducible. A FASTA
header records length, edge count, minimum and mean coverage, and whether a
path was detected as circular.

## Scope and limitations

This is an educational and portfolio-scale implementation. It loads parsed
records into memory, accepts uncompressed input, does not correct sequencing
errors, and emits unitigs rather than attempting repeat resolution or genome
finishing. It must not be used for clinical decisions.

For gzipped data, decompress to a stream or file first:

```bash
gzip -dc reads.fastq.gz > reads.fastq
```

## License

MIT
