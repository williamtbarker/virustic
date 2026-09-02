#!/usr/bin/env bash
set -euo pipefail

if ! command -v cargo >/dev/null 2>&1; then
  echo "Rust is not installed. Install it from https://rustup.rs and rerun this script." >&2
  exit 1
fi

cd "$(dirname "$0")/.."

cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features

temporary_directory="$(mktemp -d)"
trap 'rm -rf "$temporary_directory"' EXIT

cargo run --quiet -- \
  --input examples/reads.fasta \
  --output "$temporary_directory/contigs.fasta" \
  --kmer-size 7 \
  --min-coverage 2

if ! grep -q '^>contig_' "$temporary_directory/contigs.fasta"; then
  echo "The example run produced no contigs." >&2
  exit 1
fi

echo "All checks passed, including the example assembly."

