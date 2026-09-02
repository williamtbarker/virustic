# Contributing

Small, test-backed changes are welcome. Before opening a pull request, run:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
```

Please include a focused test for algorithm or parser changes and avoid adding
dependencies when the standard library keeps the implementation clear.

