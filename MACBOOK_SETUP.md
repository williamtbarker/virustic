# Validate and publish from a MacBook

## 1. Install the Rust toolchain

If `cargo --version` does not work, install Rust with the official installer:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
rustup component add clippy rustfmt
```

## 2. Unpack and enter the repository

```bash
cd ~/Documents/GPT_Hist_Review
unzip Virustic_GitHub_Candidate_v0.1.zip
cd virustic
```

Replace `YOUR-USERNAME` in `Cargo.toml`, then run the complete local gate:

```bash
./scripts/verify.sh
```

The script checks formatting, runs Clippy with warnings treated as errors,
runs all tests, and executes the included example in a temporary directory.

If formatting is the only failed check, run `cargo fmt --all`, inspect the
changes, and rerun the script.

## 3. Create the GitHub repository

Using the GitHub CLI:

```bash
git init
git add .
git commit -m "Initial release: deterministic de Bruijn graph assembler"
git branch -M main
gh repo create virustic --public --source=. --remote=origin --push
```

Or create an empty public repository named `virustic` on GitHub and follow the
push instructions it displays.

## 4. Let CI verify the public commit

Open the repository's **Actions** tab and wait for the `CI` workflow to turn
green. If it does, create the first release:

```bash
git tag -a v0.1.0 -m "Virustic 0.1.0"
git push origin v0.1.0
```

Do not commit real sequencing data, credentials, private paths, or recovered
chat transcripts. The included synthetic example is safe to publish.

