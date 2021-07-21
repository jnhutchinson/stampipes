# Rust code

This directory houses the rust source code and tools that are written specifically for our pipeline.

# Notes for users

For now, this is a source-only release. After compilation, binaries will be
located in `target/release/` or `target/x86_64-unknown-linux-musl/`, depending
on the architecture you require. (`release` matches the architecture of the
compiling machine, but we also always create a generic x86 Linux version with
[musl libc](https://www.musl-libc.org/), for portability reasons.)

## Included binaries

### `umt_tagger`

`umt_tagger` is a tool to process a BAM/CRAM file, and move a Unique Molecular Tag from the read name to an "optional tag".

## Installation

1. Install rust and cargo, probably with [rustup](https://rustup.rs/). The current required version is specified in `rust-toolchain.toml`, and should be automatically detected by your installation.
2. Run `make`, which will build the release binaries. This will take a while.
3. Move the compiled binaries to a directory on your `PATH`.

# Notes for developers

Follow the above instructions for installation, to make sure everything works.

To do a debug build (faster compilation and better error messages), type `make debug` in this directory.

## Code organization

We use the concept of
[workspaces](https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html) in
order to be able organize the code into separate packages, but still manage and
compile it all as a unit.

Compilation targets are built into `target/`.

All other directories are individual packages, which may contain either binary or library code. Current practice/design is for each package to contain at most one executable, for which the package is named, but this isn't set in stone yet.

### Top-level rust-related files

- `Cargo.toml`: Lists the packages in the workspace and compilation options
- `Cargo.lock`: All dependencies for any packages, do not edit by hand.
- `rust-toolchain.toml`: Specifies the targets to compile for and the required version of Rust.
- `.cargo/config`: Contains options specific to Cargo. TODO: See if the linker option in there is the correct place to put it.

### Package organization

- `Cargo.toml`: Contains information about the package - most importantly, the name, version, and dependencies.
- `src/*.rs`: Contains the source code for the package.
- `tests/*.rs`: Integration tests for the package.


In contrast to some other languages, it is convention in Rust to store unit tests in the same file as the code under test.



## Candidates for new binaries

The tools we intend to write

## Code style

* Try to keep code simple. Rust has a lot of ways to make things complicated, but we don't need most of them.
* Leverage existing libraries. Binding in proven C code is okay, especially for working with bioinformatic-specific formats or algorithms.
* Stick to safe mode.
* Don't spend too much time on optimization. Even a first-draft of the solution will likely be faster than the python code it would otherwise be written in.
* Write unit tests.

## Recommended libraries

This short list of libraries has proven itself useful.

* [structopt](https://docs.rs/structopt/0.3.22/structopt/) - declarative option parsing.
* [rust_htslib](https://docs.rs/rust-htslib/0.38.2/rust_htslib/) - bindings to the C htslib for working with BAM/CRAM files.
* [proptest](https://altsysrq.github.io/rustdoc/proptest/latest/proptest/) - Allows property-based testing.


### Library selection criteria

Some libraries that will be needed have not yet been chosen. Criteria in descending order of importance:

- Correct: The library must be correct, or what's the point?
- Maintained: For anything non-trivial, assurances that the library is maintained help keep the code correct.
- Clarity of use: All else equal, a package with fewer options and methods is easier to use correctly than a larger one.
- Performance: It should be straightforward to use the package in an performant manner, *if* the package is performance-sensitive. For example, for our purposes, code or data structures used in an inner-loop are performance-sensitive, network code is not.
