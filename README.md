# faier

[![faier on crates.io](https://img.shields.io/crates/v/faier.svg)](https://crates.io/crates/faier)
[![faier documentation](https://docs.rs/faier/badge.svg)](https://docs.rs/faier)

`faier` is a maintained fork of [faer](https://github.com/sarah-quinones/faer-rs), Sarah Quiñones's pure-Rust linear algebra library, carrying fixes that upstream can't take. It is used by [geode-fem](https://github.com/rjwalters/geode-fem) and related projects.

Upstream doesn't accept contributions written with AI assistance. Our development uses AI tools, so the fixes we need are kept here instead. The library is still called `faer`, so code keeps `use faer::...`.

## Using faier

Rename the package in your `Cargo.toml` and keep `faer` as the dependency key:

```toml
[dependencies]
faer = { package = "faier", version = "0.25" }
```

```rust
use faer::Mat;
```

For fixes that are on `main` but not yet released, use the git form:

```toml
[dependencies]
faer = { package = "faier", git = "https://github.com/rjwalters/faier" }
# or pin it: add `rev = "<commit>"` (or `branch = "main"`)
```

If you depend on `faer-traits` directly, rename that too; its library name is still `faer_traits`:

```toml
faer-traits = { package = "faier-traits", version = "0.24" }
# with the git form of faier, take faier-traits from git as well, or cargo
# builds two copies of the traits crate:
# faer-traits = { package = "faier-traits", git = "https://github.com/rjwalters/faier" }
```

## What's different from faer

`faier` 0.25 is based on faer 0.24.4. It adds:

- **Generalized eigenvalues (`gevd`, QZ).** Fixes for an index-overflow panic, wrong (mostly zero) eigenvalues from eigenvalues-only solves above 32 rows, loss of accuracy in the real double-shift bulge chase, and the deflation window spinning on larger problems. Non-finite input and running out of iterations are now reported as `GevdError::NoConvergence` instead of returning meaningless eigenvalues.
- **Givens rotations.** `make_givens` uses LAPACK's safe scaling (`zlartg`/`dlartg`), so subnormal and extreme-magnitude inputs no longer overflow or lose the complex phase.
- **Tridiagonal divide and conquer.** The matrix is scaled to unit norm first, so small-norm matrices no longer merge distinct eigenvalues.
- **Sparse LU column ordering.** `LuColOrdering::{Colamd, Identity, Custom}` lets you supply your own fill-reducing ordering.

Two API differences come with the LU change:

- `factorize_symbolic_lu` takes a column-ordering argument. Pass `LuColOrdering::Colamd` for faer's behaviour.
- Sparse QR's `column_counts_ata` takes the matrix `A` instead of its transpose.

Everything else matches faer 0.24.4. [`FORK.md`](FORK.md) lists every change with its regression test and the issue it resolves.

## Versions and releases

`faier` follows upstream's minor version and makes patch releases for fork fixes. The 0.25 line is a one-time exception, a minor ahead of faer 0.24 because of the LU API change, and `faier` returns to upstream's numbering from faer 0.26. Releases are published to crates.io by [release-plz](https://release-plz.dev) through crates.io Trusted Publishing. The details are in [`FORK.md`](FORK.md#releases); per-release notes are in [`faer/CHANGELOG.md`](faer/CHANGELOG.md).

The minimum supported Rust version is 1.84.0.

## Documentation

API documentation for this fork is on [docs.rs/faier](https://docs.rs/faier). For guides, examples and benchmarks, see upstream's [website](https://faer.veganb.tw) and [benchmark page](https://faer.veganb.tw/benchmarks/). They describe faer, which `faier` matches apart from the differences above.

## Contributing

Report bugs and send pull requests to [rjwalters/faier](https://github.com/rjwalters/faier/issues), not upstream, including bugs that also exist in faer. AI-assisted contributions are welcome here. Each fix should come with a regression test and a row in [`FORK.md`](FORK.md).

## License and credit

`faier` is MIT-licensed, like faer. The original copyright, license and third-party notices (Eigen, LAPACK, SuiteSparse) are kept, and all credit for faer itself belongs to its author. To cite the library, cite upstream faer-rs: [`CITATION.cff`](CITATION.cff) and [`paper.md`](paper.md) describe it. This fork is a derivative work.
