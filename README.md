# faier — faer, with AI-assisted fixes

**About this fork.** `faier` is a maintained fork of [faer](https://github.com/sarah-quinones/faer-rs) by Sarah Quiñones, used by [geode-fem](https://github.com/rjwalters/geode-fem) and related projects. Upstream does not accept contributions written with AI assistance. Our development uses AI tools, so fixes we need cannot go upstream, and we carry them here instead.

- **Same API.** The package is named `faier`, but the library name is still `faer`, so existing code keeps `use faer::...`. See [Using faier](#using-faier) below.
- **Tracks upstream.** The fork follows upstream releases and adds a small set of fixes, each listed in [`FORK.md`](FORK.md) with the issue it resolves.
- **License and credit.** faer is MIT-licensed. The original copyright, license and third-party notices are retained, and all credit for faer itself belongs to its author.
- **Report fork bugs here,** at [rjwalters/faier](https://github.com/rjwalters/faier/issues), not upstream.

## Using faier

Rename the package in your `Cargo.toml` and keep the `faer` dependency key, so your code is unchanged:

```toml
[dependencies]
faer = { package = "faier", version = "0.24" }
```

```rust
use faer::Mat;
```

> **Not yet on crates.io.** The crates.io form above works only once `faier` has been published (tracked in [#3](https://github.com/rjwalters/faier/issues/3)). Until then, and for fixes that are on `main` but not yet released, use the git form:

```toml
[dependencies]
faer = { package = "faier", git = "https://github.com/rjwalters/faier" }
# or pin it: add `rev = "<commit>"` (or `branch = "main"`)
```

If you depend on `faer-traits` directly, rename that too; its library name is still `faer_traits`:

```toml
faer-traits = { package = "faier-traits", version = "0.24" }
```

The list of changes this fork carries on top of upstream, each with the issue it resolves, is in [`FORK.md`](FORK.md).

**Citation.** [`CITATION.cff`](CITATION.cff) and [`paper.md`](paper.md) describe upstream faer-rs and credit its author; please cite upstream. This fork is a derivative work.

The upstream README follows. Its badges and documentation links point at the upstream `faer` crate, which does not include this fork's fixes.

---

<p align="center">
  <img src="https://faer.veganb.tw/faer-logo-color.png" alt="faer logo"/ width="25%">
</p>

# faer

[![upstream faer documentation](https://docs.rs/faer/badge.svg)](https://docs.rs/faer)
[![upstream faer crate](https://img.shields.io/crates/v/faer.svg)](https://crates.io/crates/faer)
(upstream `faer`)

`faer` is a rust crate that implements low level linear algebra routines and a high level wrapper for ease of use, in pure rust.
the aim is to provide a fully featured library for linear algebra with focus on portability, correctness, and performance.

see the upstream [official website](https://faer.veganb.tw) and the upstream [docs.rs](https://docs.rs/faer/latest/faer) documentation for code examples and usage instructions.

questions about using the library, contributing, and future directions can be discussed in the [zulip server](https://faer.zulipchat.com).

# contributing

if you'd like to contribute to `faer`, check out the list of "good first issue"
issues. these are all (or should be) issues that are suitable for getting
started, and they generally include a detailed set of instructions for what to
do. please ask questions on the zulip server or the issue itself if anything
is unclear!

# minimum supported rust version

the current msrv is rust 1.84.0.

# benchmarks

see [the benchmark page](https://faer.veganb.tw/benchmarks/) on the main website.

