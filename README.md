# faier — faer, with AI-assisted fixes

**About this fork.** `faier` is a maintained fork of [faer](https://github.com/sarah-quinones/faer-rs) by Sarah Quiñones, used by [geode-fem](https://github.com/rjwalters/geode-fem) and related projects. Upstream does not accept contributions written with AI assistance. Our development uses AI tools, so fixes we need cannot go upstream, and we carry them here instead.

- **Same API.** The crate is published as `faier`, but the library name is still `faer`, so existing code keeps `use faer::...`. Depend on it with a package rename:

  ```toml
  faer = { package = "faier", version = "0.24.4" }
  ```

- **Tracks upstream.** The fork follows upstream releases and adds a small set of fixes, each listed in [`FORK.md`](FORK.md) with the issue it resolves.
- **License and credit.** faer is MIT-licensed. The original copyright, license and third-party notices are retained, and all credit for faer itself belongs to its author.
- **Report fork bugs here,** at [rjwalters/faier](https://github.com/rjwalters/faier/issues), not upstream.

The upstream README follows.

---

<p align="center">
  <img src="https://faer.veganb.tw/faer-logo-color.png" alt="faer logo"/ width="25%">
</p>

# faer

[![documentation](https://docs.rs/faer/badge.svg)](https://docs.rs/faer)
[![crate](https://img.shields.io/crates/v/faer.svg)](https://crates.io/crates/faer)

`faer` is a rust crate that implements low level linear algebra routines and a high level wrapper for ease of use, in pure rust.
the aim is to provide a fully featured library for linear algebra with focus on portability, correctness, and performance.

see the [official website](https://faer.veganb.tw) and the [docs.rs](https://docs.rs/faer/latest/faer) documentation for code examples and usage instructions.

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

