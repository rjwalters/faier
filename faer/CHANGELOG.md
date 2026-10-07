# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This file covers `faier` releases only and is maintained by release-plz (see
FORK.md, "Releases"). Upstream faer's history is in the repository's root
`CHANGELOG.md`; the fork's changes on top of upstream are listed in `FORK.md`.

## [Unreleased]

## [0.24.5]

First crates.io release of `faier`, based on upstream `faer-v0.24.4`.

### Fixed

- QZ (generalized eigendecomposition): aggressive early deflation index overflow, blocked generalized Hessenberg without `Z`, `make_givens` scaling for subnormal and extreme inputs, fail-fast on non-finite iterates (`GevdError::NoConvergence`), real double-shift bulge-chase accuracy, and blocked deflation-window spin.
- Tridiagonal divide and conquer: relative deflation tolerance for small-norm matrices.
- `gevd_scratch` for `n = 2` with eigenvectors (cherry-picked upstream `fix.gevd-313`).

See the `0.24.5` rows in FORK.md for details and regression tests.
