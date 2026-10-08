# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This file covers `faier` releases only and is maintained by release-plz (see
FORK.md, "Releases"). Upstream faer's history is in the repository's root
`CHANGELOG.md`; the fork's changes on top of upstream are listed in `FORK.md`.

## [Unreleased]

## [0.25.2](https://github.com/rjwalters/faier/compare/v0.25.1...v0.25.2) - 2026-10-08

### Other

- *(gevd)* port upstream ddbd2ae's test attribution
- *(readme)* absolute links for crates.io, name gevd-313 and the internal Givens fix
- *(readme)* replace the embedded upstream README with a single fork README

## [0.25.1](https://github.com/rjwalters/faier/compare/v0.25.0...v0.25.1) - 2026-10-08

### Other

- *(readme)* crates.io and docs.rs badges for faier; drop the pre-release note

## [0.25.0]

First crates.io release of `faier`, based on upstream `faer-v0.24.4`. It is a
minor ahead of upstream because it contains breaking API changes (below); no
`0.24.x` version of `faier` exists. Code written against `faer 0.24` needs only
the two sparse-factorization changes listed under "Changed".

### Added

- Sparse LU: `LuColOrdering::{Colamd, Identity, Custom}` selects the fill-reducing column ordering for `factorize_symbolic_lu`, so a caller-supplied ordering (e.g. nested dissection or AMD) can replace the built-in COLAMD (`Colamd` remains the default). `Custom` takes the elimination order in new-to-old orientation.

### Changed

- **Breaking:** `factorize_symbolic_lu` takes a `LuColOrdering` argument (pass `LuColOrdering::Colamd` for the previous behaviour), and `qr::column_counts_ata` takes the matrix `A` (`SymbolicSparseRowMatRef`) instead of its transpose.

### Fixed

- QZ (generalized eigendecomposition): aggressive early deflation index overflow, blocked generalized Hessenberg without `Z`, `make_givens` scaling for subnormal and extreme inputs, fail-fast on non-finite iterates anywhere in the active block (`GevdError::NoConvergence`), real double-shift bulge-chase accuracy, blocked deflation-window spin, an out-of-bounds panic on non-finite input to the real unblocked QZ, and `maxit` exhaustion on finite data now reported as `NoConvergence` instead of returning `alpha = beta = 0` (undefined `0 / 0` eigenvalues).
- Tridiagonal divide and conquer: relative deflation tolerance for small-norm matrices.
- `gevd_scratch` for `n = 2` with eigenvectors, and the eigenvalues-only real QZ now writes both slots of a conjugate pair (cherry-picked upstream `fix.gevd-313`).

See the `0.25.0` rows in FORK.md for details and regression tests.
