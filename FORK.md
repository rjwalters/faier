# Fork changes

`faier` tracks upstream faer releases. Every change carried on top of upstream is listed here.

| Version | Change | Resolves |
|---|---|---|
| 0.24.4 | Crate renamed `faer` → `faier` and `faer-traits` → `faier-traits`; library names unchanged (`faer`, `faer_traits`). No code changes. | — |
| unreleased | `qz_real` aggressive early deflation: `ihi - kwbot` wrapped (`kwbot = kwtop - 1` at row 0) and panicked under overflow checks; now `wrapping_sub`. Test: `tests/fork_qz_real_aed_index_overflow.rs`. | geode-fem [#908](https://github.com/rjwalters/geode-fem/issues/908) / [#867](https://github.com/rjwalters/geode-fem/issues/867) (QZ overflow-panic suppression) |
| unreleased | Blocked generalized Hessenberg without `Z` cleared the whole rectangle `[jcol+2.., jcol..n-jcol-2]` instead of the panel's lower trapezoid, wiping the reduced pair: eigenvalues-only `gevd_*` above 32 rows returned mostly `alpha = beta = 0`. Test: `tests/fork_gen_hessenberg_without_qz.rs`. | geode-fem #908 (found reproducing #867) |
| unreleased | `make_givens` scaling guard: `1 / h` and `1 / c` overflowed for subnormal `h = hypot(\|f\|, \|g\|)` or `\|f\| << \|g\|`, turning complex QZ iterates into NaN on degenerate-cluster shifts (the complex QZ "hang"). Test: `make_givens_tests` (unit). | geode-fem #908 / #867 item 1 |
| unreleased | `make_givens` falls back to LAPACK `zlartg`/`dlartg` (Anderson's safe scaling as in LAPACK 3.10+, which scales `f` on its own when `\|f\| / \|g\|` is below `sqrt(min_positive)`) wherever the upstream formula is unsafe. The first guard took `sgn(f)` from `f / max(\|f\|, \|g\|)`, which lost the complex phase below `\|f\| / \|g\|` ~ 2e-308 (unitarity error 1e-8 at `\|g\| = 1e15`) and gave NaN below ~5e-324. Its retry-by-recursion also overflowed the stack on NaN input. Upstream results are unchanged, bit for bit, wherever upstream is safe. Test: `make_givens_tests` (unit; sweeps subnormal to 8e307 magnitudes, zero inputs and complex phases). | geode-fem #908 / #867 item 1 (review of rjwalters/faier#1) |
| unreleased | QZ (real and complex, blocked and unblocked) stops at the first non-finite iterate instead of running out `30 n` sweeps; `gevd_*` report it as `GevdError::NoConvergence`. Test: `tests/fork_qz_nonfinite_fails_fast.rs`. | geode-fem #908 / #867 item 1 |
| unreleased | `qz_real::chase_bulge_2x2`: the 2x3 work-block row rotation hit all three columns, so each double-shift chase perturbed `B` by O(\|B\|): 1e-2 eigenvalue errors and spurious complex pairs on definite pencils from ~600 rows. Test: `tests/fork_qz_real_bulge_chase_accuracy.rs`. | geode-fem #908 / #867 item 2 |
| unreleased | Blocked QZ (real and complex) deflation-window spin: small-block whole-block window no longer capped at `(n-3)/3`, and the window QZ drops to the unblocked algorithm at recursion depth 2 (LAPACK `xlaqz0`). Several times faster from ~590 rows. Test: `tests/fork_qz_aed_window_spin.rs`. | geode-fem #908 / #867 items 1, 2 |
| unreleased | Tridiagonal divide and conquer scales `T` by a power of two to unit norm first, so the `8 eps` deflation tolerance is relative (LAPACK `dstedc`); small-norm matrices no longer merge distinct eigenvalues. Test: `tests/fork_tridiag_dc_small_norm.rs`. | geode-fem #908 / #867 item 3 |
| unreleased | Cherry-picked upstream's unreleased `fix.gevd-313` (`7628d92`, sarah quiñones): `gevd_scratch` under-allocated for `n = 2` with eigenvectors, and the eigenvalues-only real QZ now writes both slots of a conjugate pair. | upstream `fix.gevd-313` |

## Upstream base

- `faer-v0.24.4` (`0539947`), from https://github.com/sarah-quinones/faer-rs
