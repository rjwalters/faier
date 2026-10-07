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
| unreleased | Blocked QZ (real and complex) deflation-window spin: small-block whole-block window no longer capped at `(n-3)/3`, and the window QZ drops to the unblocked algorithm at recursion depth 2 (LAPACK `xlaqz0`). Several times faster from ~590 rows. Tests: `qz_cplx::aed_window_spin_tests` (unit; counts blocked-loop iterations, 155 vs 20259 before, so it does not depend on host load) and `tests/fork_qz_aed_window_spin.rs` (accuracy). | geode-fem #908 / #867 items 1, 2 |
| unreleased | Tridiagonal divide and conquer scales `T` by a power of two to unit norm first, so the `8 eps` deflation tolerance is relative (LAPACK `dstedc`); small-norm matrices no longer merge distinct eigenvalues. Test: `tests/fork_tridiag_dc_small_norm.rs`. | geode-fem #908 / #867 item 3 |
| unreleased | Cherry-picked upstream's unreleased `fix.gevd-313` (`7628d92`, sarah quiñones): `gevd_scratch` under-allocated for `n = 2` with eigenvectors, and the eigenvalues-only real QZ now writes both slots of a conjugate pair. | upstream `fix.gevd-313` |

## Upstream base

- `faer-v0.24.4` (`0539947`), from https://github.com/sarah-quinones/faer-rs

## Syncing upstream

How to bring a new upstream release (`faer-vX.Y.Z`) into `faier` without losing the fork's changes. Upstream does not take AI-assisted contributions, so each row in the table above is carried until upstream fixes the same bug on its own.

### Decisions (defaults; change them here if the policy changes)

- **Merge, not rebase.** `main` is published history: merged PRs, released crates, and commit SHAs cited above and in geode-fem issues. Merging an upstream tag keeps those SHAs valid, keeps the tag an ancestor of `main` (`git merge-base --is-ancestor faer-vX.Y.Z main` checks the base), and resolves each conflict once, in a reviewable merge commit.
- **Codeberg is canonical; the GitHub mirror is the fetch path.** Upstream lives at https://codeberg.org/sarah-quinones/faer. The GitHub mirror (`upstream` remote, https://github.com/sarah-quinones/faer-rs) carries the same `faer-v*` tags and works from GitHub Actions without extra setup. Before merging, check that the tag points to the same commit on both, which catches mirror lag or divergence. Both remotes are read-only to us: never push to them, and never open issues or PRs against them.
- **Sync from release tags only. `origin/upstream-main` is retired.** That branch tracks upstream's unreleased `main`, which this flow does not use, and keeping it around makes it easy to merge unreleased commits by accident. Take an unreleased upstream fix only by explicit cherry-pick (as with `fix.gevd-313`), and add a row for it in the table. The branch is not deleted yet; the operator will remove it (`git push origin --delete upstream-main`).

### Fork regression tests

Every fork change has a test: integration tests named `faer/tests/fork_*.rs`, plus two unit modules inside the crate. To run only those, from the repository root:

```sh
cargo test -p faier --test 'fork_*' && cargo test -p faier --lib -- make_givens_tests aed_window_spin_tests
```

The first command runs every `fork_*` integration test (cargo expands the glob). The second runs the unit modules `make_givens_tests` (`faer/src/linalg/gevd/gen_hessenberg/mod.rs`) and `qz_cplx::aed_window_spin_tests` (`faer/src/linalg/gevd/qz_cplx/mod.rs`). The second command should report 3 passed; check this, because a test filter that matches nothing still passes.

### Procedure

Example values: `X.Y.Z` is the new upstream version; `NEW` stands for the tag `faer-vX.Y.Z`.

1. **Remotes.** Run once per clone:

   ```sh
   git remote get-url upstream || git remote add upstream https://github.com/sarah-quinones/faer-rs.git
   git remote get-url upstream-codeberg || git remote add upstream-codeberg https://codeberg.org/sarah-quinones/faer.git
   git remote set-url --push upstream DISABLED
   git remote set-url --push upstream-codeberg DISABLED
   ```

   The `codeberg` remote (`codeberg.org/rjwalters/faer`) is the fork's own Codeberg copy, not upstream.

2. **Fetch the tag and cross-check it.**

   ```sh
   git fetch upstream --tags
   git ls-remote --tags upstream 'faer-v*' | sed 's#.*refs/tags/##; s#\^{}##' | sort -uV | tail -3   # newest tags
   git rev-parse 'faer-vX.Y.Z^{commit}'
   git ls-remote upstream-codeberg 'refs/tags/faer-vX.Y.Z^{}'   # must print the same SHA
   ```

   If the SHAs differ, stop and find out why before merging.

3. **Branch.** `git switch -c sync/faer-vX.Y.Z origin/main`

4. **Merge.** `git merge --no-ff faer-vX.Y.Z -m "Merge upstream faer-vX.Y.Z"`. Resolve conflicts with these points in mind:
   - Keep the renames in each `Cargo.toml`. Package names stay `faier` / `faier-traits`, and the library names stay `faer` / `faer_traits` through the `[lib] name = ...` sections. Path dependencies use `package = "faier"` / `package = "faier-traits"` (in `faer/`, `faer-ffi/` and `faer-no-std-test/`). `description` and `repository` stay the fork's. Take upstream's new `version` numbers.
   - In code covered by a table row, keep the fork's behaviour unless upstream fixed the same bug. If it did, take upstream's code and handle the row in step 6.

5. **Test.** Run the fork regression tests above, then the usual checks: `cargo fmt --all -- --check`, `cargo clippy --workspace`, and `cargo test -p faier` (or `cargo nextest run`, which CI uses). Every fork test must still pass.

6. **Update this file.**
   - Under **Upstream base**, record the new tag and its commit SHA (`git rev-parse --short 'faer-vX.Y.Z^{commit}'`).
   - Change the version on the crate-rename row, and on every `unreleased` row that ships with this sync, to the new version.
   - **Drop rows that upstream has fixed.** First check that the upstream code fixes the bug: with the fork's change reverted to upstream's code, the row's test must still pass. Then delete the row. Delete the fork's version of the code, but keep its `fork_*` test as a regression guard. Upstream may fix a bug differently, so a missing conflict does not prove the fix is in.
   - Cherry-picked rows: run `git merge-base --is-ancestor 7628d92 faer-vX.Y.Z && echo contained`. If the tag contains the commit, the row for `fix.gevd-313` is no longer a fork change, so drop it.

7. **Open a PR.** Push `sync/faer-vX.Y.Z` to `origin` and open a PR against `main`, with the fork-test output in the description. **Land it with a merge commit, not squash or rebase**, or the tag stops being an ancestor of `main`. After it lands, `git merge-base --is-ancestor faer-vX.Y.Z origin/main` must succeed.

### Adding a fork change

Every change carried on top of upstream needs both of these:

1. A row in the table above, with version `unreleased` until it ships, saying what changed and what it resolves.
2. A regression test: `faer/tests/fork_<topic>.rs` if the code can be reached through the public API. Otherwise use a unit test module inside the crate, and add its filter to the fork-test command above.

### Watching for new releases

`.github/workflows/upstream-tag-watch.yml` runs weekly and on demand (`workflow_dispatch`). It compares the newest `faer-v*` tag on the GitHub mirror with the tag under **Upstream base**. If upstream is newer, it opens an issue in this repo, unless an open issue for that tag already exists. It only reads from upstream.
