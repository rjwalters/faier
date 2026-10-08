# Fork changes

`faier` tracks upstream faer releases. Every change carried on top of upstream is listed here. Changes that don't affect library behaviour are not listed: lint allowances (e.g. #7's `clippy::eq_op` / `out_of_bounds_indexing` allows), formatting, and repository tooling (Loom, CI workflows, release configuration). Their history is in `git log`.

| Version | Change | Resolves |
|---|---|---|
| 0.24.4 | Crate renamed `faer` → `faier` and `faer-traits` → `faier-traits`; library names unchanged (`faer`, `faer_traits`). No code changes. | — |
| 0.25.0 | `qz_real` aggressive early deflation: `ihi - kwbot` wrapped (`kwbot = kwtop - 1` at row 0) and panicked under overflow checks; now `wrapping_sub`. Test: `tests/fork_qz_real_aed_index_overflow.rs`. | geode-fem [#908](https://github.com/rjwalters/geode-fem/issues/908) / [#867](https://github.com/rjwalters/geode-fem/issues/867) (QZ overflow-panic suppression) |
| 0.25.0 | Blocked generalized Hessenberg without `Z` cleared the whole rectangle `[jcol+2.., jcol..n-jcol-2]` instead of the panel's lower trapezoid, wiping the reduced pair: eigenvalues-only `gevd_*` above 32 rows returned mostly `alpha = beta = 0`. Test: `tests/fork_gen_hessenberg_without_qz.rs`. | geode-fem #908 (found reproducing #867) |
| 0.25.0 | `make_givens` scaling guard: `1 / h` and `1 / c` overflowed for subnormal `h = hypot(\|f\|, \|g\|)` or `\|f\| << \|g\|`, turning complex QZ iterates into NaN on degenerate-cluster shifts (the complex QZ "hang"). Test: `make_givens_tests` (unit). | geode-fem #908 / #867 item 1 |
| 0.25.0 | `make_givens` falls back to LAPACK `zlartg`/`dlartg` (Anderson's safe scaling as in LAPACK 3.10+, which scales `f` on its own when `\|f\| / \|g\|` is below `sqrt(min_positive)`) wherever the upstream formula is unsafe. The first guard took `sgn(f)` from `f / max(\|f\|, \|g\|)`, which lost the complex phase below `\|f\| / \|g\|` ~ 2e-308 (unitarity error 1e-8 at `\|g\| = 1e15`) and gave NaN below ~5e-324. Its retry-by-recursion also overflowed the stack on NaN input. Upstream results are unchanged, bit for bit, wherever upstream is safe. Test: `make_givens_tests` (unit; sweeps subnormal to 8e307 magnitudes, zero inputs and complex phases). | geode-fem #908 / #867 item 1 (review of rjwalters/faier#1) |
| 0.25.0 | QZ (real and complex, blocked and unblocked) stops at the first non-finite iterate instead of running out `30 n` sweeps; `gevd_*` report it as `GevdError::NoConvergence`. Test: `tests/fork_qz_nonfinite_fails_fast.rs`. | geode-fem #908 / #867 item 1 |
| 0.25.0 | `qz_real` eigenvalue-pair fix-up read the fail-fast's NaN `alphai` as the start of a conjugate pair, stepped past the last slot of odd-sized blocks and wrote `alphar[ihi + 1]`: an out-of-bounds panic in `qz_real::hessenberg_to_qz` for NaN/Inf input at n = 3, 5, 7. Now steps by one when `alphai[j]` is non-finite or `j` is the last slot (`07348a5`). Test: `real_unblocked_qz_nonfinite_does_not_overrun_small_odd_blocks` in `tests/fork_qz_nonfinite_fails_fast.rs`. | geode-fem #908 / #867 item 1 (found hardening the fail-fast; review of rjwalters/faier#1) |
| 0.25.0 | Unblocked QZ (real and complex) fail-fast checked only the last diagonal entry, so a NaN/Inf elsewhere (e.g. at (0, 0)) ran out `maxit` and returned all-zero eigenvalues; in the complex path an off-diagonal Inf made the deflation tolerance infinite and returned finite garbage. Now checks the whole input block up front and the active block each iteration. Tests: `real_unblocked_qz_fails_fast_on_nonfinite_anywhere_in_the_block`, `complex_unblocked_qz_fails_fast_on_nonfinite_anywhere_in_the_block` and companions in `tests/fork_qz_nonfinite_fails_fast.rs`. | rjwalters/faier#6 |
| 0.25.0 | `qz_real::chase_bulge_2x2`: the 2x3 work-block row rotation hit all three columns, so each double-shift chase perturbed `B` by O(\|B\|): 1e-2 eigenvalue errors and spurious complex pairs on definite pencils from ~600 rows. Test: `tests/fork_qz_real_bulge_chase_accuracy.rs`. | geode-fem #908 / #867 item 2 |
| 0.25.0 | Blocked QZ (real and complex) deflation-window spin: small-block whole-block window no longer capped at `(n-3)/3`, and the window QZ drops to the unblocked algorithm at recursion depth 2 (LAPACK `xlaqz0`). Several times faster from ~590 rows. Tests: `qz_cplx::aed_window_spin_tests` (unit; counts blocked-loop iterations, 155 vs 20259 before, so it does not depend on host load) and `tests/fork_qz_aed_window_spin.rs` (accuracy). | geode-fem #908 / #867 items 1, 2 |
| 0.25.0 | Tridiagonal divide and conquer scales `T` by a power of two to unit norm first, so the `8 eps` deflation tolerance is relative (LAPACK `dstedc`); small-norm matrices no longer merge distinct eigenvalues. Test: `tests/fork_tridiag_dc_small_norm.rs`. | geode-fem #908 / #867 item 3 |
| 0.25.0 | Cherry-picked upstream's unreleased `fix.gevd-313` (`7628d92`, sarah quiñones): `gevd_scratch` under-allocated for `n = 2` with eigenvectors, and the eigenvalues-only real QZ now writes both slots of a conjugate pair. Since 0.25.2 it also carries the test attribution from upstream's follow-up `ddbd2ae` (`// tests by @sjoelund` on the gevd test module); that commit's warning fixes were already handled differently here, where the counters feed the assertion messages. | upstream `fix.gevd-313` |
| 0.25.0 | QZ (real and complex) running out `maxit` on finite data now NaN-fills the unconverged block (LAPACK `xhgeqz` `INFO`), so `gevd_*` return `GevdError::NoConvergence` instead of `Ok` with `0 / 0` eigenvalues. A failed deflation-window QZ restores the window and keeps the eigenvalues it did converge as shifts (LAPACK `dlaqz3` / `zlaqz2`), writing them to the window's own slots. Tests: `qz_real::maxit_exhaustion_tests`, `qz_cplx::maxit_exhaustion_tests` (unit; a test-only `maxit` override forces the exhaustion). | rjwalters/faier#11 |
| 0.25.0 | Removed upstream files nothing in this fork uses: `.woodpecker/` (upstream's Codeberg Woodpecker CI; GitHub Actions runs CI here) and `faer-ffi/main.cpp` (an unreferenced C++ example). When an upstream sync modifies them, resolve the modify/delete conflict by keeping the deletion (`git rm`). | `/repo:all` hygiene pass, 2026-10-07 |
| 0.25.0 | `faer/Cargo.toml` requires `private-gemm-x86 = "0.1.22"` (upstream: `"0.1.20"`): every 0.1.x release up to 0.1.20 is yanked on crates.io and 0.1.21 does not exist. When an upstream sync touches this line, keep `0.1.22` (or the higher of the two), or `Cargo.lock` resolves a yanked version again. | rjwalters/faier#19 |
| 0.25.0 | **Breaking public API change.** Sparse LU takes a user-supplied fill-reducing column ordering: `faer::sparse::linalg::lu::factorize_symbolic_lu(A, ord, params)` gains a positional `ord: LuColOrdering` argument, with `LuColOrdering::{Colamd (default, the previous hardcoded behaviour), Identity, Custom(&[I])}`. `Custom(order)` is the elimination order in **new → old** orientation (`order[k]` is the original column eliminated at step `k`, i.e. the forward column permutation `SymbolicLu::col_perm()` returns); the inverse is built internally, and a wrong-length, out-of-range or duplicate order panics with a clear message. Also **breaking**: `faer::sparse::linalg::qr::column_counts_ata` now takes `A: SymbolicSparseRowMatRef` instead of `AT: SymbolicSparseColMatRef` (pass `A` row-major, no longer `A.transpose()`), and its doc now says the etree is of $A^\top A$ (was $A A^\top$). Operator decision (rjwalters/faier#26): keep the break, so the next `faier` release is **0.25.0**, deliberately departing from the "track upstream's minor" rule in [Versioning](#versioning) for this release. Upstream sync note: upstream callers of the two functions need the new arguments. Tests: `tests/fork_lu_col_ordering.rs` (`every_lu_col_ordering_solves_on_simplicial_and_supernodal_paths`, `default_lu_col_ordering_is_colamd`, `custom_lu_col_ordering_rejects_{duplicate,wrong_length,out_of_range}`) and `test_solver_lu_custom_col_ordering` (unit, `faer/src/sparse/linalg/lu.rs`). | codeberg [sarah-quinones/faer#307](https://codeberg.org/sarah-quinones/faer/issues/307) |

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
   git fetch origin
   git fetch upstream --tags
   git ls-remote --tags upstream 'faer-v*' | sed 's#.*refs/tags/##; s#\^{}##' | sort -uV | tail -3   # newest tags
   git rev-parse 'faer-vX.Y.Z^{commit}'
   git ls-remote upstream-codeberg 'refs/tags/faer-vX.Y.Z' 'refs/tags/faer-vX.Y.Z^{}' \
     | awk '{c = $1} /\^\{\}$/ {p = $1} END {print (p != "" ? p : c)}'   # must print the same SHA
   ```

   The `git fetch origin` makes step 3 branch from the current `main`. The cross-check compares commits, not tag objects: an annotated tag has a `^{}` line with the commit it points to, a lightweight tag has only the plain line, which is already the commit, and the `awk` prints the `^{}` SHA if there is one, else the plain one. If the SHAs differ, or the `ls-remote` command prints nothing (the tag is missing on Codeberg), stop and find out why before merging.

3. **Branch.** `git switch -c sync/faer-vX.Y.Z origin/main` (after the `git fetch origin` in step 2).

4. **Merge.** `git merge --no-ff faer-vX.Y.Z -m "Merge upstream faer-vX.Y.Z"`. Resolve conflicts with these points in mind:
   - Keep the renames in each `Cargo.toml`. Package names stay `faier` / `faier-traits`, and the library names stay `faer` / `faer_traits` through the `[lib] name = ...` sections. Path dependencies use `package = "faier"` / `package = "faier-traits"` (in `faer/`, `faer-ffi/` and `faer-no-std-test/`). `description` and `repository` stay the fork's. Set each `version` by the version rule below.
   - In code covered by a table row, keep the fork's behaviour unless upstream fixed the same bug. If it did, take upstream's code and handle the row in step 6.
   - **Files the fork owns outright: keep ours.** Upstream changes to these conflict on every sync, and the resolution is always the fork's side:
     - `README.md` is the fork's own README (it is also the crates.io page for both crates). Resolve with `git checkout --ours README.md`. If upstream's README gained something faier users need (a new MSRV, a changed feature), add it to the fork README by hand.
     - `.woodpecker/` and `faer-ffi/main.cpp` were removed (see their row). Keep the deletion: `git rm -r .woodpecker faer-ffi/main.cpp`.
     - `.github/workflows/run-tests.yml` and `code-quality.yml` are the fork's rewrites, with SHA-pinned actions, a pinned nightly, least-privilege `permissions` and no `faer-entity` steps. Keep ours (`git checkout --ours`), then port any real change from upstream's version by hand, such as a new MSRV or a new job, keeping the pins. `draft-pdf.yml` was deleted here; keep the deletion. The fork-only workflows (`release-plz.yml`, `upstream-tag-watch.yml`) and `.github/dependabot.yml` have no upstream counterpart.
     - The `private-gemm-x86` requirement in `faer/Cargo.toml` stays at `"0.1.22"` (or upstream's, if higher). See its row.
     - The gevd regression-test module (`mod t` in `faer/src/linalg/gevd/mod.rs`): once upstream releases its `ddbd2ae`, the `// tests by @sjoelund` line merges cleanly, but its other edits to the same test lines conflict with the fork's versions (the fork's counters feed the assertion messages). Keep the fork's side of those hunks.
   - Upstream's root `CHANGELOG.md` is upstream's: take theirs. The fork's per-crate changelogs (`faer/CHANGELOG.md`, `faer-traits/CHANGELOG.md`) have no upstream counterpart and never conflict.
   - **Version rule.** Apply it separately to `faier` (`faer/Cargo.toml`) and `faier-traits` (`faer-traits/Cargo.toml`). Let F be the fork's current manifest version and U the version in upstream's tag. The new version is the higher of F and U by semver:
     - U has a higher minor than F (a new upstream minor, e.g. `0.25.0`): take U. For `faier` this resets the patch.
     - Same minor, U's patch is higher than F's: take U.
     - Same minor, U's patch is equal to or lower than F's: keep F. Never move a version backwards; crates.io cannot reuse a published version.

     `faier-traits` versions independently (see "Releases"), so its F and U usually differ from `faier`'s; apply the rule to its own pair. F may be a version that is not on crates.io yet: release-plz keeps an unpublished manifest version and never lowers one, so the rule and the next release PR agree. Examples for `faier`:

     | Fork `faier` (F) | Upstream tag (U) | New `faier` version |
     |---|---|---|
     | `0.25.2` | `faer-v0.25.1` | `0.25.2` |
     | `0.25.2` | `faer-v0.26.0` | `0.26.0` |
     | `0.25.0` | `faer-v0.25.0` | `0.25.0` (already published, so the next release PR bumps it to `0.25.1`) |

     After the merge, check that the other manifests agree with the chosen versions: the `faier` dependency in `faer-ffi/Cargo.toml` (`version = "..."`), the `faier` dependency in `faer-no-std-test/Cargo.toml` (path only today; if upstream adds a `version`, it must match), and the `faier-traits` `version = "..."` requirement in `faer/Cargo.toml`. Then `cargo metadata --format-version 1 >/dev/null` must succeed, and `grep -A1 'name = "faier' Cargo.lock` must show only the chosen versions. Let cargo regenerate `Cargo.lock`; do not hand-edit it.

5. **Test.** Run the fork regression tests above, then the usual checks: `cargo +nightly-2026-06-09 fmt --all -- --check` (the toolchain `code-quality.yml` pins; `rustfmt.toml` needs a nightly rustfmt, and CI fails on any diff), `cargo clippy --workspace`, and `cargo test -p faier` (or `cargo nextest run`, which CI uses). Every fork test must still pass.

6. **Update this file.**
   - Under **Upstream base**, record the new tag and its commit SHA (`git rev-parse --short 'faer-vX.Y.Z^{commit}'`).
   - Change the version on the crate-rename row to the new upstream version `X.Y.Z`: the rename is reapplied on each base. Leave the version on every other row alone. Each one names the `faier` release that shipped it, or stays `unreleased` until a release ships it (see "Adding a fork change"); a sync does not change that.
   - **Drop rows that upstream has fixed.** First check that the upstream code fixes the bug: with the fork's change reverted to upstream's code, the row's test must still pass. Then delete the row. Delete the fork's version of the code, but keep its `fork_*` test as a regression guard. Upstream may fix a bug differently, so a missing conflict does not prove the fix is in.
   - Cherry-picked rows: run `git cherry faer-vX.Y.Z 7628d92 7628d92^`. It compares patches, so it also finds the fix if upstream squashed or rebased it before tagging. No output (the tag contains `7628d92` itself) or a line starting with `-` (the tag has an equivalent change) means the row for `fix.gevd-313` is no longer a fork change, so drop it. A line starting with `+` means no identical patch is in the tag; upstream may still have fixed the bug with different code, so apply the "drop rows that upstream has fixed" check above before keeping the row.

7. **Open a PR.** Push `sync/faer-vX.Y.Z` to `origin` and open a PR against `main`, with the fork-test output in the description. **Land it with a merge commit, not squash or rebase**, or the tag stops being an ancestor of `main`. After it lands, `git merge-base --is-ancestor faer-vX.Y.Z origin/main` must succeed.

### Adding a fork change

Every change carried on top of upstream needs both of these:

1. A row in the table above, with version `unreleased` until it ships, saying what changed and what it resolves.
2. A regression test: `faer/tests/fork_<topic>.rs` if the code can be reached through the public API. Otherwise use a unit test module inside the crate, and add its filter to the fork-test command above.

Whole-tree reformat commits are listed in `.git-blame-ignore-revs`. GitHub's blame view skips them; for local `git blame`, run `git config blame.ignoreRevsFile .git-blame-ignore-revs` once per clone.

### Watching for new releases

`.github/workflows/upstream-tag-watch.yml` runs weekly and on demand (`workflow_dispatch`). It compares the newest `faer-v*` tag on the GitHub mirror with the tag under **Upstream base**. If upstream is newer, it opens an issue titled `Sync upstream faer-vX.Y.Z` in this repo, unless an issue with exactly that title already exists, open or closed. A closed issue counts as handled: to skip a release deliberately, close its issue and the watcher will not reopen or recreate it. It only reads from upstream.

## Releases

`faier` and `faier-traits` are published to crates.io by [release-plz](https://release-plz.dev), using crates.io [Trusted Publishing](https://crates.io/docs/trusted-publishing) (GitHub OIDC) instead of a stored token. Nothing else in the workspace is published: `faer-macros`, `faer-ffi`, `faer-no-std-test` and `eigen-bench-setup` have `publish = false`, and `release-plz.toml` opts in only the two packages, by package name.

### Versioning

- **`faier` tracks upstream's minor and bumps the patch for fork releases.** Further fork releases on the same minor are patch releases (`0.25.1`, `0.25.2`, …).
- **The first release is `0.25.0`, one minor ahead of upstream.** It is based on upstream `faer-v0.24.4` plus the fixes in the table above, and it includes breaking public API changes (`LuColOrdering` and the `column_counts_ata` argument, see their row), so it cannot be a `0.24.x` patch. This is a deliberate one-time exception (operator decision, 2026-10-07; rjwalters/faier#26). No `0.24.x` version of `faier` exists. Syncs of upstream `faer-v0.24.x` keep `faier`'s `0.25.y`. Syncs of upstream `faer-v0.25.z` follow the version rule in [Syncing upstream](#syncing-upstream), step 4, as usual: `faier` keeps `0.25.y` if it is the higher version and takes `0.25.z` otherwise. From upstream's next minor (`faer-v0.26.0`) on, `faier` tracks upstream's minor again.
- **A new upstream minor resets the patch.** When a sync brings in a new upstream minor (e.g. `faer-v0.26.0`), `faier` takes that version (by the version rule in [Syncing upstream](#syncing-upstream), step 4), and later fork releases are `0.26.1`, ….
- **Upstream patch releases do not set `faier`'s patch.** After a sync of an upstream patch release (say `faer-v0.25.1` while `faier` is already at `0.25.2`), keep the higher of the two versions and let the next release bump it, because crates.io cannot reuse a published version.
- **`faier-traits` versions independently** (`0.24.0` now) and is bumped only when its own files change. `faier`'s `faer-traits = { …, version = "…" }` requirement follows it; release-plz updates that line.
- Semver build metadata (`+fork.1`) is not used: cargo ignores it when resolving versions.

### Tags and changelogs

| Package | Git tag / GitHub Release | Changelog |
|---|---|---|
| `faier` | `v{version}`, e.g. `v0.25.0` | `faer/CHANGELOG.md` |
| `faier-traits` | `faier-traits-v{version}`, e.g. `faier-traits-v0.24.0` | `faer-traits/CHANGELOG.md` |

The tags cannot collide with upstream's (`faer-v*`, `faer-traits-v*`, `faer-macros-v*`), which arrive with each sync. The changelogs are inside each package directory, so they ship in the crates, upstream has no file at those paths to conflict with, and upstream's hand-written root `CHANGELOG.md` is never rewritten. This table of fork changes stays the authoritative per-fix record; the package changelogs are the per-release summary.

### How a release happens

The workflow is `.github/workflows/release-plz.yml`. Both jobs are skipped unless the repository variable `RELEASE_PLZ_ENABLED` is `true` (see the checklist below).

1. **Release PR.** On every push to `main`, the `release-plz-pr` job opens or updates a PR from a `release-plz-*` branch. The PR bumps the versions in `faer/Cargo.toml` / `faer-traits/Cargo.toml` (and `Cargo.lock`) and adds the new commits to the package changelogs. It only covers packages whose packaged files changed since the version on crates.io (`cargo package --list -p <name>` shows those files).
2. **Review it like any PR.** Check the proposed version against the rule above and edit the changelog (it lists every commit that touched the package, including upstream's commits after a sync; trim it to what users need). cargo-semver-checks is off for `faier` (see `release-plz.toml`), so also check the public API by hand: `git diff v<last faier version> -- faer/src` must not change or remove any public item unless the release is a minor bump. A breaking change that wasn't marked `feat!:` or `BREAKING CHANGE:` would otherwise ship as a patch.
3. **Merge it.** The push of the merged release PR runs the `release-plz-release` job, which publishes `faier-traits` (if bumped) and then `faier` to crates.io, pushes the tags and creates the GitHub Releases from the changelog entries. Ordinary pushes to `main` publish nothing (`release_always = false`): release-plz only releases when the pushed commit belongs to a merged PR from a `release-plz-*` branch.

**Version bumps (release-plz 0.3.169, checked against its docs and source).** release-plz derives the bump from [conventional commits](https://www.conventionalcommits.org/) via the `next_version` crate. On `0.x`, `fix:` and `feat:` both propose a patch (`0.25.0` → `0.25.1`; `features_always_increment_minor = false` is set explicitly). A breaking change (`feat!:`, a `BREAKING CHANGE:` footer, or, for `faier-traits` only, a cargo-semver-checks failure; the check is disabled for `faier` in `release-plz.toml` because it reports false breaks on the `PermRef` methods) proposes the next minor (`0.26.0`). That conflicts with the rule above, so when it happens either fix the API break or override the version before merging: `release-plz set-version faier@0.25.1` on the release PR branch, or edit `version` in `faer/Cargo.toml` (and the `faer-ffi` path dependency). If `Cargo.toml` already holds a version that is not on crates.io (for example `0.25.0` after an upstream-minor sync), release-plz keeps it and only updates the changelog.

**CI on the release PR.** The release PR is opened with the workflow's `GITHUB_TOKEN`, and GitHub does not start workflows for events created by that token, so the test workflows do not run on it automatically. Before merging, close and reopen the PR (a `reopened` event from a person triggers `pull_request` workflows) and wait for green checks. No PAT or GitHub App is used, to keep publishing free of long-lived credentials.

**Dry runs.** `cargo publish --dry-run -p faier -p faier-traits` (cargo 1.90 or newer) must pass before any release. `cargo publish --dry-run -p faier` on its own fails until `faier-traits` is on crates.io ("no matching package named `faier-traits` found"), because cargo resolves `faier`'s dependency from the registry; publishing both packages in one command resolves it locally.

### Operator checklist

The repository's automation never publishes, tags or changes settings; these steps are the operator's. crates.io publishes are permanent (a version can be yanked, never replaced or deleted).

**Status: completed 2026-10-07/08.** `faier-traits` 0.24.0 and `faier` 0.25.0 were published by hand from `3222f12` (tags `faier-traits-v0.24.0`, `v0.25.0`). The publish token was revoked afterwards, trusted publishers were added on both crates, and `RELEASE_PLZ_ENABLED` is `true`. The first automated release, `faier` 0.25.1 and `faier-traits` 0.24.1 from release PR #37, was published through Trusted Publishing. The steps below are kept as a record and in case the setup ever has to be redone (e.g. a new crate).

**(a) crates.io ownership.** Log in to crates.io with the GitHub account that will own the crates (`rjwalters`) and verify the email address (crates.io requires a verified email to publish). Check that `faier` and `faier-traits` are still unclaimed: `https://crates.io/api/v1/crates/faier` should return 404.

**(b) First publish, with a short-lived token.** Trusted Publishing cannot create a crate: crates.io only lets you add a trusted publisher to a crate that already exists, and its docs state that the initial publish requires an API token (release-plz's docs say the same). So the first versions are published by hand, once:

1. Create a token at <https://crates.io/settings/tokens/new> with the endpoint scope `publish-new` only, crate scopes `faier` and `faier-traits`, and the shortest expiry offered (or a custom date of tomorrow).
2. On a clean checkout of `main` that contains this release setup (with `faer/Cargo.toml` at `0.25.0`), run the dry run above, then `cargo login` (paste the token at the prompt, not on the command line) and `cargo publish -p faier-traits -p faier`. Cargo publishes `faier-traits` first and waits for it to appear in the index before publishing `faier`. With cargo older than 1.90, run `cargo publish -p faier-traits`, then `cargo publish -p faier`.
3. Run `cargo logout` and revoke the token on <https://crates.io/settings/tokens>.

**(c) Add the trusted publishers.** For each crate: crates.io → the crate → Settings → Trusted Publishing → Add → GitHub, with repository owner `rjwalters`, repository name `faier`, workflow filename `release-plz.yml`, and no environment. If you add a GitHub environment as an extra gate, enter its name here and add `environment: <name>` to the `release-plz-release` job.

**(d) Allow the release PR.** GitHub → Settings → Actions → General → Workflow permissions: tick "Allow GitHub Actions to create and approve pull requests" (enabled on this repository on 2026-10-07; default workflow token permissions stay read-only). Without it the `release-plz-pr` job cannot open the release PR.

**(e) Tag the hand-published versions and enable the workflow.**

```sh
SHA=<commit the crates were published from>
git tag -a faier-traits-v0.24.0 "$SHA" -m "faier-traits 0.24.0"
git tag -a v0.25.0 "$SHA" -m "faier 0.25.0"
git push origin faier-traits-v0.24.0 v0.25.0
gh release create v0.25.0 -R rjwalters/faier --verify-tag --title v0.25.0 --notes "See faer/CHANGELOG.md and FORK.md."
gh release create faier-traits-v0.24.0 -R rjwalters/faier --verify-tag --title faier-traits-v0.24.0 --notes "See faer-traits/CHANGELOG.md."
gh variable set RELEASE_PLZ_ENABLED -R rjwalters/faier --body true
```

The next push to `main` then opens the first release PR, if anything packaged has changed since `0.25.0`. `release-plz release` skips versions that are already on crates.io, so enabling the workflow after the manual publish does not publish them again.

**(f) Badges.** After the first publish, add crates.io and docs.rs badges for `faier` to `README.md` (follow-up, not part of the release setup).
