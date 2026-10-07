//! faier fork regression test (geode-fem #908 / #867, items 1 and 2).
//!
//! The blocked QZ clamped its deflation window to `(n - 3) / 3` even when the
//! active block had shrunk below the blocking threshold and the window was
//! meant to cover the whole block. The window then never reached the top of
//! the block, no sweep ran (the block is below the threshold), and the loop
//! spun on deflation windows until it ran out its `30 n` iterations, and it
//! did so recursively inside each deflation window, which also had no
//! recursion limit. The work was finished by the final unblocked QZ, so the
//! answers were right but the blocked QZ was several times slower than the
//! unblocked one from about 590 rows (where the shift count doubles). The
//! window is now uncapped for small blocks and the deflation-window QZ drops
//! to the unblocked algorithm at recursion depth 2, as in lapack `xlaqz0`.
//!
//! The spin itself is caught deterministically by the unit test
//! `qz_cplx::aed_window_spin_tests`, which counts the blocked-loop iterations
//! on this pencil (155 with the fix, 20259 without). This test
//! checks the end-to-end accuracy of the blocked solve. It used to also assert
//! that the blocked solve took at most twice as long as the unblocked one,
//! but a wall-clock ratio of two back-to-back solves can flake on a loaded
//! host.

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::gevd::{
	ComputeEigenvectors, GevdParams, gevd_cplx, gevd_scratch,
};
use faer::{Mat, Par, c64};

struct Lcg(u64);
impl Lcg {
	fn next(&mut self) -> f64 {
		self.0 = self
			.0
			.wrapping_mul(6364136223846793005)
			.wrapping_add(1442695040888963407);
		((self.0 >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
	}
}

const N: usize = 600;

/// complex-symmetric pencil `(S^T D_K S, S^T D_M S)` with spectrum exactly
/// `D_K / D_M` and an exact null cluster of `n / 9` eigenvalues, as produced
/// by a Nedelec discretization with absorbing layers
fn null_cluster_pencil() -> (Mat<c64>, Mat<c64>, Vec<c64>, usize) {
	let n = N;
	let n_null = n / 9;
	let mut rng = Lcg(0x796);
	let s = Mat::<f64>::from_fn(n, n, |i, j| {
		(if i == j { 1.0 } else { 0.0 }) + 0.5 * rng.next() / (n as f64).sqrt()
	});
	let dk: Vec<c64> = (0..n)
		.map(|i| {
			c64::new(
				if i < n_null {
					0.0
				} else {
					1.0 + 50.0 * (rng.next() + 0.5)
				},
				0.0,
			)
		})
		.collect();
	let dm: Vec<c64> = (0..n)
		.map(|_| c64::new(1.0 + rng.next(), 0.3 * (rng.next() + 0.5)))
		.collect();
	let congruence = |d: &[c64]| {
		let ds = Mat::<c64>::from_fn(n, n, |i, j| d[i] * s[(i, j)]);
		let st = Mat::<c64>::from_fn(n, n, |i, j| c64::new(s[(j, i)], 0.0));
		&st * &ds
	};
	let want = dk[n_null..]
		.iter()
		.zip(&dm[n_null..])
		.map(|(k, m)| k / m)
		.collect();
	(congruence(&dk), congruence(&dm), want, n_null)
}

/// eigenvalues from `gevd_cplx` with the default (blocked) parameters
fn solve(a: &Mat<c64>, b: &Mat<c64>) -> Vec<c64> {
	let n = a.nrows();
	let (mut a, mut b) = (a.clone(), b.clone());
	let params: GevdParams = <GevdParams as faer::Auto<c64>>::auto();
	let mut buf = MemBuffer::new(gevd_scratch::<c64>(
		n,
		ComputeEigenvectors::No,
		ComputeEigenvectors::Yes,
		Par::Seq,
		params.into(),
	));
	let mut alpha = faer::diag::Diag::<c64>::zeros(n);
	let mut beta = faer::diag::Diag::<c64>::zeros(n);
	let mut u = Mat::<c64>::zeros(n, n);
	gevd_cplx(
		a.as_mut(),
		b.as_mut(),
		alpha.as_mut(),
		beta.as_mut(),
		None,
		Some(u.as_mut()),
		Par::Seq,
		MemStack::new(&mut buf),
		params.into(),
	)
	.unwrap();
	(0..n).map(|i| alpha[i] / beta[i]).collect()
}

#[test]
fn blocked_complex_qz_is_accurate_on_a_null_cluster_pencil() {
	let (a, b, mut want, n_null) = null_cluster_pencil();
	let got = solve(&a, &b);

	assert!(got.iter().all(|l| l.re.is_finite() && l.im.is_finite()));
	let null = got.iter().filter(|l| l.norm() < 1e-7).count();
	assert_eq!(null, n_null, "null cluster size");
	let mut phys: Vec<c64> =
		got.into_iter().filter(|l| l.norm() >= 1e-7).collect();
	phys.sort_by(|x, y| x.norm().total_cmp(&y.norm()));
	want.sort_by(|x, y| x.norm().total_cmp(&y.norm()));
	let err = phys
		.iter()
		.zip(&want)
		.map(|(g, w)| (g - w).norm() / w.norm())
		.fold(0.0, f64::max);
	assert!(err <= 1e-10, "max relative error vs exact spectrum {err:e}");
}
