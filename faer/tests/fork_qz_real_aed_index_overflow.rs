//! faier fork regression test (geode-fem #908 / #867).
//!
//! The real blocked QZ's aggressive early deflation computed the number of
//! deflated eigenvalues as `ihi - kwbot`, where `kwbot` is `kwtop - 1` and
//! wraps to `usize::MAX` when the deflation window starts at row 0. Release
//! builds get the right answer from the wrapping arithmetic, but with
//! overflow checks on (the default `test` and `dev` profiles) the subtraction
//! panics with "attempt to subtract with overflow".

use faer::Mat;

/// 1-D Dirichlet Laplacian stiffness `tridiag(-1, 2, -1)` with the consistent
/// mass `tridiag(1, 4, 1) / 6`. Its generalized eigenvalues are known in
/// closed form: `(2 - 2 cos t) / ((4 + 2 cos t) / 6)`, `t = k pi / (n + 1)`.
fn laplacian_pencil(n: usize) -> (Mat<f64>, Mat<f64>) {
	let k = Mat::<f64>::from_fn(n, n, |i, j| {
		if i == j {
			2.0
		} else if i.abs_diff(j) == 1 {
			-1.0
		} else {
			0.0
		}
	});
	let m = Mat::<f64>::from_fn(n, n, |i, j| {
		if i == j {
			4.0 / 6.0
		} else if i.abs_diff(j) == 1 {
			1.0 / 6.0
		} else {
			0.0
		}
	});
	(k, m)
}

#[test]
fn real_qz_aed_window_at_row_zero_does_not_overflow() {
	let n = 120;
	let (k, m) = laplacian_pencil(n);
	let evd = k.generalized_eigen(&m).unwrap();
	let (sa, sb) = (evd.S_a().column_vector(), evd.S_b().column_vector());
	let mut got: Vec<f64> = (0..n)
		.map(|i| {
			let l = sa[i] / sb[i];
			assert!(l.im.abs() <= 1e-10 * l.re.abs(), "eigenvalue {l:?}");
			l.re
		})
		.collect();
	got.sort_by(f64::total_cmp);
	let mut want: Vec<f64> = (1..=n)
		.map(|k| {
			let c = (k as f64 * core::f64::consts::PI / (n as f64 + 1.0)).cos();
			(2.0 - 2.0 * c) / ((4.0 + 2.0 * c) / 6.0)
		})
		.collect();
	want.sort_by(f64::total_cmp);
	for (g, w) in got.iter().zip(&want) {
		assert!((g - w).abs() <= 1e-10 * w.abs(), "{g} vs {w}");
	}
}
