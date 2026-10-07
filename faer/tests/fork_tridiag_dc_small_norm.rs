//! faier fork regression test (geode-fem #908 / #867, item 3).
//!
//! The tridiagonal divide-and-conquer eigensolver (used by
//! `self_adjoint_eigen` from 128 rows up) deflates with
//! `8 eps max(max|d|, max|z|)`, where `z` is built from unit eigenvectors. For
//! a small-norm matrix that tolerance is effectively the absolute `8 eps`, so
//! distinct eigenvalues were merged: a tridiagonal scaled by 1e-13 came back
//! with relative eigenvalue errors around 1e-3. lapack's `dstedc` scales the
//! matrix to unit norm first; the fork now does the same with an exact power
//! of two.

use faer::{Mat, Side};

fn tridiagonal(alpha: &[f64], beta: &[f64], c: f64) -> Mat<f64> {
	let k = alpha.len();
	Mat::<f64>::from_fn(k, k, |i, j| {
		if i == j {
			c * alpha[i]
		} else if i + 1 == j {
			c * beta[i]
		} else if j + 1 == i {
			c * beta[j]
		} else {
			0.0
		}
	})
}

#[test]
fn divide_and_conquer_is_scale_invariant() {
	// above the divide-and-conquer threshold, with a tight cluster of
	// near-equal diagonal entries and weak coupling
	let k = 160;
	let alpha: Vec<f64> = (0..k)
		.map(|i| {
			if i % 4 == 0 {
				-1.0 + 1e-6 * i as f64
			} else {
				0.2 + 0.01 * ((i as f64) * 0.731).sin()
			}
		})
		.collect();
	let beta: Vec<f64> = (0..k - 1)
		.map(|i| 1e-3 * (1.0 + ((i as f64) * 0.377).cos()))
		.collect();
	let reference = tridiagonal(&alpha, &beta, 1.0)
		.self_adjoint_eigen(Side::Lower)
		.unwrap();
	let mu_ref = reference.S().column_vector();
	for c in [1e-13, 1e-19, 1e-150, 1e8, 1e150] {
		let t = tridiagonal(&alpha, &beta, c);
		let evd = t.self_adjoint_eigen(Side::Lower).unwrap();
		let mu = evd.S().column_vector();
		let worst = (0..k)
			.map(|i| (mu[i] / c - mu_ref[i]).abs())
			.fold(0.0, f64::max);
		assert!(worst <= 1e-13, "scale {c:e}: eigenvalue drift {worst:e}");
		// the eigenvectors still diagonalize the scaled matrix
		let u = evd.U();
		let r = (&t * u - u * evd.S()).norm_max();
		assert!(r <= 1e-13 * c, "scale {c:e}: residual {r:e}");
	}
}
