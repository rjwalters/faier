//! faier fork regression test (geode-fem #908 / #867, item 2).
//!
//! The real blocked QZ's double-shift bulge chase (`chase_bulge_2x2`, lapack
//! `dlaqz2`) applied its first row rotation of the 2 x 3 work block to all
//! three columns instead of the last two, overwriting the entry it had just
//! annihilated. The right rotations computed from that block then failed to
//! zero `B(k + 1, k)` and `B(k + 2, k)`, which were cleared anyway, so every
//! chase step perturbed `B` by up to O(|B|). On symmetric-definite FEM-like
//! pencils this gave relative eigenvalue errors of 1e-2 and more, and spurious
//! complex-conjugate pairs, from about 600 rows up.

use faer::{Mat, Side};

/// 7-point 3-D Laplacian stiffness with a nearest-neighbour mass on an
/// `nx^3` grid: a symmetric-definite pencil with a well-conditioned,
/// strictly positive spectrum.
fn laplacian_3d(nx: usize) -> (Mat<f64>, Mat<f64>) {
	let n = nx * nx * nx;
	let idx = |i: usize, j: usize, k: usize| i + nx * (j + nx * k);
	let mut kk = Mat::<f64>::zeros(n, n);
	let mut mm = Mat::<f64>::zeros(n, n);
	for k in 0..nx {
		for j in 0..nx {
			for i in 0..nx {
				let p = idx(i, j, k);
				kk[(p, p)] = 6.0;
				mm[(p, p)] = 1.0;
				let mut nb = Vec::new();
				if i + 1 < nx {
					nb.push(idx(i + 1, j, k));
				}
				if j + 1 < nx {
					nb.push(idx(i, j + 1, k));
				}
				if k + 1 < nx {
					nb.push(idx(i, j, k + 1));
				}
				for q in nb {
					kk[(p, q)] = -1.0;
					kk[(q, p)] = -1.0;
					mm[(p, q)] = 0.1;
					mm[(q, p)] = 0.1;
				}
			}
		}
	}
	(kk, mm)
}

/// spectrum of the symmetric-definite pencil via `L^-1 K L^-T`, `M = L L^T`
fn symmetric_reference(k: &Mat<f64>, m: &Mat<f64>) -> Vec<f64> {
	let n = k.nrows();
	let l = m.llt(Side::Lower).unwrap().L().to_owned();
	let mut x = k.clone();
	faer::linalg::triangular_solve::solve_lower_triangular_in_place(
		l.as_ref(),
		x.as_mut(),
		faer::Par::Seq,
	);
	let mut y = x.transpose().to_owned();
	faer::linalg::triangular_solve::solve_lower_triangular_in_place(
		l.as_ref(),
		y.as_mut(),
		faer::Par::Seq,
	);
	let c = Mat::<f64>::from_fn(n, n, |i, j| 0.5 * (y[(i, j)] + y[(j, i)]));
	c.self_adjoint_eigenvalues(Side::Lower).unwrap()
}

#[test]
fn real_qz_matches_symmetric_reference_on_definite_pencil() {
	let (k, m) = laplacian_3d(9);
	let n = k.nrows();
	let want = symmetric_reference(&k, &m);
	let evd = k.generalized_eigen(&m).unwrap();
	let (sa, sb) = (evd.S_a().column_vector(), evd.S_b().column_vector());
	let mut got = Vec::with_capacity(n);
	for i in 0..n {
		let l = sa[i] / sb[i];
		assert!(
			l.im.abs() <= 1e-9 * l.re.abs(),
			"spurious complex eigenvalue {l:?} of a symmetric-definite pencil"
		);
		got.push(l.re);
	}
	got.sort_by(f64::total_cmp);
	let err = got
		.iter()
		.zip(&want)
		.map(|(g, w)| (g - w).abs() / w.abs())
		.fold(0.0, f64::max);
	assert!(err <= 1e-10, "max relative error {err:e}");
}
