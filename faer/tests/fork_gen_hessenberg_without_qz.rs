//! faier fork regression test (geode-fem #908 / #867).
//!
//! The blocked generalized Hessenberg reduction (`n > block_size`), when
//! called without `Z`, cleared the whole rectangle `[jcol + 2.., jcol..n -
//! jcol - 2]` instead of only the rotations stored below the subdiagonal of
//! the current panel. That wiped live entries of the reduced pair, so every
//! eigenvalues-only `gevd_cplx` / `gevd_real` above 32 rows returned mostly
//! `alpha = beta = 0` (NaN eigenvalues).

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::gevd::gen_hessenberg::{
	GeneralizedHessenbergParams, generalized_hessenberg,
	generalized_hessenberg_scratch,
};
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
	fn cplx(&mut self) -> c64 {
		c64::new(self.next(), self.next())
	}
}

fn random_pair(n: usize, seed: u64) -> (Mat<c64>, Mat<c64>) {
	let mut rng = Lcg(seed);
	let a = Mat::<c64>::from_fn(n, n, |_, _| rng.cplx());
	let b = Mat::<c64>::from_fn(n, n, |i, j| {
		if i == j {
			c64::new(2.0, 0.0) + rng.cplx()
		} else if i < j {
			rng.cplx()
		} else {
			c64::new(0.0, 0.0)
		}
	});
	(a, b)
}

#[test]
fn blocked_reduction_without_q_and_z_matches_reduction_with_them() {
	let params: GeneralizedHessenbergParams =
		<GeneralizedHessenbergParams as faer::Auto<c64>>::auto();
	for n in [33usize, 50, 100] {
		let (a, b) = random_pair(n, n as u64);
		let mut buf =
			MemBuffer::new(generalized_hessenberg_scratch::<c64>(n, params));

		let (mut a1, mut b1) = (a.clone(), b.clone());
		generalized_hessenberg(
			a1.as_mut(),
			b1.as_mut(),
			None,
			None,
			Par::Seq,
			MemStack::new(&mut buf),
			params,
		);
		let (mut a2, mut b2) = (a.clone(), b.clone());
		let mut q = Mat::<c64>::identity(n, n);
		let mut z = Mat::<c64>::identity(n, n);
		generalized_hessenberg(
			a2.as_mut(),
			b2.as_mut(),
			Some(q.as_mut()),
			Some(z.as_mut()),
			Par::Seq,
			MemStack::new(&mut buf),
			params,
		);
		let da = (&a1 - &a2).norm_max() / a2.norm_max();
		let db = (&b1 - &b2).norm_max() / b2.norm_max();
		assert!(da <= 1e-14 && db <= 1e-14, "n = {n}: {da:e} {db:e}");
		// and the reduction is exact: `Q (H, T) Z^H = (A, B)`
		let ra = (&q * &a2 * z.adjoint() - &a).norm_max() / a.norm_max();
		let rb = (&q * &b2 * z.adjoint() - &b).norm_max() / b.norm_max();
		assert!(ra <= 1e-13 && rb <= 1e-13, "n = {n}: {ra:e} {rb:e}");
	}
}

fn gevd_cplx_eigenvalues(
	a: &Mat<c64>,
	b: &Mat<c64>,
	vectors: bool,
) -> Vec<c64> {
	let n = a.nrows();
	let (mut a, mut b) = (a.clone(), b.clone());
	let params: GevdParams = <GevdParams as faer::Auto<c64>>::auto();
	let ce = if vectors {
		ComputeEigenvectors::Yes
	} else {
		ComputeEigenvectors::No
	};
	let mut buf = MemBuffer::new(gevd_scratch::<c64>(
		n,
		ComputeEigenvectors::No,
		ce,
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
		vectors.then(|| u.as_mut()),
		Par::Seq,
		MemStack::new(&mut buf),
		params.into(),
	)
	.unwrap();
	(0..n).map(|i| alpha[i] / beta[i]).collect()
}

#[test]
fn eigenvalues_only_gevd_matches_gevd_with_eigenvectors() {
	for n in [50usize, 100] {
		let (a, b) = random_pair(n, 1000 + n as u64);
		let only = gevd_cplx_eigenvalues(&a, &b, false);
		let with = gevd_cplx_eigenvalues(&a, &b, true);
		// greedy one-to-one match
		let mut used = vec![false; n];
		for x in &only {
			assert!(x.re.is_finite() && x.im.is_finite(), "n = {n}: {x:?}");
			let (j, d) = with
				.iter()
				.enumerate()
				.filter(|(j, _)| !used[*j])
				.map(|(j, y)| (j, (x - y).norm() / y.norm().max(1.0)))
				.min_by(|p, q| p.1.total_cmp(&q.1))
				.unwrap();
			used[j] = true;
			assert!(d <= 1e-9, "n = {n}: {x:?} vs {:?}", with[j]);
		}
	}
}
