//! faier fork regression test (geode-fem #908 / #867).
//!
//! A non-finite QZ iterate never satisfies a deflation test, and neither the
//! blocked nor the unblocked QZ loop checked for one, so a NaN that entered
//! the iteration (for example from the `make_givens` underflow, #867 item 1)
//! made the loop run out its whole `30 n` sweep budget before returning
//! meaningless eigenvalues: 14.8 h on a 3300 x 3300 pencil. The loops now stop
//! at the first non-finite iterate and report NaN eigenvalues, which
//! `gevd_real` / `gevd_cplx` turn into `GevdError::NoConvergence`.

use faer::diag::Diag;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::evd::ComputeEigenvectors;
use faer::linalg::gevd::{
	GeneralizedSchurParams, GevdError, GevdParams, gevd_cplx, gevd_real,
	gevd_scratch, qz_cplx, qz_real,
};
use faer::{Col, Mat, Par, c64};
use std::sync::mpsc;
use std::time::{Duration, Instant};

/// The fixed loop returns at once; the unfixed one runs `30 n` sweeps
const BUDGET: Duration = Duration::from_secs(30);

fn with_watchdog<T: Send + 'static>(
	what: &str,
	f: impl FnOnce() -> T + Send + 'static,
) -> T {
	let (tx, rx) = mpsc::channel();
	let start = Instant::now();
	std::thread::spawn(move || {
		let _ = tx.send(f());
	});
	match rx.recv_timeout(BUDGET) {
		Ok(v) => {
			eprintln!("{what}: {:.3} s", start.elapsed().as_secs_f64());
			v
		},
		Err(_) => panic!("{what} did not return within {} s", BUDGET.as_secs()),
	}
}

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

const N: usize = 400;

/// an upper Hessenberg / upper triangular pair with one NaN on the diagonal
fn poisoned_pair_cplx(rng: &mut Lcg) -> (Mat<c64>, Mat<c64>) {
	let mut mk = || c64::new(rng.next(), rng.next());
	let a = Mat::<c64>::from_fn(N, N, |i, j| {
		if i == N / 2 && j == N / 2 {
			c64::new(f64::NAN, 0.0)
		} else if i <= j + 1 {
			mk()
		} else {
			c64::new(0.0, 0.0)
		}
	});
	let b = Mat::<c64>::from_fn(N, N, |i, j| {
		if i == j {
			c64::new(2.0, 0.0) + mk()
		} else if i < j {
			mk()
		} else {
			c64::new(0.0, 0.0)
		}
	});
	(a, b)
}

fn poisoned_pair_real(rng: &mut Lcg) -> (Mat<f64>, Mat<f64>) {
	let a = Mat::<f64>::from_fn(N, N, |i, j| {
		if i == N / 2 && j == N / 2 {
			f64::NAN
		} else if i <= j + 1 {
			rng.next()
		} else {
			0.0
		}
	});
	let b = Mat::<f64>::from_fn(N, N, |i, j| {
		if i == j {
			2.0 + rng.next()
		} else if i < j {
			rng.next()
		} else {
			0.0
		}
	});
	(a, b)
}

#[test]
fn complex_qz_stops_on_a_nonfinite_iterate() {
	let mut rng = Lcg(0x908);
	let (mut a, mut b) = poisoned_pair_cplx(&mut rng);
	let (alpha, beta) =
		with_watchdog("complex qz", move || {
			let params: GeneralizedSchurParams =
				<GeneralizedSchurParams as faer::Auto<c64>>::auto();
			let mut buf = MemBuffer::new(qz_cplx::hessenberg_to_qz_scratch::<
				c64,
			>(N, Par::Seq, params));
			let mut alpha = Col::<c64>::zeros(N);
			let mut beta = Col::<c64>::zeros(N);
			qz_cplx::hessenberg_to_qz(
				a.as_mut(),
				b.as_mut(),
				None,
				None,
				alpha.as_mut(),
				beta.as_mut(),
				ComputeEigenvectors::No,
				Par::Seq,
				params,
				MemStack::new(&mut buf),
			);
			(alpha, beta)
		});
	assert!(
		(0..N).any(|i| !(alpha[i].re.is_finite() && beta[i].re.is_finite())),
		"a breakdown must be reported through non-finite eigenvalues"
	);
}

#[test]
fn real_qz_stops_on_a_nonfinite_iterate() {
	let mut rng = Lcg(0x909);
	let (mut a, mut b) = poisoned_pair_real(&mut rng);
	let (alphar, beta) =
		with_watchdog("real qz", move || {
			let params: GeneralizedSchurParams =
				<GeneralizedSchurParams as faer::Auto<f64>>::auto();
			let mut buf = MemBuffer::new(qz_real::hessenberg_to_qz_scratch::<
				f64,
			>(N, Par::Seq, params));
			let mut alphar = Col::<f64>::zeros(N);
			let mut alphai = Col::<f64>::zeros(N);
			let mut beta = Col::<f64>::zeros(N);
			qz_real::hessenberg_to_qz(
				a.as_mut(),
				b.as_mut(),
				None,
				None,
				alphar.as_mut(),
				alphai.as_mut(),
				beta.as_mut(),
				ComputeEigenvectors::No,
				Par::Seq,
				params,
				MemStack::new(&mut buf),
			);
			(alphar, beta)
		});
	assert!(
		(0..N).any(|i| !(alphar[i].is_finite() && beta[i].is_finite())),
		"a breakdown must be reported through non-finite eigenvalues"
	);
}

/// an `n x n` upper Hessenberg / upper triangular pair with `poison` at
/// `A[at]` (if any) and a zero at `B[(k, k)]` for `k` in `b_zero` (if any)
fn small_pair_real(
	n: usize,
	poison: Option<((usize, usize), f64)>,
	b_zero: Option<usize>,
	seed: u64,
) -> (Mat<f64>, Mat<f64>) {
	let mut rng = Lcg(seed);
	let a = Mat::<f64>::from_fn(n, n, |i, j| {
		let x = if i <= j + 1 { rng.next() } else { 0.0 };
		match poison {
			Some((at, poison)) if at == (i, j) => poison,
			_ => x,
		}
	});
	let b = Mat::<f64>::from_fn(n, n, |i, j| {
		if i == j {
			let x = 2.0 + rng.next();
			if b_zero == Some(i) { 0.0 } else { x }
		} else if i < j {
			rng.next()
		} else {
			0.0
		}
	});
	(a, b)
}

fn small_pair_cplx(
	n: usize,
	poison: Option<((usize, usize), f64)>,
	b_zero: Option<usize>,
	seed: u64,
) -> (Mat<c64>, Mat<c64>) {
	let mut rng = Lcg(seed);
	let mut mk = || c64::new(rng.next(), rng.next());
	let a = Mat::<c64>::from_fn(n, n, |i, j| {
		let x = if i <= j + 1 { mk() } else { c64::new(0.0, 0.0) };
		match poison {
			Some((at, poison)) if at == (i, j) => c64::new(poison, 0.0),
			_ => x,
		}
	});
	let b = Mat::<c64>::from_fn(n, n, |i, j| {
		if i == j {
			let x = c64::new(2.0, 0.0) + mk();
			if b_zero == Some(i) {
				c64::new(0.0, 0.0)
			} else {
				x
			}
		} else if i < j {
			mk()
		} else {
			c64::new(0.0, 0.0)
		}
	});
	(a, b)
}

/// default qz parameters, or a blocking threshold low enough that the
/// blocked qz runs from `nmin = 15` rows
fn qz_params<T: faer::traits::ComplexField>(
	blocked_from_15: bool,
) -> GeneralizedSchurParams {
	let mut params: GeneralizedSchurParams =
		<GeneralizedSchurParams as faer::Auto<T>>::auto();
	if blocked_from_15 {
		params.blocking_threshold = 15;
	}
	params
}

/// runs the real QZ (eigenvalues only) on `(a, b)`, returning
/// `(alphar, alphai, beta)`
fn real_qz(
	mut a: Mat<f64>,
	mut b: Mat<f64>,
	params: GeneralizedSchurParams,
) -> (Col<f64>, Col<f64>, Col<f64>) {
	let n = a.nrows();
	let mut buf = MemBuffer::new(qz_real::hessenberg_to_qz_scratch::<f64>(
		n,
		Par::Seq,
		params,
	));
	let mut alphar = Col::<f64>::zeros(n);
	let mut alphai = Col::<f64>::zeros(n);
	let mut beta = Col::<f64>::zeros(n);
	qz_real::hessenberg_to_qz(
		a.as_mut(),
		b.as_mut(),
		None,
		None,
		alphar.as_mut(),
		alphai.as_mut(),
		beta.as_mut(),
		ComputeEigenvectors::No,
		Par::Seq,
		params,
		MemStack::new(&mut buf),
	);
	(alphar, alphai, beta)
}

/// runs the complex QZ (eigenvalues only) on `(a, b)`, returning
/// `(alpha, beta)`
fn cplx_qz(
	mut a: Mat<c64>,
	mut b: Mat<c64>,
	params: GeneralizedSchurParams,
) -> (Col<c64>, Col<c64>) {
	let n = a.nrows();
	let mut buf = MemBuffer::new(qz_cplx::hessenberg_to_qz_scratch::<c64>(
		n,
		Par::Seq,
		params,
	));
	let mut alpha = Col::<c64>::zeros(n);
	let mut beta = Col::<c64>::zeros(n);
	qz_cplx::hessenberg_to_qz(
		a.as_mut(),
		b.as_mut(),
		None,
		None,
		alpha.as_mut(),
		beta.as_mut(),
		ComputeEigenvectors::No,
		Par::Seq,
		params,
		MemStack::new(&mut buf),
	);
	(alpha, beta)
}

/// runs the real QZ on an `n x n` Hessenberg / triangular pair whose `(k, k)`
/// entry of `A` is `poison`, returning `(alphar, alphai, beta)`
fn real_qz_small(
	n: usize,
	k: usize,
	poison: f64,
	seed: u64,
) -> (Col<f64>, Col<f64>, Col<f64>) {
	let (a, b) = small_pair_real(n, Some(((k, k), poison)), None, seed);
	real_qz(a, b, qz_params::<f64>(false))
}

fn real_has_nonfinite(
	alphar: &Col<f64>,
	alphai: &Col<f64>,
	beta: &Col<f64>,
) -> bool {
	(0..alphar.nrows()).any(|i| {
		!(alphar[i].is_finite() && alphai[i].is_finite() && beta[i].is_finite())
	})
}

fn cplx_has_nonfinite(alpha: &Col<c64>, beta: &Col<c64>) -> bool {
	(0..alpha.nrows()).any(|i| {
		!(alpha[i].re.is_finite()
			&& alpha[i].im.is_finite()
			&& beta[i].re.is_finite()
			&& beta[i].im.is_finite())
	})
}

const SMALL_SIZES: [usize; 5] = [3, 4, 5, 7, 8];
const POISONS: [f64; 3] = [f64::INFINITY, f64::NEG_INFINITY, f64::NAN];

/// `A` positions to poison: the top-left and an interior diagonal entry, a
/// subdiagonal entry and an above-diagonal entry, none on the trailing
/// diagonal the unblocked fail-fast used to inspect
fn small_positions(n: usize) -> [(usize, usize); 4] {
	[(0, 0), (n / 2, n / 2), (n / 2, n / 2 - 1), (0, n - 1)]
}

/// below the blocked threshold the unblocked QZ runs, and its NaN fail-fast
/// fills `alphai` with NaN; the conjugate-pair fix-up used to read
/// `NaN != 0` as a pair start and step past the end of the block for odd
/// sizes (index out of bounds). even sizes are the control.
#[test]
fn real_unblocked_qz_nonfinite_does_not_overrun_small_odd_blocks() {
	for n in SMALL_SIZES {
		for (k, poison) in [
			(0, f64::NAN),
			(0, f64::INFINITY),
			(n - 1, f64::NAN),
			(n - 1, f64::INFINITY),
		] {
			let (alphar, alphai, beta) =
				real_qz_small(n, k, poison, 0x90A + n as u64);
			assert!(
				real_has_nonfinite(&alphar, &alphai, &beta),
				"n = {n}, poison {poison} at ({k}, {k}): a breakdown must be \
				 reported through non-finite eigenvalues"
			);
		}
	}
}

/// the unblocked fail-fast used to inspect only the trailing diagonal entry
/// `H[(ilast, ilast)]` / `T[(ilast, ilast)]`. a non-finite entry anywhere
/// else never reaches it, so nothing deflated, the loop ran out `maxit` and
/// returned with the eigenvalues untouched (all `0.0` here) and no sign of
/// failure (rjwalters/faier#6). it now checks the whole active block
#[test]
fn real_unblocked_qz_fails_fast_on_nonfinite_anywhere_in_the_block() {
	let mut missed = Vec::new();
	for n in SMALL_SIZES {
		for at in small_positions(n) {
			for poison in POISONS {
				let (a, b) = small_pair_real(
					n,
					Some((at, poison)),
					None,
					0x60 + n as u64,
				);
				let (alphar, alphai, beta) =
					real_qz(a, b, qz_params::<f64>(false));
				if !real_has_nonfinite(&alphar, &alphai, &beta) {
					missed.push(format!(
						"n = {n}, {poison} at {at:?}: alphar = {:?}, beta = \
						 {:?}",
						&alphar, &beta
					));
				}
			}
		}
	}
	assert!(
		missed.is_empty(),
		"{} breakdown(s) not reported through non-finite eigenvalues:\n{}",
		missed.len(),
		missed.join("\n")
	);
}

/// complex twin of the test above (`qz_cplx` had the identical gap)
#[test]
fn complex_unblocked_qz_fails_fast_on_nonfinite_anywhere_in_the_block() {
	let mut missed = Vec::new();
	for n in SMALL_SIZES {
		for at in small_positions(n) {
			for poison in POISONS {
				let (a, b) = small_pair_cplx(
					n,
					Some((at, poison)),
					None,
					0x61 + n as u64,
				);
				let (alpha, beta) = cplx_qz(a, b, qz_params::<c64>(false));
				if !cplx_has_nonfinite(&alpha, &beta) {
					missed.push(format!(
						"n = {n}, {poison} at {at:?}: alpha = {:?}, beta = \
						 {:?}",
						&alpha, &beta
					));
				}
			}
		}
	}
	assert!(
		missed.is_empty(),
		"{} breakdown(s) not reported through non-finite eigenvalues:\n{}",
		missed.len(),
		missed.join("\n")
	);
}

/// regression guard: the blocked QZ (here from 15 rows, and at 400 rows with
/// the default threshold) already checks its active block and must keep
/// failing fast for poison at the top-left and in the interior
#[test]
fn blocked_qz_fails_fast_on_nonfinite_top_left_and_interior() {
	for (n, blocked_from_15) in [(20usize, true), (N, false)] {
		for at in [(0, 0), (n / 2, n / 2)] {
			for poison in [f64::INFINITY, f64::NAN] {
				let what = format!("n = {n}, poison {poison} at {at:?}");
				let (real_bad, cplx_bad) =
					with_watchdog(&what.clone(), move || {
						let (a, b) =
							small_pair_real(n, Some((at, poison)), None, 0x62);
						let (alphar, alphai, beta) =
							real_qz(a, b, qz_params::<f64>(blocked_from_15));
						let (a, b) =
							small_pair_cplx(n, Some((at, poison)), None, 0x63);
						let (alpha, cbeta) =
							cplx_qz(a, b, qz_params::<c64>(blocked_from_15));
						(
							real_has_nonfinite(&alphar, &alphai, &beta),
							cplx_has_nonfinite(&alpha, &cbeta),
						)
					});
				assert!(real_bad, "real qz, {what}: breakdown not reported");
				assert!(cplx_bad, "complex qz, {what}: breakdown not reported");
			}
		}
	}
}

/// a genuine infinite eigenvalue (`B[(k, k)] = 0` on finite data, so
/// `beta = 0` with a finite `alpha`) is not a breakdown: every output stays
/// finite and one `beta` is zero, unblocked and blocked
#[test]
fn qz_does_not_flag_a_genuine_infinite_eigenvalue() {
	for (n, blocked_from_15) in [
		(3usize, false),
		(4, false),
		(5, false),
		(7, false),
		(8, false),
		(20, true),
	] {
		for k in [0, n / 2, n - 1] {
			let (a, b) = small_pair_real(n, None, Some(k), 0x64 + n as u64);
			let (alphar, alphai, beta) =
				real_qz(a, b, qz_params::<f64>(blocked_from_15));
			assert!(
				!real_has_nonfinite(&alphar, &alphai, &beta),
				"real, n = {n}, B[({k}, {k})] = 0: alphar = {alphar:?}, \
				 alphai = {alphai:?}, beta = {beta:?}"
			);
			assert!(
				(0..n).any(|i| beta[i].abs() <= 1e-12),
				"real, n = {n}, B[({k}, {k})] = 0: no zero beta in {beta:?}"
			);

			let (a, b) = small_pair_cplx(n, None, Some(k), 0x65 + n as u64);
			let (alpha, beta) =
				cplx_qz(a, b, qz_params::<c64>(blocked_from_15));
			assert!(
				!cplx_has_nonfinite(&alpha, &beta),
				"complex, n = {n}, B[({k}, {k})] = 0: alpha = {alpha:?}, beta \
				 = {beta:?}"
			);
			assert!(
				(0..n).any(|i| beta[i].re.hypot(beta[i].im) <= 1e-12),
				"complex, n = {n}, B[({k}, {k})] = 0: no zero beta in {beta:?}"
			);
		}
	}
}

/// the public drivers keep reporting the poisoned pencils as
/// `GevdError::NoConvergence`
#[test]
fn gevd_reports_poisoned_small_pencils_as_no_convergence() {
	for n in SMALL_SIZES {
		for at in small_positions(n) {
			let (mut a, mut b) =
				small_pair_real(n, Some((at, f64::INFINITY)), None, 0x66);
			let params: GevdParams = <GevdParams as faer::Auto<f64>>::auto();
			let mut buf = MemBuffer::new(gevd_scratch::<f64>(
				n,
				ComputeEigenvectors::No,
				ComputeEigenvectors::No,
				Par::Seq,
				params.into(),
			));
			let mut s_re = Diag::<f64>::zeros(n);
			let mut s_im = Diag::<f64>::zeros(n);
			let mut beta = Diag::<f64>::zeros(n);
			let r = gevd_real(
				a.as_mut(),
				b.as_mut(),
				s_re.as_mut(),
				s_im.as_mut(),
				beta.as_mut(),
				None,
				None,
				Par::Seq,
				MemStack::new(&mut buf),
				params.into(),
			);
			assert_eq!(
				r,
				Err(GevdError::NoConvergence),
				"real, n = {n}, {at:?}"
			);

			let (mut a, mut b) =
				small_pair_cplx(n, Some((at, f64::INFINITY)), None, 0x67);
			let params: GevdParams = <GevdParams as faer::Auto<c64>>::auto();
			let mut buf = MemBuffer::new(gevd_scratch::<c64>(
				n,
				ComputeEigenvectors::No,
				ComputeEigenvectors::No,
				Par::Seq,
				params.into(),
			));
			let mut s = Diag::<c64>::zeros(n);
			let mut beta = Diag::<c64>::zeros(n);
			let r = gevd_cplx(
				a.as_mut(),
				b.as_mut(),
				s.as_mut(),
				beta.as_mut(),
				None,
				None,
				Par::Seq,
				MemStack::new(&mut buf),
				params.into(),
			);
			assert_eq!(
				r,
				Err(GevdError::NoConvergence),
				"complex, n = {n}, {at:?}"
			);
		}
	}
}
