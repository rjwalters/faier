//! faier fork regression test (geode-fem #908 / #867).
//!
//! A non-finite QZ iterate never satisfies a deflation test, and neither the
//! blocked nor the unblocked QZ loop checked for one, so a NaN that entered
//! the iteration (for example from the `make_givens` underflow, #867 item 1)
//! made the loop run out its whole `30 n` sweep budget before returning
//! meaningless eigenvalues: 14.8 h on a 3300 x 3300 pencil. The loops now stop
//! at the first non-finite iterate and report NaN eigenvalues, which
//! `gevd_real` / `gevd_cplx` turn into `GevdError::NoConvergence`.

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::evd::ComputeEigenvectors;
use faer::linalg::gevd::{GeneralizedSchurParams, qz_cplx, qz_real};
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

/// runs the real QZ on an `n x n` Hessenberg / triangular pair whose `(k, k)`
/// entry of `A` is `poison`, returning `(alphar, alphai, beta)`
fn real_qz_small(
	n: usize,
	k: usize,
	poison: f64,
	seed: u64,
) -> (Col<f64>, Col<f64>, Col<f64>) {
	let mut rng = Lcg(seed);
	let mut a = Mat::<f64>::from_fn(n, n, |i, j| {
		if i == k && j == k {
			poison
		} else if i <= j + 1 {
			rng.next()
		} else {
			0.0
		}
	});
	let mut b = Mat::<f64>::from_fn(n, n, |i, j| {
		if i == j {
			2.0 + rng.next()
		} else if i < j {
			rng.next()
		} else {
			0.0
		}
	});
	let params: GeneralizedSchurParams =
		<GeneralizedSchurParams as faer::Auto<f64>>::auto();
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

/// below the blocked threshold the unblocked QZ runs, and its NaN fail-fast
/// fills `alphai` with NaN; the conjugate-pair fix-up used to read
/// `NaN != 0` as a pair start and step past the end of the block for odd
/// sizes (index out of bounds). even sizes are the control.
///
/// an `inf` at `(0, 0)` is left out: it never reaches the trailing diagonal
/// the unblocked fail-fast inspects (rjwalters/faier#6)
#[test]
fn real_unblocked_qz_nonfinite_does_not_overrun_small_odd_blocks() {
	for n in [3usize, 4, 5, 7, 8] {
		for (k, poison) in
			[(0, f64::NAN), (n - 1, f64::NAN), (n - 1, f64::INFINITY)]
		{
			let (alphar, alphai, beta) =
				real_qz_small(n, k, poison, 0x90a + n as u64);
			assert!(
				(0..n).any(|i| !(alphar[i].is_finite()
					&& alphai[i].is_finite()
					&& beta[i].is_finite())),
				"n = {n}, poison {poison} at ({k}, {k}): a breakdown must be \
				 reported through non-finite eigenvalues"
			);
		}
	}
}
