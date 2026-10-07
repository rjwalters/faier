//! faier fork regression test (codeberg sarah-quinones/faer#307).
//!
//! `faer::sparse::linalg::lu::factorize_symbolic_lu` takes a
//! `LuColOrdering`: `Colamd` (the previous, hardcoded behaviour), `Identity`,
//! or `Custom(order)`, where `order[k]` is the original column eliminated at
//! step `k` (new -> old). This checks, through the public API, that:
//!
//! - every ordering yields a correct `solve` and `solve_transpose`, on both the
//!   simplicial and the supernodal numeric paths (each forced through
//!   `supernodal_flop_ratio_threshold`, and checked to be the path taken);
//! - `SymbolicLu::col_perm()` returns a custom order unchanged as its forward
//!   array, with the inverse in the matching orientation (a non-involution
//!   order is used, so a flipped orientation would be caught);
//! - passing COLAMD's own order as `Custom` reproduces `Colamd` exactly;
//! - a malformed custom order panics with a clear message.

#![allow(non_snake_case)]

use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::linalg::SupernodalThreshold;
use faer::sparse::linalg::lu::{
	LuColOrdering, LuSymbolicParams, NumericLu, SymbolicLu,
	factorize_symbolic_lu,
};
use faer::sparse::{SparseColMat, Triplet};
use faer::{Conj, Mat, Par};

/// small deterministic LCG so the test needs no rng dependency
struct Lcg(u64);
impl Lcg {
	fn next_u64(&mut self) -> u64 {
		self.0 = self
			.0
			.wrapping_mul(6364136223846793005)
			.wrapping_add(1442695040888963407);
		self.0 >> 33
	}

	/// uniform in [-1, 1)
	fn next_f64(&mut self) -> f64 {
		(self.next_u64() as f64 / (1u64 << 31) as f64) * 2.0 - 1.0
	}

	fn below(&mut self, n: usize) -> usize {
		(self.next_u64() % n as u64) as usize
	}
}

/// a structurally unsymmetric, diagonally dominant `n x n` matrix: a 2D
/// 5-point grid stencil with unsymmetric values, plus random long-range
/// entries so the column elimination tree is non-trivial
fn test_matrix(nx: usize, ny: usize) -> SparseColMat<usize, f64> {
	let n = nx * ny;
	let mut rng = Lcg(0x5EED_F00D);
	let mut triplets = Vec::new();
	for y in 0..ny {
		for x in 0..nx {
			let i = y * nx + x;
			triplets.push(Triplet::new(i, i, 10.0 + rng.next_f64()));
			if x + 1 < nx {
				triplets.push(Triplet::new(i, i + 1, rng.next_f64()));
			}
			if y + 1 < ny {
				triplets.push(Triplet::new(i, i + nx, rng.next_f64()));
			}
			// only some of the transposed couplings, so the pattern is
			// unsymmetric
			if x > 0 && i % 3 != 0 {
				triplets.push(Triplet::new(i, i - 1, rng.next_f64()));
			}
			if y > 0 && i % 2 == 0 {
				triplets.push(Triplet::new(i, i - nx, rng.next_f64()));
			}
		}
	}
	for _ in 0..n {
		let i = rng.below(n);
		let j = rng.below(n);
		if i != j {
			triplets.push(Triplet::new(i, j, 0.5 * rng.next_f64()));
		}
	}
	SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
}

/// a pseudo-random permutation of `0..n` (Fisher-Yates); almost surely not an
/// involution
fn shuffled(n: usize, seed: u64) -> Vec<usize> {
	let mut rng = Lcg(seed);
	let mut p = (0..n).collect::<Vec<_>>();
	for i in (1..n).rev() {
		p.swap(i, rng.below(i + 1));
	}
	p
}

fn is_supernodal(symbolic: &SymbolicLu<usize>) -> bool {
	// the simplicial/supernodal choice is not otherwise exposed; the `Debug`
	// representation names the inner variant
	let repr = format!("{symbolic:?}");
	let supernodal = repr.contains("Supernodal");
	assert!(supernodal != repr.contains("Simplicial"));
	supernodal
}

/// factorizes `A` with the given ordering and threshold, and checks that
/// both `A x = b` and `A^T x = b` are solved to a small residual
fn check_solves(
	A: &SparseColMat<usize, f64>,
	ord: LuColOrdering<'_, usize>,
	threshold: SupernodalThreshold,
	expect_supernodal: bool,
) -> SymbolicLu<usize> {
	let n = A.nrows();
	let symbolic = factorize_symbolic_lu(
		A.symbolic(),
		ord,
		LuSymbolicParams {
			supernodal_flop_ratio_threshold: threshold,
			..Default::default()
		},
	)
	.unwrap();
	assert_eq!(is_supernodal(&symbolic), expect_supernodal);

	let mut numeric = NumericLu::<usize, f64>::new();
	let lu = symbolic
		.factorize_numeric_lu(
			&mut numeric,
			A.as_ref(),
			Par::Seq,
			MemStack::new(&mut MemBuffer::new(
				symbolic.factorize_numeric_lu_scratch::<f64>(
					Par::Seq,
					Default::default(),
				),
			)),
			Default::default(),
		)
		.unwrap();

	let mut rng = Lcg(42);
	let rhs = Mat::<f64>::from_fn(n, 3, |_, _| rng.next_f64());
	let dense = A.to_dense();

	let mut x = rhs.clone();
	lu.solve_in_place_with_conj(
		Conj::No,
		x.as_mut(),
		Par::Seq,
		MemStack::new(&mut MemBuffer::new(
			symbolic.solve_in_place_scratch::<f64>(rhs.ncols(), Par::Seq),
		)),
	);
	let residual = (&dense * &x - &rhs).norm_max();
	assert!(residual <= 1e-10, "solve residual {residual:e} ({ord:?})");

	let mut x = rhs.clone();
	lu.solve_transpose_in_place_with_conj(
		Conj::No,
		x.as_mut(),
		Par::Seq,
		MemStack::new(&mut MemBuffer::new(
			symbolic
				.solve_transpose_in_place_scratch::<f64>(rhs.ncols(), Par::Seq),
		)),
	);
	let residual = (dense.transpose() * &x - &rhs).norm_max();
	assert!(
		residual <= 1e-10,
		"solve_transpose residual {residual:e} ({ord:?})"
	);

	symbolic
}

fn assert_col_perm_is(symbolic: &SymbolicLu<usize>, new_to_old: &[usize]) {
	let (fwd, inv) = symbolic.col_perm().arrays();
	assert_eq!(fwd, new_to_old);
	for (new, &old) in new_to_old.iter().enumerate() {
		assert_eq!(inv[old], new);
	}
}

#[test]
fn every_lu_col_ordering_solves_on_simplicial_and_supernodal_paths() {
	let A = test_matrix(12, 15);
	let n = A.ncols();

	// cyclic shift and a random shuffle: neither is an involution
	let shift = (0..n).map(|k| (k + 1) % n).collect::<Vec<_>>();
	let shuffle = shuffled(n, 7);
	assert!((0..n).any(|k| shuffle[shuffle[k]] != k));

	for (threshold, supernodal) in [
		(SupernodalThreshold::FORCE_SIMPLICIAL, false),
		(SupernodalThreshold::FORCE_SUPERNODAL, true),
	] {
		let colamd =
			check_solves(&A, LuColOrdering::Colamd, threshold, supernodal);

		let identity =
			check_solves(&A, LuColOrdering::Identity, threshold, supernodal);
		assert_col_perm_is(&identity, &(0..n).collect::<Vec<_>>());

		for order in [&shift, &shuffle] {
			let custom = check_solves(
				&A,
				LuColOrdering::Custom(order),
				threshold,
				supernodal,
			);
			assert_col_perm_is(&custom, order);
		}

		// COLAMD's own order, supplied as `Custom`, reproduces `Colamd`
		let colamd_order = colamd.col_perm().arrays().0.to_vec();
		let custom = check_solves(
			&A,
			LuColOrdering::Custom(&colamd_order),
			threshold,
			supernodal,
		);
		assert_eq!(custom.col_perm().arrays(), colamd.col_perm().arrays());
	}
}

#[test]
fn default_lu_col_ordering_is_colamd() {
	let A = test_matrix(6, 7);
	let params = LuSymbolicParams::default();
	let default =
		factorize_symbolic_lu(A.symbolic(), LuColOrdering::default(), params)
			.unwrap();
	let colamd =
		factorize_symbolic_lu(A.symbolic(), LuColOrdering::Colamd, params)
			.unwrap();
	assert_eq!(default.col_perm().arrays(), colamd.col_perm().arrays());
}

#[test]
#[should_panic(expected = "appears more than once")]
fn custom_lu_col_ordering_rejects_duplicate() {
	let A = test_matrix(3, 3);
	let mut order = (0..A.ncols()).collect::<Vec<_>>();
	order[4] = order[2];
	let _ = factorize_symbolic_lu(
		A.symbolic(),
		LuColOrdering::Custom(&order),
		Default::default(),
	);
}

#[test]
#[should_panic(expected = "expected 9")]
fn custom_lu_col_ordering_rejects_wrong_length() {
	let A = test_matrix(3, 3);
	let order = (0..A.ncols() - 1).collect::<Vec<_>>();
	let _ = factorize_symbolic_lu(
		A.symbolic(),
		LuColOrdering::Custom(&order),
		Default::default(),
	);
}

#[test]
#[should_panic(expected = "out of range")]
fn custom_lu_col_ordering_rejects_out_of_range() {
	let A = test_matrix(3, 3);
	let mut order = (0..A.ncols()).collect::<Vec<_>>();
	order[0] = A.ncols();
	let _ = factorize_symbolic_lu(
		A.symbolic(),
		LuColOrdering::Custom(&order),
		Default::default(),
	);
}
