use crate::internal_prelude::*;
use crate::utils;
use equator::assert;
use linalg::matmul::{matmul, triangular as tr};
/// computes the layout of the workspace required to compute a matrix pair's
/// generalized hessenberg decomposition
pub fn generalized_hessenberg_scratch<T: ComplexField>(
	n: usize,
	params: GeneralizedHessenbergParams,
) -> StackReq {
	if params.block_size <= 1 || n <= params.block_size {
		StackReq::EMPTY
	} else {
		linalg::temp_mat_scratch::<T>(
			n,
			match params.block_size.checked_mul(6) {
				Some(n) => n,
				None => return StackReq::OVERFLOW,
			},
		)
	}
}
/// computes a plane rotation `[c, s; -conj(s), c]` with real `c >= 0` such that
/// `[c, s; -conj(s), c] [f; g] = [r; 0]`, i.e. `c = |f| / h`,
/// `s = sgn(f) conj(g) / h` and `r = sgn(f) h`, with `h = hypot(|f|, |g|)` and
/// `sgn(f) = f / |f|` (`c = 0`, `s = conj(g) / |g|`, `r = |g|` when `f = 0`)
pub(crate) fn make_givens<T: ComplexField>(f: T, g: T) -> (T::Real, T, T) {
	if g == zero() {
		(one(), zero(), f)
	} else if f == zero() {
		let c = zero::<T::Real>();
		let d = g.abs();
		let d_inv = d.recip();
		if !(d.is_finite() && d_inv.is_finite()) {
			// `|g|` is far down in the subnormal range, or overflows
			return lartg(f, g);
		}
		let s = g.conj().mul_real(&d_inv);
		let r = d.to_cplx();
		(c, s, r)
	} else {
		// the direct (upstream) form below divides by `h` twice and by `c`.
		// `1 / h` and `1 / c` overflow once `h` or `c` drops below
		// `1 / max_positive` (subnormal inputs, or `|f| << |g|`), which made
		// `s` and `r` non-finite and turned the QZ iterates into NaN. it is
		// kept, bit for bit, where `h` and `c` are safely in range; otherwise
		// the rotation is computed by the safe-scaling algorithm of lapack
		// `xlartg`
		let f1 = f.abs();
		let g1 = g.abs();
		let h1 = f1.hypot(g1);
		let rtmin = sqrt_min_positive::<T::Real>();
		let rtmax = sqrt_max_positive::<T::Real>();
		if h1 > rtmin && h1 < rtmax {
			let h1_inv = &h1.recip();
			let c = f1 * h1_inv;
			if c > rtmin {
				let r = f.mul_real(c.recip());
				let s = g.conj() * r.mul_real(h1_inv).mul_real(h1_inv);
				return (c, s, r);
			}
		}
		lartg(f, g)
	}
}
/// `max(|re(z)|, |im(z)|)`
fn abs_max<T: ComplexField>(z: &T) -> T::Real {
	z.real().abs().fmax(z.imag().abs())
}
/// `z / x` for a real `x` in `[min_positive, max_positive]`, so `1 / x` is
/// finite
fn div_real<T: ComplexField>(z: &T, x: &T::Real) -> T {
	z.mul_real(x.recip())
}
/// plane rotation with the same convention as [`make_givens`], computed by
/// the safe-scaling algorithm of lapack `zlartg` (Anderson, "Algorithm 978:
/// safe scaling in the level 1 BLAS", 2017; lapack 3.10). every intermediate
/// stays in `[min_positive, max_positive]`, and `f` and `g` are scaled
/// separately when `|f| << |g|`, so `sgn(f)` keeps full precision for any
/// finite `f` and `g`, subnormal ones included. for real `T` this reduces to
/// lapack `dlartg`
fn lartg<T: ComplexField>(f: T, g: T) -> (T::Real, T, T) {
	let safmin = min_positive::<T::Real>();
	let safmax = max_positive::<T::Real>();
	let rtmin = sqrt_min_positive::<T::Real>();
	let two = from_f64::<T::Real>(2.0);
	let four = from_f64::<T::Real>(4.0);
	let clamp = |x: T::Real| safmin.fmax(&x).fmin(&safmax);

	if g == zero() {
		return (one(), zero(), f);
	}
	if f == zero() {
		let c = zero::<T::Real>();
		let g1 = abs_max(&g);
		let rtmax = (&safmax / &two).sqrt();
		let (s, d) = if g1 > rtmin && g1 < rtmax {
			let d = g.abs2().sqrt();
			(div_real(&g.conj(), &d), d)
		} else {
			let u = clamp(g1);
			let gs = div_real(&g, &u);
			let d = gs.abs2().sqrt();
			(div_real(&gs.conj(), &d), d * &u)
		};
		return (c, s, d.to_cplx());
	}

	let f1 = abs_max(&f);
	let g1 = abs_max(&g);
	let rtmax = (&safmax / &four).sqrt();
	// `(fs, gs)` = `(f / (w u), g / u)`, with `w = 1` unless `|f| << |g|`
	let (fs, gs, u, w) = if f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax
	{
		(f, g, one::<T::Real>(), one::<T::Real>())
	} else {
		let u = clamp(f1.fmax(&g1));
		let gs = div_real(&g, &u);
		if &f1 / &u < rtmin {
			// `f` is not well scaled by `u`: scale it by its own magnitude
			let v = clamp(f1);
			let w = &v / &u;
			(div_real(&f, &v), gs, u, w)
		} else {
			(div_real(&f, &u), gs, u, one::<T::Real>())
		}
	};
	let f2 = fs.abs2();
	let g2 = gs.abs2();
	let h2 = &(&(&f2 * &w) * &w) + &g2;
	// `min_positive <= f2 <= h2 <= max_positive`
	let (c, r, s);
	if f2 >= &h2 * &safmin {
		// `min_positive <= f2 / h2 <= 1`, and `h2 / f2` is finite
		c = (&f2 / &h2).sqrt();
		r = div_real(&fs, &c);
		if f2 > rtmin && h2 < &rtmax * &two {
			// `min_positive <= sqrt(f2 h2) <= max_positive`
			s = gs.conj() * div_real(&fs, &(&f2 * &h2).sqrt());
		} else {
			s = gs.conj() * div_real(&r, &h2);
		}
	} else {
		// `f2 / h2 <= min_positive` may be subnormal, and `h2 / f2` may
		// overflow, but `sqrt(min_positive) <= sqrt(f2 h2) <= sqrt(max_positive)`
		let d = (&f2 * &h2).sqrt();
		c = &f2 / &d;
		if c >= safmin {
			r = div_real(&fs, &c);
		} else {
			r = fs.mul_real(&h2 / &d);
		}
		s = gs.conj() * div_real(&fs, &d);
	}
	(c * &w, s, r.mul_real(&u))
}
/// returns `true` if the diagonal and subdiagonal of `A` and the diagonal of
/// `B` are finite over the active block `istart..=istop`. a non-finite iterate
/// never satisfies a deflation test, so the qz loops use this to stop early
/// instead of running out their `30 n` iteration cap
pub(crate) fn active_block_is_finite<T: ComplexField>(
	A: MatRef<'_, T>,
	B: MatRef<'_, T>,
	istart: usize,
	istop: usize,
) -> bool {
	(istart..istop + 1).all(|k| {
		A[(k, k)].is_finite()
			&& B[(k, k)].is_finite()
			&& (k == istart || A[(k, k - 1)].is_finite())
	})
}
pub(crate) fn rot<T: ComplexField>(
	c: T::Real,
	s: T,
	x: RowMut<'_, T>,
	y: RowMut<'_, T>,
) {
	let (c, s) = &(c, s);
	zip!(x, y).for_each(|unzip!(x, y): Zip!(&mut _, &mut _)| {
		(*x, *y) = (x.mul_real(c) + &*y * s, y.mul_real(c) - &*x * s.conj());
	});
}
/// generalized hessenberg factorization tuning parameters
#[derive(Copy, Clone, Debug)]
pub struct GeneralizedHessenbergParams {
	/// algorithm blocking parameter
	pub block_size: usize,
	/// threshold at which blocking should be disabled
	pub blocking_threshold: usize,
	#[doc(hidden)]
	pub non_exhaustive: NonExhaustive,
}
impl<T: ComplexField> Auto<T> for GeneralizedHessenbergParams {
	fn auto() -> Self {
		Self {
			block_size: 32,
			blocking_threshold: 256 * 256,
			non_exhaustive: NonExhaustive(()),
		}
	}
}
pub(crate) fn trot<T: ComplexField>(
	c: T::Real,
	s: T,
	x: ColMut<'_, T>,
	y: ColMut<'_, T>,
) {
	rot(c, -s.conj(), x.transpose_mut(), y.transpose_mut());
}
fn generalized_hessenberg_unblocked<T: ComplexField>(
	A: MatMut<'_, T>,
	B: MatMut<'_, T>,
	Q: Option<MatMut<'_, T>>,
	Z: Option<MatMut<'_, T>>,
	Q_is_I: bool,
	Z_is_I: bool,
	par: Par,
	stack: &mut MemStack,
	params: GeneralizedHessenbergParams,
) {
	_ = params;
	_ = stack;
	_ = Q_is_I;
	_ = Z_is_I;
	_ = par;
	let mut A = A;
	let mut B = B;
	let mut Q = Q;
	let mut Z = Z;
	let n = A.nrows();
	let (Q_nrows, Q_ncols) =
		Q.rb().map(|Q| (Q.nrows(), Q.ncols())).unwrap_or((n, n));
	let (Z_nrows, Z_ncols) =
		Z.rb().map(|Z| (Z.nrows(), Z.ncols())).unwrap_or((n, n));
	assert!(all(
		A.nrows() == n,
		A.ncols() == n,
		B.nrows() == n,
		B.ncols() == n,
		Q_nrows == n,
		Q_ncols == n,
		Z_nrows == n,
		Z_ncols == n,
	));
	if n <= 2 {
		return;
	}
	for jcol in 0..n - 2 {
		for jrow in (jcol + 2..n).rev() {
			let (c, s, r) =
				make_givens(A[(jrow - 1, jcol)].copy(), A[(jrow, jcol)].copy());
			A[(jrow - 1, jcol)] = r;
			A[(jrow, jcol)] = zero();
			let (x, y) = A
				.rb_mut()
				.get_mut(.., jcol + 1..)
				.two_rows_mut(jrow - 1, jrow);
			rot(c.copy(), s.copy(), x, y);
			let (x, y) = B
				.rb_mut()
				.get_mut(.., jrow - 1..)
				.two_rows_mut(jrow - 1, jrow);
			rot(c.copy(), s.copy(), x, y);
			if let Some(mut Q) = Q.rb_mut() {
				let (x, y) = Q.rb_mut().two_cols_mut(jrow - 1, jrow);
				rot(c.copy(), s.conj(), x.transpose_mut(), y.transpose_mut());
			}
			let (c, s, r) =
				make_givens(B[(jrow, jrow)].copy(), B[(jrow, jrow - 1)].copy());
			B[(jrow, jrow)] = r;
			B[(jrow, jrow - 1)] = zero();
			let (x, y) = A.rb_mut().two_cols_mut(jrow, jrow - 1);
			rot(c.copy(), s.copy(), x.transpose_mut(), y.transpose_mut());
			let (x, y) =
				B.rb_mut().get_mut(..jrow, ..).two_cols_mut(jrow, jrow - 1);
			rot(c.copy(), s.copy(), x.transpose_mut(), y.transpose_mut());
			if let Some(mut Z) = Z.rb_mut() {
				let (x, y) = Z.rb_mut().two_cols_mut(jrow, jrow - 1);
				rot(c, s, x.transpose_mut(), y.transpose_mut());
			}
		}
	}
}
fn apply_U_from_the_right<T: ComplexField>(
	M: MatMut<'_, T>,
	U: MatRef<'_, T>,
	par: Par,
	stack: &mut MemStack,
) {
	let n = U.nrows() / 2;
	let (U11, U12, U21, U22) = U.split_at(n, n);
	let mut M = M;
	let (mut tmp, _) = unsafe {
		linalg::temp_mat_uninit::<T, _, _>(M.nrows(), M.ncols(), stack)
	};
	let mut tmp = tmp.as_mat_mut();
	let (mut T1, mut T2) = tmp.rb_mut().split_at_col_mut(n);
	let (M1, M2) = M.rb().split_at_col(n);
	utils::thread::join_raw(
		|par| {
			matmul(T1.rb_mut(), Accum::Replace, M1.rb(), U11, one(), par);
			tr::matmul(
				T1.rb_mut(),
				tr::BlockStructure::Rectangular,
				Accum::Add,
				M2.rb(),
				tr::BlockStructure::Rectangular,
				U21,
				tr::BlockStructure::TriangularUpper,
				one(),
				par,
			);
		},
		|par| {
			tr::matmul(
				T2.rb_mut(),
				tr::BlockStructure::Rectangular,
				Accum::Replace,
				M1.rb(),
				tr::BlockStructure::Rectangular,
				U12,
				tr::BlockStructure::TriangularLower,
				one(),
				par,
			);
			matmul(T2.rb_mut(), Accum::Add, M2.rb(), U22, one(), par);
		},
		par,
	);
	M.copy_from(tmp);
}
/// computes a matrix pair $(A, B)$'s generalized hessenberg decomposition such
/// that
/// - B is an upper triangular matrix
/// - $A = Q H Z^H, A = Q T Z^H$,
/// - $H$ is a hessenberg matrix stored in the upper triangular half of $A$
///   (plus the subdiagonal),
/// - $T$ is an upper triangular matrix,
/// - $Q$ and $Z$ are unitary matrices.
///
/// # warning
/// $B$ is assumed to be upper triangular on input.
///
/// $Q$ and $Z$ are postmultiplied into the input-output parameters `Q` and `Z`.
///
/// i.e.: $Q_{\text{out}} = Q_{\text{in}} * Q$ and $Z_{\text{out}} =
/// Z_{\text{in}} * Z$.
///
/// if this behavior is not desired then $Q$ and $Z$ should be overwritten by
/// the identity matrix before calling this function.
pub fn generalized_hessenberg<T: ComplexField>(
	A: MatMut<'_, T>,
	B: MatMut<'_, T>,
	Q_inout: Option<MatMut<'_, T>>,
	Z_inout: Option<MatMut<'_, T>>,
	par: Par,
	stack: &mut MemStack,
	params: GeneralizedHessenbergParams,
) {
	let n = A.nrows();
	let (Q_nrows, Q_ncols) = Q_inout
		.rb()
		.map(|Q| (Q.nrows(), Q.ncols()))
		.unwrap_or((n, n));
	let (Z_nrows, Z_ncols) = Z_inout
		.rb()
		.map(|Z| (Z.nrows(), Z.ncols()))
		.unwrap_or((n, n));
	let Q_is_I = Q_inout.rb().is_some_and(|Q| Q.is_identity());
	let Z_is_I = Z_inout.rb().is_some_and(|Z| Z.is_identity());
	assert!(all(
		A.nrows() == n,
		A.ncols() == n,
		B.nrows() == n,
		B.ncols() == n,
		Q_nrows == n,
		Q_ncols == n,
		Z_nrows == n,
		Z_ncols == n,
	));
	if n <= 2 {
		return;
	}
	if params.block_size <= 1 || A.nrows() <= params.block_size {
		return generalized_hessenberg_unblocked(
			A, B, Q_inout, Z_inout, Q_is_I, Z_is_I, par, stack, params,
		);
	}
	let mut A = A;
	let mut B = B;
	let mut Q = Q_inout;
	let mut Z = Z_inout;
	let mut jcol = 0;
	while jcol < n - 2 {
		let nnb = params.block_size.clamp(1, n - 2 - jcol);
		let top = if jcol < 2 { 0 } else { jcol + 1 };
		let n2nb = ((n - 2 - jcol) / nnb).saturating_sub(1);
		let nblst = (n - 1 - jcol) - n2nb * nnb;
		alloca!('stack: {
			let mut work0 = unsafe { mat![uninit::<T>, nblst, nblst] };
			let mut work1 =
				unsafe { mat![uninit::<T>, 2 * nnb, 2 * nnb * n2nb] };
		});
		work0.fill(zero());
		work0.rb_mut().diagonal_mut().fill(one());
		work1.fill(zero());
		for i in 0..n2nb {
			work1
				.rb_mut()
				.subcols_mut(2 * nnb * i, 2 * nnb)
				.diagonal_mut()
				.fill(one());
		}
		for j in jcol..jcol + nnb {
			for i in (j + 2..n).rev() {
				let (c, s, r) =
					make_givens(A[(i - 1, j)].copy(), A[(i, j)].copy());
				A[(i - 1, j)] = r;
				A[(i, j)] = c.to_cplx();
				B[(i, j)] = s;
			}
			let mut len = 2 + j - jcol;
			let jrow = j + n2nb * nnb + 2;
			for i in (jrow..n).rev() {
				let ref c = A[(i, j)].real();
				let ref s = B[(i, j)].copy();
				let col = nblst - (n - i) - 1;
				for jj in nblst - len..nblst {
					let ref temp0 = work0[(jj, col)].copy();
					let ref temp1 = work0[(jj, col + 1)].copy();
					work0[(jj, col + 1)] = temp1.mul_real(c) - s * temp0;
					work0[(jj, col)] = temp0.mul_real(c) + s.conj() * temp1;
				}
				len += 1;
			}
			let mut j0 = jrow;
			let mut idx = 0usize;
			while j0 > j + 2 {
				j0 -= nnb;
				let mut start = nnb;
				let mut len = 2 + j - jcol;
				let mut col = nnb + j - jcol;
				let mut U = work1.rb_mut().submatrix_mut(
					0,
					idx * 2 * nnb,
					2 * nnb,
					2 * nnb,
				);
				for i in (j0..j0 + nnb).rev() {
					let ref c = A[(i, j)].real();
					let ref s = B[(i, j)].copy();
					col -= 1;
					start -= 1;
					for jj in start..start + len {
						let ref temp0 = U[(jj, col)].copy();
						let ref temp1 = U[(jj, col + 1)].copy();
						U[(jj, col + 1)] = temp1.mul_real(c) - s * temp0;
						U[(jj, col)] = s.conj() * temp1 + temp0.mul_real(c);
					}
					len += 1;
				}
				idx += 1;
			}
			for jj in (j + 1..n).rev() {
				for i in (j + 2..Ord::min(n, jj + 2)).rev() {
					let ref c = A[(i, j)].real();
					let ref s = B[(i, j)].copy();
					let ref temp0 = B[(i - 1, jj)].copy();
					let ref temp1 = B[(i, jj)].copy();
					B[(i, jj)] = temp1.mul_real(c) - s.conj() * temp0;
					B[(i - 1, jj)] = temp0.mul_real(c) + s * temp1;
				}
				if jj + 1 < n {
					let (c, s, r) = make_givens(
						B[(jj + 1, jj + 1)].copy(),
						B[(jj + 1, jj)].copy(),
					);
					B[(jj + 1, jj + 1)] = r;
					B[(jj + 1, jj)] = zero();
					let (bjj, bjj1) = B
						.rb_mut()
						.get_mut(top..jj + 1, ..)
						.two_cols_mut(jj, jj + 1);
					rot(
						c.copy(),
						s.copy(),
						bjj1.transpose_mut(),
						bjj.transpose_mut(),
					);
					A[(jj + 1, j)] = c.to_cplx();
					B[(jj + 1, j)] = -s.conj();
				}
			}
			let jj = (n - 2 - j) % 3;
			let mut i = n - 1 - j;
			while i > jj + 1 {
				i -= 3;
				let ref c0 = A[(j + i + 1, j)].real();
				let ref s0 = -&B[(j + i + 1, j)];
				let ref c1 = A[(j + i + 2, j)].real();
				let ref s1 = -&B[(j + i + 2, j)];
				let ref c2 = A[(j + i + 3, j)].real();
				let ref s2 = -&B[(j + i + 3, j)];
				for k in top..n {
					let ref temp0 = A[(k, j + i + 0)].copy();
					let ref temp1 = A[(k, j + i + 1)].copy();
					let ref temp2 = A[(k, j + i + 2)].copy();
					let ref temp3 = A[(k, j + i + 3)].copy();
					A[(k, j + i + 3)] = temp3.mul_real(c2) + temp2 * s2.conj();
					let temp2 = temp2.mul_real(c2) - temp3 * s2;
					A[(k, j + i + 2)] = temp2.mul_real(c1) + temp1 * s1.conj();
					let temp1 = temp1.mul_real(c1) - temp2 * s1;
					A[(k, j + i + 1)] = temp1.mul_real(c0) + temp0 * s0.conj();
					let temp0 = temp0.mul_real(c0) - temp1 * s0;
					A[(k, j + i + 0)] = temp0;
				}
			}
			for i in (1..jj + 1).rev() {
				let c = A[(j + i + 1, j)].real();
				let s = B[(j + i + 1, j)].copy();
				let (aj1, aj) = A
					.rb_mut()
					.get_mut(top.., ..)
					.two_cols_mut(j + i + 1, j + i);
				trot(c, s, aj1, aj);
			}
			if j < jcol + nnb - 1 {
				let len = j - jcol + 1;
				let jrow = n - nblst;
				{
					let (mut work2, _) = unsafe {
						linalg::temp_mat_uninit::<T, _, _>(nblst, 1, stack)
					};
					let mut work2 = work2.as_mat_mut().col_mut(0);
					matmul(
						work2.rb_mut().get_mut(..len),
						Accum::Replace,
						work0.rb().get(.., ..len).adjoint(),
						A.rb().get(jrow.., j + 1),
						one(),
						par,
					);
					tr::matmul(
						work2.rb_mut().get_mut(len..),
						tr::BlockStructure::Rectangular,
						Accum::Replace,
						work0.rb().get(..nblst - len, len..).adjoint(),
						tr::BlockStructure::TriangularUpper,
						A.rb().get(jrow..jrow + nblst - len, j + 1),
						tr::BlockStructure::Rectangular,
						one(),
						par,
					);
					matmul(
						work2.rb_mut().get_mut(len..),
						Accum::Add,
						work0.rb().get(nblst - len.., len..).adjoint(),
						A.rb().get(jrow + nblst - len.., j + 1),
						one(),
						par,
					);
					for i in jrow..jrow + nblst {
						A[(i, j + 1)] = work2[i - jrow].copy();
					}
				}
				let mut j0 = jrow;
				let mut idx = 0usize;
				while j0 > jcol + nnb {
					j0 -= nnb;
					let U = work1.rb().submatrix(
						0,
						idx * 2 * nnb,
						nnb + len,
						nnb + len,
					);
					let (U11, U12, U21, U22) = U.split_at(nnb, len);
					let (mut work2, _) = unsafe {
						linalg::temp_mat_uninit::<T, _, _>(nnb + len, 1, stack)
					};
					let mut work2 = work2.as_mat_mut().col_mut(0);
					matmul(
						work2.rb_mut().get_mut(..len),
						Accum::Replace,
						U11.adjoint(),
						A.rb().get(j0..j0 + nnb, j + 1),
						one(),
						par,
					);
					tr::matmul(
						work2.rb_mut().get_mut(..len),
						tr::BlockStructure::Rectangular,
						Accum::Add,
						U21.adjoint(),
						tr::BlockStructure::TriangularLower,
						A.rb().get(j0 + nnb..j0 + nnb + len, j + 1),
						tr::BlockStructure::Rectangular,
						one(),
						par,
					);
					tr::matmul(
						work2.rb_mut().get_mut(len..),
						tr::BlockStructure::Rectangular,
						Accum::Replace,
						U12.adjoint(),
						tr::BlockStructure::TriangularUpper,
						A.rb().get(j0..j0 + nnb, j + 1),
						tr::BlockStructure::Rectangular,
						one(),
						par,
					);
					matmul(
						work2.rb_mut().get_mut(len..),
						Accum::Add,
						U22.adjoint(),
						A.rb().get(j0 + nnb..j0 + nnb + len, j + 1),
						one(),
						par,
					);
					for i in j0..j0 + len + nnb {
						A[(i, j + 1)] = work2[i - j0].copy();
					}
					idx += 1;
				}
			}
		}
		let cola = n - jcol - nnb;
		let j = n - nblst;
		{
			let (mut work2, _) = unsafe {
				linalg::temp_mat_uninit::<T, _, _>(nblst, cola, stack)
			};
			let mut work2 = work2.as_mat_mut();
			matmul(
				work2.rb_mut(),
				Accum::Replace,
				work0.rb().adjoint(),
				A.rb().get(j.., jcol + nnb..),
				one(),
				par,
			);
			A.rb_mut().get_mut(j.., jcol + nnb..).copy_from(&work2);
		}
		let mut j0 = j;
		let mut idx = 0usize;
		while j0 > jcol + nnb {
			j0 -= nnb;
			let U = work1.rb().submatrix(0, idx * 2 * nnb, 2 * nnb, 2 * nnb);
			let (U11, U12, U21, U22) = U.split_at(nnb, nnb);
			let (mut work2, _) = unsafe {
				linalg::temp_mat_uninit::<T, _, _>(2 * nnb, cola, stack)
			};
			let mut work2 = work2.as_mat_mut();
			let (mut work2_top, mut work2_bot) =
				work2.rb_mut().split_at_row_mut(nnb);
			utils::thread::join_raw(
				|par| {
					matmul(
						work2_top.rb_mut(),
						Accum::Replace,
						U11.adjoint(),
						A.rb().get(j0..j0 + nnb, jcol + nnb..),
						one(),
						par,
					);
					tr::matmul(
						work2_top.rb_mut(),
						tr::BlockStructure::Rectangular,
						Accum::Add,
						U21.adjoint(),
						tr::BlockStructure::TriangularLower,
						A.rb().get(j0 + nnb..j0 + 2 * nnb, jcol + nnb..),
						tr::BlockStructure::Rectangular,
						one(),
						par,
					);
				},
				|par| {
					tr::matmul(
						work2_bot.rb_mut(),
						tr::BlockStructure::Rectangular,
						Accum::Replace,
						U12.adjoint(),
						tr::BlockStructure::TriangularUpper,
						A.rb().get(j0..j0 + nnb, jcol + nnb..),
						tr::BlockStructure::Rectangular,
						one(),
						par,
					);
					matmul(
						work2_bot.rb_mut(),
						Accum::Add,
						U22.adjoint(),
						A.rb().get(j0 + nnb..j0 + 2 * nnb, jcol + nnb..),
						one(),
						par,
					);
				},
				par,
			);
			A.rb_mut()
				.get_mut(j0..j0 + 2 * nnb, jcol + nnb..)
				.copy_from(&work2);
			idx += 1;
		}
		if let Some(mut Q) = Q.rb_mut() {
			let topq;
			let nh;
			if Q_is_I {
				topq = Ord::max(1, j - jcol);
				nh = n - topq;
			} else {
				topq = 0;
				nh = n;
			}
			{
				let (mut work2, _) = unsafe {
					linalg::temp_mat_uninit::<T, _, _>(nh, nblst, stack)
				};
				let mut work2 = work2.as_mat_mut();
				matmul(
					work2.rb_mut(),
					Accum::Replace,
					Q.rb().get(topq.., j..),
					work0.rb(),
					one(),
					par,
				);
				Q.rb_mut().get_mut(topq.., j..).copy_from(&work2);
			}
			let mut j0 = j;
			let mut idx = 0usize;
			while j0 > jcol + nnb {
				j0 -= nnb;
				let topq;
				if Q_is_I {
					topq = Ord::max(1, j0 - jcol);
				} else {
					topq = 0;
				}
				let U =
					work1.rb().submatrix(0, idx * 2 * nnb, 2 * nnb, 2 * nnb);
				apply_U_from_the_right(
					Q.rb_mut().get_mut(topq.., j0..j0 + 2 * nnb),
					U,
					par,
					stack,
				);
				idx += 1;
			}
		}
		if Z.is_some() || top > 0 {
			work0.fill(zero());
			work0.rb_mut().diagonal_mut().fill(one());
			work1.fill(zero());
			for i in 0..n2nb {
				work1
					.rb_mut()
					.subcols_mut(2 * nnb * i, 2 * nnb)
					.diagonal_mut()
					.column_vector_mut()
					.fill(one());
			}
			for j in jcol..jcol + nnb {
				let mut len = 2 + j - jcol;
				let jrow = j + n2nb * nnb + 2;
				for i in (jrow..n).rev() {
					let ref c = A[(i, j)].real();
					let ref s = B[(i, j)].copy();
					A[(i, j)] = zero();
					B[(i, j)] = zero();
					let col = nblst - (n - i) - 1;
					for jj in nblst - len..nblst {
						let ref temp0 = work0[(jj, col)].copy();
						let ref temp1 = work0[(jj, col + 1)].copy();
						work0[(jj, col + 1)] =
							temp1.mul_real(c) - temp0 * s.conj();
						work0[(jj, col)] = temp0.mul_real(c) + temp1 * s;
					}
					len += 1;
				}
				let mut j0 = jrow;
				let mut idx = 0usize;
				while j0 > j + 2 {
					j0 -= nnb;
					let mut len = 2 + j - jcol;
					let mut col = nnb + j - jcol;
					let mut start = nnb;
					let mut U = work1.rb_mut().submatrix_mut(
						0,
						idx * 2 * nnb,
						2 * nnb,
						2 * nnb,
					);
					for i in (j0..j0 + nnb).rev() {
						let ref c = A[(i, j)].real();
						let ref s = B[(i, j)].copy();
						A[(i, j)] = zero();
						B[(i, j)] = zero();
						col -= 1;
						start -= 1;
						for jj in start..start + len {
							let ref temp0 = U[(jj, col)].copy();
							let ref temp1 = U[(jj, col + 1)].copy();
							U[(jj, col + 1)] =
								temp1.mul_real(c) - temp0 * s.conj();
							U[(jj, col)] = temp0.mul_real(c) + temp1 * s;
						}
						len += 1;
					}
					idx += 1;
				}
			}
		} else {
			// clear the rotations stored below the first subdiagonal of the
			// panel's `nnb` columns (lapack: `xlaset('lower', n - jcol - 2,
			// nnb, ..., a(jcol + 2, jcol))`). this previously cleared the whole
			// rectangle `[jcol + 2.., jcol..n - jcol - 2]`, wiping live entries
			// of the reduced pair whenever `Z` was not requested
			for k in 0..nnb {
				zip!(A.rb_mut().get_mut(jcol + 2 + k.., jcol + k))
					.for_each(|unzip!(x)| *x = zero());
				zip!(B.rb_mut().get_mut(jcol + 2 + k.., jcol + k))
					.for_each(|unzip!(x)| *x = zero());
			}
		}
		if top > 0 {
			for mut M in [A.rb_mut(), B.rb_mut()] {
				let j = n - nblst;
				{
					let (mut work2, _) = unsafe {
						linalg::temp_mat_uninit::<T, _, _>(top, nblst, stack)
					};
					let mut work2 = work2.as_mat_mut();
					matmul(
						work2.rb_mut(),
						Accum::Replace,
						M.rb().get(..top, j..),
						work0.rb(),
						one(),
						par,
					);
					M.rb_mut().get_mut(..top, j..).copy_from(&work2);
				}
				let mut j0 = j;
				let mut idx = 0usize;
				while j0 > jcol + nnb {
					j0 -= nnb;
					let U = work1.rb().submatrix(
						0,
						idx * 2 * nnb,
						2 * nnb,
						2 * nnb,
					);
					apply_U_from_the_right(
						M.rb_mut().get_mut(..top, j0..j0 + 2 * nnb),
						U,
						par,
						stack,
					);
					idx += 1;
				}
			}
		}
		if let Some(mut Z) = Z.rb_mut() {
			let topq;
			let nh;
			if Z_is_I {
				topq = Ord::max(1, j - jcol);
				nh = n - topq;
			} else {
				topq = 0;
				nh = n;
			}
			{
				let (mut work2, _) = unsafe {
					linalg::temp_mat_uninit::<T, _, _>(nh, nblst, stack)
				};
				let mut work2 = work2.as_mat_mut();
				matmul(
					work2.rb_mut(),
					Accum::Replace,
					Z.rb().get(topq.., j..),
					work0.rb(),
					one(),
					par,
				);
				Z.rb_mut().get_mut(topq.., j..).copy_from(&work2);
			}
			let mut j0 = j;
			let mut idx = 0usize;
			while j0 > jcol + nnb {
				j0 -= nnb;
				let topq;
				if Z_is_I {
					topq = Ord::max(1, j0 - jcol);
				} else {
					topq = 0;
				}
				let U =
					work1.rb().submatrix(0, idx * 2 * nnb, 2 * nnb, 2 * nnb);
				apply_U_from_the_right(
					Z.rb_mut().get_mut(topq.., j0..j0 + 2 * nnb),
					U,
					par,
					stack,
				);
				idx += 1;
			}
		}
		jcol += nnb;
	}
}
#[cfg(test)]
mod tests {
	use super::*;
	use crate::{linalg, stats};
	use dyn_stack::MemBuffer;
	use equator::assert;
	use stats::prelude::*;
	#[test]
	fn test_givens() {
		let rng = &mut StdRng::seed_from_u64(0);
		let rand = ComplexDistribution::new(StandardUniform, StandardUniform);
		let mut sample = || -> c32 { rand.sample(rng) };
		let shift = c32::new(0.5, 0.5);
		let x = sample() - shift;
		let y = sample() - shift;
		let (c, s, r) = make_givens(x, y);
		assert!((x * c + y * s - r).norm() < 1e-6);
		assert!((-x * s.conj() + y * c).norm() < 1e-6);
	}
	#[test]
	fn test_hessenberg() {
		let rng = &mut StdRng::seed_from_u64(0);
		for n in [12, 35, 128, 255] {
			let rand = stats::CwiseMatDistribution {
				nrows: n,
				ncols: n,
				dist: ComplexDistribution::new(StandardNormal, StandardNormal),
			};
			let mut sample = || -> Mat<c64> { rand.rand(rng) };
			let A = sample();
			let mut B = sample();
			zip!(&mut B).for_each_triangular_lower(
				linalg::zip::Diag::Skip,
				|unzip!(x)| {
					*x = 0.0.into();
				},
			);
			let B = B;
			let mut mem =
				MemBuffer::new(generalized_hessenberg_scratch::<c64>(
					n,
					GeneralizedHessenbergParams {
						block_size: 32,
						..auto!(c64)
					},
				));
			for gen_hessenberg in [
				|A: MatMut<'_, _>,
				 B: MatMut<'_, _>,
				 Q: Option<MatMut<'_, _>>,
				 Z: Option<MatMut<'_, _>>,
				 _,
				 _,
				 par,
				 stack: &mut MemStack,
				 params| {
					generalized_hessenberg(A, B, Q, Z, par, stack, params)
				},
				generalized_hessenberg_unblocked,
			] {
				for (Q_is_I, Z_is_I) in
					[(false, false), (true, true), (true, false), (false, true)]
				{
					let Q0 = &if Q_is_I {
						Mat::identity(n, n)
					} else {
						UnitaryMat {
							dim: n,
							standard_normal: ComplexDistribution::new(
								StandardNormal,
								StandardNormal,
							),
						}
						.rand::<Mat<c64>>(rng)
					};
					let Z0 = &if Z_is_I {
						Mat::identity(n, n)
					} else {
						UnitaryMat {
							dim: n,
							standard_normal: ComplexDistribution::new(
								StandardNormal,
								StandardNormal,
							),
						}
						.rand::<Mat<c64>>(rng)
					};
					let mut Q = Q0.to_owned();
					let mut Z = Z0.to_owned();
					let mut H = A.clone();
					let mut T = B.clone();
					gen_hessenberg(
						H.as_mut(),
						T.as_mut(),
						Some(Q.as_mut()),
						Some(Z.as_mut()),
						Q_is_I,
						Z_is_I,
						Par::Seq,
						MemStack::new(&mut mem),
						GeneralizedHessenbergParams {
							block_size: 32,
							..auto!(c64)
						},
					);
					Q = Q0.adjoint() * &Q;
					Z = Z0.adjoint() * &Z;
					assert!((&Q * &H * Z.adjoint() - &A).norm_max() < 1e-13);
					assert!((&Q * &T * Z.adjoint() - &B).norm_max() < 1e-13);
					for j in 0..n {
						for i in j + 1..n {
							assert!(T[(i, j)] == c64::ZERO);
						}
						for i in j + 2..n {
							assert!(H[(i, j)] == c64::ZERO);
						}
					}
				}
			}
		}
	}
}
#[cfg(test)]
mod make_givens_tests {
	use super::make_givens;
	use crate::c64;
	use crate::internal_prelude::*;
	use std::assert;

	const EPS: f64 = f64::EPSILON;
	/// smallest positive subnormal
	const TINY: f64 = 4.9406564584124654e-324;
	/// error bounds, in units of `EPS` (relative) and `TINY` (absolute, for
	/// subnormal intermediates)
	const K_EPS: f64 = 4.0;
	const K_TINY: f64 = 4.0;

	/// magnitudes swept for `|f|` and `|g|`: subnormal, around
	/// `sqrt(min_positive)` and `sqrt(max_positive)`, ordinary and huge
	const MAGS: [f64; 22] = [
		TINY,
		3.0 * TINY,
		1e-320,
		1e-310,
		f64::MIN_POSITIVE,
		1e-305,
		1e-300,
		1e-200,
		1e-160,
		1.4916681462400413e-154, // sqrt(min_positive)
		1e-150,
		1e-30,
		0.7,
		1.0,
		1e30,
		1e150,
		6.703903964971299e153, // sqrt(max_positive)
		1e160,
		1e200,
		1e300,
		1e307,
		8e307,
	];

	/// the upstream (faer 0.24.4) `make_givens`, verbatim, and whether the
	/// fork keeps it (bit for bit) for `(f, g)`: everywhere `1 / h` and `1 / c`
	/// (`1 / |g|` when `f = 0`) are safely in range
	fn upstream<T: ComplexField>(f: T, g: T) -> Option<(T::Real, T, T)> {
		if g == zero() {
			Some((one(), zero(), f))
		} else if f == zero() {
			let c = zero::<T::Real>();
			let d = g.abs();
			let d_inv = d.recip();
			let s = g.conj().mul_real(&d_inv);
			let r = d.to_cplx();
			(d.is_finite() && d_inv.is_finite()).then_some((c, s, r))
		} else {
			let rtmin = sqrt_min_positive::<T::Real>();
			let rtmax = sqrt_max_positive::<T::Real>();
			let f1 = f.abs();
			let g1 = g.abs();
			let h1 = f1.hypot(g1);
			let h1_inv = &h1.recip();
			let c = f1 * h1_inv;
			let r = f.mul_real(c.recip());
			let s = g.conj() * r.mul_real(h1_inv).mul_real(h1_inv);
			(h1 > rtmin && h1 < rtmax && c > rtmin).then_some((c, s, r))
		}
	}

	/// `z 2^k`, exact unless the result is subnormal
	fn scale(mut z: c64, mut k: i32) -> c64 {
		while k != 0 {
			let step = k.clamp(-1000, 1000);
			z *= 2.0f64.powi(step);
			k -= step;
		}
		z
	}

	/// `k` such that `x 2^k` is in `[1, 2)`, for finite `x > 0`
	fn exponent_shift(x: f64) -> i32 {
		-(x.log2().floor() as i32)
	}

	/// `z / |z|` to full precision, also for subnormal `z`
	fn phase(z: c64) -> c64 {
		let z = scale(z, exponent_shift(z.norm()));
		z / z.norm()
	}

	/// checks `(c, s, r) = make_givens(f, g)` against the defining relations
	/// `c = |f| / h`, `s = sgn(f) conj(g) / h`, `r = sgn(f) h` (with
	/// `sgn(0) = 1`), all evaluated with exact power-of-two scaling, and the
	/// rotation property `[c, s; -conj(s), c] [f; g] = [r; 0]`
	#[track_caller]
	fn check_cplx(f: c64, g: c64) {
		let (c, s, r) = make_givens(f, g);
		let ctx =
			format!("f = {f:e}, g = {g:e} -> c = {c:e}, s = {s:e}, r = {r:e}");
		assert!(
			c.is_finite()
				&& s.re.is_finite()
				&& s.im.is_finite()
				&& r.re.is_finite()
				&& r.im.is_finite(),
			"non-finite: {ctx}"
		);
		assert!(c >= 0.0, "{ctx}");
		let unit = c * c + s.norm_sqr();
		assert!(
			(unit - 1.0).abs() <= K_EPS * EPS,
			"not unitary ({:e}): {ctx}",
			unit - 1.0
		);

		let c0 = c64::new(0.0, 0.0);
		if g == c0 {
			assert!(c == 1.0 && s == c0 && r == f, "{ctx}");
			return;
		}
		let m = f.norm().max(g.norm());
		let k = exponent_shift(m);
		let (fk, gk) = (scale(f, k), scale(g, k));
		let h = fk.norm().hypot(gk.norm());
		let sgn = if f == c0 {
			c64::new(1.0, 0.0)
		} else {
			phase(f)
		};
		let c_ref = fk.norm() / h;
		let s_ref = sgn * gk.conj() / h;
		let r_ref = sgn * h;
		let floor = K_TINY * TINY;
		assert!((c - c_ref).abs() <= K_EPS * EPS * c_ref + floor, "c: {ctx}");
		assert!((s - s_ref).norm() <= K_EPS * EPS + floor, "s: {ctx}");
		// `r` is subnormal when `m` is: allow its spacing, scaled by `2^k`
		let r_floor = scale(c64::new(floor, 0.0), k).re;
		assert!(
			(scale(r, k) - r_ref).norm() <= K_EPS * EPS * h + r_floor,
			"r: {ctx}"
		);

		// the rotation zeroes the second component and maps `f` to `r`
		let tol = 2.0 * K_EPS * EPS * m + floor;
		assert!((f * c + s * g - r).norm() <= tol, "top: {ctx}");
		assert!((g * c - s.conj() * f).norm() <= tol, "bottom: {ctx}");

		if let Some((cu, su, ru)) = upstream(f, g) {
			assert!(
				c.to_bits() == cu.to_bits()
					&& s.re.to_bits() == su.re.to_bits()
					&& s.im.to_bits() == su.im.to_bits()
					&& r.re.to_bits() == ru.re.to_bits()
					&& r.im.to_bits() == ru.im.to_bits(),
				"differs from upstream ({cu:e}, {su:e}, {ru:e}): {ctx}"
			);
		}
	}

	/// real counterpart of [`check_cplx`] (the `dlartg` path)
	#[track_caller]
	fn check_real(f: f64, g: f64) {
		let (c, s, r) = make_givens(f, g);
		let ctx =
			format!("f = {f:e}, g = {g:e} -> c = {c:e}, s = {s:e}, r = {r:e}");
		assert!(c.is_finite() && s.is_finite() && r.is_finite(), "{ctx}");
		assert!(c >= 0.0, "{ctx}");
		let unit = c * c + s * s;
		assert!((unit - 1.0).abs() <= K_EPS * EPS, "not unitary: {ctx}");
		if g == 0.0 {
			assert!(c == 1.0 && s == 0.0 && r == f, "{ctx}");
			return;
		}
		let m = f.abs().max(g.abs());
		let k = exponent_shift(m);
		let fk = scale(c64::new(f, 0.0), k).re;
		let gk = scale(c64::new(g, 0.0), k).re;
		let h = fk.hypot(gk);
		let sgn = if f == 0.0 { 1.0 } else { f.signum() };
		let floor = K_TINY * TINY;
		let c_ref = fk.abs() / h;
		assert!((c - c_ref).abs() <= K_EPS * EPS * c_ref + floor, "c: {ctx}");
		assert!((s - sgn * gk / h).abs() <= K_EPS * EPS + floor, "s: {ctx}");
		let r_floor = scale(c64::new(floor, 0.0), k).re;
		let rk = scale(c64::new(r, 0.0), k).re;
		assert!(
			(rk - sgn * h).abs() <= K_EPS * EPS * h + r_floor,
			"r: {ctx}"
		);
		let tol = 2.0 * K_EPS * EPS * m + floor;
		assert!((f * c + s * g - r).abs() <= tol, "top: {ctx}");
		assert!((g * c - s * f).abs() <= tol, "bottom: {ctx}");
		if let Some((cu, su, ru)) = upstream(f, g) {
			assert!(
				c.to_bits() == cu.to_bits()
					&& s.to_bits() == su.to_bits()
					&& r.to_bits() == ru.to_bits(),
				"differs from upstream: {ctx}"
			);
		}
	}

	/// geode-fem #908 / #867: upstream `make_givens` formed `1 / h` with
	/// `h = hypot(|f|, |g|)` and `1 / c` with `c = |f| / h`, both of which
	/// overflow once their argument drops below `1 / max_positive` (~5.6e-309).
	/// in the complex QZ, chaining shifts from a degenerate eigenvalue cluster
	/// drives the bulge entries that far down; the iterates turned into NaN and
	/// the blocked loop ran out its `30 n` sweeps (an effective hang). the
	/// first fix computed `sgn(f)` from `f / max(|f|, |g|)`, which still lost
	/// precision (and gave `0 / 0`) once `|f| / |g|` dropped below
	/// `min_positive`; this sweep covers that regime too (lapack `zlartg`)
	#[test]
	fn givens_is_finite_and_unitary_across_the_exponent_range() {
		let phases = [
			c64::new(1.0, 0.0),
			c64::new(0.0, 1.0),
			c64::new(-1.0, 0.0),
			c64::new(0.6, 0.8),
			c64::new(-0.28, -0.96),
			c64::new(0.3, -0.4),
		];
		let c0 = c64::new(0.0, 0.0);
		for &mf in &MAGS {
			for &pf in &phases {
				let f = pf * mf;
				check_cplx(f, c0);
				check_cplx(c0, f);
				for &mg in &MAGS {
					for &pg in &phases {
						check_cplx(f, pg * mg);
					}
				}
			}
			for sf in [1.0, -1.0] {
				check_real(sf * mf, 0.0);
				check_real(0.0, sf * mf);
				for &mg in &MAGS {
					for sg in [1.0, -1.0] {
						check_real(sf * mf, sg * mg);
					}
				}
			}
		}
		check_cplx(c0, c0);
		check_real(0.0, 0.0);

		// the case from the review of rjwalters/faier#1: `|f| / |g| ~ 5e-316`
		check_cplx(c64::new(3e-301, 4e-301), c64::new(1e15, 0.0));
		check_cplx(c64::new(3e-301, 4e-301), c64::new(-2e19, 1e20));
	}

	#[test]
	fn givens_terminates_on_non_finite_input() {
		for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
			for (f, g) in [
				(c64::new(bad, 0.0), c64::new(1.0, 0.0)),
				(c64::new(1.0, 0.0), c64::new(0.0, bad)),
				(c64::new(0.0, 0.0), c64::new(bad, 1.0)),
			] {
				let _ = make_givens(f, g);
			}
			let _ = make_givens(bad, 1.0);
			let _ = make_givens(1.0, bad);
			let _ = make_givens(0.0, bad);
		}
	}
}
