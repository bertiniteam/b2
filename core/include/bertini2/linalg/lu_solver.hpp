//This file is part of Bertini 2.
//
//bertini2/linalg/lu_solver.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/linalg/lu_solver.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/linalg/lu_solver.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// silviana amethyst, university of wisconsin-eau claire

/**
 \file bertini2/linalg/lu_solver.hpp

 \brief A stateful partial-pivot LU solver that owns precision-adjustable temporary banks.

 Eigen's generic dense LU allocates one multiprecision temporary per elementary scalar
 operation -- O(N^3) heap allocations (mpc_init2 / mpc_clear, plus the per-alloc MPFR
 get_memory_functions fetch, the thread-local working-precision read, and the mimalloc
 ownership check) per factorization.  For runtime-precision mpc that allocation churn is
 a large fraction of the cost and it scales with N^3.

 This class factors in place, reusing a single scratch scalar across the whole
 elimination, so a factorization performs O(1) allocations instead of O(N^3).  It owns
 its workspace across calls and exposes a ChangePrecision() so the owner (a Newton
 corrector / predictor) can move it up and down in precision the way the AMP tracker and
 endgames re-precision their scratch -- see NewtonCorrector::ChangePrecision.
*/

#ifndef BERTINI_LINALG_LU_SOLVER_HPP
#define BERTINI_LINALG_LU_SOLVER_HPP

#include "bertini2/eigen_extensions.hpp"

namespace bertini {
namespace linalg {

/**
 \class PartialPivLU

 \brief Stateful partial-pivot LU factor/solve with preallocated, precision-adjustable banks.

 \tparam NumT The (complex) scalar type -- complex_dbl or complex_mp.

 Numerically this is the same right-looking Gaussian elimination with partial pivoting as
 Eigen's PartialPivLU (unit lower-triangular L, upper-triangular U, row pivoting), so
 solutions agree with Eigen to backward-stable accuracy.  The win is allocation behavior,
 not the arithmetic.
*/
template<typename NumT>
class PartialPivLU
{
	using RealT = typename Eigen::NumTraits<NumT>::Real;
	static constexpr bool is_mp = !std::is_same<NumT, complex_dbl>::value;

public:
	/// \brief (Re)size the workspace to hold n x n systems.  Call on a system-size change.
	/// \param n The matrix dimension.
	void ChangeSize(unsigned n)
	{
		n_ = n;
		lu_.resize(n, n);   // new mp elements come up at the thread default precision
		piv_.resize(n);
		y_.resize(n);
		if (is_mp && precision_ != 0)   // re-establish a previously-set precision on the new storage
			ApplyPrecision();
	}

	/// \brief Re-set the precision of every mp temporary bank (a no-op for complex_dbl).
	/// \param prec The new working precision (digits).
	void ChangePrecision(unsigned prec)
	{
		precision_ = prec;
		if (is_mp)
			ApplyPrecision();
	}

	/// \brief Get the current working precision (digits).
	unsigned precision() const { return precision_; }

	/// \brief Factor a copy of A with partial pivoting; A is preserved.  O(1) allocations.
	/// \param A The square matrix to factor.
	/// \return Success unless the factored diagonal looks near-singular / ill-conditioned.
	MatrixSuccessCode Factor(Mat<NumT> const& A)
	{
		lu_ = A;                 // one copy into reusable workspace (A must survive)
		return FactorInPlace();
	}

	/// \brief Factor A destructively (no copy): on return A holds the L\\U factors.
	/// \param A The square matrix; overwritten with the factorization.
	/// \return Success unless the factored diagonal looks near-singular / ill-conditioned.
	MatrixSuccessCode FactorDestructive(Mat<NumT> & A)
	{
		lu_.swap(A);             // steal A's storage; A is left holding our old workspace
		auto code = FactorInPlace();
		return code;
	}

	/// \brief Solve A x = b using the stored factors.  Reuses scratch; O(1) allocations.
	/// \tparam OutDerived The output type -- a Vec or a writable block such as K.col(stage).
	/// \param b The right-hand side (must be materialized storage, distinct from x).
	/// \param[out] x_out The solution (sized n).  A temporary expression (e.g. .col()) is fine.
	template<typename OutDerived>
	void Solve(Vec<NumT> const& b, Eigen::MatrixBase<OutDerived> const& x_out) const
	{
		// Eigen idiom for a write-to output that may be a temporary Block (e.g. K.col(stage)).
		Eigen::MatrixBase<OutDerived>& x = const_cast<Eigen::MatrixBase<OutDerived>&>(x_out);

		// y_ = P b : copy then replay the row swaps recorded during factorization.
		y_ = b;
		for (unsigned k = 0; k < n_; ++k)
			if (piv_[k] != k)
				std::swap(y_[k], y_[piv_[k]]);

		// forward substitution: L y = P b, L unit lower-triangular (implicit unit diagonal).
		for (unsigned i = 1; i < n_; ++i)
			for (unsigned j = 0; j < i; ++j)
			{
				t_ = lu_(i, j);
				t_ *= y_[j];
				y_[i] -= t_;      // y_[i] -= lu_(i,j) * y_[j]
			}

		// back substitution: U x = y, U upper-triangular.
		for (unsigned ii = n_; ii-- > 0; )
		{
			for (unsigned j = ii + 1; j < n_; ++j)
			{
				t_ = lu_(ii, j);
				t_ *= x(j);
				y_[ii] -= t_;     // y_[ii] -= lu_(ii,j) * x(j)
			}
			x(ii) = y_[ii];
			x(ii) /= lu_(ii, ii);
		}
	}

	/// \brief The factored matrix (L below the diagonal, U on and above), for the caller's
	/// LUPartialPivotDecompositionSuccessful()-style health check.
	Mat<NumT> const& Factors() const { return lu_; }

private:
	/// \brief Right-looking partial-pivot Gaussian elimination in place on lu_.
	MatrixSuccessCode FactorInPlace()
	{
		for (unsigned k = 0; k < n_; ++k)
		{
			// partial pivot: pick the row p >= k with the largest |lu_(p,k)|.
			unsigned p = k;
			amax_ = abs(lu_(k, k));
			for (unsigned i = k + 1; i < n_; ++i)
			{
				acur_ = abs(lu_(i, k));
				if (acur_ > amax_) { amax_ = acur_; p = i; }
			}
			piv_[k] = p;
			if (p != k)
				lu_.row(k).swap(lu_.row(p));

			// scale the sub-column below the pivot by the reciprocal of the pivot -- one (costly)
			// mpc division per column instead of one per sub-row, then cheap multiplies.  Pin recip_
			// to the pivot's own precision first: during AMP the working precision changes, and a
			// reciprocal left at a stale precision would contaminate the factorization (and, through
			// the solution, trip the endgame's uniform-precision checks).
			if constexpr (is_mp)
				Precision(recip_, Precision(lu_(k, k)));
			recip_ = 1;
			recip_ /= lu_(k, k);                 // recip_ = 1 / a_kk, at the pivot's precision
			for (unsigned i = k + 1; i < n_; ++i)
				lu_(i, k) *= recip_;             // l_ik = a_ik * (1 / a_kk)

			// rank-1 trailing-submatrix update, reusing t_ for every multiply-subtract so the
			// whole O(N^3) sweep allocates nothing beyond the scratch scalar.
			for (unsigned j = k + 1; j < n_; ++j)
				for (unsigned i = k + 1; i < n_; ++i)
				{
					t_ = lu_(i, k);
					t_ *= lu_(k, j);
					lu_(i, j) -= t_;   // lu_(i,j) -= lu_(i,k) * lu_(k,j)
				}
		}
		return LUPartialPivotDecompositionSuccessful(lu_);
	}

	/// \brief Set every mp bank to precision_ (called on resize / precision change).
	void ApplyPrecision()
	{
		if constexpr (is_mp)
		{
			Precision(lu_, precision_);
			Precision(y_, precision_);
			Precision(t_, precision_);
			Precision(recip_, precision_);   // re-pinned per pivot in FactorInPlace, but keep it valid
			Precision(amax_, precision_);
			Precision(acur_, precision_);
		}
	}

	unsigned n_ = 0;                 ///< The current matrix dimension.
	unsigned precision_ = 0;         ///< The current working precision (digits).
	Mat<NumT> lu_;                   ///< Factored-matrix workspace (reused across calls).
	std::vector<unsigned> piv_;      ///< Row pivots: piv_[k] is the row swapped into position k.
	mutable Vec<NumT> y_;            ///< Solve scratch (permuted rhs / forward-sub result).
	mutable NumT t_;                 ///< Reused scalar scratch for multiply-subtract.
	NumT recip_;                     ///< Reused scratch: reciprocal of the current pivot.
	RealT amax_;                     ///< Reused scratch: running max magnitude in pivot search.
	RealT acur_;                     ///< Reused scratch: current magnitude in pivot search.
};

} // namespace linalg
} // namespace bertini

#endif
