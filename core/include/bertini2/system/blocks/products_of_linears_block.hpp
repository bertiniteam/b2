//This file is part of Bertini 2.
//
//products_of_linears_block.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//products_of_linears_block.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with products_of_linears_block.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/system/blocks/products_of_linears_block.hpp

\brief A System evaluation block whose functions are each a product of linear forms.

This is the non-node representation of products of linears used by the m-homogeneous
start system and by regeneration.  Each function is

    f_i(x) = prod_r ( c_{i,r} . [x ; 1] )

i.e. a product over factors r of a linear form in the variables (the last column of
each coefficient row is the constant / homogenizing-variable coefficient, so x is
evaluated augmented with a trailing 1).  Evaluation is matrix-multiplies followed by
row-wise products; the Jacobian is the product rule, done with prefix/suffix products.
Nothing here is a function-tree node, so it never touches the SLP compiler.

Mirrors the multiprecision-coefficient pattern of bertini::LinearSlice (system/slice.hpp):
an mpfr master plus per-type working copies, with Precision() recasting the mpfr copy.
*/

#pragma once

#include <vector>
#include <tuple>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"

namespace bertini {
namespace blocks {

/**
\brief A block of system functions, each a product of linear forms.

Value-in / value-out: evaluation takes the variable vector and writes into a caller-
provided segment, so the block is self-contained and unit-testable without a System.
*/
class ProductsOfLinearsBlock
{
public:
	ProductsOfLinearsBlock() : num_vars_(0), precision_(DefaultPrecision()) {}

	/**
	\param num_vars The number of variables (n).  Each coefficient row has n+1 entries
	(the trailing entry multiplies the augmenting 1 -- the constant / homogenizing term).
	\param factors One augmented coefficient matrix per function; matrix i has shape
	(number-of-factors_i) x (num_vars + 1).
	*/
	ProductsOfLinearsBlock(size_t num_vars, std::vector<Mat<mpfr_complex>> factors)
		: num_vars_(num_vars), factors_highest_precision_(std::move(factors)), precision_(DefaultPrecision())
	{
#ifndef NDEBUG
		for (auto const& M : factors_highest_precision_)
			assert(static_cast<size_t>(M.cols()) == num_vars_ + 1 &&
			       "every products-of-linears coefficient matrix must have num_vars+1 columns");
#endif
		BuildWorking();
	}

	/// Number of functions (rows the block contributes to the system).
	size_t NumFunctions() const { return factors_highest_precision_.size(); }

	/// Number of variables the block expects in the input vector.
	size_t NumVariables() const { return num_vars_; }

	/// Products of linears do not depend on the path variable.
	bool DependsOnPathVariable() const { return false; }

	/// The Jacobian is not constant (it depends on x), unlike a slice.
	bool HasConstantJacobian() const { return false; }

	unsigned Precision() const { return precision_; }

	/// Set the working precision; recasts the mpfr working coefficients from the master.
	void Precision(unsigned new_precision) const
	{
		if (new_precision > DoublePrecision())
		{
			auto& wm = std::get<std::vector<Mat<mpfr_complex>>>(factors_working_);
			for (size_t i = 0; i < wm.size(); ++i)
				for (Eigen::Index r = 0; r < wm[i].rows(); ++r)
					for (Eigen::Index c = 0; c < wm[i].cols(); ++c)
					{
						wm[i](r, c).precision(new_precision);
						if (new_precision > precision_)
							wm[i](r, c) = factors_highest_precision_[i](r, c);
					}
		}
		precision_ = new_precision;
	}

	/**
	\brief Evaluate the block's function values into a caller-provided segment.

	\param result Length-NumFunctions() segment to write into.
	\param vars   Length-NumVariables() current variable values.
	\param path_value The path-variable value (ignored: products of linears are autonomous).
	*/
	template <typename T>
	void EvalInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& /*path_value*/) const
	{
		const auto& W = Working<T>();
		// aug = [vars ; 1].  The trailing 1 lets each coefficient row carry its constant
		// term in its last column, so a linear factor is just the dot product (row . aug).
		const Vec<T> aug = Augment<T>(vars);

		// Function i is a product of linear factors.  W[i] is its k x (n+1) coefficient
		// matrix, one row per factor, so W[i] * aug evaluates all k factors at once and the
		// function value is their product.
		for (size_t i = 0; i < W.size(); ++i)
		{
			const Vec<T> fvals = W[i] * aug;          // the k factor values L_0 .. L_{k-1}
			T val(1);
			for (Eigen::Index r = 0; r < fvals.size(); ++r)
				val *= fvals(r);                      // f_i = prod_r L_r
			result(static_cast<Eigen::Index>(i)) = val;
		}
	}

	/**
	\brief Evaluate the block's Jacobian (d f_i / d x_j) into a caller-provided block.

	\param J  A NumFunctions() x NumVariables() block to write into.
	\param vars Length-NumVariables() current variable values.
	\param path_value The path-variable value (ignored: products of linears are autonomous).
	*/
	template <typename T>
	void JacobianInPlace(Eigen::Ref<Mat<T>> J, Vec<T> const& vars, T const& /*path_value*/) const
	{
		const auto& W = Working<T>();
		const Vec<T> aug = Augment<T>(vars);          // [vars ; 1]; see EvalInPlace

		// Function i is f_i = prod_r L_r, where L_r = (row r of M) . aug is the r-th linear
		// factor's value.  By the product rule, the partial derivative w.r.t. variable c is
		//
		//     d f_i / d x_c = sum_r (d L_r / d x_c) * prod_{s != r} L_s
		//                   = sum_r      M(r,c)     * weight_r,
		//
		// since d L_r / d x_c is just the coefficient M(r,c) of x_c in factor r, and
		// weight_r := prod_{s != r} L_s is the product of every factor value except r's.
		for (size_t i = 0; i < W.size(); ++i)
		{
			const Mat<T>& M = W[i];                   // k x (n+1): rows are factors, last col is constant
			const Eigen::Index k = M.rows();

			if (k == 0)                               // a 0-factor product is the constant 1; derivative 0
			{
				J.row(static_cast<Eigen::Index>(i)).setZero();
				continue;
			}

			const Vec<T> fvals = M * aug;             // the factor values L_0 .. L_{k-1}

			// weight_r = prod_{s != r} L_s, computed in O(k) and WITHOUT division -- so it
			// stays correct even when some factor value L_s is zero (dividing the total
			// product by L_r would not).  We split the "all but r" product into the factors
			// before r and the factors after r, each built by a running accumulator:
			//   forward pass:  weight_r <- prod_{s < r} L_s            (the prefix product)
			//   backward pass: weight_r <- weight_r * prod_{s > r} L_s (times the suffix product)
			Vec<T> weight(k);
			{
				T acc(1);
				for (Eigen::Index r = 0; r < k; ++r) { weight(r) = acc; acc *= fvals(r); }      // prefix
				acc = T(1);
				for (Eigen::Index r = k - 1; r >= 0; --r) { weight(r) *= acc; acc *= fvals(r); } // suffix
			}

			// Row i of the Jacobian is (sum_r M(r,c) * weight_r) over the variable columns
			// c only; M.leftCols(num_vars_) drops the trailing constant column (not a variable).
			J.row(static_cast<Eigen::Index>(i)) = weight.transpose() * M.leftCols(static_cast<Eigen::Index>(num_vars_));
		}
	}

	/**
	\brief Time-derivative into a caller-provided segment.  Products of linears are
	autonomous (no path-variable dependence), so this is identically zero.
	*/
	template <typename T>
	void TimeDerivInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& /*vars*/, T const& /*path_value*/) const
	{
		result.setZero();
	}

private:
	template <typename T>
	const std::vector<Mat<T>>& Working() const
	{
		return std::get<std::vector<Mat<T>>>(factors_working_);
	}

	template <typename T>
	Vec<T> Augment(Vec<T> const& vars) const
	{
		Vec<T> aug(static_cast<Eigen::Index>(num_vars_ + 1));
		aug.head(static_cast<Eigen::Index>(num_vars_)) = vars;
		T one(1);
		if constexpr (!std::is_same<T, dbl>::value)
			one.precision(precision_);
		aug(static_cast<Eigen::Index>(num_vars_)) = one;
		return aug;
	}

	void BuildWorking() const
	{
		auto& wd = std::get<std::vector<Mat<dbl>>>(factors_working_);
		auto& wm = std::get<std::vector<Mat<mpfr_complex>>>(factors_working_);
		const size_t n = factors_highest_precision_.size();
		wd.resize(n);
		wm.resize(n);
		for (size_t i = 0; i < n; ++i)
		{
			const auto& M = factors_highest_precision_[i];
			wd[i].resize(M.rows(), M.cols());
			wm[i].resize(M.rows(), M.cols());
			for (Eigen::Index r = 0; r < M.rows(); ++r)
				for (Eigen::Index c = 0; c < M.cols(); ++c)
				{
					wd[i](r, c) = dbl(M(r, c));
					wm[i](r, c) = M(r, c);
				}
		}
	}

	size_t num_vars_;
	std::vector<Mat<mpfr_complex>> factors_highest_precision_; ///< master coefficients, one matrix per function
	mutable std::tuple<std::vector<Mat<dbl>>, std::vector<Mat<mpfr_complex>>> factors_working_;
	mutable unsigned precision_;
};

} // namespace blocks
} // namespace bertini
