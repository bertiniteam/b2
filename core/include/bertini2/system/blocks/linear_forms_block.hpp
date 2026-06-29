//This file is part of Bertini 2.
//
//linear_forms_block.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//linear_forms_block.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with linear_forms_block.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/system/blocks/linear_forms_block.hpp

\brief A System evaluation block whose functions are linear forms: f(x) = M x + b.

This is the linear-algebra workhorse for the block-composed System: a stack of affine
linear forms evaluated as one matrix-vector product, rather than as expanded scalar
function-tree expressions.  It is the natural home for linear slices, randomization /
linear projections, and the linear conditions emitted by the bertini.linalg layer when
the coefficients are concrete (multiprecision) numbers rather than symbolic.

Each function is

    f_i(x) = ( c_{i} . [x ; 1] )

i.e. a single linear form whose last coefficient (the column past the variables) carries
the constant term, so x is evaluated augmented with a trailing 1 -- exactly the augmented
convention of ProductsOfLinearsBlock, but with one factor per function instead of a
product.  Evaluation is a single matrix-vector multiply; the Jacobian is the constant
coefficient matrix (its variable columns).

Mirrors the multiprecision-coefficient pattern of ProductsOfLinearsBlock: an mpfr master
plus per-type working copies, with Precision() recasting the mpfr copy.  bertini::Slice
(system/slice.hpp) is a thin wrapper around this block.
*/

#pragma once

#include <tuple>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/system/blocks/describe.hpp"

namespace bertini {
namespace blocks {

/**
\brief A block of system functions, each an affine linear form f_i(x) = c_i . [x ; 1].

Value-in / value-out: evaluation takes the variable vector and writes into a caller-
provided segment, so the block is self-contained and unit-testable without a System.
*/
class LinearFormsBlock
{
public:
	LinearFormsBlock() : num_vars_(0), precision_(DefaultPrecision()) {}

	/**
	\param num_vars The number of variables (n).  The coefficient matrix has n+1 columns
	(the trailing column multiplies the augmenting 1 -- the constant term of each form).
	\param coefficients The augmented coefficient matrix: one row per function, shape
	(number-of-functions) x (num_vars + 1).
	*/
	LinearFormsBlock(size_t num_vars, Mat<complex_mp> coefficients)
		: num_vars_(num_vars), coefficients_highest_precision_(std::move(coefficients)),
		  precision_(DefaultPrecision())
	{
		assert(static_cast<size_t>(coefficients_highest_precision_.cols()) == num_vars_ + 1 &&
		       "a linear-forms coefficient matrix must have num_vars+1 columns");
		BuildWorking();
	}

	/// Number of functions (rows the block contributes to the system).
	size_t NumFunctions() const { return static_cast<size_t>(coefficients_highest_precision_.rows()); }

	/// Linear forms are degree 1.
	std::vector<int> Degrees() const { return std::vector<int>(NumFunctions(), 1); }
	std::vector<int> Degrees(VariableGroup const&) const { return Degrees(); }

	/// Linear forms are polynomial.
	bool IsPolynomial(VariableGroup const&) const { return true; }
	/// An affine form a.x + b is inhomogeneous (it carries a constant); once homogenized the
	/// constant has become a homogenizing-variable coefficient, so it is degree-1 homogeneous.
	bool IsHomogeneous(VariableGroup const&) const { return homogeneous_; }

	/// Homogenize: fold the constant column onto the homogenizing variable.  An augmented form
	/// a.x + b (last column = constant) becomes the homogeneous a.x + b*h, where the
	/// homogenizing variable h is prepended to the variable ordering -- so its column is the
	/// old constant column.  The block then evaluates M*[vars] directly (no augmenting 1).
	/// Currently supports a single affine variable group (the dominant bertini.linalg case);
	/// a second call (a second affine group) throws.
	void Homogenize(VariableGroup const& /*group*/, std::shared_ptr<node::Variable> const& /*hom_var*/)
	{
		if (homogeneous_)
			throw std::runtime_error("LinearFormsBlock::Homogenize: block is already homogenized "
				"(multiple affine variable groups are not yet supported for linear-forms blocks)");
		const auto& M = coefficients_highest_precision_;
		Mat<complex_mp> Mh(M.rows(), M.cols());                                   // same #cols: n+1
		Mh.col(0) = M.col(static_cast<Eigen::Index>(num_vars_));                    // constant -> h column (front)
		Mh.rightCols(static_cast<Eigen::Index>(num_vars_)) =
			M.leftCols(static_cast<Eigen::Index>(num_vars_));                       // original variable columns
		coefficients_highest_precision_ = Mh;
		num_vars_ += 1;          // h is now a real variable; all n+1 columns are variable columns
		homogeneous_ = true;
		BuildWorking();
	}

	/// Number of variables the block expects in the input vector.
	size_t NumVariables() const { return num_vars_; }

	/// The master coefficient matrix (one row per form).  Affine: num_vars+1 columns, the last being
	/// the constant term.  Homogeneous (post-Homogenize): num_vars columns, all variable columns.
	/// Exposed for the function-tree expansion (System::NaturalFunctionsAsNodes).
	Mat<complex_mp> const& Coefficients() const { return coefficients_highest_precision_; }
	/// Whether Homogenize has folded the constant column onto a homogenizing variable.
	bool IsHomogenized() const { return homogeneous_; }

	/// Human-facing description: each form prints as the placeholder 'f_k = c.[x, y, 1]' (structure
	/// stays legible) followed by its actual coefficient row in a 'c =' legend below -- short (4
	/// significant figures) in terse, full precision in verbose.  Terse truncates after kTerseRowCap
	/// forms so a large slice does not flood the terminal.
	void Describe(std::ostream& out, size_t& row, VariableGroup const& vars, bool verbose) const
	{
		auto const& M = coefficients_highest_precision_;
		const Eigen::Index n = M.rows();
		if (n == 0)
			return;
		const Eigen::Index cap   = static_cast<Eigen::Index>(describe_detail::kTerseRowCap);
		const Eigen::Index shown = (verbose || n <= cap) ? n : cap;

		// placeholder line per form: f_k = c.[x, y, 1]
		for (Eigen::Index r = 0; r < shown; ++r)
		{
			out << "  f_" << row++ << " = c.";
			describe_detail::PrintAugmentedVars(out, vars, num_vars_, homogeneous_);
			out << "\n";
		}

		// the actual coefficients below, as a named-expression-style legend (row r <-> f_k above)
		const int sig = describe_detail::CoeffSig(verbose);
		out << "    c =\n";
		for (Eigen::Index r = 0; r < shown; ++r)
		{
			out << "      [ ";
			for (Eigen::Index c = 0; c < M.cols(); ++c)
			{
				if (c) out << ", ";
				describe_detail::PrintCoeff(out, M(r, c), sig);
			}
			out << " ]\n";
		}

		if (shown < n)
		{
			row += static_cast<size_t>(n - shown);     // keep the global row index correct
			out << "    ... (" << (n - shown) << " more form" << (n - shown == 1 ? "" : "s")
			    << "; describe(verbose=True) for all)\n";
		}
	}

	/// Linear forms do not depend on the path variable.
	bool DependsOnPathVariable() const { return false; }

	/// The Jacobian is the (constant) coefficient matrix, independent of x and t.
	bool HasConstantJacobian() const { return true; }

	/// Analytic block: nothing symbolic to differentiate.
	void Differentiate() const {}

	unsigned Precision() const { return precision_; }

	/// Set the working precision; recasts the mpfr working coefficients from the master.
	void Precision(unsigned new_precision) const
	{
		if (new_precision > DoublePrecision())
		{
			auto& wm = std::get<Mat<complex_mp>>(coefficients_working_);
			for (Eigen::Index r = 0; r < wm.rows(); ++r)
				for (Eigen::Index c = 0; c < wm.cols(); ++c)
				{
					wm(r, c).precision(new_precision);
					if (new_precision > precision_)
						wm(r, c) = coefficients_highest_precision_(r, c);
				}
		}
		precision_ = new_precision;
	}

	/**
	\brief Evaluate the block's function values into a caller-provided segment.

	The path variable is ignored (linear forms are autonomous).

	\param result Length-NumFunctions() segment to write into.
	\param vars   Length-NumVariables() current variable values.
	*/
	template <typename T>
	void EvalInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& /*path_value*/) const
	{
		// affine: f(x) = W * [x ; 1] (the trailing 1 carries each row's constant in the last
		// column).  homogeneous (post-Homogenize): every column is a variable column, so it is
		// just W * vars (the old constant is now the homogenizing variable's coefficient).
		if (homogeneous_)
			result.noalias() = Working<T>() * vars;
		else
			result.noalias() = Working<T>() * Augment<T>(vars);
	}

	/**
	\brief Evaluate the block's Jacobian (d f_i / d x_j) into a caller-provided block.

	The Jacobian of f(x) = M x + b is simply M (its variable columns) -- constant in x.
	The variable values and path variable are unused (constant Jacobian; autonomous).

	\param J  A NumFunctions() x NumVariables() block to write into.
	*/
	template <typename T>
	void JacobianInPlace(Eigen::Ref<Mat<T>> J, Vec<T> const& /*vars*/, T const& /*path_value*/) const
	{
		// homogeneous: every column is d f / d x.  affine: drop the trailing constant column.
		if (homogeneous_)
			J = Working<T>();
		else
			J = Working<T>().leftCols(static_cast<Eigen::Index>(num_vars_));
	}

	/**
	\brief Time-derivative into a caller-provided segment.  Linear forms are autonomous
	(no path-variable dependence), so this is identically zero.
	*/
	template <typename T>
	void TimeDerivInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& /*vars*/, T const& /*path_value*/) const
	{
		result.setZero();
	}

private:
	template <typename T>
	const Mat<T>& Working() const
	{
		return std::get<Mat<T>>(coefficients_working_);
	}

	template <typename T>
	Vec<T> Augment(Vec<T> const& vars) const
	{
		Vec<T> aug(static_cast<Eigen::Index>(num_vars_ + 1));
		aug.head(static_cast<Eigen::Index>(num_vars_)) = vars;
		T one(1);
		if constexpr (!std::is_same<T, complex_dbl>::value)
			one.precision(precision_);
		aug(static_cast<Eigen::Index>(num_vars_)) = one;
		return aug;
	}

	void BuildWorking() const
	{
		const auto& M = coefficients_highest_precision_;
		auto& wd = std::get<Mat<complex_dbl>>(coefficients_working_);
		auto& wm = std::get<Mat<complex_mp>>(coefficients_working_);
		wd.resize(M.rows(), M.cols());
		wm.resize(M.rows(), M.cols());
		for (Eigen::Index r = 0; r < M.rows(); ++r)
			for (Eigen::Index c = 0; c < M.cols(); ++c)
			{
				wd(r, c) = complex_dbl(M(r, c));
				wm(r, c) = M(r, c);
			}
	}

	size_t num_vars_;
	bool homogeneous_ = false; ///< false: augmented affine (M*[x;1]); true: post-Homogenize (M*x)
	Mat<complex_mp> coefficients_highest_precision_; ///< master: rows = functions, cols = num_vars (homogeneous) or num_vars+1 (affine)
	mutable std::tuple<Mat<complex_dbl>, Mat<complex_mp>> coefficients_working_;
	mutable unsigned precision_;

	friend class boost::serialization::access;

	template <typename Archive>
	void serialize(Archive& ar, const unsigned /*version*/)
	{
		ar & num_vars_;
		ar & homogeneous_;
		ar & precision_;
		ar & coefficients_highest_precision_;
		ar & std::get<0>(coefficients_working_);
		ar & std::get<1>(coefficients_working_);
	}
};

} // namespace blocks
} // namespace bertini
