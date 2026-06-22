//This file is part of Bertini 2.
//
//randomization_block.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//randomization_block.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with randomization_block.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/system/blocks/randomization_block.hpp

\brief A System evaluation block that randomizes an overdetermined operand system down to a
square one, so its isolated solutions can be found by homotopy continuation.

An overdetermined system F = (f_1, ..., f_N) of N functions in n variables (N > n) cannot be
fed to a total-degree start system (which requires squareness).  Randomization replaces F with
n generic combinations whose isolated-solution set still contains that of F; the square result
is solved and the extraneous solutions discarded.

This block holds the overdetermined system as a `shared_ptr<SystemT>` operand (mirroring
`BlendBlock<SystemT>`) and a constant coefficient matrix, and emits the n randomized rows

    g_i(x) = sum_j  c_ij * f_j(x) * prod_g  h_g^{(D_{i,g} - d_{j,g})}

where f_j is the operand's j-th natural function (homogenized to its own multidegree d_j),
D_i is row i's target multidegree, h_g is the homogenizing variable of variable group g, and
c_ij = 0 whenever any deficit D_{i,g} - d_{j,g} is negative.  The h-power factors are what make
a constant-coefficient combination of functions of *differing* degree homogenize correctly:
before homogenization every h_g = 1 and the formula collapses to g_i = sum_j c_ij f_j.

Construction (degrees, sorting, the [I|C] / common-target-multidegree matrix) is done by
`System::Randomize`, which hands this block the finished coefficient matrix and the per-row /
per-function multidegrees.  The block is multidegree-native so single-affine-group and
multihomogeneous randomization share one evaluator.

Templated on the system type so `SystemT` stays a dependent name in the `Block` variant (a
System contains a Block which contains a RandomizationBlock<System> holding a shared_ptr<System>
-- the template parameter breaks the recursive type), exactly as for `BlendBlock<SystemT>`.

\note Like BlendBlock, the operand is held by shared_ptr and its mutable per-evaluation state
(variable values) is shared across copies; per-thread tracking must Clone (the deep-clone of a
block operand is a shared follow-up with BlendBlock).
*/

#pragma once

#include <vector>
#include <memory>
#include <stdexcept>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/split_member.hpp>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/function_tree.hpp"
#include "bertini2/system/blocks/describe.hpp"

namespace bertini {
namespace blocks {

template <typename SystemT>
class RandomizationBlock
{
public:
	using Var = std::shared_ptr<node::Variable>;
	using OperandPtr = std::shared_ptr<SystemT>;

	RandomizationBlock() : num_groups_(0), precision_(DefaultPrecision()) {}

	/**
	\param operand              The overdetermined system being randomized (N natural functions).
	\param coefficients         The n x N constant matrix; row i, column j is c_ij.
	\param target_multidegrees  n rows, each a length-num_groups vector D_i (row i's target degree
	                            per variable group).
	\param operand_multidegrees N rows, each a length-num_groups vector d_j (operand function j's
	                            degree per variable group), measured on the affine operand.
	\param num_groups           The number of (affine) variable groups G.
	*/
	RandomizationBlock(OperandPtr operand,
	                   Mat<mpfr_complex> coefficients,
	                   std::vector<std::vector<int>> target_multidegrees,
	                   std::vector<std::vector<int>> operand_multidegrees,
	                   size_t num_groups)
		: operand_(std::move(operand)),
		  coefficients_highest_precision_(std::move(coefficients)),
		  target_multidegrees_(std::move(target_multidegrees)),
		  operand_multidegrees_(std::move(operand_multidegrees)),
		  num_groups_(num_groups),
		  precision_(DefaultPrecision())
	{
		BuildWorking();
	}

	// Memory-isolating copy (ADR-0027): the operand System is deep-copied via the
	// System copy constructor, which shares its immutable node DAG and compiled SLP Program but
	// gives it its own per-thread evaluation Memory.  Everything else (the coefficient matrices,
	// multidegrees, working buffers) is value state, copied per copy.  This lets path tracking
	// Clone a System into per-thread copies that share no mutable state.
	RandomizationBlock(RandomizationBlock const& other)
		: operand_(other.operand_ ? std::make_shared<SystemT>(*other.operand_) : nullptr),
		  coefficients_highest_precision_(other.coefficients_highest_precision_),
		  coefficients_working_(other.coefficients_working_),
		  target_multidegrees_(other.target_multidegrees_),
		  operand_multidegrees_(other.operand_multidegrees_),
		  num_groups_(other.num_groups_),
		  homogenized_(other.homogenized_),
		  hom_vars_(other.hom_vars_),
		  hom_var_index_(other.hom_var_index_),
		  precision_(other.precision_)
	{}

	RandomizationBlock& operator=(RandomizationBlock const& other)
	{
		if (this != &other)
		{
			operand_ = other.operand_ ? std::make_shared<SystemT>(*other.operand_) : nullptr;
			coefficients_highest_precision_ = other.coefficients_highest_precision_;
			coefficients_working_ = other.coefficients_working_;
			target_multidegrees_ = other.target_multidegrees_;
			operand_multidegrees_ = other.operand_multidegrees_;
			num_groups_ = other.num_groups_;
			homogenized_ = other.homogenized_;
			hom_vars_ = other.hom_vars_;
			hom_var_index_ = other.hom_var_index_;
			precision_ = other.precision_;
		}
		return *this;
	}

	RandomizationBlock(RandomizationBlock&&) = default;
	RandomizationBlock& operator=(RandomizationBlock&&) = default;

	/// Number of randomized functions (rows this block contributes) = n.
	size_t NumFunctions() const { return static_cast<size_t>(coefficients_highest_precision_.rows()); }

	/// The randomized functions' degrees: the (total) target degree of each row, i.e. the sum of
	/// its target multidegree over the variable groups.  This is the degree-reporting hook -- the
	/// path count of the squared system is the product of these.
	std::vector<int> Degrees() const
	{
		std::vector<int> d(NumFunctions(), 0);
		for (size_t i = 0; i < target_multidegrees_.size(); ++i)
			for (int e : target_multidegrees_[i])
				d[i] += e;
		return d;
	}
	/// Per-group target degree.  We match the group against the operand's variable groups by
	/// identity so the correct column of the target multidegree is returned; an unrecognized group
	/// reports the total (defensive -- the owning System asks per its own groups).
	std::vector<int> Degrees(VariableGroup const& vars) const
	{
		const int g = GroupIndex(vars);
		if (g < 0)
			return Degrees();
		std::vector<int> d(NumFunctions(), 0);
		for (size_t i = 0; i < target_multidegrees_.size(); ++i)
			d[i] = target_multidegrees_[i][static_cast<size_t>(g)];
		return d;
	}

	/// The randomized system is polynomial / homogeneous exactly when the operand is (System's
	/// IsPolynomial/IsHomogeneous are whole-system, no-argument checks -- as BlendBlock delegates).
	bool IsPolynomial(VariableGroup const& /*vars*/) const { return operand_->IsPolynomial(); }
	bool IsHomogeneous(VariableGroup const& /*vars*/) const { return homogenized_ ? operand_->IsHomogeneous() : false; }

	/**
	\brief Homogenize: thread the owning system's homogenizing variable for this group into the
	operand, so the operand's functions and this block's h-power deficit factors share one set of
	homogenizing variables.

	The owning System calls this once per affine variable group, in group order.  We accumulate the
	supplied hom vars; once we have one per group we homogenize the operand with exactly those
	variables (System::Homogenize(VariableGroup const&) reuses them rather than minting fresh ones)
	and record each hom var's position in the operand's variable ordering for the h-power factors.
	*/
	void Homogenize(VariableGroup const& /*group*/, Var const& hom_var)
	{
		if (homogenized_)
			return;
		hom_vars_.push_back(hom_var);
		if (hom_vars_.size() == num_groups_)
		{
			operand_->Homogenize(VariableGroup(hom_vars_));
			homogenized_ = true;
			// the hom var of group g now sits at a fixed position of the (shared) variable ordering;
			// the block reads h_g from there during evaluation.
			hom_var_index_.assign(num_groups_, -1);
			auto const& ordering = operand_->Variables();
			for (size_t g = 0; g < num_groups_; ++g)
				for (size_t k = 0; k < ordering.size(); ++k)
					if (ordering[k] == hom_vars_[g]) { hom_var_index_[g] = static_cast<int>(k); break; }
		}
	}

	/// The randomization itself adds no path-variable dependence; it inherits the operand's.
	bool DependsOnPathVariable() const { return operand_->HavePathVariable(); }

	/// The randomized Jacobian varies in x (operand Jacobian times constants, plus h-power terms).
	bool HasConstantJacobian() const { return false; }

	/// Analytic block: nothing symbolic to differentiate (the operand differentiates itself).
	void Differentiate() const {}

	/// The overdetermined operand, for the function-tree expansion oracle and the matrix getter.
	OperandPtr const& Operand() const { return operand_; }
	/// The randomization matrix R (master, highest precision); n x N, row i column j is c_ij.
	Mat<mpfr_complex> const& RandomizationMatrix() const { return coefficients_highest_precision_; }

	// Accessors for the function-tree expansion oracle (System::NaturalFunctionsAsNodes), which
	// rebuilds g_i = sum_j c_ij * node(f_j) * prod_g h_g^{(D_{i,g}-d_{j,g})}.
	bool IsHomogenized() const { return homogenized_; }
	std::vector<Var> const& HomVars() const { return hom_vars_; }
	std::vector<std::vector<int>> const& TargetMultidegrees() const { return target_multidegrees_; }
	std::vector<std::vector<int>> const& OperandMultidegrees() const { return operand_multidegrees_; }

	/// Human-facing description: terse shows `f_a..f_b = R . g  (R: nxN)` then the underlying
	/// functions `g_j` indented (they are the interesting part); verbose additionally prints R's
	/// entries.
	void Describe(std::ostream& out, size_t& row, VariableGroup const& /*vars*/, bool verbose) const
	{
		const size_t n = NumFunctions();
		const size_t N = static_cast<size_t>(coefficients_highest_precision_.cols());
		out << "  ";
		describe_detail::PrintRowLabel(out, row, n);
		out << " = R . g   (R: " << n << "x" << N << " randomization matrix)\n";
		row += n;
		auto g = operand_->NaturalFunctionsAsNodes();
		for (size_t j = 0; j < g.size(); ++j)
			out << "      g_" << j << " = " << g[j] << "\n";
		if (verbose)
		{
			out << "    R =\n";
			for (Eigen::Index i = 0; i < coefficients_highest_precision_.rows(); ++i)
			{
				out << "      [ ";
				for (Eigen::Index j = 0; j < coefficients_highest_precision_.cols(); ++j)
				{
					if (j) out << ", ";
					describe_detail::PrintCoeff(out, coefficients_highest_precision_(i, j));
				}
				out << " ]\n";
			}
		}
	}

	unsigned Precision() const { return precision_; }

	/// Set the working precision: recast the mpfr working coefficients from the master and bring
	/// the operand to the same precision (so a high-precision coefficient is not multiplied against
	/// a low-precision operand value).
	void Precision(unsigned new_precision) const
	{
		if (new_precision > DoublePrecision())
		{
			auto& wm = std::get<Mat<mpfr_complex>>(coefficients_working_);
			for (Eigen::Index r = 0; r < wm.rows(); ++r)
				for (Eigen::Index c = 0; c < wm.cols(); ++c)
				{
					wm(r, c).precision(new_precision);
					if (new_precision > precision_)
						wm(r, c) = coefficients_highest_precision_(r, c);
				}
		}
		operand_->precision(new_precision);
		precision_ = new_precision;
	}

	/**
	\brief Evaluate the randomized function values into a caller-provided segment.

	g_i = sum_j c_ij * f_j * prod_g h_g^{(D_{i,g} - d_{j,g})}.
	*/
	template <typename T>
	void EvalInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& /*path_value*/) const
	{
		SyncPrecision(vars);
		const size_t N = static_cast<size_t>(coefficients_highest_precision_.cols());
		const Vec<T> f = operand_->template Eval<T>(vars).head(static_cast<Eigen::Index>(N));
		const auto& C = Working<T>();

		result.setZero();
		for (size_t i = 0; i < NumFunctions(); ++i)
		{
			T acc = Zero<T>();
			for (size_t j = 0; j < N; ++j)
			{
				const T& c = C(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
				if (IsZero(c))
					continue;
				acc += c * f(static_cast<Eigen::Index>(j)) * HPower<T>(i, j, vars);
			}
			result(static_cast<Eigen::Index>(i)) = acc;
		}
	}

	/**
	\brief Evaluate the randomized Jacobian into a caller-provided block.

	By the product rule on  c_ij * f_j * W_ij(x)  (W_ij = prod_g h_g^{e_{ij,g}}):
	    d g_i / d x_l = sum_j c_ij [ (d f_j / d x_l) W_ij + f_j (d W_ij / d x_l) ],
	and dW_ij/dx_l is nonzero only when x_l is a homogenizing variable h_{g'} appearing in W_ij.
	*/
	template <typename T>
	void JacobianInPlace(Eigen::Ref<Mat<T>> J, Vec<T> const& vars, T const& /*path_value*/) const
	{
		SyncPrecision(vars);
		const size_t N = static_cast<size_t>(coefficients_highest_precision_.cols());
		const Vec<T> f  = operand_->template Eval<T>(vars).head(static_cast<Eigen::Index>(N));
		const Mat<T> Jf = operand_->template Jacobian<T>(vars).topRows(static_cast<Eigen::Index>(N));
		const auto& C = Working<T>();

		J.setZero();
		for (size_t i = 0; i < NumFunctions(); ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				const T& c = C(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
				if (IsZero(c))
					continue;

				const T W = HPower<T>(i, j, vars);                 // W_ij = prod_g h_g^{e_g}
				// operand-Jacobian term: c_ij * W_ij * (d f_j / d x), across all columns.
				J.row(static_cast<Eigen::Index>(i)) += (c * W) * Jf.row(static_cast<Eigen::Index>(j));

				// h-power terms: only the homogenizing-variable columns pick up f_j * dW/dh_{g'}.
				if (homogenized_)
				{
					for (size_t g = 0; g < num_groups_; ++g)
					{
						const int e = Deficit(i, j, g);
						if (e <= 0 || hom_var_index_[g] < 0)
							continue;                              // e<0 means c_ij==0 (already skipped); e==0 -> no dependence
						// dW/dh_g = e * h_g^{e-1} * prod_{g2 != g} h_{g2}^{e_{g2}}
						const T dW = MixedPartial<T>(i, j, g, vars);
						J(static_cast<Eigen::Index>(i), hom_var_index_[g]) += c * f(static_cast<Eigen::Index>(j)) * dW;
					}
				}
			}
		}
	}

	/// Time-derivative.  The randomization is autonomous in the path variable; if the operand
	/// depends on t (unusual for a randomized target) its rows would contribute, but the common
	/// case is an autonomous operand, so this is zero.
	template <typename T>
	void TimeDerivInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& /*vars*/, T const& /*path_value*/) const
	{
		result.setZero();
	}

private:
	template <typename T>
	const Mat<T>& Working() const { return std::get<Mat<T>>(coefficients_working_); }

	/// A precision-correct zero / one for the accumulators.
	template <typename T>
	T Zero() const
	{
		T z(0);
		if constexpr (!std::is_same<T, dbl>::value) z.precision(precision_);
		return z;
	}

	static bool IsZero(dbl const& c) { return c == dbl(0); }
	static bool IsZero(mpfr_complex const& c) { return c.real() == 0 && c.imag() == 0; }

	/// The (signed) degree deficit of operand function j against row i's target, in group g.
	int Deficit(size_t i, size_t j, size_t g) const
	{
		return target_multidegrees_[i][g] - operand_multidegrees_[j][g];
	}

	/// W_ij = prod_g h_g^{(D_{i,g} - d_{j,g})}.  Before homogenization every h_g = 1, so 1.
	template <typename T>
	T HPower(size_t i, size_t j, Vec<T> const& vars) const
	{
		T w(1);
		if constexpr (!std::is_same<T, dbl>::value) w.precision(precision_);
		if (!homogenized_)
			return w;
		for (size_t g = 0; g < num_groups_; ++g)
		{
			const int e = Deficit(i, j, g);
			if (e == 0 || hom_var_index_[g] < 0)
				continue;
			w *= IntPow<T>(vars(hom_var_index_[g]), e);
		}
		return w;
	}

	/// d W_ij / d h_{g} = e_g * h_g^{e_g - 1} * prod_{g2 != g} h_{g2}^{e_{g2}}  (e_g > 0 assumed).
	template <typename T>
	T MixedPartial(size_t i, size_t j, size_t g, Vec<T> const& vars) const
	{
		const int eg = Deficit(i, j, g);
		T d(eg);
		if constexpr (!std::is_same<T, dbl>::value) d.precision(precision_);
		d *= IntPow<T>(vars(hom_var_index_[g]), eg - 1);
		for (size_t g2 = 0; g2 < num_groups_; ++g2)
		{
			if (g2 == g)
				continue;
			const int e2 = Deficit(i, j, g2);
			if (e2 == 0 || hom_var_index_[g2] < 0)
				continue;
			d *= IntPow<T>(vars(hom_var_index_[g2]), e2);
		}
		return d;
	}

	/// base^e for integer e >= 0, by repeated multiplication (exact, no branch-cut surprises).
	template <typename T>
	static T IntPow(T const& base, int e)
	{
		T result(1);
		if constexpr (!std::is_same<T, dbl>::value) result.precision(bertini::Precision(base));
		for (int k = 0; k < e; ++k)
			result *= base;
		return result;
	}

	/// Match a variable group against the operand's affine variable groups by identity; -1 if none.
	int GroupIndex(VariableGroup const& vars) const
	{
		auto const& groups = operand_->VariableGroups();
		for (size_t g = 0; g < groups.size(); ++g)
			if (SameGroup(groups[g], vars))
				return static_cast<int>(g);
		return -1;
	}
	static bool SameGroup(VariableGroup const& a, VariableGroup const& b)
	{
		if (a.size() != b.size())
			return false;
		for (size_t k = 0; k < a.size(); ++k)
			if (a[k] != b[k])
				return false;
		return true;
	}

	template <typename T>
	void SyncPrecision(Vec<T> const& vars) const
	{
		if constexpr (!std::is_same<T, dbl>::value)
		{
			if (vars.size() > 0)
			{
				const unsigned p = bertini::Precision(vars(0));
				if (p != precision_)
					Precision(p);
			}
		}
	}

	void BuildWorking() const
	{
		const auto& M = coefficients_highest_precision_;
		auto& wd = std::get<Mat<dbl>>(coefficients_working_);
		auto& wm = std::get<Mat<mpfr_complex>>(coefficients_working_);
		wd.resize(M.rows(), M.cols());
		wm.resize(M.rows(), M.cols());
		for (Eigen::Index r = 0; r < M.rows(); ++r)
			for (Eigen::Index c = 0; c < M.cols(); ++c)
			{
				wd(r, c) = dbl(M(r, c));
				wm(r, c) = M(r, c);
			}
	}

	OperandPtr operand_;
	Mat<mpfr_complex> coefficients_highest_precision_;             ///< master n x N randomization matrix
	mutable std::tuple<Mat<dbl>, Mat<mpfr_complex>> coefficients_working_;
	std::vector<std::vector<int>> target_multidegrees_;            ///< n x G
	std::vector<std::vector<int>> operand_multidegrees_;           ///< N x G
	size_t num_groups_;                                            ///< G (affine variable groups)
	bool homogenized_ = false;
	std::vector<Var> hom_vars_;                                    ///< the (shared) homogenizing variable per group
	std::vector<int> hom_var_index_;                               ///< its position in the variable ordering
	mutable unsigned precision_;

	friend class boost::serialization::access;

	// Operand held as shared_ptr<SystemT> (non-const, since Homogenize mutates the twin); split
	// like BlendBlock so serialization never deserializes into a const object.
	template <typename Archive>
	void save(Archive& ar, const unsigned /*version*/) const
	{
		ar & operand_;
		ar & coefficients_highest_precision_;
		ar & std::get<0>(coefficients_working_);
		ar & std::get<1>(coefficients_working_);
		ar & target_multidegrees_;
		ar & operand_multidegrees_;
		ar & num_groups_;
		ar & homogenized_;
		ar & hom_vars_;
		ar & hom_var_index_;
		ar & precision_;
	}

	template <typename Archive>
	void load(Archive& ar, const unsigned /*version*/)
	{
		ar & operand_;
		ar & coefficients_highest_precision_;
		ar & std::get<0>(coefficients_working_);
		ar & std::get<1>(coefficients_working_);
		ar & target_multidegrees_;
		ar & operand_multidegrees_;
		ar & num_groups_;
		ar & homogenized_;
		ar & hom_vars_;
		ar & hom_var_index_;
		ar & precision_;
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

} // namespace blocks
} // namespace bertini
