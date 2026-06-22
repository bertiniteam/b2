//This file is part of Bertini 2.
//
//polynomial_block.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//polynomial_block.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with polynomial_block.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/system/blocks/polynomial_block.hpp

\brief The evaluation block holding the classic polynomial path: function-tree functions
(and their derivatives) and/or the compiled straight-line program.

This is the fold of System's historical `functions_` / `slp_` / `eval_method_` machinery
into a first-class block, so a System becomes a uniform loop over blocks ("everything is a
block").  It owns the function trees, their differentiation, and the SLP, and it carries
the system's variable ordering + path variable (the function-tree Jacobian/time-derivative
evaluate node trees through them).  It is value-in: each Eval/Jacobian/TimeDeriv sets the
variables (and path value) from its arguments, then evaluates.

The SLP it can compile itself: SLPCompiler reads exactly the variable ordering, path
variable, functions, and derivatives this block owns (see SLPCompiler::Compile).
*/

#pragma once

#include <vector>
#include <memory>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/function_tree.hpp"
#include "bertini2/system/straight_line_program.hpp"

namespace bertini {
namespace blocks {

/**
\brief The polynomial evaluation block: function trees + derivatives + SLP.
*/
class PolynomialBlock
{
public:
	using Fn  = std::shared_ptr<node::Function>;
	using Nd  = std::shared_ptr<node::Node>;
	using Var = std::shared_ptr<node::Variable>;

	PolynomialBlock() : precision_(DefaultPrecision()) {}

	// ---- construction (System forwards AddFunction / AddSubFunction here) ----
	void AddFunction(Fn const& f)    { functions_.push_back(f);    Invalidate(); }
	void AddSubFunction(Fn const& f) { subfunctions_.push_back(f); Invalidate(); }
	void AddConstant(Fn const& f)    { constant_subfunctions_.push_back(f); Invalidate(); }

	/// The variable ordering + path variable the function trees are evaluated against; the
	/// owning System keeps these in sync (they change as variable groups are added / the
	/// system is homogenized).
	void SetVariableOrdering(VariableGroup const& vars) const { variables_ = vars; Invalidate(); }
	void SetPathVariable(Var const& t) const { path_variable_ = t; Invalidate(); }
	void ClearPathVariable() const { path_variable_.reset(); Invalidate(); }

	// Mutable access for the owning System's construction-time manipulations (Homogenize walks
	// the trees in place; Reorder/Simplify reassign entries).  The System keeps the variable
	// groups / ordering; the block keeps the functions and their derivatives.
	std::vector<Fn>&       Functions()       { return functions_; }
	std::vector<Fn> const& Functions() const { return functions_; }
	std::vector<Fn> const& ConstantSubfunctions() const { return constant_subfunctions_; }
	size_t NumConstants() const { return constant_subfunctions_.size(); }

	bool IsDifferentiated() const { return is_differentiated_; }
	void Invalidate() const { is_differentiated_ = false; }

	/// Per-function degrees (total, and with respect to a variable group).
	std::vector<int> Degrees() const
	{
		std::vector<int> d; d.reserve(functions_.size());
		for (auto const& f : functions_) d.push_back(f->Degree());
		return d;
	}
	std::vector<int> Degrees(VariableGroup const& vars) const
	{
		std::vector<int> d; d.reserve(functions_.size());
		for (auto const& f : functions_) d.push_back(f->Degree(vars));
		return d;
	}

	/// Homogenize each function w.r.t. the group + its homogenizing var, functionally:
	/// each function is rebound to a freshly homogenized copy, so any external holder of the
	/// original function node never observes it change (shared variables are preserved).
	void Homogenize(VariableGroup const& group, Var const& hom_var)
	{
		for (auto& f : functions_)
			f = std::static_pointer_cast<node::Function>(f->Homogenized(group, hom_var));
		Invalidate();
	}
	bool IsHomogeneous(VariableGroup const& vars) const
	{
		for (auto const& f : functions_) if (!f->IsHomogeneous(vars)) return false;
		return true;
	}
	bool IsPolynomial(VariableGroup const& vars) const
	{
		for (auto const& f : functions_) if (!f->IsPolynomial(vars)) return false;
		return true;
	}

	/// Reset the cached values in the function / derivative trees (function-tree eval path).
	void Reset() const
	{
		for (auto const& f : functions_)             f->Reset();
		for (auto const& f : subfunctions_)          f->Reset();
		for (auto const& f : constant_subfunctions_) f->Reset();
		for (auto const& n : space_derivatives_) n->Reset();
		for (auto const& n : time_derivatives_)  n->Reset();
	}

	// ---- block contract: metadata ----
	size_t NumFunctions() const { return functions_.size(); }
	bool DependsOnPathVariable() const { return static_cast<bool>(path_variable_); }

	/// Human-facing description: one line per function, `f_k = <expression>`.  Polynomials are the
	/// content, so terse and verbose are the same (the expression is shown either way).
	void Describe(std::ostream& out, size_t& row, VariableGroup const& /*vars*/, bool /*verbose*/) const
	{
		for (auto const& f : functions_)
			out << "  f_" << row++ << " = " << f->EntryNode() << "\n";
	}

	unsigned Precision() const { return precision_; }
	void Precision(unsigned new_precision) const
	{
		for (auto const& f : functions_)            f->precision(new_precision);
		for (auto const& f : subfunctions_)         f->precision(new_precision);
		for (auto const& f : constant_subfunctions_) f->precision(new_precision);
		if (path_variable_) path_variable_->precision(new_precision);
		for (auto const& v : variables_) v->precision(new_precision);
		if (is_differentiated_)
		{
			for (auto const& n : space_derivatives_) n->precision(new_precision);
			for (auto const& n : time_derivatives_)  n->precision(new_precision);
		}
		slp_.precision(new_precision);
		precision_ = new_precision;
	}

	// ---- block contract: evaluation (value-in) ----
	// The compiled SLP is the sole evaluator; the function/derivative trees survive only as the
	// thing the SLP is compiled from (and that Simplify/Differentiate operate on).
	template <typename T>
	void EvalInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& path_value) const
	{
		EnsureDifferentiated();
		SetValues<T>(vars, path_value);
		slp_.template GetFuncValsInPlace<T>(result);
	}

	template <typename T>
	void JacobianInPlace(Eigen::Ref<Mat<T>> J, Vec<T> const& vars, T const& path_value) const
	{
		EnsureDifferentiated();
		SetValues<T>(vars, path_value);
		slp_.template GetJacobianInPlace<T>(J);
	}

	template <typename T>
	void TimeDerivInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& path_value) const
	{
		if (!path_variable_) { result.setZero(); return; }
		EnsureDifferentiated();
		SetValues<T>(vars, path_value);
		slp_.template GetTimeDerivInPlace<T>(result);
	}

	// ---- accessors the SLP compiler reads (mirror the System names) ----
	VariableGroup const& VariableOrdering() const { return variables_; }
	bool HavePathVariable() const { return static_cast<bool>(path_variable_); }
	Var GetPathVariable() const { return path_variable_; }
	size_t NumNaturalFunctions() const { return functions_.size(); }
	std::vector<Fn> const& GetNaturalFunctions() const { return functions_; }
	// The SLP is built from the explicit per-variable derivative trees (space/time derivatives).
	std::vector<Nd> const& GetSpaceDerivatives() const
	{
		if (space_derivatives_.empty())
			DifferentiateUsingDerivatives();
		return space_derivatives_;
	}
	std::vector<Nd> const& GetTimeDerivatives() const
	{
		if (path_variable_ && time_derivatives_.empty())
			DifferentiateUsingDerivatives();
		return time_derivatives_;
	}

	void SetAutoSimplify(bool b) const { auto_simplify_ = b; Invalidate(); }

	/// Build the symbolic derivatives (and, for SLP eval, compile the SLP from this block).
	void Differentiate() const
	{
		if (is_differentiated_) return;
		DifferentiateUsingDerivatives();
		is_differentiated_ = true;  // set before SimplifyDerivatives/Compile, which read the deriv state
		if (auto_simplify_)
			SimplifyDerivatives();
		slp_ = SLPCompiler().Compile(*this);
	}

	/// Simplify the function trees (and invalidate the derivatives, which must be rebuilt).
	void SimplifyFunctions() const
	{
		using bertini::Simplify;
		// functional (non-mutating) simplify: rebind each function to its simplified form.
		// (Handle::Simplified() currently returns self -- a Function is an opaque boundary --
		// so this is inert for top-level Function wrappers, as it has always been; the
		// machinery is now functional and ready for when that changes.)
		for (auto& f : functions_)
			f = std::static_pointer_cast<node::Function>(Simplify(f));
		Invalidate();
	}

	/// Simplify the derivative trees in place (evaluated at a random point, per the System's
	/// historical logic, to drive the simplifier).  Requires the derivatives to already exist.
	void SimplifyDerivatives() const
	{
		using bertini::Simplify;

		const size_t num_vars = variables_.size();
		std::vector<dbl> old_vals(num_vars); dbl old_path_var_val{};
		for (size_t ii = 0; ii < num_vars; ++ii)
		{
			old_vals[ii] = variables_[ii]->template Eval<dbl>();
			variables_[ii]->template SetToRandUnit<dbl>();
		}
		if (path_variable_)
		{
			old_path_var_val = path_variable_->template Eval<dbl>();
			path_variable_->template SetToRandUnit<dbl>();
		}

		for (auto const& n : space_derivatives_) n->Reset();
		for (auto const& n : time_derivatives_)  n->Reset();

		for (auto& n : space_derivatives_) n = Simplify(n);
		for (auto& n : time_derivatives_)  n = Simplify(n);

		for (size_t ii = 0; ii < num_vars; ++ii)
			variables_[ii]->template set_current_value<dbl>(old_vals[ii]);
		if (path_variable_)
			path_variable_->template set_current_value<dbl>(old_path_var_val);

		for (auto const& n : space_derivatives_) n->Reset();
		for (auto const& n : time_derivatives_)  n->Reset();
	}

private:
	void EnsureDifferentiated() const { if (!is_differentiated_) Differentiate(); }

	template <typename T>
	void SetValues(Vec<T> const& vars, T const& path_value) const
	{
		slp_.SetVariableValues(vars);
		if (path_variable_) { path_variable_->template set_current_value<T>(path_value); slp_.template SetPathVariable<T>(path_value); }
	}

	void DifferentiateUsingDerivatives() const
	{
		const size_t n = functions_.size();
		space_derivatives_.resize(n * variables_.size());
		for (size_t jj = 0; jj < variables_.size(); ++jj)
			for (size_t ii = 0; ii < n; ++ii)
				space_derivatives_[ii + jj * n] = functions_[ii]->Differentiate(variables_[jj]);
		if (path_variable_)
		{
			time_derivatives_.resize(n);
			for (size_t ii = 0; ii < n; ++ii)
				time_derivatives_[ii] = functions_[ii]->Differentiate(path_variable_);
		}
	}

	mutable std::vector<Fn> functions_;
	mutable std::vector<Fn> subfunctions_;
	mutable std::vector<Fn> constant_subfunctions_;

	mutable std::vector<Nd>  space_derivatives_;
	mutable std::vector<Nd>  time_derivatives_;
	mutable StraightLineProgram slp_;

	mutable VariableGroup variables_;   ///< the system's variable ordering (kept in sync by the owning System)
	mutable Var path_variable_;         ///< the path variable, or null

	mutable bool auto_simplify_ = false;  ///< kept in sync with the owning System's auto_simplify_
	mutable bool is_differentiated_ = false;
	mutable unsigned precision_;

	friend class boost::serialization::access;

	template <typename Archive>
	void serialize(Archive& ar, const unsigned /*version*/)
	{
		ar & constant_subfunctions_;
		ar & subfunctions_;
		ar & functions_;
		ar & is_differentiated_;
		ar & space_derivatives_;
		ar & time_derivatives_;
		ar & slp_;
		ar & variables_;
		ar & path_variable_;
		ar & precision_;
	}
};

} // namespace blocks
} // namespace bertini
