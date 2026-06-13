//This file is part of Bertini 2.
//
//blend_block.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//blend_block.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with blend_block.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/system/blocks/blend_block.hpp

\brief A System evaluation block that linearly combines operands with t-coefficients.

This is how a coupled / start-target homotopy is represented inside the block-composed
System: a blend

    H(x,t) = sum_i  c_i(t) * f_i(x,t)

of operand evaluables f_i (e.g. an SLP polynomial target and a products-of-linears
start) weighted by coefficient functions c_i(t).  The coefficients are ordinary
function-tree nodes in the path variable t -- so (1-t), gamma*t, and later t^lifting
all work, and the time-derivative is exact via Node::Differentiate(t):

    dH/dt = sum_i [ c_i'(t) * f_i(x,t) + c_i(t) * df_i/dt ]

(the second term vanishes for autonomous operands like products-of-linears).  The
Jacobian is sum_i c_i(t) * (df_i/dx).

The block is templated on the operand type (duck-typed to the evaluation-block
contract); for the coupled homotopy the operands are Systems.  Nothing here is type-
erased -- it is a plain class template combined by the std::variant of the System.

\note Evaluating the coefficient nodes sets the shared path-variable node's value, which
is mutable node state; per-thread tracking must Clone the homotopy, as it already must
for Systems (see the system-shallow-copy-vs-threads note).  For mpfr evaluation the
caller's working precision (ThreadPrecision/DefaultPrecision) must match the operands'.
*/

#pragma once

#include <vector>
#include <memory>
#include <cassert>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/function_tree.hpp"

namespace bertini {
namespace blocks {

template <typename Operand>
class BlendBlock
{
public:
	using Nd = std::shared_ptr<node::Node>;
	using Var = std::shared_ptr<node::Variable>;

	BlendBlock() : precision_(DefaultPrecision()) {}

	/**
	\param path_variable The shared path variable t the coefficients are functions of.
	\param coefficients  One coefficient node c_i(t) per operand.
	\param operands      The evaluables being blended; all must share NumFunctions().
	*/
	BlendBlock(Var path_variable, std::vector<Nd> coefficients, std::vector<Operand> operands)
		: path_variable_(std::move(path_variable)),
		  coefficients_(std::move(coefficients)),
		  operands_(std::move(operands)),
		  precision_(DefaultPrecision())
	{
		assert(!operands_.empty() && "BlendBlock needs at least one operand");
		assert(coefficients_.size() == operands_.size() && "one coefficient per operand");
#ifndef NDEBUG
		for (auto const& op : operands_)
			assert(op.NumFunctions() == operands_.front().NumFunctions() &&
			       "all blend operands must produce the same number of functions");
#endif
		derivative_coefficients_.reserve(coefficients_.size());
		for (auto const& c : coefficients_)
			derivative_coefficients_.push_back(c->Differentiate(path_variable_));
	}

	size_t NumFunctions() const { return operands_.empty() ? 0 : operands_.front().NumFunctions(); }

	bool DependsOnPathVariable() const { return true; }

	bool HasConstantJacobian() const { return false; }

	unsigned Precision() const { return precision_; }

	void Precision(unsigned new_precision) const
	{
		for (auto const& op : operands_)
			op.Precision(new_precision);
		precision_ = new_precision;
	}

	/// H(x,t) = sum_i c_i(t) * f_i(x,t)
	template <typename T>
	void EvalInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& path_value) const
	{
		result.setZero();
		Vec<T> scratch(static_cast<Eigen::Index>(NumFunctions()));
		for (size_t i = 0; i < operands_.size(); ++i)
		{
			operands_[i].template EvalInPlace<T>(scratch, vars, path_value);
			result += EvalNode<T>(coefficients_[i], path_value) * scratch;
		}
	}

	/// dH/dx = sum_i c_i(t) * (df_i/dx)
	template <typename T>
	void JacobianInPlace(Eigen::Ref<Mat<T>> J, Vec<T> const& vars, T const& path_value) const
	{
		J.setZero();
		Mat<T> scratch(static_cast<Eigen::Index>(NumFunctions()), J.cols());
		for (size_t i = 0; i < operands_.size(); ++i)
		{
			operands_[i].template JacobianInPlace<T>(scratch, vars, path_value);
			J += EvalNode<T>(coefficients_[i], path_value) * scratch;
		}
	}

	/// dH/dt = sum_i [ c_i'(t) * f_i(x,t) + c_i(t) * df_i/dt ]
	template <typename T>
	void TimeDerivInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& path_value) const
	{
		result.setZero();
		Vec<T> scratch(static_cast<Eigen::Index>(NumFunctions()));
		for (size_t i = 0; i < operands_.size(); ++i)
		{
			operands_[i].template EvalInPlace<T>(scratch, vars, path_value);
			result += EvalNode<T>(derivative_coefficients_[i], path_value) * scratch;

			if (operands_[i].DependsOnPathVariable())
			{
				operands_[i].template TimeDerivInPlace<T>(scratch, vars, path_value);
				result += EvalNode<T>(coefficients_[i], path_value) * scratch;
			}
		}
	}

private:
	/// Evaluate a coefficient (or derivative) node at the given path-variable value.
	template <typename T>
	T EvalNode(Nd const& n, T const& path_value) const
	{
		path_variable_->template set_current_value<T>(path_value);
		n->Reset();
		return n->template Eval<T>();
	}

	Var path_variable_;
	std::vector<Nd> coefficients_;
	std::vector<Nd> derivative_coefficients_;
	std::vector<Operand> operands_;
	mutable unsigned precision_;
};

} // namespace blocks
} // namespace bertini
