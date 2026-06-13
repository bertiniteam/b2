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

\brief A System evaluation block that linearly combines whole Systems with t-coefficients.

This is how a coupled / start-target homotopy is represented inside the block-composed
System: a blend

    H(x,t) = sum_i  c_i(t) * f_i(x)

of operand Systems f_i (e.g. an SLP polynomial target and a products-of-linears MHom
start) weighted by coefficient functions c_i(t) that are ordinary function-tree nodes in
the path variable t -- so (1-t), gamma*t, and later t^lifting all work, and the
time-derivative is exact via Node::Differentiate(t):

    dH/dt = sum_i c_i'(t) * f_i(x)

(the operand Systems are autonomous in t, so df_i/dt = 0).  The Jacobian is
sum_i c_i(t) * (df_i/dx).

Only the operands' NATURAL functions are blended; the homotopy System that owns this
block carries the (shared) patch, appended once after the block -- so patches are not
double-counted or scaled.

Templated on the system type so that `System` can remain a forward declaration in the
Block variant (a System contains a Block which contains a BlendBlock<System> which holds
shared_ptr<const System> -- the template parameter makes `System` a dependent name, so
its member functions are only required complete at instantiation, breaking the cycle).
Operands are held by shared_ptr, never by value, for the same reason.

\note Evaluating a coefficient node sets the shared path-variable node's value (mutable
node state); per-thread tracking must Clone the homotopy, as it already must for Systems.
The operand Systems are shared_ptr<const System>; deep-cloning them per thread is a
follow-up (block-operand Clone semantics).
*/

#pragma once

#include <vector>
#include <memory>
#include <cassert>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/split_member.hpp>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/function_tree.hpp"

namespace bertini {
namespace blocks {

template <typename SystemT>
class BlendBlock
{
public:
	using Nd = std::shared_ptr<node::Node>;
	using Var = std::shared_ptr<node::Variable>;
	using OperandPtr = std::shared_ptr<const SystemT>;

	BlendBlock() : precision_(DefaultPrecision()) {}

	/**
	\param path_variable The shared path variable t the coefficients are functions of.
	\param coefficients  One coefficient node c_i(t) per operand.
	\param operands      The Systems being blended; all must share NumNaturalFunctions().
	*/
	BlendBlock(Var path_variable, std::vector<Nd> coefficients, std::vector<OperandPtr> operands)
		: path_variable_(std::move(path_variable)),
		  coefficients_(std::move(coefficients)),
		  operands_(std::move(operands)),
		  precision_(DefaultPrecision())
	{
		assert(!operands_.empty() && "BlendBlock needs at least one operand");
		assert(coefficients_.size() == operands_.size() && "one coefficient per operand");
		derivative_coefficients_.reserve(coefficients_.size());
		for (auto const& c : coefficients_)
			derivative_coefficients_.push_back(c->Differentiate(path_variable_));
	}

	/// The number of (natural) functions the blend contributes; the owning System adds any patch.
	size_t NumFunctions() const
	{
		return operands_.empty() ? 0 : operands_.front()->NumNaturalFunctions();
	}

	bool DependsOnPathVariable() const { return true; }

	bool HasConstantJacobian() const { return false; }

	unsigned Precision() const { return precision_; }

	void Precision(unsigned new_precision) const
	{
		for (auto const& op : operands_)
			op->precision(new_precision);
		// The coefficient nodes (and the shared path variable) must move too, or a blend
		// of a low-precision operand value with a high-precision coefficient yields a
		// high-precision result that the tracker then carries as the path point, mismatching
		// the system's working precision.
		if (path_variable_)
			path_variable_->precision(new_precision);
		for (auto const& c : coefficients_)
			c->precision(new_precision);
		for (auto const& c : derivative_coefficients_)
			c->precision(new_precision);
		precision_ = new_precision;
	}

	/// H(x,t) = sum_i c_i(t) * f_i(x)   (natural functions only)
	template <typename T>
	void EvalInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& path_value) const
	{
		SyncPrecision(vars);
		result.setZero();
		const Eigen::Index k = static_cast<Eigen::Index>(NumFunctions());
		for (size_t i = 0; i < operands_.size(); ++i)
		{
			const Vec<T> fi = operands_[i]->template Eval<T>(vars);
			result += EvalNode<T>(coefficients_[i], path_value) * fi.head(k);
		}
	}

	/// dH/dx = sum_i c_i(t) * (df_i/dx)   (natural rows only)
	template <typename T>
	void JacobianInPlace(Eigen::Ref<Mat<T>> J, Vec<T> const& vars, T const& path_value) const
	{
		SyncPrecision(vars);
		J.setZero();
		const Eigen::Index k = static_cast<Eigen::Index>(NumFunctions());
		for (size_t i = 0; i < operands_.size(); ++i)
		{
			const Mat<T> Ji = operands_[i]->template Jacobian<T>(vars);
			J += EvalNode<T>(coefficients_[i], path_value) * Ji.topRows(k);
		}
	}

	/// dH/dt = sum_i c_i'(t) * f_i(x)   (operands are autonomous in t)
	template <typename T>
	void TimeDerivInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& path_value) const
	{
		SyncPrecision(vars);
		result.setZero();
		const Eigen::Index k = static_cast<Eigen::Index>(NumFunctions());
		for (size_t i = 0; i < operands_.size(); ++i)
		{
			const Vec<T> fi = operands_[i]->template Eval<T>(vars);
			result += EvalNode<T>(derivative_coefficients_[i], path_value) * fi.head(k);
		}
	}

private:
	/// Bring the operand systems to the precision of the evaluation point, so their
	/// SetVariables precision checks pass as the adaptive tracker changes precision.
	/// (No-op for double.)
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
	std::vector<OperandPtr> operands_;
	mutable unsigned precision_;

	friend class boost::serialization::access;

	// Operands are shared_ptr<const SystemT>; serialize them as non-const pointers so we
	// never ask boost to deserialize into a const object (const-ness is a local concern).
	template <typename Archive>
	void save(Archive& ar, const unsigned /*version*/) const
	{
		ar & precision_;
		ar & path_variable_;
		ar & coefficients_;
		ar & derivative_coefficients_;
		std::vector<std::shared_ptr<SystemT>> ops;
		ops.reserve(operands_.size());
		for (auto const& o : operands_)
			ops.push_back(std::const_pointer_cast<SystemT>(o));
		ar & ops;
	}

	template <typename Archive>
	void load(Archive& ar, const unsigned /*version*/)
	{
		ar & precision_;
		ar & path_variable_;
		ar & coefficients_;
		ar & derivative_coefficients_;
		std::vector<std::shared_ptr<SystemT>> ops;
		ar & ops;
		operands_.assign(ops.begin(), ops.end());
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

} // namespace blocks
} // namespace bertini
