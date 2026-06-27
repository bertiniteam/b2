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
#include "bertini2/system/blocks/describe.hpp"

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

	// Memory-isolating copy (ADR-0027): each operand System is deep-copied via the
	// System copy constructor, which shares its immutable node DAG and compiled SLP Program but
	// gives it its own per-thread evaluation Memory.  The coefficient-system cache is reset (rebuilt
	// lazily per copy).  The coefficient nodes / path variable are shared --- they are never written
	// during evaluation.  This lets path tracking Clone a System into per-thread copies that share
	// no mutable state.
	BlendBlock(BlendBlock const& other)
		: path_variable_(other.path_variable_),
		  coefficients_(other.coefficients_),
		  derivative_coefficients_(other.derivative_coefficients_),
		  precision_(other.precision_)
	{
		operands_.reserve(other.operands_.size());
		for (auto const& op : other.operands_)
			operands_.push_back(std::make_shared<const SystemT>(*op));
		// coefficient_system_ deliberately left null: rebuilt lazily per copy
	}

	BlendBlock& operator=(BlendBlock const& other)
	{
		if (this != &other)
		{
			path_variable_ = other.path_variable_;
			coefficients_ = other.coefficients_;
			derivative_coefficients_ = other.derivative_coefficients_;
			precision_ = other.precision_;
			operands_.clear();
			operands_.reserve(other.operands_.size());
			for (auto const& op : other.operands_)
				operands_.push_back(std::make_shared<const SystemT>(*op));
			coefficient_system_.reset();
		}
		return *this;
	}

	BlendBlock(BlendBlock&&) = default;
	BlendBlock& operator=(BlendBlock&&) = default;

	/// The number of (natural) functions the blend contributes; the owning System adds any patch.
	size_t NumFunctions() const
	{
		return operands_.empty() ? 0 : operands_.front()->NumNaturalFunctions();
	}

	/// Accessors for the function-tree expansion (System::ExpandToFunctionTree): the blend is
	/// H = sum_i coefficients_[i](t) * operands_[i], so the expansion needs the coefficient nodes,
	/// the operand systems (each itself expanded), and the shared path variable.
	std::vector<Nd> const& Coefficients() const { return coefficients_; }
	std::vector<OperandPtr> const& Operands() const { return operands_; }
	Var const& PathVariable() const { return path_variable_; }

	/// Human-facing description: terse shows the blend `f_a..f_b = c_0(t)*A + c_1(t)*B`, with the
	/// coefficient nodes shown symbolically and the operands labelled A, B, ...; verbose lists each
	/// operand's functions indented.
	void Describe(std::ostream& out, size_t& row, VariableGroup const& /*vars*/, bool verbose) const
	{
		const size_t k = NumFunctions();
		out << "  ";
		describe_detail::PrintRowLabel(out, row, k);
		out << " = ";
		for (size_t i = 0; i < coefficients_.size(); ++i)
			out << (i ? " + " : "") << "(" << coefficients_[i] << ")*"
			    << static_cast<char>('A' + (i < 26 ? i : 25));
		out << "   (blend of " << operands_.size() << " systems)\n";
		row += k;
		if (verbose)
			for (size_t i = 0; i < operands_.size(); ++i)
			{
				const char label = static_cast<char>('A' + (i < 26 ? i : 25));
				out << "    " << label << ":\n";
				auto f = operands_[i]->NaturalFunctionsAsNodes();
				for (size_t j = 0; j < f.size(); ++j)
					out << "      " << label << "_" << j << " = " << f[j] << "\n";
			}
	}

	/// The blend c_0(t)f_0 + c_1(t)f_1 + ... has, per function, the max degree of its operands
	/// in the space variables (the t-coefficients are constant in space).
	std::vector<int> Degrees() const
	{
		std::vector<int> d(NumFunctions(), 0);
		for (auto const& op : operands_)
		{
			auto od = op->Degrees();
			for (size_t i = 0; i < d.size() && i < od.size(); ++i) d[i] = std::max(d[i], od[i]);
		}
		return d;
	}
	std::vector<int> Degrees(VariableGroup const& vars) const
	{
		std::vector<int> d(NumFunctions(), 0);
		for (auto const& op : operands_)
		{
			auto od = op->Degrees(vars);
			for (size_t i = 0; i < d.size() && i < od.size(); ++i) d[i] = std::max(d[i], od[i]);
		}
		return d;
	}

	// A blend is the coupling homotopy, built from already-prepared (homogenized) operands; its
	// operands are shared_ptr<const System> and cannot be mutated, so Homogenize is a no-op.  It
	// is homogeneous/polynomial iff all its operands are.
	void Homogenize(VariableGroup const&, std::shared_ptr<node::Variable> const&) {}
	bool IsHomogeneous(VariableGroup const&) const
	{
		for (auto const& op : operands_) if (!op->IsHomogeneous()) return false;
		return true;
	}
	bool IsPolynomial(VariableGroup const&) const
	{
		for (auto const& op : operands_) if (!op->IsPolynomial()) return false;
		return true;
	}

	bool DependsOnPathVariable() const { return true; }

	bool HasConstantJacobian() const { return false; }

	/// Analytic block: nothing symbolic to differentiate.
	void Differentiate() const {}

	unsigned Precision() const { return precision_; }

	void Precision(unsigned new_precision) const
	{
		for (auto const& op : operands_)
			op->precision(new_precision);
		// The coefficients are evaluated through the coefficient sub-system's SLP (which carries
		// its own precision; see EvalCoefficients), so the coefficient nodes and the shared path
		// variable are no longer evaluated during tracking.  Their precision is vestigial and
		// left untouched, keeping the shared node DAG read-only across threads (ADR-0027).
		if (coefficient_system_)
			coefficient_system_->precision(new_precision);
		precision_ = new_precision;
	}

	/// H(x,t) = sum_i c_i(t) * f_i(x)   (natural functions only)
	template <typename T>
	void EvalInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& path_value) const
	{
		SyncPrecision(vars);
		result.setZero();
		const Eigen::Index k = static_cast<Eigen::Index>(NumFunctions());
		const Vec<T> coeffs = EvalCoefficients<T>(path_value);
		for (size_t i = 0; i < operands_.size(); ++i)
		{
			const Vec<T> fi = operands_[i]->template Eval<T>(vars);
			result += coeffs(static_cast<Eigen::Index>(i)) * fi.head(k);
		}
	}

	/// dH/dx = sum_i c_i(t) * (df_i/dx)   (natural rows only)
	template <typename T>
	void JacobianInPlace(Eigen::Ref<Mat<T>> J, Vec<T> const& vars, T const& path_value) const
	{
		SyncPrecision(vars);
		J.setZero();
		const Eigen::Index k = static_cast<Eigen::Index>(NumFunctions());
		const Vec<T> coeffs = EvalCoefficients<T>(path_value);
		for (size_t i = 0; i < operands_.size(); ++i)
		{
			const Mat<T> Ji = operands_[i]->template Jacobian<T>(vars);
			J += coeffs(static_cast<Eigen::Index>(i)) * Ji.topRows(k);
		}
	}

	/// dH/dt = sum_i c_i'(t) * f_i(x)   (operands are autonomous in t)
	template <typename T>
	void TimeDerivInPlace(Eigen::Ref<Vec<T>> result, Vec<T> const& vars, T const& path_value) const
	{
		SyncPrecision(vars);
		result.setZero();
		const Eigen::Index k = static_cast<Eigen::Index>(NumFunctions());
		const Vec<T> coeffs = EvalCoefficients<T>(path_value);
		const size_t n = operands_.size();
		for (size_t i = 0; i < n; ++i)
		{
			const Vec<T> fi = operands_[i]->template Eval<T>(vars);
			result += coeffs(static_cast<Eigen::Index>(n + i)) * fi.head(k);
		}
	}

private:
	/// Bring the operand systems to the precision of the evaluation point, so their
	/// SetVariables precision checks pass as the adaptive tracker changes precision.
	/// (No-op for double.)
	template <typename T>
	void SyncPrecision(Vec<T> const& vars) const
	{
		if constexpr (!std::is_same<T, complex_dbl>::value)
		{
			if (vars.size() > 0)
			{
				const unsigned p = bertini::Precision(vars(0));
				if (p != precision_)
					Precision(p);
			}
		}
	}

	/// Build (once) a small System whose functions are the coefficients c_i(t) followed by
	/// the derivative coefficients c_i'(t), with the path variable t as their sole variable.
	/// Evaluating it through its SLP yields every c_i(t) and c_i'(t) in one run --- so the
	/// blend carries no node-level tree evaluation, and the coefficients' program is compiled
	/// once and reused across tracker steps rather than rebuilt per evaluation.
	SystemT& EnsureCoefficientSystem() const
	{
		if (!coefficient_system_)
		{
			auto sys = std::make_shared<SystemT>();
			for (auto const& c : coefficients_)
				sys->AddFunction(c);
			for (auto const& c : derivative_coefficients_)
				sys->AddFunction(c);
			sys->AddVariableGroup(VariableGroup{path_variable_});
			sys->precision(precision_);
			coefficient_system_ = sys;
		}
		return *coefficient_system_;
	}

	/// Evaluate [c_0(t) .. c_{n-1}(t), c_0'(t) .. c_{n-1}'(t)] at the given path-variable value.
	template <typename T>
	Vec<T> EvalCoefficients(T const& path_value) const
	{
		auto& cs = EnsureCoefficientSystem();
		if constexpr (!std::is_same<T, complex_dbl>::value)
			cs.precision(bertini::Precision(path_value));
		Vec<T> t_point(1);
		t_point(0) = path_value;
		return cs.template Eval<T>(t_point);
	}

	Var path_variable_;
	std::vector<Nd> coefficients_;
	std::vector<Nd> derivative_coefficients_;
	std::vector<OperandPtr> operands_;
	mutable unsigned precision_;

	/// Cached program for the coefficients (lazily built, see EnsureCoefficientSystem).  Not
	/// serialized: it is a pure cache, rebuilt on first evaluation after a load.
	mutable std::shared_ptr<SystemT> coefficient_system_;

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
