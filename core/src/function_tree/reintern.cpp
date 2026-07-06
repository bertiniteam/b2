//This file is part of Bertini 2.
//
//reintern.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//reintern.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with reintern.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file reintern.cpp

\brief The re-interning rebuilder (ADR-0042): a visitor over the closed node universe
(ADR-0014) that reconstructs each loaded node through its ordinary Make factory with
already-reinterned children, so the result lands in (or comes from) the live intern
tables.
*/

#include "bertini2/function_tree/reintern.hpp"

#include "bertini2/function_tree.hpp"
#include "bertini2/detail/visitor.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace bertini {
namespace node {

namespace {

	using Nd = std::shared_ptr<Node>;

	// Rebuilds one node (children first, via the shared memo) and leaves the interned
	// replacement in result_.  Mirrors the canonical encoder's node universe.
	class ReinternVisitor : public VisitorBase,
			// symbols and roots
			public Visitor<node::Variable>,
			public Visitor<node::Integer>,
			public Visitor<node::Complex>,
			public Visitor<node::Rational>,
			public Visitor<node::NamedExpression>,
			public Visitor<node::Differential>,

			// arithmetic
			public Visitor<node::SumOperator>,
			public Visitor<node::MultOperator>,
			public Visitor<node::IntegerPowerOperator>,
			public Visitor<node::PowerOperator>,
			public Visitor<node::ExpOperator>,
			public Visitor<node::LogOperator>,
			public Visitor<node::NegateOperator>,
			public Visitor<node::SqrtOperator>,

			// the trig operators
			public Visitor<node::SinOperator>,
			public Visitor<node::ArcSinOperator>,
			public Visitor<node::CosOperator>,
			public Visitor<node::ArcCosOperator>,
			public Visitor<node::TanOperator>,
			public Visitor<node::ArcTanOperator>,

			public Visitor<node::special_number::Pi>,
			public Visitor<node::special_number::E>
	{
	public:
		explicit ReinternVisitor(ReinternMemo& memo) : memo_(memo)
		{}

		Nd Rebuild(Nd const& n)
		{
			auto it = memo_.find(n.get());
			if (it != memo_.end())
				return it->second;
			n->Accept(*this);
			memo_.emplace(n.get(), result_);
			return result_;
		}

		virtual void Visit(node::Variable const& n) override
		{
			result_ = Variable::Make(n.name());       // unifies with the live canonical variable
		}

		virtual void Visit(node::Differential const& n) override
		{
			auto v = std::static_pointer_cast<Variable>(
				Rebuild(std::const_pointer_cast<Node>(
					std::static_pointer_cast<Node const>(n.GetVariable()))));
			result_ = Differential::Make(v, n.name());
		}

		virtual void Visit(node::Integer const& n) override
		{
			result_ = Integer::Make(n.GetValue());
		}

		virtual void Visit(node::Rational const& n) override
		{
			result_ = Rational::Make(n.GetValueReal(), n.GetValueImag());
		}

		virtual void Visit(node::Complex const& n) override
		{
			result_ = Complex::Make(n.GetValue());
		}

		virtual void Visit(node::special_number::Pi const& /*n*/) override
		{
			result_ = Pi();
		}

		virtual void Visit(node::special_number::E const& /*n*/) override
		{
			result_ = E();
		}

		virtual void Visit(node::NamedExpression const& n) override
		{
			result_ = NamedExpression::Make(Rebuild(n.EntryNode()), n.name());
		}

		virtual void Visit(node::SumOperator const& n) override
		{
			std::vector<std::pair<Nd, bool>> terms;
			auto const& operands = n.Operands();
			auto const& signs = n.GetSigns();
			terms.reserve(operands.size());
			for (size_t ii = 0; ii < operands.size(); ++ii)
				terms.emplace_back(Rebuild(operands[ii]), static_cast<bool>(signs[ii]));
			result_ = SumOperator::Make(terms);
		}

		virtual void Visit(node::MultOperator const& n) override
		{
			std::vector<std::pair<Nd, bool>> factors;
			auto const& operands = n.Operands();
			auto const& mult_or_div = n.GetMultOrDiv();
			factors.reserve(operands.size());
			for (size_t ii = 0; ii < operands.size(); ++ii)
				factors.emplace_back(Rebuild(operands[ii]), static_cast<bool>(mult_or_div[ii]));
			result_ = MultOperator::Make(factors);
		}

		virtual void Visit(node::PowerOperator const& n) override
		{
			auto base = Rebuild(n.GetBase());
			auto exponent = Rebuild(n.GetExponent());
			result_ = PowerOperator::Make(base, exponent);
		}

		virtual void Visit(node::IntegerPowerOperator const& n) override
		{
			result_ = IntegerPowerOperator::Make(Rebuild(n.Operand()), n.exponent());
		}

		virtual void Visit(node::NegateOperator const& n) override { result_ = NegateOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::SqrtOperator const& n) override { result_ = SqrtOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::ExpOperator const& n) override { result_ = ExpOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::LogOperator const& n) override { result_ = LogOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::SinOperator const& n) override { result_ = SinOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::ArcSinOperator const& n) override { result_ = ArcSinOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::CosOperator const& n) override { result_ = CosOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::ArcCosOperator const& n) override { result_ = ArcCosOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::TanOperator const& n) override { result_ = TanOperator::Make(Rebuild(n.Operand())); }
		virtual void Visit(node::ArcTanOperator const& n) override { result_ = ArcTanOperator::Make(Rebuild(n.Operand())); }

	private:
		ReinternMemo& memo_;
		Nd result_;
	};

} // unnamed namespace


std::shared_ptr<Node> Reintern(std::shared_ptr<Node> const& root, ReinternMemo& memo)
{
	if (!root)
		throw std::invalid_argument("Reintern: null root");
	ReinternVisitor visitor(memo);
	return visitor.Rebuild(root);
}

std::shared_ptr<Node> Reintern(std::shared_ptr<Node> const& root)
{
	ReinternMemo memo;
	return Reintern(root, memo);
}

} // namespace node
} // namespace bertini
