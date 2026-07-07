//This file is part of Bertini 2.
//
//canonical_encoding.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//canonical_encoding.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with canonical_encoding.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file canonical_encoding.cpp

\brief The canonical-encoding visitor (ADR-0042).  Every emitted byte here is part of
the persistent digest contract: changes are digest-breaking and require a version bump
of the system-encoding header plus a golden-fixture update in the same commit.
*/

#include "bertini2/function_tree/canonical_encoding.hpp"

#include "bertini2/function_tree.hpp"
#include "bertini2/detail/visitor.hpp"

#include <sstream>
#include <stdexcept>

namespace bertini {
namespace node {

namespace {

	// A name as a netstring: "<byte-length>:<bytes>".  Unambiguous for any UTF-8 name
	// (emoji included), with no escaping rules to get wrong.
	void EmitName(std::ostream& out, std::string const& name)
	{
		out << name.size() << ':' << name;
	}

	// An mpfr value at full stored precision, in a pinned format (scientific, all digits).
	template<typename RealT>
	std::string MpfrDigits(RealT const& v)
	{
		return v.str(0, std::ios::scientific);
	}


	// The encoding traversal.  Mirrors the SLPCompiler visitor's node universe: if a node
	// kind is compilable, it is encodable -- and an unknown kind throws (fail loud; a
	// silent skip would mint wrong digests).
	class CanonicalEncoder : public VisitorBase,
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
		CanonicalEncoder(std::ostream& out, EncodingContext& ctx) : out_(out), ctx_(ctx)
		{}

		// Encode one node: back-reference if already seen, else index it and dispatch.
		void Encode(std::shared_ptr<Node> const& n)
		{
			auto it = ctx_.seen.find(n.get());
			if (it != ctx_.seen.end())
			{
				out_ << '#' << it->second;
				return;
			}
			ctx_.seen.emplace(n.get(), ctx_.next_index++);
			n->Accept(*this);
		}

		virtual void Visit(node::Variable const& n) override
		{
			out_ << "(var ";
			EmitName(out_, n.name());
			out_ << ')';
		}

		virtual void Visit(node::Differential const& n) override
		{
			out_ << "(diff ";
			EmitName(out_, n.GetVariable()->name());
			out_ << ')';
		}

		virtual void Visit(node::Integer const& n) override
		{
			out_ << "(int " << n.GetValue().str() << ')';
		}

		virtual void Visit(node::Rational const& n) override
		{
			out_ << "(rat " << n.GetValueReal().str() << ' ' << n.GetValueImag().str() << ')';
		}

		virtual void Visit(node::Complex const& n) override
		{
			auto const& v = n.GetValue();
			out_ << "(cplx " << v.precision() << ' '
			     << MpfrDigits(v.real()) << ' ' << MpfrDigits(v.imag()) << ')';
		}

		virtual void Visit(node::special_number::Pi const& /*n*/) override
		{
			out_ << "(pi)";
		}

		virtual void Visit(node::special_number::E const& /*n*/) override
		{
			out_ << "(e)";
		}

		virtual void Visit(node::NamedExpression const& n) override
		{
			out_ << "(named ";
			EmitName(out_, n.name());
			out_ << ' ';
			Encode(n.EntryNode());
			out_ << ')';
		}

		virtual void Visit(node::SumOperator const& n) override
		{
			out_ << "(sum";
			auto const& operands = n.Operands();
			auto const& signs = n.GetSigns();
			for (size_t ii = 0; ii < operands.size(); ++ii)
			{
				out_ << ' ' << (signs[ii] ? '+' : '-');
				Encode(operands[ii]);
			}
			out_ << ')';
		}

		virtual void Visit(node::MultOperator const& n) override
		{
			out_ << "(mul";
			auto const& operands = n.Operands();
			auto const& mult_or_div = n.GetMultOrDiv();
			for (size_t ii = 0; ii < operands.size(); ++ii)
			{
				out_ << ' ' << (mult_or_div[ii] ? '*' : '/');
				Encode(operands[ii]);
			}
			out_ << ')';
		}

		virtual void Visit(node::PowerOperator const& n) override
		{
			out_ << "(pow ";
			Encode(n.GetBase());
			out_ << ' ';
			Encode(n.GetExponent());
			out_ << ')';
		}

		virtual void Visit(node::IntegerPowerOperator const& n) override
		{
			out_ << "(ipow " << n.exponent() << ' ';
			Encode(n.Operand());
			out_ << ')';
		}

		virtual void Visit(node::NegateOperator const& n) override { EncodeUnary("neg", n); }
		virtual void Visit(node::SqrtOperator const& n) override { EncodeUnary("sqrt", n); }
		virtual void Visit(node::ExpOperator const& n) override { EncodeUnary("exp", n); }
		virtual void Visit(node::LogOperator const& n) override { EncodeUnary("log", n); }
		virtual void Visit(node::SinOperator const& n) override { EncodeUnary("sin", n); }
		virtual void Visit(node::ArcSinOperator const& n) override { EncodeUnary("asin", n); }
		virtual void Visit(node::CosOperator const& n) override { EncodeUnary("cos", n); }
		virtual void Visit(node::ArcCosOperator const& n) override { EncodeUnary("acos", n); }
		virtual void Visit(node::TanOperator const& n) override { EncodeUnary("tan", n); }
		virtual void Visit(node::ArcTanOperator const& n) override { EncodeUnary("atan", n); }

	private:
		void EncodeUnary(char const* tag, node::UnaryOperator const& n)
		{
			out_ << '(' << tag << ' ';
			Encode(n.Operand());
			out_ << ')';
		}

		std::ostream& out_;
		EncodingContext& ctx_;
	};

} // unnamed namespace


void EncodeCanonical(std::shared_ptr<Node> const& root, std::ostream& out, EncodingContext& ctx)
{
	if (!root)
		throw std::invalid_argument("EncodeCanonical: null root");
	CanonicalEncoder encoder(out, ctx);
	encoder.Encode(root);
}

std::string CanonicalEncoding(std::shared_ptr<Node> const& root)
{
	std::ostringstream out;
	EncodingContext ctx;
	EncodeCanonical(root, out, ctx);
	return out.str();
}

} // namespace node
} // namespace bertini
