//This file is part of Bertini 2.
//
//content_identity.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//content_identity.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with content_identity.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file content_identity.cpp

\brief System::CanonicalEncodingText / ContentDigest / Hash / IsSame (ADR-0042).

The System-level canonical encoding: format version `b2sysenc/1`.  Every byte emitted
here is part of the persistent digest contract -- a change to this encoding must bump
the version header and the golden fixture (system_identity_test) in the same commit.

Included (everything evaluation-relevant; randomness IS identity): session
canonicalization settings, variable-group structure + time order, homogenizing
variables, path variable, parameters, every block's exact payload (poly functions;
linear-forms / products-of-linears / randomization coefficient matrices; blend
coefficients incl. gamma; operand systems recursively), the patch, and the
pre-homogenization function snapshot.  Excluded (transient/derived): working
precision, current variable/path values, differentiation caches, SLPs, and the
derivable variable-ordering cache.
*/

#include "bertini2/system/system.hpp"

#include "bertini2/function_tree/canonical.hpp"
#include "bertini2/function_tree/canonical_encoding.hpp"

#include <cstring>
#include <sstream>

namespace bertini {

namespace {

	// A name as a netstring, matching the node encoder's convention.
	void EmitName(std::ostream& out, std::string const& name)
	{
		out << name.size() << ':' << name;
	}

	// One exact complex_mp entry: stored precision + full digits (pinned scientific format).
	// Value-equal entries at different precisions encode differently by design (ADR-0042).
	void EmitComplexEntry(std::ostream& out, complex_mp const& z)
	{
		out << z.precision() << ' '
		    << z.real().str(0, std::ios::scientific) << ' '
		    << z.imag().str(0, std::ios::scientific);
	}

	// An exact matrix: dims then entries in row-major order.
	void EmitComplexMatrix(std::ostream& out, Mat<complex_mp> const& m)
	{
		out << m.rows() << 'x' << m.cols();
		for (Eigen::Index r = 0; r < m.rows(); ++r)
			for (Eigen::Index c = 0; c < m.cols(); ++c)
			{
				out << ' ';
				EmitComplexEntry(out, m(r, c));
			}
	}

	// A group of variables, by name, in stored (FIFO-semantic) order.
	void EmitVariableGroup(std::ostream& out, VariableGroup const& group)
	{
		out << group.size();
		for (auto const& v : group)
		{
			out << ' ';
			EmitName(out, v->name());
		}
	}

	char GroupTag(VariableGroupType t)
	{
		switch (t)
		{
			case VariableGroupType::Affine: return 'A';
			case VariableGroupType::Homogeneous: return 'H';
			case VariableGroupType::Ungrouped: return 'U';
		}
		throw std::runtime_error("unencodable VariableGroupType");
	}

	char const* OrderName(node::MonomialOrder o)
	{
		switch (o)
		{
			case node::MonomialOrder::Lex: return "Lex";
			case node::MonomialOrder::RevLex: return "RevLex";
			case node::MonomialOrder::GrevLex: return "GrevLex";
		}
		throw std::runtime_error("unencodable MonomialOrder");
	}

	// Per-row multidegree tables (randomization bookkeeping).
	void EmitMultidegrees(std::ostream& out, std::vector<std::vector<int>> const& mds)
	{
		out << mds.size();
		for (auto const& row : mds)
		{
			out << " [" << row.size();
			for (auto d : row)
				out << ' ' << d;
			out << ']';
		}
	}

} // unnamed namespace


void System::EncodeCanonical(std::ostream& out, node::EncodingContext& ctx) const
{
	// 1. format version + the session-global, identity-affecting canonicalization settings
	out << "b2sysenc/1 order=" << OrderName(node::CurrentMonomialOrder())
	    << " canon=" << (node::CanonicalizeByDefault() ? 1 : 0)
	    << " powerfold=" << (node::PowerFoldByDefault() ? 1 : 0) << '\n';

	// 2. variable structure (FIFO semantics: stored order is authored order)
	out << "timeorder";
	for (auto t : time_order_of_variable_groups_)
		out << ' ' << GroupTag(t);
	out << '\n';

	out << "affine " << variable_groups_.size();
	for (auto const& g : variable_groups_)
	{
		out << ' ';
		EmitVariableGroup(out, g);
	}
	out << '\n';

	out << "hom " << hom_variable_groups_.size();
	for (auto const& g : hom_variable_groups_)
	{
		out << ' ';
		EmitVariableGroup(out, g);
	}
	out << '\n';

	out << "ungrouped ";
	EmitVariableGroup(out, ungrouped_variables_);
	out << '\n';

	out << "homvars ";
	EmitVariableGroup(out, homogenizing_variables_);
	out << '\n';

	out << "pathvar ";
	if (have_path_variable_ && path_variable_)
		EmitName(out, path_variable_->name());
	else
		out << '-';
	out << '\n';

	// 3. parameters
	out << "implicit ";
	EmitVariableGroup(out, implicit_parameters_);
	out << '\n';

	out << "explicit " << explicit_parameters_.size();
	for (auto const& p : explicit_parameters_)
	{
		out << ' ';
		node::EncodeCanonical(p, out, ctx);
	}
	out << '\n';

	// 4. blocks, in stored order, each with its exact payload
	out << "blocks " << blocks_.size() << '\n';
	for (auto const& blk : blocks_)
	{
		std::visit([&](auto const& b) {
			using BlockT = std::decay_t<decltype(b)>;
			if constexpr (std::is_same_v<BlockT, blocks::PolynomialBlock>)
			{
				auto const& functions = b.Functions();
				out << "(block poly " << functions.size();
				for (auto const& f : functions)
				{
					out << ' ';
					node::EncodeCanonical(f, out, ctx);
				}
				out << ")\n";
			}
			else if constexpr (std::is_same_v<BlockT, blocks::LinearFormsBlock>)
			{
				out << "(block linforms " << b.NumVariables()
				    << ' ' << (b.IsHomogenized() ? 1 : 0) << ' ';
				EmitComplexMatrix(out, b.Coefficients());
				out << ")\n";
			}
			else if constexpr (std::is_same_v<BlockT, blocks::ProductsOfLinearsBlock>)
			{
				auto const& factors = b.Factors();
				out << "(block prodlin " << b.NumVariables() << ' ' << factors.size();
				for (auto const& m : factors)
				{
					out << ' ';
					EmitComplexMatrix(out, m);
				}
				out << ")\n";
			}
			else if constexpr (std::is_same_v<BlockT, blocks::RandomizationBlock<System>>)
			{
				out << "(block randomization homogenized=" << (b.IsHomogenized() ? 1 : 0)
				    << " groups=" << b.NumGroups() << " R=";
				EmitComplexMatrix(out, b.RandomizationMatrix());
				out << " target=";
				EmitMultidegrees(out, b.TargetMultidegrees());
				out << " operandmd=";
				EmitMultidegrees(out, b.OperandMultidegrees());
				out << " homvars=";
				out << b.HomVars().size();
				for (auto const& v : b.HomVars())
				{
					out << ' ';
					EmitName(out, v->name());
				}
				out << " operand=(\n";
				b.Operand()->EncodeCanonical(out, ctx);  // recursion, same context
				out << "))\n";
			}
			else if constexpr (std::is_same_v<BlockT, blocks::BlendBlock<System>>)
			{
				out << "(block blend pathvar ";
				if (b.PathVariable())
					EmitName(out, b.PathVariable()->name());
				else
					out << '-';
				out << " coefficients " << b.Coefficients().size();
				for (auto const& c : b.Coefficients())
				{
					out << ' ';
					node::EncodeCanonical(c, out, ctx);
				}
				out << " operands " << b.Operands().size();
				for (auto const& op : b.Operands())
				{
					out << " (\n";
					op->EncodeCanonical(out, ctx);  // recursion, same context
					out << ')';
				}
				out << ")\n";
			}
			else
			{
				static_assert(!sizeof(BlockT), "unencodable block kind: extend the canonical "
					"encoding (and bump b2sysenc + golden fixture) when adding block types");
			}
		}, blk);
	}

	// 5. the patch (its coefficients are random -- and randomness is identity)
	out << "patched " << (is_patched_ ? 1 : 0);
	if (is_patched_)
	{
		out << " sizes " << patch_.VariableGroupSizes().size();
		for (auto s : patch_.VariableGroupSizes())
			out << ' ' << s;
		out << " coefficients " << patch_.Coefficients().size();
		for (auto const& vec : patch_.Coefficients())
		{
			out << " [" << vec.size();
			for (Eigen::Index ii = 0; ii < vec.size(); ++ii)
			{
				out << ' ';
				EmitComplexEntry(out, vec(ii));
			}
			out << ']';
		}
	}
	out << '\n';

	// 6. the pre-homogenization snapshot (observable via user-coordinate SymbolicJacobian:
	// authored-homogeneous and homogenized-from-affine are distinct identities)
	out << "prehom " << pre_homogenization_functions_.size();
	for (auto const& f : pre_homogenization_functions_)
	{
		out << ' ';
		node::EncodeCanonical(f, out, ctx);
	}
	out << '\n';
}


std::string System::CanonicalEncodingText() const
{
	std::ostringstream out;
	node::EncodingContext ctx;
	EncodeCanonical(out, ctx);
	return out.str();
}

detail::Digest256 System::ContentDigest() const
{
	return detail::Sha256(CanonicalEncodingText());
}

std::size_t System::Hash() const
{
	auto const digest = ContentDigest();
	std::size_t h = 0;
	static_assert(sizeof(h) <= 32, "Digest256 must cover std::size_t");
	std::memcpy(&h, digest.bytes.data(), sizeof(h));
	return h;
}

bool System::IsSame(System const& other) const
{
	return ContentDigest() == other.ContentDigest();
}

} // namespace bertini
