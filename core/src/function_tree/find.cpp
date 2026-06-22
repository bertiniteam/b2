//This file is part of Bertini 2.
//
//src/function_tree/find.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/find.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/find.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire


#include "bertini2/function_tree/find.hpp"
#include "bertini2/function_tree.hpp"

#include <algorithm>
#include <set>

namespace bertini {
namespace node {

	namespace {

		// Recursively descend the tree, collecting distinct nodes of kind T.  A node that is a T
		// is collected AND descended into (so nested Ts -- e.g. a NamedExpression whose entry
		// contains another NamedExpression -- are all found).  Descent uses only public child
		// accessors, so every operator subtype is handled via its base class.
		template<typename T>
		void FindImpl(std::shared_ptr<const Node> const& n,
		              std::vector<std::shared_ptr<T>>& out,
		              std::set<T const*>& seen,
		              std::set<Node const*>& visited)
		{
			if (!n)
				return;

			if (!visited.insert(n.get()).second) // already processed this (possibly shared) node
				return;

			if (auto t = std::dynamic_pointer_cast<const T>(n))
				if (seen.insert(t.get()).second)
					out.push_back(std::const_pointer_cast<T>(t));

			// NamedExpression -- descend into the wrapped expression
			if (auto h = std::dynamic_pointer_cast<const NamedExpression>(n))
				FindImpl<T>(h->EntryNode(), out, seen, visited);
			// Sum, Mult, ... -- any number of operands
			else if (auto nary = std::dynamic_pointer_cast<const NaryOperator>(n))
				for (auto const& child : nary->Operands())
					FindImpl<T>(child, out, seen, visited);
			// Negate, Sqrt, Exp, Log, IntegerPower, trig, ... -- a single operand
			else if (auto u = std::dynamic_pointer_cast<const UnaryOperator>(n))
				FindImpl<T>(u->Operand(), out, seen, visited);
			// PowerOperator derives from Operator directly (base and exponent)
			else if (auto p = std::dynamic_pointer_cast<const PowerOperator>(n))
			{
				FindImpl<T>(p->GetBase(), out, seen, visited);
				FindImpl<T>(p->GetExponent(), out, seen, visited);
			}
			// numbers, special numbers, variables, differentials are leaves with no children.
		}

		template<typename T>
		std::vector<std::shared_ptr<T>>& SortByName(std::vector<std::shared_ptr<T>>& v)
		{
			std::sort(v.begin(), v.end(),
			          [](std::shared_ptr<T> const& a, std::shared_ptr<T> const& b)
			          { return a->name() < b->name(); });
			return v;
		}

	} // unnamed namespace


	template<typename T>
	std::vector<std::shared_ptr<T>> Find(std::shared_ptr<const Node> const& n)
	{
		std::vector<std::shared_ptr<T>> out;
		std::set<T const*> seen;
		std::set<Node const*> visited;
		FindImpl<T>(n, out, seen, visited);
		return SortByName(out);
	}

	template<typename T>
	std::vector<std::shared_ptr<T>> Find(std::vector<std::shared_ptr<const Node>> const& roots)
	{
		std::vector<std::shared_ptr<T>> out;
		std::set<T const*> seen;
		std::set<Node const*> visited;
		for (auto const& r : roots)
			FindImpl<T>(r, out, seen, visited);
		return SortByName(out);
	}

	// One explicit instantiation per discoverable kind (keeps the traversal in this TU).
	template std::vector<std::shared_ptr<Variable>> Find<Variable>(std::shared_ptr<const Node> const&);
	template std::vector<std::shared_ptr<Variable>> Find<Variable>(std::vector<std::shared_ptr<const Node>> const&);
	template std::vector<std::shared_ptr<NamedExpression>> Find<NamedExpression>(std::shared_ptr<const Node> const&);
	template std::vector<std::shared_ptr<NamedExpression>> Find<NamedExpression>(std::vector<std::shared_ptr<const Node>> const&);

} // namespace node
} // namespace bertini
