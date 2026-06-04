//This file is part of Bertini 2.
//
//src/function_tree/gather.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/gather.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/gather.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire


#include "bertini2/function_tree/gather.hpp"
#include "bertini2/function_tree.hpp"

#include <algorithm>
#include <set>

namespace bertini {
namespace node {

	namespace {

		// Recursively descend the tree, collecting distinct Variables.  Descent uses
		// only public child accessors, so every operator subtype is handled by its
		// base class (UnaryOperator / NaryOperator) without enumerating concrete types.
		void GatherImpl(std::shared_ptr<const Node> const& n,
		                std::vector<std::shared_ptr<Variable>>& ordered,
		                std::set<Variable const*>& seen_vars,
		                std::set<Node const*>& visited)
		{
			if (!n)
				return;

			if (!visited.insert(n.get()).second) // already processed this (possibly shared) node
				return;

			if (auto v = std::dynamic_pointer_cast<const Variable>(n))
			{
				if (seen_vars.insert(v.get()).second)
					ordered.push_back(std::const_pointer_cast<Variable>(v));
				return;
			}

			// Function / Jacobian -- descend into the entry (root) node
			if (auto h = std::dynamic_pointer_cast<const Handle>(n))
			{
				GatherImpl(h->EntryNode(), ordered, seen_vars, visited);
				return;
			}

			// Sum, Mult, ... -- any number of operands
			if (auto nary = std::dynamic_pointer_cast<const NaryOperator>(n))
			{
				for (auto const& child : nary->Operands())
					GatherImpl(child, ordered, seen_vars, visited);
				return;
			}

			// Negate, Sqrt, Exp, Log, IntegerPower, trig, ... -- a single operand
			if (auto u = std::dynamic_pointer_cast<const UnaryOperator>(n))
			{
				GatherImpl(u->Operand(), ordered, seen_vars, visited);
				return;
			}

			// PowerOperator derives from Operator directly (base and exponent)
			if (auto p = std::dynamic_pointer_cast<const PowerOperator>(n))
			{
				GatherImpl(p->GetBase(), ordered, seen_vars, visited);
				GatherImpl(p->GetExponent(), ordered, seen_vars, visited);
				return;
			}

			// LinearProduct stores its variables internally
			if (auto lp = std::dynamic_pointer_cast<const LinearProduct>(n))
			{
				VariableGroup lp_vars;
				lp->GetVariables(lp_vars);
				for (auto const& v : lp_vars)
					if (v && seen_vars.insert(v.get()).second)
						ordered.push_back(v);
				return;
			}

			// numbers, special numbers, differentials, etc. are leaves with no
			// solve-variables to contribute.
		}

		VariableGroup SortedByName(std::vector<std::shared_ptr<Variable>>& vars)
		{
			std::sort(vars.begin(), vars.end(),
			          [](std::shared_ptr<Variable> const& a, std::shared_ptr<Variable> const& b)
			          { return a->name() < b->name(); });
			return VariableGroup(vars.begin(), vars.end());
		}

	} // unnamed namespace


	VariableGroup GatherVariables(std::shared_ptr<const Node> const& n)
	{
		std::vector<std::shared_ptr<Variable>> ordered;
		std::set<Variable const*> seen_vars;
		std::set<Node const*> visited;

		GatherImpl(n, ordered, seen_vars, visited);

		return SortedByName(ordered);
	}


	VariableGroup GatherVariables(std::vector<std::shared_ptr<Function>> const& functions)
	{
		std::vector<std::shared_ptr<Variable>> ordered;
		std::set<Variable const*> seen_vars;
		std::set<Node const*> visited;

		for (auto const& f : functions)
			GatherImpl(f, ordered, seen_vars, visited);

		return SortedByName(ordered);
	}

} // namespace node
} // namespace bertini
