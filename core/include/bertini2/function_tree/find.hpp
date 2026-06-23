//This file is part of Bertini 2.
//
//include/bertini2/function_tree/find.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/function_tree/find.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/function_tree/find.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire

/**
\file include/bertini2/function_tree/find.hpp

\brief Discover named nodes of a given kind in a function tree (sympy's `find`).

`Find<T>(root)` returns the distinct nodes of kind `T` appearing anywhere in the tree,
de-duplicated by identity and sorted by name.  It is the general form of the old
`GatherVariables`: `Find<Variable>` are the solve-variables, `Find<NamedExpression>` the
named subexpressions (`a = x^2+y^2`), `Find<...>` the core's named symbols.  Descent uses
only public child accessors, so every operator subtype is handled via its base, and it
descends through Handle/NamedExpression into the wrapped expression.

The definitions live in find.cpp with explicit instantiations (one per kind we discover),
to keep this traversal out of every translation unit.
*/

#pragma once

#include <memory>
#include <vector>

#include "bertini2/function_tree/forward_declares.hpp"

namespace bertini {
namespace node {

	/// All distinct nodes of kind T in the subtree rooted at n, sorted by name.
	template<typename T>
	std::vector<std::shared_ptr<T>> Find(std::shared_ptr<const Node> const& n);

	/// The union of Find<T> across several roots (one shared traversal), sorted by name.
	template<typename T>
	std::vector<std::shared_ptr<T>> Find(std::vector<std::shared_ptr<const Node>> const& roots);

} // namespace node
} // namespace bertini
