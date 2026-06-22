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
#include "bertini2/function_tree/find.hpp"
#include "bertini2/function_tree.hpp"

namespace bertini {
namespace node {

	// GatherVariables is the variable-specific case of the general Find (sympy's find): the
	// one traversal lives in find.cpp, and these delegate to Find<Variable>.

	VariableGroup GatherVariables(std::shared_ptr<const Node> const& n)
	{
		return Find<Variable>(n);
	}

	VariableGroup GatherVariables(std::vector<std::shared_ptr<Node>> const& functions)
	{
		std::vector<std::shared_ptr<const Node>> roots(functions.begin(), functions.end());
		return Find<Variable>(roots);
	}

} // namespace node
} // namespace bertini
