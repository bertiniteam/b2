//This file is part of Bertini 2.
//
//include/bertini2/function_tree/gather.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/function_tree/gather.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/function_tree/gather.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire


/**
\file include/bertini2/function_tree/gather.hpp

\brief Free functions for discovering the variables appearing in a function tree.
*/

#pragma once

#include <memory>
#include <vector>

#include "bertini2/function_tree/node.hpp"

namespace bertini {
namespace node {

	/**
	\brief Collect the distinct Variables appearing in a function-tree subtree.

	Traverses the tree non-intrusively, via the public child accessors of the
	operator and root node types.  The returned group is de-duplicated (by node
	identity) and ordered alphabetically by Variable name.

	\param n The root of the subtree to scan.  May be null, in which case an empty
	         group is returned.
	\return A VariableGroup of the distinct variables found, sorted by name.
	*/
	VariableGroup GatherVariables(std::shared_ptr<const Node> const& n);

	/**
	\brief Collect the distinct Variables appearing across a collection of functions.

	The union of the variables of each function, de-duplicated (by node identity)
	and ordered alphabetically by Variable name.  This is what the
	function-list System constructor uses to auto-build a single variable group.

	\param functions The functions to scan.
	\return A VariableGroup of the distinct variables found, sorted by name.
	*/
	VariableGroup GatherVariables(std::vector<std::shared_ptr<Function>> const& functions);

} // namespace node
} // namespace bertini
