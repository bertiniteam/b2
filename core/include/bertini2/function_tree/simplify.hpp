//This file is part of Bertini 2.
//
//b2/core/include/bertini2/function_tree/simplify.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//b2/core/include/bertini2/function_tree/simplify.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with b2/core/include/bertini2/function_tree/simplify.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  Dani Brake
//  University of Wisconsin Eau Claire
//  ACMS
//  Fall 2017

/**
\file General methods exploting polymorphism, for simplifying a Node, and all subnodes
*/


#ifndef BERTINI_FUNCTION_TREE_SIMPLIFY_HPP
#define BERTINI_FUNCTION_TREE_SIMPLIFY_HPP

#pragma once

#include "bertini2/function_tree/node.hpp"

namespace bertini {

unsigned Simplify(std::shared_ptr<bertini::node::Node> const& n);

} // namespace bertini


#endif //include guards
