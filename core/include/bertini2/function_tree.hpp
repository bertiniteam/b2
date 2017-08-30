//This file is part of Bertini 2.
//
//function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function_tree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function_tree.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, University of Wisconsin Eau Claire
// jeb collins, west texas a&m

/**
\file function_tree.hpp 

\brief Collects the various header files which define the Bertini2 function tree node types.
*/


#ifndef BERTINI_FUNCTION_TREE_HPP
#define BERTINI_FUNCTION_TREE_HPP

#include "bertini2/function_tree/node.hpp"

#include "bertini2/function_tree/operators/operator.hpp"

#include "bertini2/function_tree/operators/arithmetic.hpp"
#include "bertini2/function_tree/operators/trig.hpp"


#include "bertini2/function_tree/symbols/symbol.hpp"
#include "bertini2/function_tree/symbols/variable.hpp"
#include "bertini2/function_tree/symbols/number.hpp"
#include "bertini2/function_tree/symbols/special_number.hpp"

#include "bertini2/function_tree/roots/function.hpp"
#include "bertini2/function_tree/roots/jacobian.hpp"

#include "bertini2/function_tree/factory.hpp"

#include "bertini2/function_tree/simplify.hpp"

#endif




