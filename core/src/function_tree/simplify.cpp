//This file is part of Bertini 2.
//
//b2/core/src/bertini2/function_tree/simplify.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//b2/core/src/bertini2/function_tree/simplify.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with b2/core/src/bertini2/function_tree/simplify.cpp.  If not, see <http://www.gnu.org/licenses/>.
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





#include "bertini2/function_tree/simplify.hpp"

namespace bertini {

unsigned Simplify(std::shared_ptr<bertini::node::Node> const& n)
{
	unsigned num_reductions = 0;
	unsigned num_ones = 0;
	unsigned num_zeros = 0;

	num_reductions = n->ReduceDepth();
	num_ones = n->EliminateOnes();
	num_zeros = n->EliminateZeros();

	unsigned num_rounds{0};
	while (num_reductions || num_ones || num_zeros)
	{
		++num_rounds;
		num_reductions = n->ReduceDepth();
		num_ones = n->EliminateOnes();
		num_zeros = n->EliminateZeros();
	}
	return num_rounds;
}

} // namespace bertini

