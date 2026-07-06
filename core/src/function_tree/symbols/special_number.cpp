//This file is part of Bertini 2.
//
//src/function_tree/symbols/special_number.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/symbols/special_number.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/symbols/special_number.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame
// Jeb Collins, West Texas A&M


#include "bertini2/function_tree/symbols/special_number.hpp"

#include <boost/math/constants/constants.hpp>




namespace bertini{
	namespace node{
		namespace special_number{
using ::pow;

// Return value of constant.
//
// pi is computed from the authoritative source -- Boost.Math's `pi` constant, which for the mpfr
// backend defers to MPFR's `mpfr_const_pi` (correctly rounded to the working precision, and cached)
// -- rather than `acos(-1)`, whose result is not guaranteed correctly rounded and varies with the
// inverse-cosine implementation.  This keeps pi identical across ranks/runs at a given precision,
// which matters because the Cauchy endgame's roots of unity and the total-degree start points are
// built from pi.  See https://github.com/bertiniteam/b2/issues/156.












			}// special number namespace


std::shared_ptr<Node> Pi()
{
	return special_number::Pi::Make();
}

std::shared_ptr<Node> E()
{
	return special_number::E::Make();
}

std::shared_ptr<Node> I()
{
	return Complex::Make(0,1);
}


std::shared_ptr<Node> Two()
{
	return Integer::Make(2);
}

std::shared_ptr<Node> One()
{
	return Integer::Make(1);
}

std::shared_ptr<Node> Zero()
{
	return Integer::Make(0);
}
	} // re: namespace node
} // re: namespace bertini
