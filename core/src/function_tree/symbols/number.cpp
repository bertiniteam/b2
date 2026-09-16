//This file is part of Bertini 2.
//
//src/function_tree/symbols/number.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/symbols/number.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/symbols/number.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame
// Jeb Collins, West Texas A&M


#include "bertini2/function_tree/symbols/number.hpp"




namespace bertini{
	namespace node{
		using ::pow;
		
std::shared_ptr<Node> Number::Differentiate(std::shared_ptr<Variable> const& /*v*/) const
{
	return Integer::Make(0);
}

///////////////////
//
//  INTEGERS
//
////////////////////////


void Integer::print(std::ostream & target) const
{
	target << true_value_;
}


/////////////
//
//  Floats
//
////////////////


void Complex::print(std::ostream & target) const
{
	// real-valued floats print bare; the complex pair form is reserved for
	// genuinely complex values.  str(0) prints EVERY digit the stored value
	// carries: streaming the number directly would use the ostream's default
	// precision (6 significant digits), silently truncating every printed
	// system -- a coefficient must round-trip exactly through classic input.
	if (highest_precision_value_.imag() == 0)
		target << highest_precision_value_.real().str(0);
	else
		target << "(" << highest_precision_value_.real().str(0) << ","
		       << highest_precision_value_.imag().str(0) << ")";
}


//
//  Rational
//



void Rational::print(std::ostream & target) const
{
	// real-valued rationals print bare (parseable as a rational literal); the
	// complex pair form is reserved for genuinely complex values
	if (true_value_imag_ == 0)
		target << true_value_real_;
	else
		target << "(" << true_value_real_ << "," << true_value_imag_ << ")";
}


	} // re: namespace node
} // re: namespace bertini
