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
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame
// Jeb Collins, West Texas A&M


#include "function_tree/symbols/special_number.hpp"




namespace bertini{
	namespace node{
		namespace special_number{
using ::pow;

// Return value of constant
dbl Pi::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return acos(-1.0);
}

void Pi::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = acos(-1.0);
}


mpfr_complex Pi::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return mpfr_complex(mpfr_float(acos(mpfr_float(-1))));
}

void Pi::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = mpfr_complex(mpfr_float(acos(mpfr_float(-1))));
}




// Return value of constant
dbl E::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const 
{
	return dbl(exp(1.0),0.0);
}

void E::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const 
{
	evaluation_value = dbl(exp(1.0),0.0);
}


mpfr_complex E::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const 
{
	return mpfr_complex(mpfr_float(exp(mpfr_float(1))));
}

void E::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const 
{
	evaluation_value = mpfr_complex(mpfr_float(exp(mpfr_float(1))));
}
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
	return Float::Make(0,1);
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
