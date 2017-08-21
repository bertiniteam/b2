//This file is part of Bertini 2.
//
//src/function_tree/symbols/symbol.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/symbols/symbol.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/symbols/symbol.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame
// Jeb Collins, West Texas A&M


#include "function_tree/symbols/symbol.hpp"




namespace bertini{
	namespace node{

unsigned Symbol::EliminateZeros()
{
	return 0;
}
unsigned Symbol::EliminateOnes()
{
	return 0;
}


const std::string& NamedSymbol::name() const
{
	return name_;
};

void NamedSymbol::name(const std::string & new_name)
{name_ = new_name;}


NamedSymbol::NamedSymbol(const std::string & new_name) : name_(new_name)
{}


void NamedSymbol::print(std::ostream& target) const
{
	target << name();
}

	} // re: namespace node
} // re: namespace bertini
