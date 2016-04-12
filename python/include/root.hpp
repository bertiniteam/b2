//This file is part of Bertini 2.
//
//python/root.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/root.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/root.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Jeb Collins
//  West Texas A&M University
//  Mathematics
//  Fall 2015
//
//
//  python/root.hpp:  the header file for the python interface for root classes.

#ifndef BERTINI_PYTHON_ROOT_HPP
#define BERTINI_PYTHON_ROOT_HPP

#include "function_tree.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		
		// Function class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(funcDeg1Overloads, node::Function::Degree, 0, 1)
		int (node::Function::*funcDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::Function::Degree;
		int (node::Function::*funcDeg2)(VariableGroup const&) const  = &bertini::node::Function::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(funcIsHom1Overloads, node::Function::IsHomogeneous, 0, 1)
		bool (node::Function::*funcIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::Function::IsHomogeneous;
		bool (node::Function::*funcIsHom2)(VariableGroup const& vars) const= &bertini::node::Function::IsHomogeneous;

		
		
	} // re: python
} // re: bertini

#endif
