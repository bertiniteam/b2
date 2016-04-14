//This file is part of Bertini 2.
//
//python/symbol.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/symbol.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/symbol.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/symbol.hpp:  the header file for the python interface for symbol classes.

#ifndef BERTINI_PYTHON_SYMBOL_HPP
#define BERTINI_PYTHON_SYMBOL_HPP

#include "function_tree.hpp"


namespace bertini{
	namespace python{
		
		using namespace boost::python;
		// Number class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(numberDeg1Overloads, node::Number::Degree, 0, 1)
		int (bertini::node::Number::*numberDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::Number::Degree;
		int (bertini::node::Number::*numberDeg2)(VariableGroup const&) const  = &bertini::node::Number::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(numberIsHom1Overloads, node::Number::IsHomogeneous, 0, 1)
		bool (bertini::node::Number::*numberIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::Number::IsHomogeneous;
		bool (bertini::node::Number::*numberIsHom2)(VariableGroup const& vars) const= &bertini::node::Number::IsHomogeneous;
		
		// Pi class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(piDeg1Overloads, node::special_number::Pi::Degree, 0, 1)
		int (node::special_number::Pi::*piDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::special_number::Pi::Degree;
		int (node::special_number::Pi::*piDeg2)(VariableGroup const&) const  = &bertini::node::special_number::Pi::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(piIsHom1Overloads, node::special_number::Pi::IsHomogeneous, 0, 1)
		bool (node::special_number::Pi::*piIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::special_number::Pi::IsHomogeneous;
		bool (node::special_number::Pi::*piIsHom2)(VariableGroup const& vars) const= &bertini::node::special_number::Pi::IsHomogeneous;
		
		// E class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(eDeg1Overloads, node::special_number::E::Degree, 0, 1)
		int (node::special_number::E::*eDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::special_number::E::Degree;
		int (node::special_number::E::*eDeg2)(VariableGroup const&) const  = &bertini::node::special_number::E::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(eIsHom1Overloads, node::special_number::E::IsHomogeneous, 0, 1)
		bool (node::special_number::E::*eIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::special_number::E::IsHomogeneous;
		bool (node::special_number::E::*eIsHom2)(VariableGroup const& vars) const= &bertini::node::special_number::E::IsHomogeneous;
		
		// Variable class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(varDeg1Overloads, node::Variable::Degree, 0, 1)
		int (node::Variable::*varDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::Variable::Degree;
		int (node::Variable::*varDeg2)(VariableGroup const&) const  = &bertini::node::Variable::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(varIsHom1Overloads, node::Variable::IsHomogeneous, 0, 1)
		bool (node::Variable::*varIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::Variable::IsHomogeneous;
		bool (node::Variable::*varIsHom2)(VariableGroup const& vars) const= &bertini::node::Variable::IsHomogeneous;
		
		// Differential class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(diffDeg1Overloads, node::Differential::Degree, 0, 1)
		int (node::Differential::*diffDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::Differential::Degree;
		int (node::Differential::*diffDeg2)(VariableGroup const&) const  = &bertini::node::Differential::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(diffIsHom1Overloads, node::Differential::IsHomogeneous, 0, 1)
		bool (node::Differential::*diffIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::Differential::IsHomogeneous;
		bool (node::Differential::*diffIsHom2)(VariableGroup const& vars) const= &bertini::node::Differential::IsHomogeneous;
		

		
		
		
	} // re: python
} // re: bertini


#endif
