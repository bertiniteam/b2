// python/function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/function_tree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/function_tree.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Jeb Collins
//  West Texas A&M University
//  Mathematics
//  Fall 2015
//
//
//  python/symbol.cpp:  the source file for the python interface for symbol classes.

#include <stdio.h>

#include "symbol.hpp"
#include "boost_mp.hpp"


using namespace boost::python;


namespace bertini{
	namespace python{

		
		
		void ExportSymbols()
		{
			
			// Symbol class(abstract)
			class_<bertini::node::Symbol, bases<bertini::node::Node>,std::shared_ptr<node::Symbol>, boost::noncopyable>("Symbol", no_init);
			
			
			// NamedSymbol class(abstract)
			class_<bertini::node::NamedSymbol, bases<bertini::node::Symbol>,std::shared_ptr<node::NamedSymbol>, boost::noncopyable >("NamedSymbol", no_init)
			.def("print", &bertini::node::NamedSymbol::print)
			;
			
			
			// Number class(abstract)
			class_<node::Number, bases<bertini::node::Symbol>, boost::noncopyable >("Number", no_init)
			.def("degree", numberDeg1, numberDeg1Overloads())
			.def("degree", numberDeg2)
			.def("multidegree", &node::Number::MultiDegree)
			.def("is_homogeneous", numberIsHom1, numberIsHom1Overloads())
			.def("is_homogeneous", numberIsHom2)
			.def("precision", &node::Number::precision)
			.def("homogenize", &node::Number::Homogenize)
			.def("reset", &node::Number::Reset)
			;
			
			
			// Float class
			class_<bertini::node::Float, bases<bertini::node::Number>, std::shared_ptr<node::Float> >("Float", init< optional<double, double> >())
			.def(init<dbl>())
			.def(init<mpfr>())
			.def(init< bmp, bmp >())
			.def(init< optional<std::string, std::string> >())
			
			.def("differentiate", &node::Float::Differentiate)
			;
			
			
			// Pi class
			class_<node::special_number::Pi, bases<node::NamedSymbol> >("Pi", init<>())
			.def("degree", piDeg1, piDeg1Overloads())
			.def("degree", piDeg2)
			.def("multidegree", &node::special_number::Pi::MultiDegree)
			.def("is_homogeneous", piIsHom1, piIsHom1Overloads())
			.def("is_homogeneous", piIsHom2)
			.def("precision", &node::special_number::Pi::precision)
			.def("homogenize", &node::special_number::Pi::Homogenize)
			.def("reset", &node::special_number::Pi::Reset)
			.def("differentiate", &node::special_number::Pi::Differentiate)
			;
			
			
			// E class
			class_<node::special_number::E, bases<node::NamedSymbol> >("E", init<>())
			.def("degree", eDeg1, eDeg1Overloads())
			.def("degree", eDeg2)
			.def("multidegree", &node::special_number::E::MultiDegree)
			.def("is_homogeneous", eIsHom1, eIsHom1Overloads())
			.def("is_homogeneous", eIsHom2)
			.def("precision", &node::special_number::E::precision)
			.def("homogenize", &node::special_number::E::Homogenize)
			.def("reset", &node::special_number::E::Reset)
			.def("differentiate", &node::special_number::E::Differentiate)
			;
			
			
			// Variable class
			class_<node::Variable, bases<node::NamedSymbol>, std::shared_ptr<node::Variable> >("Variable", init< optional <std::string> >())
			.def("degree", varDeg1, varDeg1Overloads())
			.def("degree", varDeg2)
			.def("multidegree", &node::Variable::MultiDegree)
			.def("is_homogeneous", varIsHom1, varIsHom1Overloads())
			.def("is_homogeneous", varIsHom2)
			.def("precision", &node::Variable::precision)
			.def("homogenize", &node::Variable::Homogenize)
			.def("reset", &node::Variable::Reset)
			.def("differentiate", &node::Variable::Differentiate)
			.def("set_current_value", &node::Variable::set_current_value<dbl>)
			;
			
			
			// Differential class
			class_<node::Differential, bases<node::NamedSymbol>,std::shared_ptr<node::Differential> >
			("Differential", init<>())
			.def(init<std::shared_ptr<node::Variable>,std::string>())
			.def("get_variable", &node::Differential::GetVariable)
			.def("degree", diffDeg1, diffDeg1Overloads())
			.def("degree", diffDeg2)
			.def("multidegree", &node::Differential::MultiDegree)
			.def("is_homogeneous", diffIsHom1, diffIsHom1Overloads())
			.def("is_homogeneous", diffIsHom2)
			.def("precision", &node::Differential::precision)
			.def("homogenize", &node::Differential::Homogenize)
			.def("differentiate", &node::Differential::Differentiate)
			;
			
		}
		
	} // re: python
} // re: bertini
