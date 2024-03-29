//This file is part of Bertini 2.
//
//python/symbol_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/symbol_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/symbol_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//
//  python/symbol_export.cpp:  Source file for exposing symbol nodes to python.




#include <stdio.h>

#include "symbol_export.hpp"



namespace bertini{
	namespace python{
		
		template<typename NodeBaseT>
		template<class PyClass>
		void NamedSymbolVisitor<NodeBaseT>::visit(PyClass& cl) const
		{			
			cl
			.add_property("name", getname, setname)
			;
		}

		
		template<typename NodeBaseT>
		template<class PyClass>
		void RationalVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("rand", 		&NodeBaseT::template Rand<16>)
			.def("rand_real", 	&NodeBaseT::template RandReal<16>)
			.staticmethod("rand")
			.staticmethod("rand_real")
			;
		}

		
		
		
		template<typename NodeBaseT>
		template<class PyClass>
		void VariableVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("set_current_value", &NodeBaseT::template set_current_value<dbl>)
			.def("set_current_value", &NodeBaseT::template set_current_value<mpfr_complex>)
			;
		}

		
		template<typename NodeBaseT>
		template<class PyClass>
		void DifferentialVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("get_variable", &NodeBaseT::GetVariable,return_value_policy<reference_existing_object>())
			;
		}

		
		
		
		void ExportSymbols()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".symbol");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("symbol") = new_submodule;

			scope new_submodule_scope = new_submodule;


			// Symbol class
			class_<Symbol, boost::noncopyable, bases<Node>, std::shared_ptr<Symbol> >("AbstractSymbol", no_init)
			;
			
			// NamedSymbol class
			class_<NamedSymbol, boost::noncopyable, bases<Symbol>, std::shared_ptr<NamedSymbol> >("AbstractNamedSymbol", no_init)
			;
			
			// Number class
			class_<Number, boost::noncopyable, bases<Symbol>, std::shared_ptr<Number> >("AbstractNumber", no_init)
			;
			
			// Float class
			class_<Float, bases<Number>, std::shared_ptr<Float> >("Float", no_init)
			.def("__init__", make_constructor(&Float::template Make<mpfr_float const&, mpfr_float const&>))
			.def("__init__", make_constructor(&Float::template Make<std::string const&>))
			.def("__init__", make_constructor(&Float::template Make<std::string const&, std::string const&>))
			.def("__init__", make_constructor(&Float::template Make<mpfr_complex const&>))
			;
			
			
			// Pi class
			class_<special_number::Pi, bases<NamedSymbol>, std::shared_ptr<special_number::Pi> >("Pi", no_init)
			.def("__init__", make_constructor(&special_number::Pi::template Make<>))
			;


			// E class
			class_<special_number::E, bases<NamedSymbol>, std::shared_ptr<special_number::E> >("E", no_init)
			.def("__init__", make_constructor(&special_number::E::template Make<>))
			;

			
			// Integer class
			class_<Integer, bases<Number>, std::shared_ptr<Integer> >("Integer",no_init)
			.def("__init__", make_constructor(&Integer::template Make<int const&>))
			.def("__init__", make_constructor(&Integer::template Make<mpz_int const&>))
			.def("__init__", make_constructor(&Integer::template Make<std::string const&>))
			;

			
			// Rational class
			class_<Rational, bases<Number>, std::shared_ptr<Rational> >("Rational", no_init)
			.def("__init__", make_constructor(&Rational::template Make<int const&>))
			.def("__init__", make_constructor(&Rational::template Make<int const&, int const&, int const&, int const&>))
			.def("__init__", make_constructor(&Rational::template Make<std::string const&>))
			.def("__init__", make_constructor(&Rational::template Make<std::string const&, std::string const&>))
			.def("__init__", make_constructor(&Rational::template Make<mpq_rational const&, mpq_rational const&>))
			.def(RationalVisitor<Rational>())
			;

			
			// Variable class
			class_<Variable, bases<NamedSymbol>, std::shared_ptr<Variable> >("Variable", no_init)
			.def("__init__", make_constructor(&Variable::template Make< std::string const& >))
			.def(VariableVisitor<Variable>())
			;
			
			
			// Differential class
			class_<Differential, bases<NamedSymbol>,std::shared_ptr<node::Differential> >("Differential", no_init)
			.def("__init__", make_constructor(Differential::template Make<std::shared_ptr<Variable> const&,std::string const&>))
			.def(DifferentialVisitor<Differential>())
			;
			
			
			def("make_pi", &Pi);
			def("make_i", &I);
			def("make_e", &E);
			
		};
		
		
	} //namespace python
} // namespace bertini
