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
// Copyright(C) Bertini2 Development Team
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
			.add_property("name",
				make_function(getname, return_value_policy<copy_const_reference>()),
				setname)
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
		void VariableVisitor<NodeBaseT>::visit(PyClass& /*cl*/) const
		{
			// Variables no longer carry a value: evaluation is through the SLP (System.eval /
			// Node.eval), so there is no longer a set_current_value to bind.
		}

		
		template<typename NodeBaseT>
		template<class PyClass>
		void DifferentialVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			// shared_ptr<Variable const> has no registered python class; cast away const
			// and return by value so the conversion (and downcast) machinery applies
			.def("get_variable", +[](NodeBaseT const& d) { return std::const_pointer_cast<Variable>(d.GetVariable()); }, (arg("self")), "the variable this differential is with respect to")
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
			.def(NamedSymbolVisitor<NamedSymbol>())
			;
			
			// Number class
			class_<Number, boost::noncopyable, bases<Symbol>, std::shared_ptr<Number> >("AbstractNumber", no_init)
			;
			
			// Complex class -- a complex-number literal node (NOT a real "float"): it stores a full
			// arbitrary-precision complex value; a real literal is the zero-imaginary special case.
			// Construct from one argument (the real part, or a multiprec value) or two (real, imag).
			// Prefer Integer/Rational for exact coefficients -- they are faster and need no stored sample.
			class_<Complex, bases<Number>, std::shared_ptr<Complex> >("Complex",
				"A complex-number literal node in an expression tree.  Holds a full arbitrary-precision "
				"complex value (a real literal is just zero imaginary part).  Build it as Complex(real) "
				"or Complex(real, imag) from exact strings/multiprec values; prefer Integer or Rational "
				"for exact non-irrational coefficients.", no_init)
			.def("__init__", make_constructor(&Complex::template Make<real_mp const&, real_mp const&>))
			.def("__init__", make_constructor(&Complex::template Make<std::string const&>))
			.def("__init__", make_constructor(&Complex::template Make<std::string const&, std::string const&>))
			.def("__init__", make_constructor(&Complex::template Make<complex_mp const&>))
			.def("value", &Complex::GetValue, return_value_policy<copy_const_reference>(), "the literal value this node represents, at its stored (highest) precision")
			;
			
			
			// Pi class
			// no longer a NamedSymbol in C++ (it's a Number + the Named capability),
			// so .name is bound directly rather than inherited from AbstractNamedSymbol
			class_<special_number::Pi, bases<Number>, std::shared_ptr<special_number::Pi> >("Pi", no_init)
			.def("__init__", make_constructor(&special_number::Pi::template Make<>))
			.def(NamedSymbolVisitor<special_number::Pi>())
			;


			// E class
			class_<special_number::E, bases<Number>, std::shared_ptr<special_number::E> >("E", no_init)
			.def("__init__", make_constructor(&special_number::E::template Make<>))
			.def(NamedSymbolVisitor<special_number::E>())
			;

			
			// Integer class
			class_<Integer, bases<Number>, std::shared_ptr<Integer> >("Integer",no_init)
			.def("__init__", make_constructor(&Integer::template Make<int const&>))
			.def("__init__", make_constructor(&Integer::template Make<mpz_int const&>))
			.def("__init__", make_constructor(&Integer::template Make<std::string const&>))
			.def("value", &Integer::GetValue, return_value_policy<copy_const_reference>(), "the literal value this node represents")
			;

			
			// Rational class
			class_<Rational, bases<Number>, std::shared_ptr<Rational> >("Rational", no_init)
			.def("__init__", make_constructor(&Rational::template Make<int const&>))
			.def("__init__", make_constructor(&Rational::template Make<int const&, int const&, int const&, int const&>))
			.def("__init__", make_constructor(&Rational::template Make<std::string const&>))
			.def("__init__", make_constructor(&Rational::template Make<std::string const&, std::string const&>))
			.def("__init__", make_constructor(&Rational::template Make<mpq_rational const&, mpq_rational const&>))
			.def(RationalVisitor<Rational>())
			.def("value_real", &Rational::GetValueReal, return_value_policy<copy_const_reference>(), "the real part of the literal value this node represents")
			.def("value_imag", &Rational::GetValueImag, return_value_policy<copy_const_reference>(), "the imaginary part of the literal value this node represents")
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

		};
		
		
	} //namespace python
} // namespace bertini
