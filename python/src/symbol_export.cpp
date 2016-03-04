//
//  symbol_export.cpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/2/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

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
			.def("rand", &NodeBaseT::Rand)
			.def("rand_real", &NodeBaseT::RandReal)
			;
		}

		
		
		
		template<typename NodeBaseT>
		template<class PyClass>
		void VariableVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("set_current_value", &NodeBaseT::template set_current_value<dbl>)
			.def("set_current_value", &NodeBaseT::template set_current_value<mpfr>)
			;
		}

		
		template<typename NodeBaseT>
		template<class PyClass>
		void DifferentialVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("get_variable", &NodeBaseT::GetVariable)
			;
		}

		
		
		
		void ExportSymbols()
		{
			
			// Symbol class
			class_<Symbol, boost::noncopyable, bases<Node>, std::shared_ptr<Symbol> >("Symbol", no_init)
			;
			
			// NamedSymbol class
			class_<NamedSymbol, boost::noncopyable, bases<Symbol>, std::shared_ptr<NamedSymbol> >("NamedSymbol", no_init)
			;
			
			// Number class
			class_<Number, boost::noncopyable, bases<Symbol>, std::shared_ptr<Number> >("Number", no_init)
			;
			
			// Float class
			class_<Float, bases<Number>, std::shared_ptr<Float> >("Float", init< double, double >())
			.def(init<dbl>())
			.def(init<mpfr>())
			.def(init< bmp, bmp >())
			.def(init< std::string, std::string >())
			;
			
			
			// Pi class
			class_<special_number::Pi, bases<NamedSymbol>, std::shared_ptr<special_number::Pi> >("Pi", init<>())
			;


			// E class
			class_<special_number::E, bases<NamedSymbol>, std::shared_ptr<special_number::E> >("E", init<>())
			;

			
			// Integer class
			class_<Integer, bases<Number>, std::shared_ptr<Integer> >("Integer", init< int >())
			.def(init<mpz_int>())
			.def(init<std::string const&>())
			;

			
			// Rational class
			class_<Rational, bases<Number>, std::shared_ptr<Rational> >("Rational", init< int >())
			.def(init<int, int>())
			.def(init<int, int, int, int>())
			.def(init<std::string>())
			.def(init<std::string, std::string>())
			.def(init<mpq_rational const&, mpq_rational const&>())
			
			.def(RationalVisitor<Rational>())
			;

			
			// Variable class
			class_<Variable, bases<NamedSymbol>, std::shared_ptr<Variable> >("Variable", init< std::string >())
			.def(VariableVisitor<Variable>())
			;
			
			
			// Differential class
			class_<Differential, bases<NamedSymbol>,std::shared_ptr<node::Differential> >("Differential", init<std::shared_ptr<Variable>,std::string>())
			
			.def(DifferentialVisitor<Differential>())
			;
			
			
			def("make_pi", &Pi);
			def("make_i", &I);
			def("make_e", &E);
			
		};
		
		
	} //namespace python
} // namespace bertini
