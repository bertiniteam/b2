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
			class_<Float, bases<Number>, std::shared_ptr<Float> >("Float", init< optional<double, double> >())
			.def(init<dbl>())
			.def(init<mpfr>())
			.def(init< bmp, bmp >())
			.def(init< optional<std::string, std::string> >())
			;
			
			
			// Pi class
			class_<special_number::Pi, bases<NamedSymbol>, std::shared_ptr<special_number::Pi> >("Pi", init<>())
			;


			// E class
			class_<special_number::E, bases<NamedSymbol>, std::shared_ptr<special_number::E> >("E", init<>())
			;

			
			// Variable class
			class_<Variable, bases<NamedSymbol>, std::shared_ptr<Variable> >("Variable", init< optional <std::string> >())
			.def(VariableVisitor<Variable>())
			;
			
			
			// Differential class
			class_<Differential, bases<NamedSymbol>,std::shared_ptr<node::Differential> >("Differential", init<>())
			.def(init<std::shared_ptr<Variable>,std::string>())
			
			.def(DifferentialVisitor<Differential>())
			;
			
			
			def("Pi", &Pi);
			def("I", &I);
			def("E", &E);
			
		};
		
		
	} //namespace python
} // namespace bertini
