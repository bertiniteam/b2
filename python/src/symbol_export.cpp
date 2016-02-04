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
			
			
			//			// Pi class
			//			class_<special_number::Pi, bases<NamedSymbol>, std::shared_ptr<special_number::Pi> >("Pi", init<>())
			//			.def(PiVisitor<special_number::Pi>())
			//			;
			//
			//
			//			// E class
			//			class_<special_number::E, bases<NamedSymbol>, std::shared_ptr<special_number::E> >("E", init<>())
			//			.def(EVisitor<special_number::E>())
			//			;
			//
			
			// Variable class
			class_<Variable, bases<NamedSymbol>, std::shared_ptr<Variable> >("Variable", init< optional <std::string> >())
			;
			
			
			//			// Differential class
			//			class_<Differential, bases<NamedSymbol>,std::shared_ptr<node::Differential> >("Differential", init<>())
			//			.def(init<std::shared_ptr<Variable>,std::string>())
			//
			//			.def(DifferentialVisitor<Differential>())
			//			;
			
		};
		
		
	} //namespace python
} // namespace bertini
