//
//  symbol_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 1/30/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_symbol_export_hpp
#define Xcode_b2_symbol_export_hpp

#include "symbol_visitors.hpp"



namespace bertini{
	namespace python{
		
		void ExportSymbols()
		{
			using bertini::node;
			
			// Symbol class
			class_<Symbol, boost::noncopyable, bases<Node>, std::shared_Ptr<Symbol> >("Symbol", no_init)
			.def(SymbolVisitor<Symbol>())
			;

			// NamedSymbol class
			class_<NamedSymbol, boost::noncopyable, bases<Symbol>, std::shared_Ptr<Symbol> >("NamedSymbol", no_init)
			.def(NamedSymbolVisitor<NamedSymbol>())
			;

			// Number class
			class_<Number, boost::noncopyable, bases<NamedSymbol>, std::shared_Ptr<Number> >("Number", no_init)
			.def(NumberVisitor<Number>())
			;

			// Float class
			class_<Float, bases<Number>, std::shared_ptr<Float> >("Float", init< optional<double, double> >())
			.def(init<dbl>())
			.def(init<mpfr>())
			.def(init< bmp, bmp >())
			.def(init< optional<std::string, std::string> >())
			
			.def(FloatVisitor<Float>())
			;

			
			// Pi class
			class_<special_number::Pi, bases<NamedSymbol>, std::shared_ptr<special_number::Pi> >("Pi", init<>())
			.def(PiVisitor<special_number::Pi>())
			;
			
			
			// E class
			class_<special_number::E, bases<NamedSymbol>, std::shared_ptr<special_number::E> >("E", init<>())
			.def(EVisitor<special_number::E>())
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

		};
		
		
	} //namespace python
} // namespace bertini

#endif
