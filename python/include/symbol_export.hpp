//
//  symbol_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 1/30/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_symbol_visitors_hpp
#define Xcode_b2_symbol_visitors_hpp
#include <bertini2/function_tree/symbols/symbol.hpp>
#include <bertini2/function_tree/symbols/number.hpp>
#include <bertini2/function_tree/symbols/special_number.hpp>
#include <bertini2/function_tree/symbols/variable.hpp>
#include <bertini2/function_tree/symbols/differential.hpp>



#include "python_common.hpp"


namespace bertini{
	namespace python{
		
		using namespace bertini::node;
		using mpz_int = boost::multiprecision::mpz_int;
		using mpq_rational = boost::multiprecision::mpq_rational;

		void ExportSymbols();

		

		
		///////// NamedSymbol class(abstract) ////////////////
		template<typename NodeBaseT>
		class NamedSymbolVisitor: public def_visitor<NamedSymbolVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:
			std::string (NodeBaseT::*getname)() const = &NodeBaseT::name;
			void (NodeBaseT::*setname)(const std::string &) = &NodeBaseT::name;

		};

		

		///////// Integer class ////////////////
		template<typename NodeBaseT>
		class IntegerVisitor: public def_visitor<IntegerVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};

		
		
		
		///////// Rational class ////////////////
		template<typename NodeBaseT>
		class RationalVisitor: public def_visitor<RationalVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};


		
		
		///////// Variable class ////////////////
		template<typename NodeBaseT>
		class VariableVisitor: public def_visitor<VariableVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};

		
		
//		///////// Differential class ////////////////
		template<typename NodeBaseT>
		class DifferentialVisitor: public def_visitor<DifferentialVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};
		
	}//namespace python
} // namespace bertini

#endif
