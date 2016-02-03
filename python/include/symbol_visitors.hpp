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


#include "python_common.hpp"
#include "node_visitors.hpp"

namespace bertini{
	namespace python{
		
		using namespace bertini::node;

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

		
		


//		///////// Pi class ////////////////
//		template<typename NodeBaseT>
//		class PiVisitor: public def_visitor<PiVisitor<NodeBaseT> >
//		{
//			friend class def_visitor_access;
//			
//		public:
//			template<class PyClass>
//			void visit(PyClass& cl) const
//			{
//				NamedSymbolBaseVisitor<NodeBaseT>().visit(cl);
//			}
//		};
//
//		
//		
//		///////// E class ////////////////
//		template<typename NodeBaseT>
//		class EVisitor: public def_visitor<EVisitor<NodeBaseT> >
//		{
//			friend class def_visitor_access;
//			
//		public:
//			template<class PyClass>
//			void visit(PyClass& cl) const
//			{
//				NamedSymbolBaseVisitor<NodeBaseT>().visit(cl);
//			}
//		};

		
		
		///////// Variable class ////////////////
//		template<typename NodeBaseT>
//		class VariableVisitor: public def_visitor<VariableVisitor<NodeBaseT> >
//		{
//			friend class def_visitor_access;
//			
//		public:
//			template<class PyClass>
//			void visit(PyClass& cl) const
//			{
//				NodeBaseVisitor<NodeBaseT>().visit(cl);
//			}
//		};

		
		
//		///////// Differential class ////////////////
//		template<typename NodeBaseT>
//		class DifferentialVisitor: public def_visitor<DifferentialVisitor<NodeBaseT> >
//		{
//			friend class def_visitor_access;
//			
//		public:
//			template<class PyClass>
//			void visit(PyClass& cl) const
//			{
//				NamedSymbolBaseVisitor<NodeBaseT>().visit(cl);
//				
//				cl
//				.def("get_variable", &NodeBaseT::GetVariable)
//				;
//			}
//		};
		
	}//namespace python
} // namespace bertini

#endif
