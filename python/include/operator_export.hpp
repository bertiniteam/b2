//
//  operator_visitors.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 1/30/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_operator_visitors_hpp
#define Xcode_b2_operator_visitors_hpp

#include <bertini2/operator.hpp>
#include <bertini2/arithmetic.hpp>

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		using namespace bertini::node;
		
		
		template<typename NodeBaseT>
		class UnaryOpBaseVisitor: public def_visitor<UnaryOpBaseVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				NodeBaseVisitor<NodeBaseT>().visit(cl);
				
				cl
				.def("set_child", &NodeBaseT::SetChild)
				.def("first_child", &NodeBaseT::first_child)
				;
			}
		}

	
		template<typename NodeBaseT>
		class BinaryOpBaseVisitor: public def_visitor<BinaryOpBaseVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				NodeBaseVisitor<NodeBaseT>().visit(cl);
			}
		}

		
		
		template<typename NodeBaseT>
		class NaryOpBaseVisitor: public def_visitor<NaryOpBaseVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				NodeBaseVisitor<NodeBaseT>().visit(cl);
				
				cl
				.def("add_child", addChild1)
				.def("first_child", &NodeBaseT::first_child)
				.def("children_size", &NodeBaseT::children_size)
				;
			}
			
			
		private:
			void (NodeBaseT::*addChild1)(std::shared_ptr<Node> child) = &NodeBaseT::AddChild;

		}

		
		
		
		
		///////// SumOperator class ////////////////
		template<typename NodeBaseT>
		class SumOpVisitor: public def_visitor<SumOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				NaryOpBaseVisitor<NodeBaseT>().visit(cl);
				
				cl
				.def("add_child", addChild2)
				;
			}
			
			
			
			
		private:
			void (NodeBaseT::*addChild2)(std::shared_ptr<Node> child, bool) = &NodeBaseT::AddChild;

		}

		
		
		
		
		///////// NegateOperator class ////////////////
		template<typename NodeBaseT>
		class NegateOpVisitor: public def_visitor<NegateOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				UnaryOpBaseVisitor<NodeBaseT>().visit(cl);
			}
		}

		
		
		
		
		///////// MultOperator class ////////////////
		template<typename NodeBaseT>
		class MultOpVisitor: public def_visitor<MultOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				NaryOpBaseVisitor<NodeBaseT>().visit(cl);
				
				cl
				.def("add_child", addChild2)
				;
			}
			
			
			
			
		private:
			void (NodeBaseT::*addChild2)(std::shared_ptr<Node> child, bool) = &NodeBaseT::AddChild;
			
		}

		

		///////// PowerOperator class ////////////////
		template<typename NodeBaseT>
		class PowerOpVisitor: public def_visitor<PowerOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				BinaryOpBaseVisitor<NodeBaseT>().visit(cl);
				
				cl
				.def("set_base", &NodeBaseT::SetBase)
				.def("set_exponent", &NodeBaseT::SetExponent)
				;
			}
			
		}

	
		
		
		///////// IntegerPowerOperator class ////////////////
		template<typename NodeBaseT>
		class IntPowerOpVisitor: public def_visitor<IntPowerOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				UnaryOpBaseVisitor<NodeBaseT>().visit(cl);
				
				cl
				.add_property("exponent", &NodeBaseT::exponent, &NodeBaseT::set_exponent)
				;
			}
			
		}

	
		
		///////// SqrtOperator class ////////////////
		template<typename NodeBaseT>
		class SqrtOpVisitor: public def_visitor<SqrtOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				UnaryOpBaseVisitor<NodeBaseT>().visit(cl);
			}
			
		}

		
		
		///////// ExpOperator class ////////////////
		template<typename NodeBaseT>
		class ExpOpVisitor: public def_visitor<ExpOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				UnaryOpBaseVisitor<NodeBaseT>().visit(cl);
			}
			
		}

		
		
		///////// LogOperator class ////////////////
		template<typename NodeBaseT>
		class LogOpVisitor: public def_visitor<LogOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				UnaryOpBaseVisitor<NodeBaseT>().visit(cl);
			}
			
		}

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
#endif
