//
//  operator_visitors.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 1/30/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_operator_visitors_hpp
#define Xcode_b2_operator_visitors_hpp

#include <bertini2/function_tree/operators/operator.hpp>
#include <bertini2/function_tree/operators/arithmetic.hpp>
#include <bertini2/function_tree/operators/trig.hpp>

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		using namespace bertini::node;

		using Node = Node;
		using Nodeptr = std::shared_ptr<Node>;
		
		void ExportOperators();
		
		template<typename NodeBaseT>
		class UnaryOpVisitor: public def_visitor<UnaryOpVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};

	
		
		template<typename NodeBaseT>
		class NaryOpVisitor: public def_visitor<NaryOpVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};

		
		
		
		
		///////// SumOperator & MultOperator class ////////////////
		template<typename NodeBaseT>
		class SumMultOpVisitor: public def_visitor<SumMultOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
			
			
		private:
			void (NodeBaseT::*addChild2)(std::shared_ptr<Node> child, bool) = &NodeBaseT::AddChild;

		};

		
		///////// PowerOperator class ////////////////
		template<typename NodeBaseT>
		class PowerOpVisitor: public def_visitor<PowerOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};

		
		///////// IntegerPowerOperator class ////////////////
		template<typename NodeBaseT>
		class IntPowOpVisitor: public def_visitor<IntPowOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:
			std::string (NodeBaseT::*getexp)() const = &NodeBaseT::exponent;
			void (NodeBaseT::*setexp)(const std::string &) = &NodeBaseT::set_exponent;
			
		};

		
		
	} //re: namespace python
}//re: namespace bertini

		
		
		
		
		
		
		
		
		
		
#endif
