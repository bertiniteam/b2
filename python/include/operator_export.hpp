//This file is part of Bertini 2.
//
//python/operator_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/operator_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/operator_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  Danielle Brake
//  UWEC
//  Spring 2018
//
//
//  python/operator_export.hpp:  Header file for exposing operator nodes to python.



#pragma once
#ifndef BERTINI_PYTHON_OPERATOR_EXPORT_HPP
#define BERTINI_PYTHON_OPERATOR_EXPORT_HPP

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
		
		
		
		
		
		/**
		 UnaryOperator class(abstract)
		 */
		template<typename NodeBaseT>
		class UnaryOpVisitor: public def_visitor<UnaryOpVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};

	
		
		
		/** NaryOperator class(abstract)
		 */
		template<typename NodeBaseT>
		class NaryOpVisitor: public def_visitor<NaryOpVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};

		
		
		
		
		/**
		 SumOperator and MultOperator classes
		 */
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

		
		
		
		/**
		 PowerOperator class
		 */
		template<typename NodeBaseT>
		class PowerOpVisitor: public def_visitor<PowerOpVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};

		
		/**
		 IntegerPowerOperator class 
		 */
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
