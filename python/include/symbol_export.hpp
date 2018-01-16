//This file is part of Bertini 2.
//
//python/symbol_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/symbol_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/symbol_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//
//  python/symbol_export.hpp:  Header file for exposing symbol nodes to python.


#pragma once

#ifndef BERTINI_PYTHON_SYMBOLS_EXPORT_HPP
#define BERTINI_PYTHON_SYMBOLS_EXPORT_HPP
#include <bertini2/function_tree/symbols/symbol.hpp>
#include <bertini2/function_tree/symbols/number.hpp>
#include <bertini2/function_tree/symbols/special_number.hpp>
#include <bertini2/function_tree/symbols/variable.hpp>
#include <bertini2/function_tree/symbols/differential.hpp>



#include "python_common.hpp"


namespace bertini{
	namespace python{
		
		using namespace bertini::node;


		
		// the main function for exporting symbols into python.
		void ExportSymbols();

		

		
		/**
		 NamedSymbol class(abstract)
		 */
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

		

		/**
		 Integer class 
		 */
		template<typename NodeBaseT>
		class IntegerVisitor: public def_visitor<IntegerVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};

		
		
		
		/** 
		 Rational class
		 */
		template<typename NodeBaseT>
		class RationalVisitor: public def_visitor<RationalVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};


		
		
		/**
		 Variable class 
		 */
		template<typename NodeBaseT>
		class VariableVisitor: public def_visitor<VariableVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};

		
		
		/**
		 Differential class
		*/
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
