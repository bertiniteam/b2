//This file is part of Bertini 2.
//
//python/root_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/root_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/root_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/root_export.hpp:  Header file for exposing root nodes to python.


#pragma once

#ifndef BERTINI_PYTHON_ROOT_EXPORT_HPP
#define BERTINI_PYTHON_ROOT_EXPORT_HPP
#include <bertini2/function_tree/roots/function.hpp>
#include <bertini2/function_tree/roots/jacobian.hpp>

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		using namespace bertini::node;
		
		using Node = Node;
		using Nodeptr = std::shared_ptr<Node>;
		
		void ExportRoots();
		
		
		/**
		 Function class
		 */
		template<typename NodeBaseT>
		class FunctionVisitor: public def_visitor<FunctionVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};

		
		
		/** 
		 Jacobian class
		 */
		template<typename NodeBaseT>
		class JacobianVisitor: public def_visitor<JacobianVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};

		
	}
}

#endif
