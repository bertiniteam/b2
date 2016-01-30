//This file is part of Bertini 2.0.
//
// python/function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/function_tree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/function_tree.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
//
//  python/function_tree.hpp:  the header file for the python interface for function_tree class.

#ifndef BERTINI_PYTHON_FUNCTION_TREE_HPP
#define BERTINI_PYTHON_FUNCTION_TREE_HPP

#include <bertini2/mpfr_complex.hpp>
#include <bertini2/function_tree/node.hpp>
#include <bertini2/function_tree.hpp>
#include <bertini2/system.hpp>


#include <boost/python.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <boost/python/wrapper.hpp>

#include "export_mpfr.hpp"


namespace bertini{
	namespace python{

		
		using namespace boost::python;
		using Node = node::Node;
		using Nd = std::shared_ptr<node::Node>;
		
		void SetupFunctionTree();
		void ExportNode();
		void ExportSymbols();
		void ExportOperators();
		void ExportRoots();
		void ExportSystem();
		
		
	}
}


#endif
