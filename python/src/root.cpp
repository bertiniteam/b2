//This file is part of Bertini 2.
//
//python/root.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/root.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/root.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 by Bertini2 Development Team
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
//
//  python/root.cpp:  the source file for the python interface for root classes.

#include <stdio.h>
#include "root.hpp"

using namespace boost::python;
using namespace bertini::node;


namespace bertini{
	namespace python{

		
		
		void ExportRoots()
		{
			using Nodeptr = std::shared_ptr<node::Node>;
			
			
			// Function class
			class_< Function, bases<NamedSymbol>, std::shared_ptr<Function> >("Function", init< optional< std::string > >())
			.def(init< const Nodeptr& >())
			.def("entry_node", &Function::entry_node)
			.def("set_root", &Function::SetRoot)
			.def("reset", &Function::Reset)
			.def("ensure_not_empty", &Function::EnsureNotEmpty)
			.def("differentiate", &Function::Differentiate)
			.def("degree", funcDeg1, funcDeg1Overloads())
			.def("degree", funcDeg2)
			.def("multidegree", &node::Function::MultiDegree)
			.def("is_homogeneous", funcIsHom1, funcIsHom1Overloads())
			.def("is_homogeneous", funcIsHom2)
			.def("homogenize", &node::Function::Homogenize)
			.def("precision", &Function::precision)
			;
			
			
			// Jacobian class
			class_< Jacobian, bases<Function>, std::shared_ptr<Jacobian> >("Jacobian", init< optional< const Nodeptr&> >())
			.def("EvalJ", &Jacobian::EvalJ<dbl>)
			;
			
		}

	} // re: python
} // re: bertini
