//This file is part of Bertini 2.
//
//python/root_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/root_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/root_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/root_export.cpp:  Source file for exposing root nodes to python.




#include <stdio.h>
#include "root_export.hpp"


namespace bertini{
	namespace python{
		
		
		
		template<typename NodeBaseT>
		template<class PyClass>
		void FunctionVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("root", &Function::entry_node,return_value_policy<reference_existing_object>())
			.def("root", &Function::SetRoot)
			.def("ensure_not_empy", &Function::EnsureNotEmpty)
			;
		}

		
		template<typename NodeBaseT>
		template<class PyClass>
		void JacobianVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("evalJd", &Jacobian::template EvalJ<dbl>)
			.def("evalJmp", &Jacobian::template EvalJ<mpfr>)
			;
		}

		
		
		void ExportRoots()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".root");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("root") = new_submodule;

			scope new_submodule_scope = new_submodule;

			// Function class
			class_<Function, bases<NamedSymbol>, std::shared_ptr<Function> >("Function", init<std::string>())
			.def(init<const Nodeptr&>() )
			
			.def(FunctionVisitor<Function>())
			;

			
			// Jacobian class
			class_<Jacobian, bases<Function>, std::shared_ptr<Jacobian> >("Jacobian", init<const Nodeptr&>())
			
			.def(JacobianVisitor<Jacobian>())
			;
			
		}

	}
}

