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
// Copyright(C) Bertini2 Development Team
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
//  silviana amethyst
//  UWEC
//  Spring 2018, Fall 2023
//
//
//
//  python/root_export.cpp:  Source file for exposing root nodes to python.




#include <stdio.h>
#include "root_export.hpp"


namespace bertini{
	namespace python{
		
		

		void ExportRoots()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".root");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("root") = new_submodule;

			scope new_submodule_scope = new_submodule;


			// NamedExpression: Named(expr, "a") -- a user-named subexpression that prints as its
			// name and evaluates to its expression.  This is the sole surviving handle node
			// (Function and Handle were deleted; it inherits NamedSymbol directly).
			class_<NamedExpression, bases<NamedSymbol>, std::shared_ptr<NamedExpression> >("NamedExpression", no_init)
			.def("__init__",make_constructor(&NamedExpression::template Make<const std::shared_ptr<Node>&, std::string const&>))
			.def("root", +[](NamedExpression const& h) { return h.EntryNode(); }, (arg("self")), "the defining expression this name stands for")
			;

		}

	}
}

