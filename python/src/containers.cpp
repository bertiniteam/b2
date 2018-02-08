//This file is part of Bertini 2.
//
//python/src/containers.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/src/containers.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/src/containers.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017-2018 by Bertini2 Development Team
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
//  python/src/containers.cpp:  source file for exposing trackers to python.


#include "containers_export.hpp"

namespace bertini{
	namespace python{


void ExportContainers()
{
	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".list");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("list") = new_submodule;

	scope new_submodule_scope = new_submodule;
	new_submodule_scope.attr("__doc__") = "List types for PyBertini";


	// std::vector of Rational Node ptrs
	class_< std::vector<std::shared_ptr< bertini::node::Rational > > >("Rational")
	.def(vector_indexing_suite< std::vector<std::shared_ptr< bertini::node::Rational> > , true >())
	;

	// The VariableGroup deque container
	class_< bertini::VariableGroup >("VariableGroup")
	.def(vector_indexing_suite< bertini::VariableGroup, true >())
	;
	
	// std::vector of ints
	class_< std::vector<int> >("int")
	.def(vector_indexing_suite< std::vector<int> , true >())
	;

	// std::vector of VariableGroups
	class_< std::vector<bertini::VariableGroup> >("VariableGroup")
	.def(vector_indexing_suite< std::vector<bertini::VariableGroup> , true >())
	;

	// std::vector of Function Node ptrs
	class_< std::vector<std::shared_ptr< bertini::node::Function > > >("Function")
	.def(vector_indexing_suite< std::vector<std::shared_ptr< bertini::node::Function> > , true >())
	;

	// std::vector of Jacobian Node ptrs
	class_< std::vector<std::shared_ptr< bertini::node::Jacobian > > >("Jacobian")
	.def(vector_indexing_suite< std::vector<std::shared_ptr< bertini::node::Jacobian> > , true >())
	;

};

	}
}