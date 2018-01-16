//This file is part of Bertini 2.
//
//python/containers_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/containers_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/containers_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/containers_export.hpp:  Exports all needed containers from Bertini 2.0 to python.

#pragma once
#ifndef BERTINI_PYTHON_CONTAINERS_EXPORT_HPP
#define BERTINI_PYTHON_CONTAINERS_EXPORT_HPP


#include <bertini2/function_tree.hpp>

#include "python_common.hpp"




void ExportContainers()
{
	// std::vector of Rational Node ptrs
	class_< std::vector<std::shared_ptr< bertini::node::Rational > > >("ListRational")
	.def(vector_indexing_suite< std::vector<std::shared_ptr< bertini::node::Rational> > , true >())
	;

	// The VariableGroup deque container
	class_< bertini::VariableGroup >("VariableGroup")
	.def(vector_indexing_suite< bertini::VariableGroup, true >())
	;
	
	// std::vector of ints
	class_< std::vector<int> >("ListInt")
	.def(vector_indexing_suite< std::vector<int> , true >())
	;

	// std::vector of VariableGroups
	class_< std::vector<bertini::VariableGroup> >("ListVariableGroup")
	.def(vector_indexing_suite< std::vector<bertini::VariableGroup> , true >())
	;

	// std::vector of Function Node ptrs
	class_< std::vector<std::shared_ptr< bertini::node::Function > > >("ListFunction")
	.def(vector_indexing_suite< std::vector<std::shared_ptr< bertini::node::Function> > , true >())
	;

	// std::vector of Jacobian Node ptrs
	class_< std::vector<std::shared_ptr< bertini::node::Jacobian > > >("ListJacobian")
	.def(vector_indexing_suite< std::vector<std::shared_ptr< bertini::node::Jacobian> > , true >())
	;

};








#endif
