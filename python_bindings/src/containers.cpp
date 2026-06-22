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
//  Spring 2018
//
//
//  python/src/containers.cpp:  source file for exposing trackers to python.


#include "containers_export.hpp"

#include <boost/python/iterator.hpp>

namespace bertini{
	namespace python{

template<typename T>
template<typename PyClass>
void ListVisitor<T>::visit(PyClass& cl) const
{
	cl

	.def(vector_indexing_suite< T , true >())
	// By default indexed elements are returned by proxy. This can be
    // disabled by supplying *true* in the NoProxy template parameter.

	// vector_indexing_suite only provides the __getitem__/__len__ sequence protocol; add a real
	// __iter__ so `for s in solver.solutions(): ...` (and any other list container) iterates
	// directly rather than relying on the index-fallback.
	.def("__iter__", boost::python::iterator<T>())

	.def("__str__", &ListVisitor::__str__)
	.def("__repr__", &ListVisitor::__repr__)
	;
}


void ExportContainers()
{
	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".container");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("container") = new_submodule;

	scope new_submodule_scope = new_submodule;
	new_submodule_scope.attr("__doc__") = "Various container types";
 

	boost::python::converter::registry::push_back(&pylist_converter<bertini::VariableGroup>::convertible
	    , &pylist_converter<bertini::VariableGroup>::construct
	    , boost::python::type_id<bertini::VariableGroup>());

	// allow a Python list of Functions to convert to std::vector<Function ptr>
	// (used by System([f0,f1,...]) and add_functions)
	using VecFn = std::vector<std::shared_ptr<bertini::node::Function>>;
	boost::python::converter::registry::push_back(&pylist_converter<VecFn>::convertible
	    , &pylist_converter<VecFn>::construct
	    , boost::python::type_id<VecFn>());

	// allow a Python list of VariableGroups to convert to std::vector<VariableGroup>
	// (used by System::set_variable_groups)
	using VecVarGroup = std::vector<bertini::VariableGroup>;
	boost::python::converter::registry::push_back(&pylist_converter<VecVarGroup>::convertible
	    , &pylist_converter<VecVarGroup>::construct
	    , boost::python::type_id<VecVarGroup>());



	// std::vector of Rational Node ptrs
	using T1 = std::vector<std::shared_ptr< bertini::node::Rational > >;
	class_< T1 >("ListOfRational")
	.def(ListVisitor<T1>())
	;

	// The VariableGroup vector container
	using T2 = bertini::VariableGroup;
	class_< T2 >("VariableGroup")
	.def(ListVisitor<T2>())
	.def("__init__", boost::python::make_constructor(&create_MyClass<T2>))
	;
	
	// std::vector of ints
	using T3 = std::vector<int>;
	class_< T3 >("ListOfInt")
	.def(ListVisitor<T3>())
	;


	// std::vector of VariableGroups
	using T4 = std::vector<bertini::VariableGroup>;
	class_< T4 >("ListOfVariableGroup")
	.def(ListVisitor<T4>())
	;


	// std::vector of Function Node ptrs
	using T5 = std::vector<std::shared_ptr< bertini::node::Function > >;
	class_< T5 >("ListOfFunction")
	.def(ListVisitor<T5>())
	;




	// std::vector of Eigen::matrix
	using T7 = std::vector<bertini::Vec<dbl_complex>>;
	class_< T7 >("ListOfVectorComplexDoublePrecision")
	.def(ListVisitor<T7>())
	;

	// std::vector of Eigen::matrix
	using T8 = std::vector<bertini::Vec<mpfr_complex>>;
	class_< T8 >("ListOfVectorComplexVariablePrecision")
	.def(ListVisitor<T8>())
	;

	using T9 = std::vector<bertini::algorithm::SolutionMetaData<dbl_complex>>;
	class_< T9 >("ListOfSolutionMetaData_DoublePrec")
	.def(ListVisitor<T9>())
	;

	using T10 = std::vector<bertini::algorithm::SolutionMetaData<mpfr_complex>>;
	class_< T10 >("ListOfSolutionMetaData_MultiPrec")
	.def(ListVisitor<T10>())
	;

	using T11 = std::vector<bertini::algorithm::EGBoundaryMetaData<dbl_complex>>;
	class_< T11 >("ListOfEGBoundaryMetaData_DoublePrec")
	.def(ListVisitor<T11>())
	;

	using T12 = std::vector<bertini::algorithm::EGBoundaryMetaData<mpfr_complex>>;
	class_< T12 >("ListOfEGBoundaryMetaData_MultiPrec")
	.def(ListVisitor<T12>())
	;
}; // export containers

	}
}