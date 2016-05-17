//This file is part of Bertini 2.
//
//python/system.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/system.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/system.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/system.cpp:  the source file for the python interface for System classes.

#include <stdio.h>
#include "system.hpp"

using namespace boost::python;
using namespace bertini;


namespace bertini{
	namespace python{
		
		
		void ExportSystem()
		{
			class_< bertini::System, std::shared_ptr<bertini::System> >("System", init<>())
			.def("precision", &System::precision)
			.def("differentiate", &System::Differentiate)
			.def("eval", sysEval1<dbl>)
			.def("eval", sysEval2<dbl>)
			.def("jacobian", sysJac1<dbl>)
			.def("jacobian", sysJac2<dbl>)
			.def("homogenize", &System::Homogenize)
			.def("is_homogenous", &System::IsHomogeneous)
			.def("is_polynomial", &System::IsPolynomial)
			
			.def("num_functions", &System::NumFunctions)
			.def("num_variables", &System::NumVariables)
			.def("num_hom_variables", &System::NumHomVariables)
			.def("num_variable_groups", &System::NumVariableGroups)
			.def("num_ungrouped_variables", &System::NumUngroupedVariables)
			.def("num_hom_variable_groups", &System::NumHomVariableGroups)
			.def("num_constants", &System::NumConstants)
			.def("num_parameters", &System::NumParameters)
			.def("num_implicit_parameters", &System::NumImplicitParameters)
			
			.def("set_variables_dbl", &System::SetVariables<dbl>)
			.def("set_path_variable_dbl", &System::SetPathVariable<dbl>)
			.def("set_implicit_parameters_dbl", &System::SetImplicitParameters<dbl>)
			
			.def("add_variable_group", &System::AddVariableGroup)
			.def("add_hom_variable_group", &System::AddHomVariableGroup)
			.def("add_ungrouped_variable", &System::AddUngroupedVariable)
			.def("add_ungrouped_variables", &System::AddUngroupedVariables)
			.def("add_implicit_parameter", &System::AddImplicitParameter)
			.def("add_implicit_parameters", &System::AddImplicitParameters)
			.def("add_parameter", &System::AddParameter)
			.def("add_parameters", &System::AddParameters)
			.def("add_subfunction", &System::AddSubfunction)
			.def("add_subfunctions", &System::AddSubfunctions)
			.def("add_function", sysAddFunc1)
			.def("add_function", sysAddFunc2)
			.def("add_functions", &System::AddFunctions)
			.def("add_constant", &System::AddConstant)
			.def("add_constants", &System::AddConstants)
			.def("add_path_variable", &System::AddPathVariable)
			.def("have_path_variable", &System::HavePathVariable)
			
			.def("function", &System::Function)
			.def("variable_groups", &System::VariableGroups)
			.def("hom_variable_groups", &System::HomVariableGroups)
			.def("degrees", sysDeg1)
			.def("degrees", sysDeg2)
			.def("reorder_functions_by_degree_decreasing", &System::ReorderFunctionsByDegreeDecreasing)
			.def("reorder_functions_by_degree_increasing", &System::ReorderFunctionsByDegreeIncreasing)
			.def("clear_variables", &System::ClearVariables)
			.def("copy_variable_structure", &System::CopyVariableStructure)
			
			
			.def(self_ns::str(self_ns::self))
			.def(self_ns::repr(self_ns::self))
			.def(self += self)
			.def(self + self) //Test this, operator may not have correct signature
			.def(self *= Nd())
			.def(self * Nd())
			.def(Nd() * self)
			;
		}
		
		
	}// re: python
}// re: bertini
