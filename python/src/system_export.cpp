//
//  system_export.cpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/9/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#include <stdio.h>
#include "system_export.hpp"



namespace bertini{
	namespace python{
		
		template<typename SystemBaseT>
		template<class PyClass>
		void SystemVisitor<SystemBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("precision", &SystemBaseT::precision)
			.def("differentiate", &SystemBaseT::Differentiate)
			.def("eval", return_Eval1_ptr<dbl>() )
			.def("eval", return_Eval1_ptr<mpfr>() )
			.def("eval", return_Eval2_ptr<dbl>() )
			.def("eval", return_Eval2_ptr<mpfr>() )
			.def("jacobian", return_Jac1_ptr<dbl>() )
			.def("jacobian", return_Jac1_ptr<mpfr>() )
			.def("jacobian", return_Jac2_ptr<dbl>() )
			.def("jacobian", return_Jac2_ptr<mpfr>() )
			.def("homogenize", &SystemBaseT::Homogenize)
			.def("is_homogenous", &SystemBaseT::IsHomogeneous)
			.def("is_polynomial", &SystemBaseT::IsPolynomial)
			
			.def("num_functions", &SystemBaseT::NumFunctions)
			.def("num_variables", &SystemBaseT::NumVariables)
			.def("num_hom_variables", &SystemBaseT::NumHomVariables)
			.def("num_variable_groups", &SystemBaseT::NumVariableGroups)
			.def("num_ungrouped_variables", &SystemBaseT::NumUngroupedVariables)
			.def("num_hom_variable_groups", &SystemBaseT::NumHomVariableGroups)
			.def("num_constants", &SystemBaseT::NumConstants)
			.def("num_parameters", &SystemBaseT::NumParameters)
			.def("num_implicit_parameters", &SystemBaseT::NumImplicitParameters)
			
			.def("set_variables", &SystemBaseT::template SetVariables<dbl>)
			.def("set_variables", &SystemBaseT::template SetVariables<mpfr>)
			.def("set_path_variable", &SystemBaseT::template SetPathVariable<dbl>)
			.def("set_path_variable", &SystemBaseT::template SetPathVariable<mpfr>)
			.def("set_implicit_parameters", &SystemBaseT::template SetImplicitParameters<dbl>)
			.def("set_implicit_parameters", &SystemBaseT::template SetImplicitParameters<mpfr>)
			
			.def("add_variable_group", &SystemBaseT::AddVariableGroup)
			.def("add_hom_variable_group", &SystemBaseT::AddHomVariableGroup)
			.def("add_ungrouped_variable", &SystemBaseT::AddUngroupedVariable)
			.def("add_ungrouped_variables", &SystemBaseT::AddUngroupedVariables)
			.def("add_implicit_parameter", &SystemBaseT::AddImplicitParameter)
			.def("add_implicit_parameters", &SystemBaseT::AddImplicitParameters)
			.def("add_parameter", &SystemBaseT::AddParameter)
			.def("add_parameters", &SystemBaseT::AddParameters)
			.def("add_subfunction", &SystemBaseT::AddSubfunction)
			.def("add_subfunctions", &SystemBaseT::AddSubfunctions)
			.def("add_function", sysAddFunc1)
			.def("add_function", sysAddFunc2)
			.def("add_functions", &SystemBaseT::AddFunctions)
			.def("add_constant", &SystemBaseT::AddConstant)
			.def("add_constants", &SystemBaseT::AddConstants)
			.def("add_path_variable", &SystemBaseT::AddPathVariable)
			.def("have_path_variable", &SystemBaseT::HavePathVariable)
			
			.def("function", &SystemBaseT::Function)
			.def("variable_groups", &SystemBaseT::VariableGroups)
			.def("hom_variable_groups", &SystemBaseT::HomVariableGroups)
			.def("degrees", sysDeg1)
			.def("degrees", sysDeg2)
			.def("reorder_functions_by_degree_decreasing", &SystemBaseT::ReorderFunctionsByDegreeDecreasing)
			.def("reorder_functions_by_degree_increasing", &SystemBaseT::ReorderFunctionsByDegreeIncreasing)
			.def("clear_variables", &SystemBaseT::ClearVariables)
			.def("copy_variable_structure", &SystemBaseT::CopyVariableStructure)
			
			
			.def(self_ns::str(self_ns::self))
			.def(self_ns::repr(self_ns::self))
			.def(self += self)
			.def(self + self) //Test this, operator may not have correct signature
			.def(self *= std::shared_ptr<node::Node>())
			.def(self * std::shared_ptr<node::Node>())
			.def(std::shared_ptr<node::Node>() * self)
			;
		}

		
		void ExportSystem()
		{
			
			// System class
			class_<System, std::shared_ptr<bertini::System> >("System", init<>())
			.def(SystemVisitor<bertini::System>())
			;
			
		}


	}
}