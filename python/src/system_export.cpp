//This file is part of Bertini 2.
//
//python/system_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/system_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/system_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/system_export.cpp:  Source file for exposing systems to python, including start systems.

#include <stdio.h>
#include "system_export.hpp"



namespace bertini{
	namespace python{
		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

		
		struct StartSystemWrap : start_system::StartSystem, wrapper<start_system::StartSystem>
		{
			unsigned long long NumStartPoints() const {return this->get_override("NumStartPoints")(); }
		}; // re: StartSystemWrap

		
		
		
		
		template<typename SystemBaseT>
		template<class PyClass>
		void SystemVisitor<SystemBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("precision", get_prec_)
			.def("precision", set_prec_)
			.def("differentiate", &SystemBaseT::Differentiate)

			.def("eval", return_Eval0_ptr<dbl>() ,"evaluate the system in double precision, using already-set variable values.")
			.def("eval", return_Eval0_ptr<mpfr>() ,"evaluate the system in multiple precision, using already-set variable values.")
			.def("eval", return_Eval1_ptr<dbl>() ,"evaluate the system in double precision, using space variable values passed into this function.")
			.def("eval", return_Eval1_ptr<mpfr>() ,"evaluate the system in multiple precision, using space variable values passed into this function.")
			.def("eval", return_Eval2_ptr<dbl>() ,"evaluate the system in double precision using space and time values passed into this function")
			.def("eval", return_Eval2_ptr<mpfr>() ,"evaluate the system in multiple precision using space and time values passed into this function")
			
			.def("jacobian", return_Jac0_ptr<dbl>() ,"evaluate the jacobian of the system, using already-set time and space value.")
			.def("jacobian", return_Jac0_ptr<mpfr>() ,"evaluate the jacobian of the system, using already-set time and space value.")
			.def("jacobian", return_Jac1_ptr<dbl>() ,"evaluate the jacobian of the system, using space values you pass in to this function")
			.def("jacobian", return_Jac1_ptr<mpfr>() ,"evaluate the jacobian of the system, using space values you pass in to this function")
			.def("jacobian", return_Jac2_ptr<dbl>() , "evaluate the jacobian of the system, using time and space values passed into this function")
			.def("jacobian", return_Jac2_ptr<mpfr>() , "evaluate the jacobian of the system, using time and space values passed into this function")

			.def("homogenize", &SystemBaseT::Homogenize)
			.def("is_homogeneous", &SystemBaseT::IsHomogeneous)
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
			
			.def("auto_patch",&SystemBaseT::AutoPatch,"Apply a patch to the system, given its current variable group structure.")
			.def("copy_patches",&SystemBaseT::CopyPatches,"Copy the patch from another system into this one.")
			.def("get_patch",&SystemBaseT::GetPatch,"Get the patch from the system.")
			.def("is_patched",&SystemBaseT::IsPatched,"Check whether the system is patched.")

			.def("rescale_point_to_fit_patch",&SystemBaseT::template RescalePointToFitPatch<dbl>,"Return a rescaled version of the input point, which fits the patch for the system.")
			.def("rescale_point_to_fit_patch",&SystemBaseT::template RescalePointToFitPatch<mpfr>,"Return a rescaled version of the input point, which fits the patch for the system.")

			.def("rescale_point_to_fit_patch_in_place",&SystemBaseT::template RescalePointToFitPatchInPlace<dbl>,"Re-scale the input point, in place, to fit the patch for the system.  This assumes you have properly set the variable groups and auto-patched the system.")

			.def("rescale_point_to_fit_patch_in_place",&SystemBaseT::template RescalePointToFitPatchInPlace<mpfr>,"Re-scale the input point, in place, to fit the patch for the system.  This assumes you have properly set the variable groups and auto-patched the system.")

			.def("dehomogenize_point",&SystemBaseT::template DehomogenizePoint<dbl>)
			.def("dehomogenize_point",&SystemBaseT::template DehomogenizePoint<mpfr>)

			.def(self_ns::str(self_ns::self))
			.def(self_ns::repr(self_ns::self))
			.def(self += self)
			.def(self + self) 
			.def(self *= std::shared_ptr<node::Node>())
			.def(self * std::shared_ptr<node::Node>())
			.def(std::shared_ptr<node::Node>() * self)
			;
		}
		
		
		
		
		
		template<typename SystemBaseT>
		template<class PyClass>
		void StartSystemVisitor<SystemBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("num_start_points", &SystemBaseT::NumStartPoints)
			.def("start_pointd", return_GenStart_ptr<dbl>() )
			.def("start_pointmp", return_GenStart_ptr<mpfr>() )
			;


		};
		
		void ExportSystem()
		{
			
			// System class
			class_<System, std::shared_ptr<System> >("System", init<>())
			.def(SystemVisitor<System>())
			;
			
			// StartSystem class
			class_<StartSystemWrap, boost::noncopyable, bases<System>, std::shared_ptr<start_system::StartSystem> >("StartSystem", no_init)
			.def(StartSystemVisitor<start_system::StartSystem>())
			;

			// TotalDegree class
			class_<start_system::TotalDegree, bases<start_system::StartSystem>, std::shared_ptr<start_system::TotalDegree> >("TotalDegree", init<System const&>())
			.def("random_value", &start_system::TotalDegree::RandomValue<dbl>)
			.def("random_value", &start_system::TotalDegree::RandomValue<mpfr>)
			.def("random_values", &start_system::TotalDegree::RandomValues, return_value_policy<copy_const_reference>())
			;

			
		}


	}
}
