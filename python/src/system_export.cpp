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
//  Danielle Brake
//  UWEC
//  Spring 2018
//
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
			.def("precision", get_prec_,"Get the current precision of the system.  Returns a postive number, representing the number of digits (not bits) at which the system is currently represented.  (there is a reference-level precision stored, so you can change this up / down mostly fearlessly)")
			.def("precision", set_prec_,"Set / change the precision of the system.  Feed in a positive number, representing the digits (not bits) of the precision.  Double precision is 16, but that only effects the multi-precision precision...  you can eval in double precision without changing the precision to 16.")
			.def("differentiate", &SystemBaseT::Differentiate)

			.def("eval", return_Eval0_ptr<dbl>() ,"Evaluate the system in double precision, using already-set variable values.")
			.def("eval", return_Eval0_ptr<mpfr>() ,"Evaluate the system in multiple precision, using already-set variable values.")
			.def("eval", return_Eval1_ptr<dbl>() ,"Evaluate the system in double precision, using space variable values passed into this function.")
			.def("eval", return_Eval1_ptr<mpfr>() ,"Evaluate the system in multiple precision, using space variable values passed into this function.")
			.def("eval", return_Eval2_ptr<dbl>() ,"Evaluate the system in double precision using space and time values passed into this function.  Throws if doesn't use a time variable")
			.def("eval", return_Eval2_ptr<mpfr>() ,"Evaluate the system in multiple precision using space and time values passed into this function.  Throws if doesn't use a time variable")
			
			.def("jacobian", return_Jac0_ptr<dbl>() ,"Evaluate the Jacobian (martix of partial derivatives) of the system, using already-set time and space value.")
			.def("jacobian", return_Jac0_ptr<mpfr>() ,"Evaluate the Jacobian (martix of partial derivatives) of the system, using already-set time and space value.")
			.def("jacobian", return_Jac1_ptr<dbl>() ,"Evaluate the Jacobian (martix of partial derivatives) of the system, using space values you pass in to this function")
			.def("jacobian", return_Jac1_ptr<mpfr>() ,"Evaluate the Jacobian (martix of partial derivatives) of the system, using space values you pass in to this function")
			.def("jacobian", return_Jac2_ptr<dbl>() , "Evaluate the Jacobian (martix of partial derivatives) of the system, using time and space values passed into this function.  Throws if doesn't use a time variable")
			.def("jacobian", return_Jac2_ptr<mpfr>() , "Evaluate the Jacobian (martix of partial derivatives) of the system, using time and space values passed into this function.  Throws if doesn't use a time variable")

			.def("homogenize", &SystemBaseT::Homogenize,"Homogenize the system, adding new homogenizing variables if necessary.  This may change your polynomials; that is, it has side effects.")
			.def("is_homogeneous", &SystemBaseT::IsHomogeneous, "Determines whether all polynomials in the system have the same degree.  Non-polynomial functions are not homogeneous.")
			.def("is_polynomial", &SystemBaseT::IsPolynomial, "Determines whether all polynomials are polynomial.  Transcendental functions, e.g., are non-polynomial.  Returns a bool.")
			
			.def("num_functions", &SystemBaseT::NumFunctions,"The total number of functions in the system.  Does not include patches.")
			.def("num_variables", &SystemBaseT::NumVariables,"the *total* number of variables in the system.  Includes homogenizing variables")
			.def("num_hom_variables", &SystemBaseT::NumHomVariables, "The number of homogenizing variables defined in the system.  Should be equal to the number of homvargroups")
			.def("num_variable_groups", &SystemBaseT::NumVariableGroups,"The number of affine variable groups.  This should probably be renamed to num_affine_variable_groups")
			.def("num_ungrouped_variables", &SystemBaseT::NumUngroupedVariables,"The number of variables, not grouped into an affine or projective space")
			.def("num_hom_variable_groups", &SystemBaseT::NumHomVariableGroups,"The number of homogeneous or projective variable groups.  The number of homogenizing variables should eventually equal this.")
			// .def("num_constants", &SystemBaseT::NumConstants,"Has no impact on anything.  The number of constants in the system.")
			// .def("num_parameters", &SystemBaseT::NumParameters,"Has no impact on anything.  The number of 'parameters' in the system.")
			// .def("num_implicit_parameters", &SystemBaseT::NumImplicitParameters,"Has no impact on anything.  The number of 'implicit parameters' in the system.") // commented out until implemented
			
			.def("set_variables", &SystemBaseT::template SetVariables<dbl>,"Set the values of the variables. Expects a vector of doubles")
			.def("set_variables", &SystemBaseT::template SetVariables<mpfr>,"Set the values of the variables. Expects a vector of complex mpfr's")
			.def("set_path_variable", &SystemBaseT::template SetPathVariable<dbl>,"Set the value of the path variable.  This one's double-precision.  Throws if path variable not defined.")
			.def("set_path_variable", &SystemBaseT::template SetPathVariable<mpfr>,"Set the value of the path variable.  This one's variable-precision.  Throws if path variable not defined.")
			// .def("set_implicit_parameters", &SystemBaseT::template SetImplicitParameters<dbl>,"Doesn't do anything.  Sets the values of algebraically constrained parameters")
			// .def("set_implicit_parameters", &SystemBaseT::template SetImplicitParameters<mpfr>,"Doesn't do anything.  Sets the values of algebraically constrained parameters")
			
			.def("add_variable_group", &SystemBaseT::AddVariableGroup,"Add a (affine) variable group to the System")
			.def("add_hom_variable_group", &SystemBaseT::AddHomVariableGroup,"Add a projective or homogeneous variable group to the System")
			// .def("add_ungrouped_variable", &SystemBaseT::AddUngroupedVariable,"Add an ungrouped variable to the system.  I honestly don't know why you'd do that.  This should be removed, and is a holdover from Bertini 1")
			// .def("add_ungrouped_variables", &SystemBaseT::AddUngroupedVariables,"Add some ungrouped variables to the system.  I honestly don't know why you'd do that.  This should be removed, and is a holdover from Bertini 1")
			// .def("add_implicit_parameter", &SystemBaseT::AddImplicitParameter)
			// .def("add_implicit_parameters", &SystemBaseT::AddImplicitParameters)
			// .def("add_parameter", &SystemBaseT::AddParameter)
			// .def("add_parameters", &SystemBaseT::AddParameters)
			// .def("add_subfunction", &SystemBaseT::AddSubfunction)
			// .def("add_subfunctions", &SystemBaseT::AddSubfunctions)
			.def("add_function", sysAddFunc1,"Add a function to the System")
			.def("add_function", sysAddFunc2,"Add a function to the System")
			.def("add_functions", &SystemBaseT::AddFunctions,"Add some functions to the System.  Expects a list of functions")
			// .def("add_constant", &SystemBaseT::AddConstant)
			// .def("add_constants", &SystemBaseT::AddConstants)
			.def("add_path_variable", &SystemBaseT::AddPathVariable,"Add a path variable to the System")
			.def("have_path_variable", &SystemBaseT::HavePathVariable,"Asks whether the System has a path variable defined")
			
			.def("function", &SystemBaseT::Function,"Get a function with a given index.  Problems ensue if out of range -- uses un-rangechecked version of underlying getter")
			.def("variable_groups", &SystemBaseT::VariableGroups, "Get the list of (affine) variable_groups from the system")
			.def("hom_variable_groups", &SystemBaseT::HomVariableGroups, "Get the list of projective / homogeneous variable_groups from the system")
			.def("degrees", sysDeg1, "Get a list of the degrees of the functions in the system, with respect to all variables in all groups (and in fact overall)")
			.def("degrees", sysDeg2, "Get a list of the degrees of the functions in the system, with respect to a variable_group passed in to this function.  Negative numbers indicate non-polynomial")
			.def("reorder_functions_by_degree_decreasing", &SystemBaseT::ReorderFunctionsByDegreeDecreasing,"Change the order of the functions to be in decreasing order")
			.def("reorder_functions_by_degree_increasing", &SystemBaseT::ReorderFunctionsByDegreeIncreasing,"Change the order of the functions to be in decreasing order")
			.def("clear_variables", &SystemBaseT::ClearVariables, "Remove the variable structure from the system")
			.def("copy_variable_structure", &SystemBaseT::CopyVariableStructure, "Copy the variable structure from another System")
			
			.def("auto_patch",&SystemBaseT::AutoPatch,"Apply a patch to the system, given its current variable group structure.")
			.def("copy_patches",&SystemBaseT::CopyPatches,"Copy the patches from another system into this one.")
			.def("get_patch",&SystemBaseT::GetPatch,"Get (a reference to) the patches from the system.")
			.def("is_patched",&SystemBaseT::IsPatched,"Check whether the system is patched.")

			.def("rescale_point_to_fit_patch",&SystemBaseT::template RescalePointToFitPatch<dbl>,"Return a rescaled version of the input point, which fits the patch for the system.")
			.def("rescale_point_to_fit_patch",&SystemBaseT::template RescalePointToFitPatch<mpfr>,"Return a rescaled version of the input point, which fits the patch for the system.")

			.def("rescale_point_to_fit_patch_in_place",&SystemBaseT::template RescalePointToFitPatchInPlace<dbl>,"Re-scale the input point, in place, to fit the patch for the system.  This assumes you have properly set the variable groups and auto-patched the system.")

			.def("rescale_point_to_fit_patch_in_place",&SystemBaseT::template RescalePointToFitPatchInPlace<mpfr>,"Re-scale the input point, in place, to fit the patch for the system.  This assumes you have properly set the variable groups and auto-patched the system.")

			.def("dehomogenize_point",&SystemBaseT::template DehomogenizePoint<dbl>, "Dehomogenize a vector of doubles (complex), using the variable structure in this System")
			.def("dehomogenize_point",&SystemBaseT::template DehomogenizePoint<mpfr>, "Dehomogenize a vector of mpfr's (complex), using the variable structure in this System")

			.def(self_ns::str(self_ns::self))//, "String representation of the system"
			.def(self_ns::repr(self_ns::self))//, "Round-trippable representation of the system.  Probably not functional"
			.def(self += self)
			.def(self + self) 
			.def(self *= std::shared_ptr<node::Node>())//, "'Scalar-multiply' a system"
			.def(self * std::shared_ptr<node::Node>())//, "'Scalar-multiply' a system"
			.def(std::shared_ptr<node::Node>() * self)//, "'Scalar-multiply' a system"
			;
		}
		
		
		
		
		
		template<typename SystemBaseT>
		template<class PyClass>
		void StartSystemVisitor<SystemBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("num_start_points", &SystemBaseT::NumStartPoints, "Get the number of start points that would be required by the system.  Non-negative, unsigned")
			.def("start_point_d", return_GenStart_ptr<dbl>(),"Get the k-th start point in double precision")
			.def("start_point_mp", return_GenStart_ptr<mpfr>(),"Get the k-th start point in current multiple precision")
			;


		};
		
		void ExportAllSystems()
		{

			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".system");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("system") = new_submodule;
			

			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "Systems of functions, for tracking &c.";

			ExportSystem();
			ExportStartSystems();
		}

		void call_simplify(object obj){
				System& sys=extract<System&>(obj)();
				Simplify(sys);
			};

		void ExportSystem()
		{
			
			// System class
			class_<System, std::shared_ptr<System> >("System", init<>())
			.def(SystemVisitor<System>())
			;

			// free functions
			def("concatenate", &Concatenate, "concatenate two Systems to produce a new one.  Appends the second onto what was the first.");
			def("clone", &Clone, "Make a complete clone of a System.  Includes all functions, variables, etc.  Truly and genuinely distinct.");

			


			def("simplify", &call_simplify, "Perform all possible simplifications.  Has side effects of modifying your functions, if held separately.  Shared nodes between multiple systems may have adverse effects");
			
		}

		void ExportStartSystems()
		{

			{ // enter a scope for config types
				scope current_scope;
				std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
				new_submodule_name.append(".start_system");
				object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
				current_scope.attr("start_system") = new_submodule;

				scope new_submodule_scope = new_submodule;

				ExportStartSystemBase();
				ExportTotalDegree();
			}
		}



		void ExportStartSystemBase()
		{
			//
			// StartSystem class
			class_<StartSystemWrap, boost::noncopyable, bases<System>, std::shared_ptr<start_system::StartSystem> >("AbstractStartSystem", no_init)
			.def(StartSystemVisitor<start_system::StartSystem>())
			;
		}



		void ExportTotalDegree()
		{
			class_<start_system::TotalDegree, bases<start_system::StartSystem>, std::shared_ptr<start_system::TotalDegree> >("TotalDegree",init<System const&>())//,"Only constructor for a TotalDegree start system, requires a system.  You cannot construct one without.  If this is a problem, please contact the authors for help.")
			.def("random_value", &start_system::TotalDegree::RandomValue<dbl>, "Get the k-th random value, in double precision")
			.def("random_value", &start_system::TotalDegree::RandomValue<mpfr>, "Get the k-th random value, in current multiple precision")
			.def("random_values", &start_system::TotalDegree::RandomValues, return_value_policy<copy_const_reference>(), "Get (a reference to) the random values for the start system, as Nodes")
			;

			
		}
	}
}
