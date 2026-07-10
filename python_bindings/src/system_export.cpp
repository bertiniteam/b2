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
//
//
//  python/system_export.cpp:  Source file for exposing systems to python, including start systems.

#include <stdio.h>
#include <sstream>
#include <boost/python/stl_iterator.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "system_export.hpp"
#include <bertini2/system/slice.hpp>   // for System::Slices() -> list of Slice
#include <bertini2/io/classic_writer.hpp>



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
			.def("precision", get_prec_, (arg("self")), "Get the current precision of the system.  Returns a postive number, representing the number of digits (not bits) at which the system is currently represented.  (there is a reference-level precision stored, so you can change this up / down mostly fearlessly)")
			.def("precision", set_prec_, (arg("self"), arg("precision")),"Set / change the precision of the system.  Feed in a positive number, representing the digits (not bits) of the precision.  Double precision is 16, but that only effects the multi-precision precision...  you can eval in double precision without changing the precision to 16.")
			.def("differentiate", &SystemBaseT::Differentiate, (arg("self")), "differentiate the system with respect to the declared variable groups")

			.def("to_classic_input",
				+[](SystemBaseT const& self, int mptype, int odepredictor,
				    double tracktolbeforeeg, double tracktolduringeg, double finaltol,
				    double maxstepsize, double stepsuccessfactor, double stepfailfactor,
				    unsigned stepsforincrease, unsigned long maxnumbersteps, unsigned maxnewtonits,
				    unsigned maxcrossedpathresolves){
					bertini::classic::ClassicWriteOptions opt;
					opt.mptype = mptype;                       opt.odepredictor = odepredictor;
					opt.tracktolbeforeeg = tracktolbeforeeg;   opt.tracktolduringeg = tracktolduringeg;
					opt.finaltol = finaltol;                   opt.maxstepsize = maxstepsize;
					opt.stepsuccessfactor = stepsuccessfactor; opt.stepfailfactor = stepfailfactor;
					opt.stepsforincrease = stepsforincrease;   opt.maxnumbersteps = maxnumbersteps;
					opt.maxnewtonits = maxnewtonits;           opt.maxcrossedpathresolves = maxcrossedpathresolves;
					return bertini::classic::SystemToClassicFile(self, opt);
				},
				(arg("self"), arg("mptype") = 2, arg("odepredictor") = 5,
				 arg("tracktolbeforeeg") = 1e-5, arg("tracktolduringeg") = 1e-6, arg("finaltol") = 1e-11,
				 arg("maxstepsize") = 0.1, arg("stepsuccessfactor") = 2.0, arg("stepfailfactor") = 0.5,
				 arg("stepsforincrease") = 5u, arg("maxnumbersteps") = 100000ul, arg("maxnewtonits") = 2u,
				 arg("maxcrossedpathresolves") = 2u),
				"Emit this system as a Bertini 1 classic input file (a CONFIG + INPUT string) so the same problem can be solved in Bertini 1 for cross-validation.  Every tracking knob that governs path resolution -- predictor, tolerances, and the FULL step-size cadence (maxstepsize / stepsuccessfactor / stepfailfactor / stepsforincrease) plus maxnewtonits -- is settable, so the emitted file is fully controlled against a Bertini 2 solve (defaults mirror Bertini 2's).  mptype: 0 double, 1 fixed-multiple, 2 adaptive.  odepredictor: 0 Euler, 2 RK4, 5 RKF45, 6 Cash-Karp.  AMP coeff/degree bounds are derived from the system.  Run `bertini1` on the result in a SCRATCH dir (it writes many files into its CWD).")

			// Register mpfr overloads first so complex_dbl overloads have highest priority
			// (boost::python resolves in LIFO order). Without this, a numpy int64
			// array causes eigenpy to probe Vec<mpfr> construction first; with
			// thread_default_precision=0 on Boost>=1.87 that aborts in mpfr_init2.
			.def("eval", return_Eval0_ptr<mpfr>(), (arg("self")) ,"Evaluate the system in multiple precision, using already-set variable values.")
			.def("eval", return_Eval0_ptr<complex_dbl>(), (arg("self")) ,"Evaluate the system in double precision, using already-set variable values.")
			.def("eval", return_Eval1_ptr<mpfr>(), (arg("self")) ,"Evaluate the system in multiple precision, using space variable values passed into this function.")
			.def("eval", return_Eval1_ptr<complex_dbl>(), (arg("self")) ,"Evaluate the system in double precision, using space variable values passed into this function.")
			.def("eval", return_Eval2_ptr<mpfr>(), (arg("self")) ,"Evaluate the system in multiple precision using space and time values passed into this function.  Throws if doesn't use a time variable")
			.def("eval", return_Eval2_ptr<complex_dbl>(), (arg("self")) ,"Evaluate the system in double precision using space and time values passed into this function.  Throws if doesn't use a time variable")
			
			// these two commented out because i don't need in-place Eigen::Ref wrapping here. 
			// but if you did, you'd use two lines like this, ha.
			// .def("eval", &eval_wrap_1<mpfr>)
			// .def("eval", &eval_wrap_1<complex_dbl>)

			.def("eval_jacobian", return_Jac0_ptr<complex_dbl>(), (arg("self")) ,"Evaluate the Jacobian (martix of partial derivatives) of the system, using already-set time and space value.")
			.def("eval_jacobian", return_Jac0_ptr<mpfr>(), (arg("self")) ,"Evaluate the Jacobian (martix of partial derivatives) of the system, using already-set time and space value.")
			.def("eval_jacobian", return_Jac1_ptr<complex_dbl>(), (arg("self")) ,"Evaluate the Jacobian (martix of partial derivatives) of the system, using space values you pass in to this function")
			.def("eval_jacobian", return_Jac1_ptr<mpfr>(), (arg("self")) ,"Evaluate the Jacobian (martix of partial derivatives) of the system, using space values you pass in to this function")
			.def("eval_jacobian", return_Jac2_ptr<complex_dbl>(), (arg("self")) , "Evaluate the Jacobian (martix of partial derivatives) of the system, using time and space values passed into this function.  Throws if doesn't use a time variable")
			.def("eval_jacobian", return_Jac2_ptr<mpfr>(), (arg("self")) , "Evaluate the Jacobian (martix of partial derivatives) of the system, using time and space values passed into this function.  Throws if doesn't use a time variable")

			.def("eval_time_derivative",
				+[](SystemBaseT const& self, bertini::Vec<mpfr> const& v, mpfr const& t) { return self.TimeDerivative(v, t); },
				(arg("self"), arg("space"), arg("time")), "Evaluate dH/dt (the time derivative) in multiple precision at the given space and time values.  Rows of t-independent blocks are zero.")
			.def("eval_time_derivative",
				+[](SystemBaseT const& self, bertini::Vec<complex_dbl> const& v, complex_dbl const& t) { return self.TimeDerivative(v, t); },
				(arg("self"), arg("space"), arg("time")), "Evaluate dH/dt (the time derivative) in double precision at the given space and time values.  Rows of t-independent blocks are zero.")

			.def("homogenize", static_cast<void (SystemBaseT::*)()>(&SystemBaseT::Homogenize), (arg("self")),"Homogenize the system, adding new homogenizing variables if necessary.  This may change your polynomials; that is, it has side effects.")
			.def("is_homogeneous", &SystemBaseT::IsHomogeneous, (arg("self")), "Determines whether all polynomials in the system have the same degree.  Non-polynomial functions are not homogeneous.")
			.def("is_polynomial", &SystemBaseT::IsPolynomial, (arg("self")), "Determines whether all polynomials are polynomial.  Transcendental functions, e.g., are non-polynomial.  Returns a bool.")

			.def("content_digest",
				+[](SystemBaseT const& self){ return self.ContentDigest().Hex(); },
				(arg("self")),
				"The persistent content digest of the system: SHA-256 of its canonical exact encoding, as 64 lowercase hex characters.  Stable across sessions, machines, and versions of the encoding format (a format change bumps the version inside the encoding, changing digests loudly rather than silently).  Everything evaluation-relevant is identity -- functions, variable groups and ordering, path variable, patch and randomization coefficients, gamma; randomness included.  Transient state (precision, current variable values, differentiation) is not.  This is the key a database of solutions references systems and homotopies by.")
			.def("is_same",
				+[](SystemBaseT const& self, System const& other){ return self.IsSame(other); },
				(arg("self"), arg("other")),
				"Content equality: True iff the two systems have equal content digests (identical canonical encodings).  Independently built systems with the same mathematical content compare equal; systems differing in any identity-bearing way (functions, grouping, patch, gamma, ...) do not.  Does NOT change ==/hash semantics of the Python object.")
			.def("seal",
				+[](SystemBaseT& self){ self.Seal(); },
				(arg("self")),
				"Seal the system: memoize its content digest and forbid structural mutation (hashcons-on-freeze).  After sealing, structural mutators (add_function, homogenize, auto_patch, ...) raise; evaluation, precision changes, and differentiation still work.  Copying (clone / deepcopy) yields an unsealed copy.  Idempotent.")
			.def("is_sealed",
				+[](SystemBaseT const& self){ return self.IsSealed(); },
				(arg("self")),
				"Whether the system has been sealed against structural mutation.")
			
			.def("num_functions", &SystemBaseT::NumTotalFunctions, (arg("self")),"The total number of functions in the system.  Does not include patches.")
			.def("num_variables", &SystemBaseT::NumVariables, (arg("self")),"the *total* number of variables in the system.  Includes homogenizing variables")
			.def("num_hom_variables", &SystemBaseT::NumHomVariables, (arg("self")), "The number of homogenizing variables defined in the system.  Should be equal to the number of homvargroups")
			.def("num_variable_groups", &SystemBaseT::NumVariableGroups, (arg("self")),"The number of affine variable groups.  This should probably be renamed to num_affine_variable_groups")
			.def("num_ungrouped_variables", &SystemBaseT::NumUngroupedVariables, (arg("self")),"The number of variables, not grouped into an affine or projective space")
			.def("num_hom_variable_groups", &SystemBaseT::NumHomVariableGroups, (arg("self")),"The number of homogeneous or projective variable groups.  The number of homogenizing variables should eventually equal this.")
			// .def("num_constants", &SystemBaseT::NumConstants,"Has no impact on anything.  The number of constants in the system.")
			// .def("num_parameters", &SystemBaseT::NumParameters,"Has no impact on anything.  The number of 'parameters' in the system.")
			// .def("num_implicit_parameters", &SystemBaseT::NumImplicitParameters,"Has no impact on anything.  The number of 'implicit parameters' in the system.") // commented out until implemented
			
			.def("set_variables", &SystemBaseT::template SetVariables<complex_dbl>, (arg("self"), arg("values")), "Set the values of the variables. Expects a vector of doubles")
			.def("set_variables", &SystemBaseT::template SetVariables<mpfr>, (arg("self"), arg("values")), "Set the values of the variables. Expects a vector of complex mpfr's")
			.def("set_path_variable", &SystemBaseT::template SetPathVariable<complex_dbl>, (arg("self"), arg("values")), "Set the value of the path variable.  This one's double-precision.  Throws if path variable not defined.")
			.def("set_path_variable", &SystemBaseT::template SetPathVariable<mpfr>, (arg("self"), arg("values")), "Set the value of the path variable.  This one's variable-precision.  Throws if path variable not defined.")
			// .def("set_implicit_parameters", &SystemBaseT::template SetImplicitParameters<complex_dbl>,"Doesn't do anything.  Sets the values of algebraically constrained parameters")
			// .def("set_implicit_parameters", &SystemBaseT::template SetImplicitParameters<mpfr>,"Doesn't do anything.  Sets the values of algebraically constrained parameters")
			
			.def("add_variable_group", &SystemBaseT::AddVariableGroup, (arg("self"), arg("group")), "Add a (affine) variable group to the System")
			.def("set_variable_groups", &System::SetVariableGroups, (arg("self"), arg("groups")), "Replace the entire variable-group structure of the System with the given list of (affine) variable groups.  Clears existing groups but preserves the path variable.")
			.def("add_hom_variable_group", &SystemBaseT::AddHomVariableGroup, (arg("self"), arg("group")), "Add a projective or homogeneous variable group to the System")
			.def("add_linear_forms_block",
				+[](SystemBaseT& self, std::size_t num_vars, bertini::Mat<mpfr> const& coefficients) {
					self.AddBlock(bertini::blocks::LinearFormsBlock(num_vars, coefficients));
				},
				(arg("self"), arg("num_vars"), arg("coefficients")),
				"Add a block of affine linear forms f(x) = M [x;1] to the System, evaluated as a single matrix-vector product rather than as scalar expressions.  coefficients is an complex_mp matrix with one row per function and num_vars+1 columns; the trailing column carries each form's constant term.")
			.def("add_products_of_linears_block",
				+[](SystemBaseT& self, std::size_t num_vars, boost::python::list const& factors) {
					std::vector<bertini::Mat<mpfr>> v{
						boost::python::stl_input_iterator<bertini::Mat<mpfr>>(factors),
						boost::python::stl_input_iterator<bertini::Mat<mpfr>>() };
					self.AddBlock(bertini::blocks::ProductsOfLinearsBlock(num_vars, std::move(v)));
				},
				(arg("self"), arg("num_vars"), arg("factors")),
				"Add a block of products-of-linear-forms f_i(x) = prod_r ( c_{i,r} . [x;1] ) to the System, evaluated as matrix-multiplies-then-row-products rather than as scalar expressions.  factors is a list with one entry per function; entry i is an complex_mp matrix with one row per linear factor and num_vars+1 columns (the trailing column carries each factor's constant term).  Each function's degree is its number of factors.")
			.def("slices",
				+[](SystemBaseT const& self) {
					boost::python::list out;
					for (auto const& s : self.Slices())
						out.append(s);
					return out;
				},
				(arg("self")),
				"The linear-form slices embedded in this system, one per linear-forms block (an empty list if none).  Lets you back out the slice structure of a system that was built from a slice -- the inverse of Slice.as_system().")
			// .def("add_ungrouped_variable", &SystemBaseT::AddUngroupedVariable,"Add an ungrouped variable to the system.  I honestly don't know why you'd do that.  This should be removed, and is a holdover from Bertini 1")
			// .def("add_ungrouped_variables", &SystemBaseT::AddUngroupedVariables,"Add some ungrouped variables to the system.  I honestly don't know why you'd do that.  This should be removed, and is a holdover from Bertini 1")
			// .def("add_implicit_parameter", &SystemBaseT::AddImplicitParameter)
			// .def("add_implicit_parameters", &SystemBaseT::AddImplicitParameters)
			// .def("add_parameter", &SystemBaseT::AddParameter)
			// .def("add_parameters", &SystemBaseT::AddParameters)
			// .def("add_subfunction", &SystemBaseT::AddSubfunction)
			// .def("add_subfunctions", &SystemBaseT::AddSubfunctions)
			.def("add_function", AddJustFn, (arg("self"), arg("f")), "Add a function (a bare expression) to the System")
			
			.def("add_functions", &SystemBaseT::AddFunctions, (arg("self"), arg("functions")), "Add some functions to the System.  Expects a list of functions")
			// .def("add_constant", &SystemBaseT::AddConstant)
			// .def("add_constants", &SystemBaseT::AddConstants)
			.def("add_path_variable", &SystemBaseT::AddPathVariable, (arg("self"), arg("pathvar")), "Add a path variable to the System")
			.def("have_path_variable", &SystemBaseT::HavePathVariable, (arg("self")), "Asks whether the System has a path variable defined")
			
			.def("function", &SystemBaseT::Function, (arg("self"), arg("index")), "Get a function with a given index.  Problems ensue if out of range -- uses un-rangechecked version of underlying getter")
			.def("functions",
				+[](SystemBaseT const& self) {
					boost::python::list out;
					for (auto const& f : self.GetNaturalFunctions()) out.append(f);
					return out;
				},
				(arg("self")),
				"The system's functions, as a list of function-tree nodes (issue #297; structured blocks are expanded).  So `critpt_sys.add_functions(sys.functions())` copies them all in.")
			.def("copy_functions",
				+[](SystemBaseT& self, SystemBaseT const& other) -> SystemBaseT& {
					self.CopyFunctions(other);
					return self;
				},
				return_internal_reference<>(),
				(arg("self"), arg("other")),
				"Append another system's functions to this one and return self (issue #297; sugar for add_functions(other.functions())).")
			.def("coordinates_of",
				+[](SystemBaseT const& self, Vec<mpfr> point, boost::python::object group) -> boost::python::object {
					boost::python::extract<unsigned> as_idx(group);
					unsigned idx = as_idx.check() ? as_idx()
						: self.FIFOIndexOfGroup(boost::python::extract<VariableGroup const&>(group)());
					return boost::python::object(self.CoordinatesOfGroup(point, idx));
				},
				(arg("self"), arg("point"), arg("group")),
				"Project a user-coordinate point onto one variable group: return just that group's coordinates.  "
				"`group` is either the VariableGroup object or its 0-based FIFO index.  Affine groups return their "
				"affine coordinates; projective groups are returned as-is (not dehomogenized).  Handy for an "
				"augmented system (e.g. a critical-point system) where you only care about one group.")
			.def("coordinates_of",
				+[](SystemBaseT const& self, Vec<complex_dbl> point, boost::python::object group) -> boost::python::object {
					boost::python::extract<unsigned> as_idx(group);
					unsigned idx = as_idx.check() ? as_idx()
						: self.FIFOIndexOfGroup(boost::python::extract<VariableGroup const&>(group)());
					return boost::python::object(self.CoordinatesOfGroup(point, idx));
				},
				(arg("self"), arg("point"), arg("group")))
			.def("symbolic_jacobian",
				+[](SystemBaseT const& self, bool usercoordinates) {
					auto J = self.SymbolicJacobian(usercoordinates);
					boost::python::list rows;
					for (std::size_t i = 0; i < J.rows; ++i) {
						boost::python::list row;
						for (std::size_t j = 0; j < J.cols; ++j)
							row.append(J.entries[i*J.cols + j]);
						rows.append(row);
					}
					return rows;
				},
				(arg("self"), arg("usercoordinates")=true),
				"The symbolic Jacobian of the system, as a list of rows of expression nodes (NOT numeric -- contrast eval_jacobian).  usercoordinates=True (default): differentiate the natural (pre-homogenization) functions w.r.t. the user-declared affine/projective variable groups -- homogenizing variables never appear, patches omitted.  usercoordinates=False: differentiate the current (possibly homogenized) functions w.r.t. the full internal variable ordering (homogenizing variables included), with the patch's rows appended when patched.  Prefer bertini.System.jacobian(...), which returns a 2-D numpy object array.")
			.def("variable_groups", &SystemBaseT::VariableGroups, (arg("self")), "Get the list of (affine) variable_groups from the system")
			.def("hom_variable_groups", &SystemBaseT::HomVariableGroups, (arg("self")), "Get the list of projective / homogeneous variable_groups from the system")
			.def("degrees", sysDeg1, (arg("self")), "Get a list of the degrees of the functions in the system, with respect to all variables in all groups (and in fact overall)")
			.def("degrees", sysDeg2, (arg("self"), arg("group")), "Get a list of the degrees of the functions in the system, with respect to a variable_group passed in to this function.  Negative numbers indicate non-polynomial")
			.def("randomize",
				+[](SystemBaseT const& self) { return self.Randomize(); },
				(arg("self")),
				"Randomize an overdetermined system (N functions, n variables, N>n) down to a square one, returning a NEW system; this one is left untouched.  The square result has n generic combinations of the original functions, whose isolated solutions still contain this system's -- solve it, then discard the extraneous solutions by re-evaluating this system.  For a single affine variable group the functions are sorted by descending degree and R=[I|C], giving the optimal total-degree path count.")
			.def("randomize",
				+[](SystemBaseT const& self, int codimension) { return self.Randomize(codimension); },
				(arg("self"), arg("codimension")),
				"Randomize down to `codimension` functions (generic combinations of the natural functions), returning a NEW system; this one is left untouched.  Like the no-arg form, but reduces to a chosen number of functions rather than squaring -- `codimension=1` yields a single function, the general \"randomize down to codimension c\" used for lower-dimensional components.  Raises if codimension < 1 or >= the number of natural functions (randomization must reduce the count).")
			.def("randomize",
				+[](SystemBaseT const& self, bertini::Mat<mpfr> const& R) { return self.Randomize(R); },
				(arg("self"), arg("matrix")),
				"Randomize using a supplied coefficient matrix R (one row per randomized function, one column per natural function of this system); the functions are kept in their current order.  Returns a NEW system.")
			.def("randomization_matrix",
				+[](SystemBaseT const& self) { return self.RandomizationMatrix(); },
				(arg("self")),
				"The randomization matrix R (codimension x N, complex_mp) of a system produced by randomize().  Raises if the system has no randomization block.")
			.def("reorder_functions_by_degree_decreasing", &SystemBaseT::ReorderFunctionsByDegreeDecreasing, (arg("self")),"Change the order of the functions to be in decreasing order")
			.def("reorder_functions_by_degree_increasing", &SystemBaseT::ReorderFunctionsByDegreeIncreasing, (arg("self")),"Change the order of the functions to be in decreasing order")
			.def("clear_variables", &SystemBaseT::ClearVariables, (arg("self")), "Remove the variable structure from the system")
			.def("copy_variable_structure", &SystemBaseT::CopyVariableStructure, (arg("self"), arg("other")), "Copy the variable structure from another System")
			
			.def("auto_patch",&SystemBaseT::AutoPatch, (arg("self")),"Apply a patch to the system, given its current variable group structure.")
			.def("copy_patches",&SystemBaseT::CopyPatches, (arg("self"), arg("other")),"Copy the patches from another system into this one.")
			.def("get_patch",&SystemBaseT::GetPatch, (arg("self")),"Get (a reference to) the patches from the system.")
			.def("is_patched",&SystemBaseT::IsPatched, (arg("self")),"Check whether the system is patched.")

			.def("rescale_point_to_fit_patch",&SystemBaseT::template RescalePointToFitPatch<complex_dbl>,(arg("self"), arg("point")),"Return a rescaled version of the input point, which fits the patch for the system.")
			.def("rescale_point_to_fit_patch",&SystemBaseT::template RescalePointToFitPatch<mpfr>,(arg("self"), arg("point")),"Return a rescaled version of the input point, which fits the patch for the system.")

			.def("rescale_point_to_fit_patch_in_place",&SystemBaseT::template RescalePointToFitPatchInPlace<complex_dbl>,(arg("self"), arg("point")),"Re-scale the input point, in place, to fit the patch for the system.  This assumes you have properly set the variable groups and auto-patched the system.")

			// .def("rescale_point_to_fit_patch_in_place",&SystemBaseT::template RescalePointToFitPatchInPlace<mpfr>,"Re-scale the input point, in place, to fit the patch for the system.  This assumes you have properly set the variable groups and auto-patched the system.")
			.def("rescale_point_to_fit_patch_in_place",&rescale_wrap_inplace_mpfr,(arg("self"), arg("point")),"Re-scale the input point, in place, to fit the patch for the system.  This assumes you have properly set the variable groups and auto-patched the system.")

			.def("dehomogenize_point",&SystemBaseT::template DehomogenizePoint<complex_dbl>,(arg("self"), arg("point")), "Dehomogenize a vector of doubles (complex), using the variable structure in this System")
			.def("dehomogenize_point",&SystemBaseT::template DehomogenizePoint<mpfr>,(arg("self"), arg("point")), "Dehomogenize a vector of mpfr's (complex), using the variable structure in this System")

			.def("homogenize_point",&SystemBaseT::template HomogenizePoint<complex_dbl>,(arg("self"), arg("point")), "Take a point in user (dehomogenized) coordinates to this system's internal coordinates: inserts the homogenizing coordinate for each affine variable group, then rescales onto the system's patch if patched.  Inverse of dehomogenize_point.")
			.def("homogenize_point",&SystemBaseT::template HomogenizePoint<mpfr>,(arg("self"), arg("point")), "Take a point in user (dehomogenized) coordinates to this system's internal coordinates: inserts the homogenizing coordinate for each affine variable group, then rescales onto the system's patch if patched.  Inverse of dehomogenize_point.")

			.def("variable_ordering",&SystemBaseT::VariableOrdering,(arg("self")), "The ordering of variables saying what each coordinate of a point in THIS system's coordinates means.  On your original system these are your variables; on a solver's target_system() the homogenizing variables appear too.")

			.def("describe",
				+[](SystemBaseT const& self, bool verbose) { std::ostringstream ss; self.Describe(ss, verbose); return ss.str(); },
				(arg("self"), arg("verbose") = false),
				"A human-facing description of the system, block by block (the same as str(system) when verbose=False).  verbose=True reveals the actual coefficients/matrices and the underlying functions of randomization / blend blocks.  For reading, not re-parsing.")
			.def(self_ns::str(self_ns::self))//, "String representation of the system (terse; structured blocks shown with placeholder symbols)
			.def(self_ns::repr(self_ns::self))//, "String representation of the system
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
			.def("num_start_points", &SystemBaseT::NumStartPoints,(arg("self")), "Get the number of start points that would be required by the system.  Non-negative, unsigned")
			.def("start_point_d", return_GenStart_ptr<complex_dbl>(),(arg("self"), arg("index")),"Get the k-th start point in double precision")
			.def("start_point_mp", return_GenStart_ptr<mpfr>(),(arg("self"), arg("index")),"Get the k-th start point in current multiple precision")
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

		// Pickle support for System, backed by the same Boost text-archive serialization that
		// bertini::Clone and the MPI system-broadcast use.  setstate applies the same
		// post-deserialize fixups as Clone: rebuild the SLP's derivatives (the archived SLP's
		// derivative outputs do not round-trip faithfully) and normalize precision across the tree.
		// This makes copy.copy / copy.deepcopy work and lets Systems cross process boundaries
		// (multiprocessing).  As with clone, the result's variables are distinct node objects from
		// the original's.
		struct SystemPickleSuite : boost::python::pickle_suite
		{
			static boost::python::tuple getinitargs(System const&)
			{
				return boost::python::make_tuple();
			}

			static boost::python::object getstate(System const& sys)
			{
				std::ostringstream oss;
				{
					boost::archive::text_oarchive oa(oss);
					oa << sys;
				}
				return boost::python::str(oss.str());
			}

			static void setstate(System& sys, boost::python::object state)
			{
				std::string s = boost::python::extract<std::string>(state)();
				std::istringstream iss(s);
				{
					boost::archive::text_iarchive ia(iss);
					ia >> sys;
				}
				sys.Differentiate();
				sys.precision(sys.precision());
			}
		};

		void ExportSystem()
		{
			
			// System class
		class_<System, std::shared_ptr<System> >("System", "The type in Bertini for systems of simultaneous equations.  Add functions and variable groups via member functions.", init<>())
			.def(init< std::vector<std::shared_ptr<node::Node>> >((arg("functions")), "Construct a System from a list of functions (bare expressions).  The variables are auto-discovered from the functions and placed into a single affine variable group, ordered alphabetically by name."))
			.def(SystemVisitor<System>())
			.def_pickle(SystemPickleSuite())
			;

			// free functions
			def("concatenate", &Concatenate,(arg("self"), arg("other")), "concatenate two Systems to produce a new one.  Appends the second's functions onto a copy of the first.  The two must share variable ordering (cloning one from the other, or just reusing the same variables, guarantees this -- variables are canonical by name).  If exactly one is patched, the result takes that patch.");
			def("clone", &Clone,(arg("self")), "Copy a System.  The copy shares the immutable node DAG (variables, functions, subexpressions) with the original but gets its own evaluation memory, so it is safe to evaluate concurrently AND its variables line up with the original's -- which is what lets you clone a set-up system, give the clone different functions, and concatenate the two (issue #256).  Adding/removing functions on one does not affect the other.  For a fully serialized deep copy use copy.deepcopy or pickle.");
			def("make_homotopy", &MakeHomotopy,
				(arg("target"), arg("start"), arg("path_variable")="t", arg("gamma")=std::shared_ptr<node::Node>()),
				"Form the gamma-trick straight-line homotopy H = (1-t)*target + gamma*t*start, with the path variable added.  At t=1 the homotopy is gamma*start (so start's solutions are its roots) and at t=0 it is target.  When start carries a structured block (e.g. a products-of-linears start system) the two systems are combined with a blend block; otherwise node arithmetic is used.  gamma=None generates a random rational gamma.  Pair with nag_algorithm.user_homotopy to solve.");
			def("make_moving_homotopy", &MakeMovingHomotopy,
				(arg("fixed"), arg("start_moving"), arg("end_moving"), arg("path_variable")="t", arg("gamma")=std::shared_ptr<node::Node>()),
				"Form a homotopy that moves ONLY the moving rows, leaving the fixed system evaluated once.  H = [ fixed's blocks ; (1-t)*end_moving + gamma*t*start_moving ]: the fixed equations (polynomial system + any static slices) stay as their own blocks (evaluated once, contributing zero to dH/dt) while only the moving rows slide.  At t=1 the moving rows are gamma*start_moving, at t=0 they are end_moving.  start_moving and end_moving hold just the moving rows and share fixed's variable structure.  The fixed rows come first, then the moving rows; build the matching target as fixed concatenated with end_moving.  gamma=None generates a random rational gamma.  Pair with nag_algorithm.user_homotopy to solve.");

			


			def("simplify", &call_simplify,(arg("self")), "Perform all possible simplifications.  Has side effects of modifying your functions, if held separately.  Shared nodes between multiple systems may have adverse effects");

			def("intern_system",
				+[](boost::python::object system_obj){
					// Extract the HOLDER's shared_ptr, not a converter temporary: boost.python's
					// from-python shared_ptr conversion mints an ephemeral control block that dies
					// at end of call, which would immediately expire the intern table's weak_ptr.
					// The holder's control block lives exactly as long as the Python object.
					std::shared_ptr<System>& held = boost::python::extract<std::shared_ptr<System>&>(system_obj)();
					// the representative is sealed, so its structural mutators raise; handing
					// Python a non-const pointer is safe in the same sense the C++ const is
					return std::const_pointer_cast<System>(InternSystem(held));
				},
				(arg("system")),
				"Hash-cons a System: return the live SEALED representative with an equal content_digest() if one exists, else seal and register this one.  Two equal systems interned in one session come back as the SAME object (s1 is s2).  The representative is sealed -- structural mutators raise; clone() it to get a mutable copy.  The intern table holds systems weakly, so it never keeps them alive on its own.");

			def("symbolic_jacobian",
				+[](boost::python::object functions, boost::python::object variables) {
					boost::python::stl_input_iterator<std::shared_ptr<node::Node>> fbegin(functions), fend;
					std::vector<std::shared_ptr<node::Node>> funcs(fbegin, fend);
					boost::python::stl_input_iterator<std::shared_ptr<node::Variable>> vbegin(variables), vend;
					bertini::VariableGroup vars(vbegin, vend);
					auto J = bertini::SymbolicJacobian(funcs, vars);
					boost::python::list rows;
					for (std::size_t i = 0; i < J.rows; ++i) {
						boost::python::list row;
						for (std::size_t j = 0; j < J.cols; ++j)
							row.append(J.entries[i*J.cols + j]);
						rows.append(row);
					}
					return rows;
				},
				(arg("functions"), arg("variables")),
				"The symbolic Jacobian of a list of functions with respect to a list of variables: J[i][j] = d functions[i] / d variables[j], as a list of rows of expression nodes.  Purely symbolic (no evaluation).  Prefer bertini.jacobian(...), which returns a 2-D numpy object array.");

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
				ExportRootsOfUnity();
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
			// The total-degree start system is now built from random linear products (generic
			// position); it has no per-variable "random value" (that was the roots-of-unity start,
			// now TotalDegreeBinomial).  It evaluates through a products-of-linears block.
			class_<start_system::TotalDegreeLinearProduct, bases<start_system::StartSystem>, std::shared_ptr<start_system::TotalDegreeLinearProduct> >("TotalDegreeLinearProduct",init<System const&>())
			;
		}

		void ExportRootsOfUnity()
		{
			// The roots-of-unity start system: x_i^{d_i} - r_i (formerly misnamed "TotalDegreeLinearProduct").
			// A structured, teaching/novelty start; kept reachable but not the default.
			class_<start_system::TotalDegreeBinomial, bases<start_system::StartSystem>, std::shared_ptr<start_system::TotalDegreeBinomial> >("TotalDegreeBinomial",init<System const&>())
			.def("random_value", &start_system::TotalDegreeBinomial::RandomValue<complex_dbl>,(arg("self"), arg("index")), "Get the k-th random value, in double precision")
			.def("random_value", &start_system::TotalDegreeBinomial::RandomValue<mpfr>,(arg("self"), arg("index")), "Get the k-th random value, in current multiple precision")
			.def("random_values", &start_system::TotalDegreeBinomial::RandomValues,(arg("self")), return_value_policy<copy_const_reference>(), "Get (a reference to) the random values for the start system, as Nodes")
			;
		}
	}
}
