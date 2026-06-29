//This file is part of Bertini 2.
//
//python/zero_dim_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/zero_dim_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/zero_dim_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Wisconsin-Eau Claire
//  2023
//
//
//  python/zero_dim_export.hpp:  Header file for exposing the zero dim solve algorithm to Python


#ifndef BERTINI2_PYBERTINI_NAG_ALGORITHMS_ZERO_DIM
#define BERTINI2_PYBERTINI_NAG_ALGORITHMS_ZERO_DIM

#pragma once

#include "python_common.hpp"
#include "configured_visitor.hpp"
#include "generic_observable.hpp"
#include "generic_observer.hpp"

#include <bertini2/endgames.hpp>
#include <bertini2/nag_algorithms/zero_dim_solve.hpp>
#include <bertini2/nag_algorithms/events.hpp>
#include <bertini2/system/start_systems.hpp>
#include <bertini2/blackbox/config.hpp>   // blackbox::type::Start (start-system selector enum)
#include <bertini2/io/classic_writer.hpp>

#include <boost/python/stl_iterator.hpp>

#ifdef BERTINI2_HAVE_MPI
#include "bertini2/parallel/mpi_include.hpp"
#endif





namespace bertini{
	namespace python{
		
		using namespace bertini;




void ExportZeroDim();

// Registers AnyZeroDim + the nag observers submodule (CustomObserver + lifecycle
// events). Must run before the ZeroDim classes (which declare bases<AnyZeroDim>).
void ExportNagObservers();

// sub-functions (defined in zero_dim_{configs,double,mp,amp}_export.cpp)
void ExportZDConfigs();
void ExportZDDouble();
void ExportZDMP();
void ExportZDAMP();


// Template helpers — in header so zero_dim_configs_export.cpp can instantiate them.
template<typename NumT>
void ExposeSolutionMetaData(std::string const& class_name){
	using namespace bertini::algorithm;
	using MDT = SolutionMetaData<NumT>;
	class_<MDT>(class_name.c_str(),init<>())
	.def_readwrite("path_index",&MDT::path_index,
		"Index of the start path that produced this solution.")
	.def_readwrite("solution_index",&MDT::solution_index,
		"Index of this solution in the solution list.")
	.def_readwrite("precision_changed",&MDT::precision_changed,
		"Whether precision was increased while tracking this path (adaptive precision only).")
	.def_readwrite("time_of_first_prec_increase",&MDT::time_of_first_prec_increase,
		"The time value at which precision first increased on this path (adaptive precision only).")
	.def_readwrite("max_precision_used",&MDT::max_precision_used,
		"The highest precision (in digits) used while tracking this path (adaptive precision only).")
	.def_readwrite("pre_endgame_success",&MDT::pre_endgame_success,
		"The SuccessCode from tracking this path up to the endgame boundary. 0 means Success.")
	.def_readwrite("condition_number",&MDT::condition_number,
		"The latest estimate of the condition number (spectral norm) near the endpoint. Used, "
		"together with multiplicity, to classify the endpoint as singular.")
	.def_readwrite("newton_residual",&MDT::newton_residual,
		"The latest Newton step norm near the endpoint.")
	.def_readwrite("final_time_used",&MDT::final_time_used,
		"The final time value tracked to.")
	.def_readwrite("accuracy_estimate",&MDT::accuracy_estimate,
		"Accuracy estimate from the endgame, the difference between successive extrapolations.")
	.def_readwrite("accuracy_estimate_user_coords",&MDT::accuracy_estimate_user_coords,
		"Accuracy estimate in natural (dehomogenized) coordinates.")
	.def_readwrite("cycle_num",&MDT::cycle_num,
		"The cycle number used by the endgame's extrapolation.")
	.def_readwrite("endgame_success",&MDT::endgame_success,
		"The SuccessCode from the endgame. 0 means Success; anything else means the path did not "
		"converge to a finite solution (e.g. GoingToInfinity, SecurityMaxNormReached).")
	.def_readwrite("function_residual",&MDT::function_residual,
		"Infinity norm of the target system evaluated at the endpoint.")
	.def_readwrite("multiplicity",&MDT::multiplicity,
		"How many paths ended at this same point (1 for a simple solution). Computed by comparing "
		"dehomogenized endpoints with the infinity norm against final_tolerance * "
		"same_point_tolerance_multiplier.")
	.def_readwrite("multiplicity_representative",&MDT::multiplicity_representative,
		"For a multiplicity-m solution the solver returns m coincident endpoints; exactly one of "
		"them is the chosen representative (True) and the other m-1 are duplicates (False).  Use it "
		"to collapse a multiple solution to a single row -- which is what ZeroDim.to_dataframe() does "
		"by default (merge_multiplicities=True).  Simple solutions and at-infinity/failed endpoints "
		"are each their own representative (True).")
	.def_readwrite("is_real",&MDT::is_real,
		"Whether the (dehomogenized) endpoint is real, i.e. the infinity norm of its coordinates' "
		"imaginary parts is below PostProcessingConfig.real_threshold. Only meaningful for finite, "
		"successful endpoints.")
	.def_readwrite("is_finite",&MDT::is_finite,
		"Whether the endpoint is finite (not at infinity): the infinity norm of its dehomogenized "
		"coordinates is at most PostProcessingConfig.endpoint_finite_threshold. False also for paths "
		"the endgame flagged as diverging.")
	.def_readwrite("is_singular",&MDT::is_singular,
		"Whether the endpoint is singular: multiplicity > 1, or the condition-number estimate exceeds "
		"PostProcessingConfig.condition_number_threshold. Only meaningful for successful endpoints.")
	.def_readwrite("is_nonsolution",&MDT::is_nonsolution,
		"Whether the endpoint is a NONSOLUTION: a finite, successful point that is not a solution of "
		"the target system -- the extraneous nonsolutions introduced when ZeroDimSolver squares up an "
		"over-determined system.  Orthogonal to is_finite; excluded from the finite/real/singular "
		"solution accessors and surfaced by nonsolutions().  Load-bearing for regeneration cascades.")
	;
}

template<typename NumT>
void ExposeEndgameBoundaryMetaData(std::string const& class_name){
	using namespace bertini::algorithm;
	using MDT = EGBoundaryMetaData<NumT>;
	class_<MDT>(class_name.c_str(),init<>())
	.def_readwrite("path_point",&MDT::path_point)
	.def_readwrite("success_code",&MDT::success_code)
	.def_readwrite("last_used_stepsize",&MDT::last_used_stepsize)
	;
}


template<typename AlgoT>
class ZDVisitor: public def_visitor<ZDVisitor<AlgoT> >
{
	friend class ::boost::python::def_visitor_access;

	public:
		template<class PyClass>
		void visit(PyClass& cl) const;

	private:

		using MutableTrackerGetter = typename AlgoT::TrackerT& (AlgoT::*)();
		static MutableTrackerGetter GetTrackerMutable()
		{
			return &AlgoT::GetTracker;
		};

		using MutableEndgameGetter = typename AlgoT::EndgameT& (AlgoT::*)();
		static MutableEndgameGetter GetEndgameMutable()
		{
			return &AlgoT::GetEndgame;
		};
};


// Visitor body — template member function must be in header so each split TU can instantiate it.
template<typename AlgoT>
template<class PyClass>
void ZDVisitor<AlgoT>::visit(PyClass& cl) const
{
	cl
	.def(ConfiguredVisitor<AlgoT>())
	.def(ObservableVisitor<AlgoT>())
	.def("solve",
		+[](AlgoT& self, boost::python::object comm) -> void {
			// The GIL must be HELD while touching Python objects (e.g. extracting the mpi4py
			// communicator), and RELEASED only around the actual C++ solve so the worker threads
			// run truly in parallel on stock CPython (no free-threaded build needed).  Python
			// observer callbacks re-acquire the GIL in the observer trampoline (generic_observer.hpp).
			// No communicator => a purely LOCAL solve (serial or shared-memory threaded).  We call
			// Solve() directly, NOT Run(): Run() auto-promotes to MPI manager-worker whenever the
			// process was launched under mpirun (parallel::Size() > 1), which would make a plain
			// solve() secretly drive MPI_COMM_WORLD.  That breaks the common "MPI at an outer level,
			// independent local solves per rank" pattern (e.g. parameter_sweep across ranks).  MPI is
			// opt-in: you get it only by explicitly passing communicator=.
			if (comm.is_none()) {
				bertini::python::ScopedGILRelease unlock_gil;
				self.Solve();
			}
#ifdef BERTINI2_HAVE_MPI
			else {
				// Extract the communicator with the GIL held, then release it for RunParallel.
				MPI_Comm c = MPI_Comm_f2c(boost::python::extract<MPI_Fint>(comm.attr("py2f")()));
				bertini::python::ScopedGILRelease unlock_gil;
				self.RunParallel(c);
			}
#else
			else {
				PyErr_SetString(PyExc_RuntimeError,
					"solve(communicator=...) requires a bertini2 built with MPI; this build has none.");
				boost::python::throw_error_already_set();
			}
#endif
		},
		(boost::python::arg("communicator") = boost::python::object()),
		"Run the zero-dim algorithm. Pass an mpi4py communicator for parallel execution.")
	.def("get_tracker", GetTrackerMutable(), return_internal_reference<>(), "get a mutable reference to the Tracker being used")
	.def("get_endgame", GetEndgameMutable(), return_internal_reference<>(), "get a mutable reference to the Endgame being used")
	.def("all_solutions",
		+[](AlgoT const& self, bool user_coords) -> decltype(self.SolutionsInternalCoords()) {
			return user_coords ? self.SolutionsUserCoords() : self.SolutionsInternalCoords();
		},
		(boost::python::arg("self"), boost::python::arg("user_coords") = true),
		return_internal_reference<>(),
		"get ALL the computed solutions, one per tracked path (finite, at-infinity, and failed alike).  by default they are in the coordinates of YOUR variables (dehomogenized, depatched).  pass user_coords=False to decline, getting the solver's internal coordinates instead: homogenized, lying on the target system's patch -- the representation to use for continuing work.  the container is computed at most once per solve; repeated calls and indexing do not recompute it.  for the filtered view see solutions(); for the at-infinity ones see infinite_solutions.")
	.def("solutions",
		+[](AlgoT const& self, bool singular, bool real, bool nonreal, bool nonsingular,
		    bool infinite, bool nonsolution, bool user_coords){
			boost::python::list out;
			// Each category flag toggles inclusion of one kind of endpoint.  A genuine finite
			// solution is returned iff its singular-class AND its real-class are both enabled; the
			// at-infinity and nonsolution endpoints are opt-in.  Default = every finite genuine
			// solution (real + complex, simple + multiple), no at-infinity, no nonsolutions.
			auto pred = [=](auto const& m) -> bool {
				if (m.is_nonsolution) return nonsolution;   // finite, but not a solution of the target
				if (!m.is_finite)     return infinite;      // at infinity (or a failed path's placeholder)
				bool sing_ok = m.is_singular ? singular : nonsingular;
				bool real_ok = m.is_real    ? real     : nonreal;
				return sing_ok && real_ok;
			};
			for (auto const& p : self.SolutionsWhere(pred, user_coords)) out.append(p);
			return out;
		},
		(boost::python::arg("self"),
		 boost::python::arg("singular") = true, boost::python::arg("real") = true,
		 boost::python::arg("nonreal") = true, boost::python::arg("nonsingular") = true,
		 boost::python::arg("infinite") = false, boost::python::arg("nonsolution") = false,
		 boost::python::arg("user_coords") = true),
		"the solutions, filtered by category (returns points, not metadata).  By DEFAULT every finite "
		"genuine solution -- real and complex, simple and multiple -- and nothing else.  Each keyword "
		"toggles a category: singular / nonsingular select by conditioning, real / nonreal by realness "
		"(a finite solution is returned only if BOTH its conditioning class and its realness class are "
		"enabled), infinite=True also returns the at-infinity endpoints, and nonsolution=True also "
		"returns the nonsolutions.  E.g. solutions(real=False) -> complex finite solutions only; "
		"solutions(singular=False) -> nonsingular finite solutions; solutions(infinite=True) adds the "
		"divergent paths.  user_coords=False gives the solver's internal coordinates.  See also "
		"real_solutions / nonsingular_solutions / singular_solutions / infinite_solutions / nonsolutions "
		"for the common single-category views, and all_solutions for the raw per-path list.")
	.def("finite_solutions",
		+[](AlgoT const& self, bool user_coords){
			boost::python::list out;
			for (auto const& p : self.FiniteSolutions(user_coords)) out.append(p);
			return out;
		},
		(boost::python::arg("self"), boost::python::arg("user_coords") = true),
		"the FINITE solutions: successful endpoints the solver calls finite (is_finite applies endpoint_finite_threshold).  includes singular, nonsingular, and real solutions alike.  user coordinates by default; user_coords=False for internal coordinates.")
	.def("real_solutions",
		+[](AlgoT const& self, bool user_coords){
			boost::python::list out;
			for (auto const& p : self.RealSolutions(user_coords)) out.append(p);
			return out;
		},
		(boost::python::arg("self"), boost::python::arg("user_coords") = true),
		"the REAL finite solutions (is_real applies the configured tolerance).")
	.def("nonsingular_solutions",
		+[](AlgoT const& self, bool user_coords){
			boost::python::list out;
			for (auto const& p : self.NonsingularSolutions(user_coords)) out.append(p);
			return out;
		},
		(boost::python::arg("self"), boost::python::arg("user_coords") = true),
		"the NONSINGULAR finite solutions (simple, well-conditioned roots).")
	.def("singular_solutions",
		+[](AlgoT const& self, bool user_coords){
			boost::python::list out;
			for (auto const& p : self.SingularSolutions(user_coords)) out.append(p);
			return out;
		},
		(boost::python::arg("self"), boost::python::arg("user_coords") = true),
		"the SINGULAR finite solutions (multiple or ill-conditioned roots).")
	.def("nonsolutions",
		+[](AlgoT const& self, bool user_coords){
			boost::python::list out;
			for (auto const& p : self.Nonsolutions(user_coords)) out.append(p);
			return out;
		},
		(boost::python::arg("self"), boost::python::arg("user_coords") = true),
		"the NONSOLUTIONS: finite, successful endpoints that are NOT solutions of the target system "
		"-- the extraneous nonsolutions a squared-up over-determined system introduces (empty otherwise).  "
		"Excluded from finite_solutions/real/singular; a regeneration cascade reads these to discard them.")
	.def("infinite_solutions",
		+[](AlgoT const& self, bool user_coords){
			boost::python::list out;
			for (auto const& p : self.InfiniteSolutions(user_coords)) out.append(p);
			return out;
		},
		(boost::python::arg("self"), boost::python::arg("user_coords") = true),
		"the solutions AT INFINITY: endpoints not classified finite (is_finite is False).  these are the paths the endgame resolved as diverging -- its GoingToInfinity / SecurityMaxNormReached verdict, or a successful endpoint whose dehomogenized infinity norm exceeds endpoint_finite_threshold.  the complement of finite_solutions within all_solutions.  NOTE a path that FAILED before the endgame also has is_finite False (its point is not a real solution at infinity); inspect solution_metadata()/report() to tell a true divergence from a tracking failure.")
	.def("target_system",
		+[](AlgoT& self) -> decltype(self.TargetSystem()) { return self.TargetSystem(); },
		return_internal_reference<>(),
		"get the prepared target system: the homogenized, auto-patched clone of the system you supplied.  its patch is the one internal-coordinate solutions lie on; use its dehomogenize_point/homogenize_point/variable_ordering to move between representations.")
	.def("solution_metadata", &AlgoT::FinalSolutionMetadata, return_internal_reference<>(), "get the metadata for the solutions at the target time")
	.def("endgame_boundary_solutions", &AlgoT::EndgameBoundarySolutions, return_internal_reference<>(), "get the solutions (per-path point data) at the endgame boundary, where regular tracking switches to the endgame")
	.def("endgame_boundary_metadata", &AlgoT::EndgameBoundaryMetadata, return_internal_reference<>(), "get the MidpathCheckReport from the path-crossing check at the endgame boundary: how many crossings were detected, which paths, how many re-track attempts were made, and whether the check ultimately passed")
	.def("report", &AlgoT::Report, "a concise end-of-solve diagnostic summary (a SolveReport): how every path ended up -- finite solutions, diverged, or FAILED (by named reason) -- plus singular/real counts, max condition number, the path-crossing outcome, and all_paths_resolved.  print(solver.report()) for a human-readable summary; a count alone can hide a path the tracker silently lost.")
	.def("to_classic_input",
		+[](AlgoT const& self, System const& sys) -> std::string {
			using namespace bertini::tracking;
			classic::ClassicWriteOptions opt;

			// mptype: how THIS solver tracks precision -- 2 adaptive, 0 fixed-double, 1 fixed-multiple.
			using TrackerTraitsT = TrackerTraits<typename AlgoT::TrackerT>;
			if (TrackerTraitsT::IsAdaptivePrec)
				opt.mptype = 2;
			else if (std::is_same<typename TrackerTraitsT::BaseComplexT, complex_dbl>::value)
				opt.mptype = 0;
			else
				opt.mptype = 1;

			// predictor + step-size cadence + Newton: read off the live tracker.
			auto const& tracker  = self.GetTracker();
			auto const& stepping = tracker.template Get<SteppingConfig>();
			auto const& newton   = tracker.template Get<NewtonConfig>();
			opt.odepredictor      = classic::PredictorToClassic(tracker.GetPredictor());
			opt.maxstepsize       = static_cast<double>(stepping.max_step_size);
			opt.stepsuccessfactor = static_cast<double>(stepping.step_size_success_factor);
			opt.stepfailfactor    = static_cast<double>(stepping.step_size_fail_factor);
			opt.stepsforincrease  = stepping.consecutive_successful_steps_before_stepsize_increase;
			opt.maxnumbersteps    = stepping.max_num_steps;
			opt.maxnewtonits      = newton.max_num_newton_iterations;

			// tolerances + crossed-path resolve cap: read off the algorithm's own configs.
			auto const& tol = self.template Get<algorithm::TolerancesConfig>();
			opt.tracktolbeforeeg = static_cast<double>(tol.newton_before_endgame);
			opt.tracktolduringeg = static_cast<double>(tol.newton_during_endgame);
			opt.finaltol         = static_cast<double>(tol.final_tolerance);
			opt.maxcrossedpathresolves =
				self.template Get<typename AlgoT::ZeroDimConf>().max_num_crossed_path_resolve_attempts;

			return classic::SystemToClassicFile(sys, opt);
		},
		(boost::python::arg("self"), boost::python::arg("system")),
		"emit a complete Bertini 1 classic input file (CONFIG + INPUT) for the given natural system, "
		"using THIS solver's tracking settings for the CONFIG: precision mode (mptype), ODE predictor, "
		"tolerances (before/during endgame, final), the full step-size cadence, max Newton iterations, "
		"and the crossed-path resolve cap.  The system you pass supplies INPUT -- pass the natural "
		"(un-homogenized) system you constructed the solver from, since the solver homogenizes and "
		"patches its internal copy.  Use this to re-run the exact same problem with the exact same knobs "
		"in Bertini 1 for cross-validation (the random start system aside).  For sweeping settings without "
		"a solver, see ``system.to_classic_input(**kwargs)``.")
	;
}


// Map a Start enum choice to a clone factory -- the one place a concrete start-system type is named
// on the binding side.  User is rejected: it needs the user_homotopy path (RefToGiven).
inline start_system::StartSystemFactory<System> StartFactoryFor(blackbox::type::Start which)
{
	switch (which)
	{
		case blackbox::type::Start::MHom:
			return start_system::MakeStartFactory<start_system::MHomogeneous>();
		case blackbox::type::Start::RootsOfUnity:
			return start_system::MakeStartFactory<start_system::RootsOfUnity>();
		case blackbox::type::Start::User:
			throw std::runtime_error("the User start system is not a clone-owned start; build the "
			                         "homotopy with nag_algorithm.user_homotopy(...) instead");
		case blackbox::type::Start::TotalDegree:
		default:
			return start_system::MakeStartFactory<start_system::TotalDegree>();
	}
}

// Expose the blackbox start-system selector enum to Python (call once).
inline void ExportStartSystemEnum()
{
	boost::python::enum_<blackbox::type::Start>("StartSystemType",
		"Which start system a ZeroDim solver builds: total_degree (default) or mhomogeneous. "
		"User homotopies use nag_algorithm.user_homotopy(...) instead.")
		.value("total_degree", blackbox::type::Start::TotalDegree)
		.value("roots_of_unity", blackbox::type::Start::RootsOfUnity)
		.value("mhomogeneous", blackbox::type::Start::MHom)
		;
}

// Register the opaque start-system factory type + a producer function (call once).  ZeroDim no
// longer bakes the start system into its type; a non-default start is selected by passing one of
// these factories as the solver constructor's second argument.
inline void ExportStartSystemFactory()
{
	class_<start_system::StartSystemFactory<System>>("StartSystemFactory",
		"Opaque factory that builds a start system for a ZeroDim solver.  Get one from "
		"start_system_factory(StartSystemType.X) and pass it as the solver's second constructor arg.",
		no_init);
	def("start_system_factory", &StartFactoryFor,
		(boost::python::arg("which")),
		"start_system_factory(which): the start-system factory for a StartSystemType, to pass as the "
		"second argument of a ZeroDim solver constructor (the default constructor uses total_degree).");
}

// Helper template — defined here so all split TUs can use it without duplication.  ZeroDim is no
// longer templated on the start system; the concrete start system is chosen at construction (default
// TotalDegree, or an explicit factory) and held polymorphically.
template<typename TrackerT, typename EndgameT>
void ExportZeroDimSpecific(std::string const& class_name){
	using ZeroDimT = algorithm::ZeroDimSolver<TrackerT, EndgameT, System>;
	// with_custodian_and_ward<1,2>: keep the Python System (ctor arg, index 2) alive as long as the
	// solver (self, index 1).  The solver's internal system is a SHALLOW copy that shares the
	// function-tree nodes' shared_ptrs with the user's System; those shared_ptrs carry boost.python
	// deleters tied to the Python node wrappers, so if the user passes a temporary System and lets it
	// die, the endgame's Differentiate() drops a node and the deleter touches an already-freed Python
	// object -> crash.  Tying the System's lifetime to the solver makes the documented `ZeroDim(sys)`
	// usage safe even when the caller keeps no reference to `sys`.
	//
	// Two constructors: (system) defaults to the TotalDegree start; (system, start_factory) selects
	// a non-default start (e.g. MHomogeneous) -- the factory is a value, copied into the solver.
	class_<ZeroDimT, std::shared_ptr<ZeroDimT>, bases<algorithm::AnyZeroDim> >(class_name.c_str(),
		init<System const&>()[with_custodian_and_ward<1, 2>()])
	.def(init<System const&, start_system::StartSystemFactory<System> const&>(
			(boost::python::arg("system"), boost::python::arg("start_factory")))
			[with_custodian_and_ward<1, 2>()])
	.def(ZDVisitor<ZeroDimT>())
	// ZeroDimSolver-specific (the engine HomotopySolver has no system to square up):
	.def("was_randomized", &ZeroDimT::WasRandomized,
		"True if the supplied system was over-determined and was squared up by randomization "
		"(so the extraneous solutions the squaring introduces have been filtered out of "
		"finite_solutions).  False for a square system.")
	.def("randomization_matrix", &ZeroDimT::RandomizationMatrix, return_value_policy<copy_const_reference>(),
		"The exact n x N coefficient matrix used to square up an over-determined system "
		"(empty if the system was already square).")
	;
}


// --- user-homotopy ZeroDim: run a homotopy YOU built from a list of start points YOU have ---
//
// This is the SAME ZeroDim template (same Solve / pre-endgame / midpath / endgame /
// post-processing), instantiated with start_system::User (start points come from the supplied
// list) and a user-supplied homotopy (taken as-is, not formed).  Nothing about the
// solve loop is re-implemented here -- this only registers the class + a constructor.

// Build a start_system::User from a target system + a Python list of start-point vectors
// (each element an complex_mp vector via eigenpy).  The returned User references `target`, and
// RefToGiven references all three systems, so the friendly Python wrapper keeps them all alive.
inline std::shared_ptr<start_system::User>
MakeUserStartSystem(System const& target, boost::python::list const& start_points)
{
	SampCont<complex_mp> solns{
		boost::python::stl_input_iterator<Vec<complex_mp>>(start_points),
		boost::python::stl_input_iterator<Vec<complex_mp>>() };
	return std::make_shared<start_system::User>(target, solns);
}

// Register the User start system once (it is not tracker-specific).
inline void ExportUserStartSystem(){
	class_<start_system::User, std::shared_ptr<start_system::User>, boost::noncopyable>(
		"UserStartSystem",
		"A start system that is simply a list of start points you already have (e.g. solutions "
		"from an earlier solve), to be tracked through a homotopy you constructed.  Built for you "
		"by nag_algorithm.user_homotopy(...).",
		no_init)
		.def("__init__", make_constructor(&MakeUserStartSystem),
			"UserStartSystem(target_system, start_points): start_points is a list of vectors.")
		.def("num_start_points", &start_system::User::NumStartPoints)
		;
}

// Bind one ZeroDim<...,User,RefToGiven> variant.  Ctor takes (target, start, homotopy) by
// reference (RefToGiven); the Python wrapper retains all three so the references stay valid.
template<typename TrackerT, typename EndgameT>
void ExportZeroDimUserHomotopy(std::string const& class_name){
	using ZeroDimT = algorithm::HomotopySolver<TrackerT, EndgameT, System>;
	// RefToGiven: this solver holds REFERENCES to all three given systems, so each must outlive the
	// solver.  Chain with_custodian_and_ward to keep target (arg 2), start (arg 3) and homotopy
	// (arg 4) alive as long as the solver (self, arg 1).
	class_<ZeroDimT, std::shared_ptr<ZeroDimT>, bases<algorithm::AnyZeroDim> >(class_name.c_str(),
		init<System const&, start_system::User const&, System const&>(
			(boost::python::arg("target"), boost::python::arg("start"), boost::python::arg("homotopy")))
			[with_custodian_and_ward<1, 2, with_custodian_and_ward<1, 3, with_custodian_and_ward<1, 4> > >()])
	.def(ZDVisitor<ZeroDimT>())
	;
}






}} // namespaces




#endif // the include guards
