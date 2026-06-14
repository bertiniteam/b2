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

#include <bertini2/endgames.hpp>
#include <bertini2/nag_algorithms/zero_dim_solve.hpp>
#include <bertini2/system/start_systems.hpp>

#include <boost/python/stl_iterator.hpp>

#ifdef BERTINI2_HAVE_MPI
#include <mpi.h>
#endif





namespace bertini{
	namespace python{
		
		using namespace bertini;




void ExportZeroDim();

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
	.def_readwrite("path_index",&MDT::path_index)
	.def_readwrite("solution_index",&MDT::solution_index)
	.def_readwrite("precision_changed",&MDT::precision_changed)
	.def_readwrite("time_of_first_prec_increase",&MDT::time_of_first_prec_increase)
	.def_readwrite("max_precision_used",&MDT::max_precision_used)
	.def_readwrite("pre_endgame_success",&MDT::pre_endgame_success)
	.def_readwrite("condition_number",&MDT::condition_number)
	.def_readwrite("newton_residual",&MDT::newton_residual)
	.def_readwrite("final_time_used",&MDT::final_time_used)
	.def_readwrite("accuracy_estimate",&MDT::accuracy_estimate)
	.def_readwrite("accuracy_estimate_user_coords",&MDT::accuracy_estimate_user_coords)
	.def_readwrite("cycle_num",&MDT::cycle_num)
	.def_readwrite("endgame_success",&MDT::endgame_success, "this is a SuccessCode.  0 means Success.  Anything other than 0 means something happened.")
	.def_readwrite("function_residual",&MDT::function_residual)
	.def_readwrite("multiplicity",&MDT::multiplicity)
	.def_readwrite("is_real",&MDT::is_real)
	.def_readwrite("is_finite",&MDT::is_finite)
	.def_readwrite("is_singular",&MDT::is_singular)
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
	.def("solve",
		+[](AlgoT& self, boost::python::object comm) -> void {
#ifdef BERTINI2_HAVE_MPI
			if (comm.is_none()) {
				self.Run();
			} else {
				MPI_Fint f = boost::python::extract<MPI_Fint>(comm.attr("py2f")());
				self.RunParallel(MPI_Comm_f2c(f));
			}
#else
			self.Solve();
#endif
		},
		(boost::python::arg("communicator") = boost::python::object()),
		"Run the zero-dim algorithm. Pass an mpi4py communicator for parallel execution.")
	.def("get_tracker", GetTrackerMutable(), return_internal_reference<>(), "get a mutable reference to the Tracker being used")
	.def("get_endgame", GetEndgameMutable(), return_internal_reference<>(), "get a mutable reference to the Endgame being used")
	.def("solutions",
		+[](AlgoT const& self, bool user_coords) -> decltype(self.SolutionsInternalCoords()) {
			return user_coords ? self.SolutionsUserCoords() : self.SolutionsInternalCoords();
		},
		(boost::python::arg("self"), boost::python::arg("user_coords") = true),
		return_internal_reference<>(),
		"get the computed solutions.  by default they are in the coordinates of YOUR variables (dehomogenized, depatched).  pass user_coords=False to decline, getting the solver's internal coordinates instead: homogenized, lying on the target system's patch -- the representation to use for continuing work.  the container is computed at most once per solve; repeated calls and indexing do not recompute it.")
	.def("target_system",
		+[](AlgoT& self) -> decltype(self.TargetSystem()) { return self.TargetSystem(); },
		return_internal_reference<>(),
		"get the prepared target system: the homogenized, auto-patched clone of the system you supplied.  its patch is the one internal-coordinate solutions lie on; use its dehomogenize_point/homogenize_point/variable_ordering to move between representations.")
	.def("solution_metadata", &AlgoT::FinalSolutionMetadata, return_internal_reference<>(), "get the metadata for the solutions at the target time")
	.def("endgame_boundary_data", &AlgoT::EndgameBoundaryData, return_internal_reference<>(), "get the data for the state at the endgame boundary (when we switch from regular tracking to endgame tracking")
	;
}


// Helper template — defined here so all split TUs can use it without duplication.
template<typename TrackerT, typename EndgameT, typename SystemT, typename StartSystemT>
void ExportZeroDimSpecific(std::string const& class_name){
	using ZeroDimT = algorithm::ZeroDim<TrackerT, EndgameT, SystemT, StartSystemT>;
	class_<ZeroDimT, std::shared_ptr<ZeroDimT> >(class_name.c_str(), init<SystemT>())
	.def(ZDVisitor<ZeroDimT>())
	;
}


// --- user-homotopy ZeroDim: run a homotopy YOU built from a list of start points YOU have ---
//
// This is the SAME ZeroDim template (same Solve / pre-endgame / midpath / endgame /
// post-processing), instantiated with start_system::User (start points come from the supplied
// list) and policy::RefToGiven (the homotopy is taken as-is, not formed).  Nothing about the
// solve loop is re-implemented here -- this only registers the class + a constructor.

// Build a start_system::User from a target system + a Python list of start-point vectors
// (each element an mpfr_complex vector via eigenpy).  The returned User references `target`, and
// RefToGiven references all three systems, so the friendly Python wrapper keeps them all alive.
inline std::shared_ptr<start_system::User>
MakeUserStartSystem(System const& target, boost::python::list const& start_points)
{
	SampCont<mpfr_complex> solns{
		boost::python::stl_input_iterator<Vec<mpfr_complex>>(start_points),
		boost::python::stl_input_iterator<Vec<mpfr_complex>>() };
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
	using ZeroDimT = algorithm::ZeroDim<TrackerT, EndgameT, System, start_system::User, policy::RefToGiven>;
	class_<ZeroDimT, std::shared_ptr<ZeroDimT> >(class_name.c_str(),
		init<System const&, start_system::User const&, System const&>(
			(boost::python::arg("target"), boost::python::arg("start"), boost::python::arg("homotopy"))))
	.def(ZDVisitor<ZeroDimT>())
	;
}






}} // namespaces




#endif // the include guards
