//This file is part of Bertini 2.
//
//python/endgame_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/endgame_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/endgame_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Dani Brake
//  University of Notre Dame, University of Wisconsin Eau Claire
//  Summer 2016, Fall 2017
//
//
//  python/endgame_export.cpp:  source file for exposing endgames to python.

#include "endgame_export.hpp"

namespace bertini{
	namespace python{


		template<typename EndgameT>
		template<class PyClass>
		void EndgameBaseVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;
			using BCT = typename TrackerTraits<TrackerT>::BaseComplexType;
			using BRT = typename TrackerTraits<TrackerT>::BaseRealType;

			cl
			.def("cycle_number", this->GetCycleNumberFn(),"Get the cycle number as currently computed")

			.def("get_endgame_settings",&EndgameT::EndgameSettings,return_internal_reference<>(),"Get the current non-specific endgame settings")
			.def("get_security_settings",&EndgameT::SecuritySettings,return_internal_reference<>(),"Get the 'security' settings for the endgame (path truncation near infinity)")

			.def("set_endgame_settings",&EndgameT::template Set<endgame::EndgameConfig>,"Set the values of non-specific endgame settings")
			.def("set_security_settings",&EndgameT::template Set<endgame::SecurityConfig>,"Set the values of security-level settings")

			.def("get_tracker", &EndgameT::GetTracker, return_internal_reference<>(),"Get the tracker used in this endgame.  This is the same tracker as you feed the endgame object when you make it.  This is a reference variable")
			.def("get_system",  &EndgameT::GetSystem,  return_internal_reference<>(),"Get the tracked system.  This is a reference to the internal system.")

			.def("final_approximation", &EndgameT::template FinalApproximation<BCT>, return_internal_reference<>(),"Get the current approximation of the root, in the ambient numeric type for the tracker being used")
			.def("run", RunDefaultTime<BCT>(),"Run the endgame, from start point and start time, to t=0.  Expects complex numeric type matching that of the tracker being used.")
			.def("run", RunCustomTime<BCT>(),"Run the endgame, from start point and start time, to your choice of target time t.  Expects complex numeric type matching that of the tracker being used.")
			;
		}


		template<typename EndgameT>
		template<class PyClass>
		void CauchyVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;

			cl
			.def(init<TrackerT const&, endgame::CauchyConfig const&>())
			.def(init<TrackerT const&, endgame::EndgameConfig const&>())
			.def(init<TrackerT const&, endgame::SecurityConfig const&>());

		}

 
		template<typename EndgameT>
		template<class PyClass>
		void PowerSeriesVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;

			cl
			.def(init<TrackerT const&, endgame::PowerSeriesConfig const&>())
			.def(init<TrackerT const&, endgame::EndgameConfig const&>())
			.def(init<TrackerT const&, endgame::SecurityConfig const&>());
		}




		void ExportEndgameSettings()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".config");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("config") = new_submodule;
			

			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "Endgame configuration structs.";

			class_<endgame::EndgameConfig>("Endgame","Generic endgame settings.  Number of sample points, etc.  Note that some of its configs are rational numbers",init<>())
				.def_readwrite("sample_point_refinement_factor", &endgame::EndgameConfig::sample_point_refinement_factor, "Extra amount of tolerance for refining before computing the final approximation, during endgame.")
				.def_readwrite("num_sample_points", &endgame::EndgameConfig::num_sample_points,"The number of points to use for extrapolant calculation.  In the Power Series Endgame, the is the number of geometrically spaces points on the path.  For Cauchy, this is the number of points on each circle tracked around the target time value.")
				.def_readwrite("min_track_time", &endgame::EndgameConfig::min_track_time,"The minimum distance from the target time to track to.  Decreasing this may help failing runs succeed, or maybe not, because you are, after all, tracking toward a singularity.")
				.def_readwrite("sample_factor", &endgame::EndgameConfig::sample_factor,"The factor by which to space the geometrically spaced `distance' between sample points, or sample circles for Cauchy.")
				.def_readwrite("max_num_newton_iterations", &endgame::EndgameConfig::max_num_newton_iterations,"the maximum number of newton iterations to be taken during sample point sharpening.  Increasing this can help speed convergence, at the risk of path jumping.")
				.def_readwrite("final_tolerance", &endgame::EndgameConfig::final_tolerance, "The tolerance to which to track the path, using the endgame.  Endgames require two consecutive estimates to be this close to each other under the relative infinity norm.  Default value is 1e-11.")
				;


			class_<endgame::SecurityConfig>("Security","Security settings for endgames.  Control things like truncation because estimated root is near infinity",init<>())
				.def_readwrite("level", &endgame::SecurityConfig::level,"Turns on or off truncation of paths going to infinity during the endgame.  0 is off, 1 is on.")
				.def_readwrite("max_norm", &endgame::SecurityConfig::max_norm,"If on, the norm at which to truncate a path.")
				;

			class_<endgame::PowerSeriesConfig>("PowerSeriesConfig","Settings specific to the power series endgame for computing singular endpoints",init<>())
				.def_readwrite("max_cycle_number", &endgame::PowerSeriesConfig::max_cycle_number,"The maximum cycle number to consider, when calculating the cycle number which best fits the path being tracked.")
				.def_readwrite("cycle_number_amplification", &endgame::PowerSeriesConfig::cycle_number_amplification,"The maximum number allowable iterations during endgames, for points used to approximate the final solution.")
				;


			class_<endgame::CauchyConfig>("CauchyConfig","Settings specific to the Cauchy endgame for computing singular endpoints",init<>())
				.def_readwrite("cycle_cutoff_time", &endgame::CauchyConfig::cycle_cutoff_time)
				.def_readwrite("ratio_cutoff_time", &endgame::CauchyConfig::ratio_cutoff_time)
				.def_readwrite("minimum_for_c_over_k_stabilization", &endgame::CauchyConfig::minimum_for_c_over_k_stabilization)
				.def_readwrite("maximum_cauchy_ratio", &endgame::CauchyConfig::maximum_cauchy_ratio)
				.def_readwrite("num_needed_for_stabilization", &endgame::CauchyConfig::num_needed_for_stabilization,"When running stabilization testing for the cycle number when entering the endgame, this is the number of consecutive points for which the test must pass.")
				.def_readwrite("fail_safe_maximum_cycle_number", &endgame::CauchyConfig::fail_safe_maximum_cycle_number, "max number of loops before giving up." )
				;
			
		}

		void ExportEndgames()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".endgame");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("endgame") = new_submodule;
			

			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "Endgames and associated types and functions.  For tracking around singularities.";

			ExportEndgameSettings();

			ExportAMPPSEG();
			ExportFDPSEG();
			ExportFMPSEG();

			ExportAMPCauchyEG();
			ExportFDCauchyEG();
			ExportFMCauchyEG();
		}

		void ExportAMPPSEG()
		{
			using TrackerT = AMPTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("AMPPSEG",
				"The adaptive precision implementation of the power series endgame.",
				init<TrackerT const&>())
			.def(EndgameBaseVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}


		void ExportFMPSEG()
		{
			using TrackerT = MultiplePrecisionTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("FixedMultiplePSEG",
				"The fixed but arbitrary precision implementation of the power series endgame",
				init<TrackerT const&>())
			.def(EndgameBaseVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}


		void ExportFDPSEG()
		{
			using TrackerT = DoublePrecisionTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("FixedDoublePSEG",
				"The double-precision implementation of the power series endgame",
				init<TrackerT const&>())
			.def(EndgameBaseVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}





		void ExportFDCauchyEG()
		{
			using TrackerT = DoublePrecisionTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("FixedDoubleCauchyEG",
				"The fixed double precision implementation of the Cauchy endgame",
				init<TrackerT const&>())
			.def(EndgameBaseVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}

		void ExportFMCauchyEG()
		{
			using TrackerT = MultiplePrecisionTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("FixedMultipleCauchyEG",
				"The fixed multiple precision implementation of the Cauchy endgame",
				init<TrackerT const&>())
			.def(EndgameBaseVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}


		void ExportAMPCauchyEG()
		{
			using TrackerT = AMPTracker;
			using EGT = typename endgame::EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("AMPCauchyEG",
				"The adaptive precision implementation of the Cauchy endgame",
				init<TrackerT const&>())
			.def(EndgameBaseVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}

}} // re: namespaces
