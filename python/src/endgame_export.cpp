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
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Dani Brake
//  University of Notre Dame
//  Summer 2016
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
			.def("cycle_number", this->GetCycleNumberFn())

			.def("get_endgame_settings",&EndgameT::EndgameSettings,return_internal_reference<>())
			.def("get_security_settings",&EndgameT::SecuritySettings,return_internal_reference<>())

			.def("set_endgame_settings",&EndgameT::template Set<endgame::EndgameConfig>)
			.def("set_security_settings",&EndgameT::template Set<endgame::SecurityConfig>)

			.def("get_tracker", &EndgameT::GetTracker, return_internal_reference<>(),"Get the tracker used in this endgame.  This is the same tracker as you feed the endgame object when you make it.")
			.def("get_system",  &EndgameT::GetSystem,  return_internal_reference<>(),"Get the tracked system")

			.def("final_approximation", &EndgameT::template FinalApproximation<BCT>, return_internal_reference<>(),"Get the current approximation of the root")
			.def("run", RunDefaultTime<BCT>(),"Run the endgame, from start point and start time, to t=0")
			.def("run", RunCustomTime<BCT>(),"Run the endgame, from start point and start time, to your choice of target time t")
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

			
			class_<endgame::SecurityConfig>("Security","Security settings for endgames.  Control things like truncation because estimated root is near infinity",init<>())
			.def(SecurityVisitor());

			class_<endgame::PowerSeriesConfig>("PowerSeriesConfig","Settings specific to the power series endgame for computing singular endpoints",init<>())
			.def(PowerSeriesConfigVisitor());


			class_<endgame::CauchyConfig>("CauchyConfig","Settings specific to the Cauchy endgame for computing singular endpoints",init<>())
			.def(CauchyConfigVisitor());

			
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
