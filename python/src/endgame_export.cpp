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
//  Daniel Brake
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
		void EndgameVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;
			using BCT = typename TrackerTraits<TrackerT>::BaseComplexType;

			cl
			.def("cycle_number", get_cycle_number_)

			.def("get_endgame_settings",&EndgameT::EndgameSettings,return_internal_reference<>())
			.def("get_security_settings",&EndgameT::SecuritySettings,return_internal_reference<>())
			.def("get_tolerances",&EndgameT::Tolerances, return_internal_reference<>())

			.def("set_endgame_settings",&EndgameT::SetEndgameSettings)
			.def("set_security_settings",&EndgameT::SetSecuritySettings)
			.def("set_tolerances",&EndgameT::SetToleranceSettings,"Set the tracking tolerances")

			.def("get_tracker", &EndgameT::GetTracker, return_internal_reference<>(),"Get the tracker used in this endgame.  This is the same tracker as you feed the endgame object when you make it.")
			.def("get_system",  &EndgameT::GetSystem,  return_internal_reference<>(),"Get the tracked system")

			.def("final_approximation", &EndgameT::template FinalApproximation<BCT>, return_internal_reference<>(),"Get the current approximation of the root")
			.def("run", &EndgameT::template Run<BCT>,"Run the endgame, from start point and start time, to t=0")
			;
		}


		template<typename EndgameT>
		template<class PyClass>
		void CauchyVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;

			cl
			.def(init<TrackerT const&, config::Cauchy<RT> const&>())
			.def(init<TrackerT const&, config::Endgame<RT> const&>())
			.def(init<TrackerT const&, config::Security<RT> const&>())
			.def(init<TrackerT const&, config::Tolerances<RT> const&>());

		}

 
		template<typename EndgameT>
		template<class PyClass>
		void PowerSeriesVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;

			cl
			.def(init<TrackerT const&, config::PowerSeries const&>())
			.def(init<TrackerT const&, config::Endgame<RT> const&>())
			.def(init<TrackerT const&, config::Security<RT> const&>())
			.def(init<TrackerT const&, config::Tolerances<RT> const&>());
		}




		void ExportEndgameSettings()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".config");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("config") = new_submodule;

			scope new_submodule_scope = new_submodule;

			class_<config::Security<double>>("Security_d",init<>())
			.def(SecurityVisitor<double>());

			class_<config::Security<mpfr_float>>("Security_mp",init<>())
			.def(SecurityVisitor<mpfr_float>());

			class_<config::PowerSeries>("PowerSeriesConfig",init<>())
			.def(PowerSeriesConfigVisitor());


			class_<config::Cauchy<double>>("CauchyConfig_d",init<>())
			.def(CauchyConfigVisitor<double>());

			class_<config::Cauchy<mpfr_float>>("CauchyConfig_mp",init<>())
			.def(CauchyConfigVisitor<mpfr_float>());

			
		}

		void ExportEndgames()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".endgame");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("endgame") = new_submodule;

			scope new_submodule_scope = new_submodule;

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
			using EGT = typename EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("AMPPSEG",init<TrackerT const&>())
			.def(EndgameVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}


		void ExportFMPSEG()
		{
			using TrackerT = MultiplePrecisionTracker;
			using EGT = typename EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("FixedMultiplePSEG",init<TrackerT const&>())
			.def(EndgameVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}


		void ExportFDPSEG()
		{
			using TrackerT = DoublePrecisionTracker;
			using EGT = typename EndgameSelector<TrackerT>::PSEG;

			class_<EGT>("FixedDoublePSEG",init<TrackerT const&>())
			.def(EndgameVisitor<EGT>())
			.def(PowerSeriesVisitor<EGT>());
		}





		void ExportFDCauchyEG()
		{
			using TrackerT = DoublePrecisionTracker;
			using EGT = typename EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("FDCauchyEG", init<TrackerT const&>())
			.def(EndgameVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}

		void ExportFMCauchyEG()
		{
			using TrackerT = MultiplePrecisionTracker;
			using EGT = typename EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("FMCauchyEG", init<TrackerT const&>())
			.def(EndgameVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}


		void ExportAMPCauchyEG()
		{
			using TrackerT = AMPTracker;
			using EGT = typename EndgameSelector<TrackerT>::Cauchy;

			class_<EGT>("AMPCauchyEG", init<TrackerT const&>())
			.def(EndgameVisitor<EGT>())
			.def(CauchyVisitor<EGT>())
			;
		}

}} // re: namespaces
