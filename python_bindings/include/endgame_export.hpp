//This file is part of Bertini 2.
//
//python/endgame_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/endgame_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/endgame_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Notre Dame
//  Summer 2016, Summer 2023, Fall 2024
//
//
//  python/endgame_export.hpp:  Header file for exposing endgames to python.

#pragma once

#include "python_common.hpp"
#include "generic_observable.hpp"

#include <bertini2/endgames.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace bertini{
	namespace python{

		using namespace bertini::tracking;

		/**
		 Abstract Endgame class
		 */
		template<typename EndgameT>
		class EndgameBaseVisitor: public def_visitor<EndgameBaseVisitor<EndgameT> >
		{
			friend class ::boost::python::def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:

			using BCT = typename TrackerTraits<typename EndgameT::TrackerType>::BaseComplexT;
			using BRT = typename TrackerTraits<typename EndgameT::TrackerType>::BaseRealT;
			using BaseEGT = typename EndgameT::BaseEGT;

			static
			SuccessCode WrapRun(EndgameT & self, Vec<BCT> const& s){
				return self.Run(s);
			}


			using unsigned_of_void = unsigned (BaseEGT::*)() const;
			static unsigned_of_void GetCycleNumberFn()
			{
				return &BaseEGT::CycleNumber;
			};

			template <typename T>
			static
			Vec<T> return_final_approximation(EndgameT const& self)
			{
				return self.template FinalApproximation<T>();
			}

		};// EndgameVisitor class



		/**
		 Particulars for the PowerSeries endgame.
		 */
		template<typename PowerSeriesT>
		class PowerSeriesVisitor: public def_visitor<PowerSeriesVisitor<PowerSeriesT> >
		{
			friend class ::boost::python::def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:

			using BCT = typename TrackerTraits<typename PowerSeriesT::TrackerType>::BaseComplexT;
			using BRT = typename TrackerTraits<typename PowerSeriesT::TrackerType>::BaseRealT;


		};// CauchyVisitor class




		/**
		 Particulars for the Cauchy endgame.
		 */
		template<typename CauchyT>
		class CauchyVisitor: public def_visitor<CauchyVisitor<CauchyT> >
		{
			friend class ::boost::python::def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:

			using BCT = typename TrackerTraits<typename CauchyT::TrackerType>::BaseComplexT;
			using BRT = typename TrackerTraits<typename CauchyT::TrackerType>::BaseRealT;


		};// CauchyVisitor class




		// Visitor body definitions — template member functions must be in the header
		// so each split TU can instantiate them for its own concrete type.

		template<typename EndgameT>
		template<class PyClass>
		void EndgameBaseVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;
			using BCT = typename TrackerTraits<TrackerT>::BaseComplexT;

			cl
			.def("cycle_number", this->GetCycleNumberFn(),arg("self"),"Get the cycle number as currently computed")

			.def("get_endgame_settings",&EndgameT::EndgameSettings,return_internal_reference<>(),arg("self"),"Get the current non-specific endgame settings")
			.def("get_security_settings",&EndgameT::SecuritySettings,return_internal_reference<>(),arg("self"),"Get the 'security' settings for the endgame (path truncation near infinity)")

			.def("set_endgame_settings",&EndgameT::template Set<endgame::EndgameConfig>,(arg("self"),arg("settings")),"Set the values of non-specific endgame settings")
			.def("set_security_settings",&EndgameT::template Set<endgame::SecurityConfig>,(arg("self"),arg("settings")),"Set the values of security-level settings")

			.def("get_tracker", &EndgameT::GetTracker, return_internal_reference<>(),arg("self"),"Get the tracker used in this endgame.  This is the same tracker as you feed the endgame object when you make it.  This is a reference variable")
			.def("get_system",  &EndgameT::GetSystem,  return_internal_reference<>(),arg("self"),"Get the tracked system.  This is a reference to the internal system.")

			.def("final_approximation", &return_final_approximation<BCT>,arg("self"),"Get the current approximation of the root, in the ambient numeric type for the tracker being used")

			.def("run", &EndgameBaseVisitor::WrapRun,
				 (arg("self"), "start_point"),
				 "Run the endgame from the stored boundary time to the stored target time. "
				 "Call set_boundary_time() before running.")

			.def("set_boundary_time", &EndgameT::SetBoundaryTime,
				 (arg("self"), arg("t")),
				 "Set the time at which the endgame begins (authoritative, stored at the given precision).")
			.def("set_target_time", &EndgameT::SetTargetTime,
				 (arg("self"), arg("t")),
				 "Set the time the endgame tracks toward (default 0).")
			.def("boundary_time", &EndgameT::BoundaryTime,
				 return_value_policy<copy_const_reference>(), arg("self"),
				 "Get the stored endgame boundary time.")
			.def("target_time", &EndgameT::TargetTime,
				 return_value_policy<copy_const_reference>(), arg("self"),
				 "Get the stored target time.")

			.def(ObservableVisitor<EndgameT>())
			;
		}


		template<typename EndgameT>
		template<class PyClass>
		void CauchyVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;

			cl
			.def(init<TrackerT const&, endgame::CauchyConfig const&>((arg("self"),arg("tracker"),arg("cauchyconfig"))))
			.def(init<TrackerT const&, endgame::EndgameConfig const&>((arg("self"),arg("tracker"),arg("endgameconfig"))))
			.def(init<TrackerT const&, endgame::SecurityConfig const&>((arg("self"),arg("tracker"),arg("securityconfig"))));
		}


		template<typename EndgameT>
		template<class PyClass>
		void PowerSeriesVisitor<EndgameT>::visit(PyClass& cl) const
		{
			using TrackerT = typename EndgameT::TrackerType;

			cl
			.def(init<TrackerT const&, endgame::PowerSeriesConfig const&>((arg("self"),arg("tracker"),arg("powerseriesconfig"))))
			.def(init<TrackerT const&, endgame::EndgameConfig const&>((arg("self"),arg("tracker"),arg("endgameconfig"))))
			.def(init<TrackerT const&, endgame::SecurityConfig const&>((arg("self"),arg("tracker"),arg("securityconfig"))));
		}


		// Prototypes for functions defined in the split .cpp files.

		void ExportEndgames();
		void ExportEndgameSettings();

		// per-tracker-family export functions (defined in endgame_{double,mp,amp}_export.cpp)
		void ExportAMPPSEG();
		void ExportFDPSEG();
		void ExportFMPSEG();

		void ExportAMPCauchyEG();
		void ExportFDCauchyEG();
		void ExportFMCauchyEG();


}}// re: namespaces

