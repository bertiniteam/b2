//This file is part of Bertini 2.
//
//python/tracker_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/tracker_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/tracker_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/tracker_export.cpp:  source file for exposing trackers to python.

#include "tracker.hpp"

namespace bertini{
	namespace python{

		template<typename TrackerT>
		template<class PyClass>
		void TrackerVisitor<TrackerT>::visit(PyClass& cl) const
		{
			cl
			.def("setup", &TrackerT::Setup)
			.def("track_path", &TrackerT::TrackPath)
			.def("get_system",&TrackerT::GetSystem,return_internal_reference<>())
			.def("predictor",get_predictor_,"Query the current predictor method used by the tracker.")
			.def("predictor",set_predictor_,"Set the predictor method used by the tracker.");
		}


		void ExportTrackers()
		{	
			using namespace bertini::tracking::config;

			enum_<Predictor>("Predictor")
	        .value("Constant", Predictor::Constant)
	        .value("Euler", Predictor::Euler)
	        .value("Heun", Predictor::Heun)
	        .value("RK4", Predictor::RK4)
	        .value("HeunEuler", Predictor::HeunEuler)
	        .value("RKNorsett34", Predictor::RKNorsett34)
	        .value("RKF45", Predictor::RKF45)
	        .value("RKCashKarp45", Predictor::RKCashKarp45)
	        .value("RKDormandPrince56", Predictor::RKDormandPrince56)
	        .value("RKVerner67", Predictor::RKVerner67)
	        .export_values()
	        ;


			ExportAMPTracker();
			ExportFixedTrackers();
		}

		void ExportAMPTracker()
		{
			class_<AMPTracker, std::shared_ptr<AMPTracker> >("AMPTracker", init<const System&>())
			.def(TrackerVisitor<AMPTracker>())
			;
		}

		void ExportFixedTrackers()
		{
			ExportFixedDoubleTracker();
			ExportFixedMultipleTracker();
		}

		void ExportFixedDoubleTracker()
		{

		}

		void ExportFixedMultipleTracker()
		{

		}



}} // namespaces

