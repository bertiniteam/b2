//This file is part of Bertini 2.
//
//python/tracker_observers.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/tracker_observers.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/tracker_observers.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Dani Brake
//  UWEC
//  Fall 2017
//
//
//  python/tracker_observers.cpp:  source file for exposing trackers to python.


#include "tracker_observers.hpp"


namespace bertini{
	namespace python{



void ExportTrackerObservers()
{

	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".tracking");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("tracking") = new_submodule;

	scope new_submodule_scope = new_submodule;
	new_submodule_scope.attr("__doc__") = "Tracking things for PyBertini.  Includes the three fundamental trackers, and utility functions.";

	{
	// this should be called in the tracking module namespace
		scope scope_B;
		std::string submodule_name_B(extract<const char*>(scope_B.attr("__name__")));
		submodule_name_B.append(".observers");
		object submodule_B(borrowed(PyImport_AddModule(submodule_name_B.c_str())));
		scope_B.attr("observers") = submodule_B;
		scope new_submodule_scope_B = submodule_B;

		{
			scope scope_C;
			std::string submodule_name_C(extract<const char*>(scope_C.attr("__name__")));
			submodule_name_C.append(".amp");
			object submodule_C(borrowed(PyImport_AddModule(submodule_name_C.c_str())));
			scope_C.attr("amp") = submodule_C;
			scope new_submodule_scope_C = submodule_C;

			class_<FirstPrecisionRecorder<AMPTracker>>("FirstPrecisionRecorder", init< >())
			.def(TrackingObserverVisitor<FirstPrecisionRecorder<AMPTracker>>())
			;


			class_<GoryDetailLogger<AMPTracker>>("GoryDetailLogger", init< >())
			.def(TrackingObserverVisitor<GoryDetailLogger<AMPTracker>>())
			;
		}
		
	}
}

}} // namespaces
