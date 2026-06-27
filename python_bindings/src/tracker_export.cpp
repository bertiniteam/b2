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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Notre Dame
//  Summer 2016, Spring 2018
//
//
//  python/tracker_export.cpp:  Tracker class registrations.
//  Visitor method bodies live in tracker_export.hpp (template definitions).
//  Config/enum registration lives in tracker_config_export.cpp.

#include "tracker_export.hpp"

namespace bertini{
	namespace python{

		// ExportConfigSettings is defined in tracker_config_export.cpp.
		void ExportConfigSettings();

		void ExportTrackers()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".tracking");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("tracking") = new_submodule;

			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "Tracking things.  Includes the three fundamental trackers, and utility functions.";

			// Register the common observable base so that a base-class pointer (e.g. the one
			// returned by a path event's tracker()) RTTI-downcasts to the concrete tracker in
			// Python -- the same mechanism AnyZeroDim uses for event.solver().  No methods are
			// exposed on it; add_observer lives on each concrete tracker via ObservableVisitor.
			// Must be registered before the trackers that declare bases<Observable>.
			class_<bertini::Observable, boost::noncopyable>("Observable", no_init);

			ExportConfigSettings();
			ExportAMPTracker();
			ExportFixedTrackers();
		}

		void ExportAMPTracker()
		{
			class_<AMPTracker, std::shared_ptr<AMPTracker>, bases<bertini::Observable> >("AMPTracker", "The adaptive multiple precision (AMP) tracker.  Ambient numeric type is multiple-precision (mpfr_complex).  Contruct one by feeding it a system -- cannot be constructed without feeding it a system.  Adjust its settings via configs and the `setup` function.  Then, call method `track_path`.", init<const System&>())
			.def(TrackerVisitor<AMPTracker>())
			.def(AMPTrackerVisitor<AMPTracker>())
			;
		}

		void ExportFixedTrackers()
		{
			ExportFixedDoubleTracker();
			ExportFixedMultipleTracker();
		}

		void ExportFixedDoubleTracker()
		{
			class_<DoublePrecisionTracker, std::shared_ptr<DoublePrecisionTracker>, bases<bertini::Observable> >("DoublePrecisionTracker", "The double precision tracker.  Tracks using only complex doubles.  Ambient numeric type is double.  Contruct one by feeding it a system -- cannot be constructed without feeding it a system.  Adjust its settings via configs and the `setup` function.  Then, call method `track_path`.", init<const System&>())
			.def(TrackerVisitor<DoublePrecisionTracker>())
			.def(FixedDoubleTrackerVisitor<DoublePrecisionTracker>())
			;
		}

		void ExportFixedMultipleTracker()
		{
			class_<MultiplePrecisionTracker, std::shared_ptr<MultiplePrecisionTracker>, bases<bertini::Observable> >("MultiplePrecisionTracker", "The fixed multiple precision tracker.  Ambient numeric type is multiple-precision (mpfr_complex).  Precision is the value of bertini.default_precision() at contruction.  Errors if you try to feed it things not at that precision.  Contruct one by feeding it a system -- cannot be constructed without feeding it a system.  Adjust its settings via configs and the `setup` function.  Then, call method `track_path`.", init<const System&>())
			.def(TrackerVisitor<MultiplePrecisionTracker>())
			.def(FixedMultipleTrackerVisitor<MultiplePrecisionTracker>())
			;
		}

}} // namespaces
