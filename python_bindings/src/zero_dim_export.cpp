//This file is part of Bertini 2.
//
//python/zero_dim_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/zero_dim_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/zero_dim_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/zero_dim_export.cpp:  Coordinator — sets up the nag_algorithms submodule and
//  calls per-tracker-family registration functions (defined in split files).  It also
//  registers AnyZeroDim + the nag observers (the one TU here that needs the heavy headers).

#include "zero_dim_export.hpp"

namespace bertini{
	namespace python{

		// Forward declarations — defined in zero_dim_{configs,double,mp,amp}_export.cpp.
		void ExportZDConfigs();
		void ExportZDDouble();
		void ExportZDMP();
		void ExportZDAMP();

		// Registers AnyZeroDim (the lifecycle-event emitter base + RTTI/bases anchor) and the
		// nag observers submodule: a single CustomObserver + event set serving every ZeroDim
		// variant.  Must run before the ZeroDim classes, which declare bases<AnyZeroDim>.
		void ExportNagObservers(){
			using namespace bertini::algorithm;
			using AZ = AnyZeroDim;

			// the abstract emitter base, in the nag_algorithms scope
			class_<AZ, boost::noncopyable>("AnyZeroDim", no_init);

			// observers submodule
			scope current_scope;
			std::string obs_name(extract<const char*>(current_scope.attr("__name__")));
			obs_name.append(".observers");
			object obs_module(borrowed(PyImport_AddModule(obs_name.c_str())));
			current_scope.attr("observers") = obs_module;
			scope obs_scope = obs_module;
			obs_scope.attr("__doc__") = "Observers for the zero-dim solve.  Subclass CustomObserver "
				"(or use it directly), then attach with the solver's add_observer.";

			// the subclassable observer base; shared_ptr holder so an attached observer is
			// co-owned (it can outlive the python reference without dangling)
			class_<ObserverWrapper<Observer<AZ>>, std::shared_ptr<ObserverWrapper<Observer<AZ>>>,
			       bases<AnyObserver>, boost::noncopyable>("CustomObserver", init< >())
				;

			// .solver() returns the emitting solver; RTTI resolves it to the concrete ZeroDim
			// (registered with bases<AnyZeroDim>), so python sees the full solver API.
			auto solver_getter = +[](const AlgorithmEvent<AZ>& e) -> const AZ& { return e.Get(); };
			class_<AlgorithmEvent<AZ>, bases<AnyEvent>, boost::noncopyable>("AlgorithmEvent", no_init)
				.def("solver", solver_getter, return_value_policy<reference_existing_object>());

			class_<AlgorithmStarted<AZ>,  bases<AlgorithmEvent<AZ>>, boost::noncopyable>("AlgorithmStarted",  no_init);
			class_<AlgorithmComplete<AZ>, bases<AlgorithmEvent<AZ>>, boost::noncopyable>("AlgorithmComplete", no_init);

			class_<PathBeginning<AZ>, bases<AlgorithmEvent<AZ>>, boost::noncopyable>("PathBeginning", no_init)
				.def("path_index", &PathBeginning<AZ>::PathIndex, "index of the solution path that is starting");
			class_<PathComplete<AZ>, bases<AlgorithmEvent<AZ>>, boost::noncopyable>("PathComplete", no_init)
				.def("path_index", &PathComplete<AZ>::PathIndex, "index of the solution path that finished");
		}

		void ExportZeroDim(){
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".nag_algorithms");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("nag_algorithms") = new_submodule;

			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "Algorithms for computing things, like point solutions to square systems (zerodim algorithm).";

			ExportNagObservers();    // registers AnyZeroDim before the ZeroDim classes below

			ExportZDConfigs();
			ExportZDDouble();
			ExportZDMP();
			ExportZDAMP();
		}

	}
}
