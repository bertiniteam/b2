//This file is part of Bertini 2.
//
//python/generic_observer.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/generic_observer.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/generic_observer.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  UWEC
//  Fall 2017, Spring 2018
//
//
//  python/generic_observer.hpp:  source file for exposing trackers to python.

#pragma once
#include "python_common.hpp"
#include <bertini2/detail/observer.hpp>

namespace bertini{
	namespace python{

// Wrapper struct to allow derived classes to overide methods in python
template<typename ObsT>
struct ObserverWrapper : ObsT, wrapper<ObsT>
{

	// Use boost::ref so Boost.Python's to_python_indirect path is taken,
	// which uses RTTI to find the most-derived registered event type and
	// enables isinstance() checks in Python.  The const_cast is safe: events
	// are short-lived temporaries and Python only reads them during the call.
	//
	// Lifetime contract (see AnyEvent in detail/events.hpp): the python `event`
	// object handed to Observe() is only valid for the duration of the call.
	// A python observer must read what it needs and copy the values out (e.g.
	// into numpy/lists) before returning; it must NOT stash the event object,
	// `event.tracker()`, or anything they return for use after Observe() ends.
	//
	// Return-value translation: python observers conventionally return None
	// (their Observe() just does side effects), which we map to KeepObserving so
	// existing observers keep working.  An observer that wants to self-detach can
	// `return bertini.ObserveResult.Unsubscribe`.
	ObserveResult Observe(AnyEvent const& e) override {
		// A threaded solve releases the GIL and fires events from C++ worker threads; re-acquire
		// the GIL before touching any Python object.  Safe on the main thread too (serial solve).
		ScopedGILAcquire acquire_gil;
		object result = this->get_override("Observe")(boost::ref(const_cast<AnyEvent&>(e)));
		if (result.is_none())
			return ObserveResult::KeepObserving;
		extract<ObserveResult> as_result(result);
		if (as_result.check())
			return as_result();
		return ObserveResult::KeepObserving;
	}
	
}; // re: ObserverWrapper

void ExportObserver();

}} // namespaces
