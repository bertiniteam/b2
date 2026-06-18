//This file is part of Bertini 2.
//
//python/generic_observers.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/generic_observers.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/generic_observers.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/generic_observers.cpp:  source file for exposing trackers to python.


#include "generic_observer.hpp"
#include <bertini2/detail/observable.hpp>
#include <boost/python/exception_translator.hpp>

namespace bertini{
	namespace python{


void ExportObserver()
{
	// Attaching an observer to an observable it cannot observe raises TypeError.
	register_exception_translator<bertini::IncompatibleObserver>(
		[](bertini::IncompatibleObserver const& e){ PyErr_SetString(PyExc_TypeError, e.what()); });

	class_<AnyEvent, boost::noncopyable>("AnyEvent", no_init);

	enum_<ObserveResult>("ObserveResult",
		"What an observer may return from Observe(): KeepObserving (the default if "
		"you return None) or Unsubscribe to ask the observable to drop this observer.")
		.value("KeepObserving", ObserveResult::KeepObserving)
		.value("Unsubscribe",   ObserveResult::Unsubscribe)
	;

	// shared_ptr holder: when a python observer is attached, the observable can
	// co-own it (weak_ptr), so dropping the python reference doesn't dangle.
	class_<ObserverWrapper<AnyObserver>, std::shared_ptr<ObserverWrapper<AnyObserver>>, boost::noncopyable>("AnyAbstractObserver",  init< >())
	;
}

}} // namespaces
