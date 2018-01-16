//This file is part of Bertini 2.
//
//python/generic_observable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/generic_observable.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/generic_observable.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Danielle Brake
//  UWEC
//  Fall 2017, Spring 2018
//
//
//  python/generic_observable.hpp:  source file for exposing trackers to python.

#pragma once
#include "python_common.hpp"
#include <bertini2/detail/observable.hpp>

namespace bertini{
	namespace python{

template <typename ObsT>
class ObservableVisitor : public def_visitor<ObservableVisitor<ObsT>>
{
	friend class def_visitor_access;

	static void AddObserver(object& obj, object& obs)
	{
		ObsT& self=extract<ObsT&>(obj)();
		AnyObserver& observer=extract<AnyObserver&>(obs)();
		self.AddObserver(observer);
	};

	static void RemoveObserver(object& obj, object& obs)
	{
		ObsT& self=extract<ObsT&>(obj)();
		AnyObserver& observer=extract<AnyObserver&>(obs)();
		self.RemoveObserver(observer);
	};

public:

	template<class PyClass>
	void visit(PyClass& cl) const{
		cl
		.def("add_observer", 	&ObservableVisitor::AddObserver)
		.def("remove_observer", &ObservableVisitor::RemoveObserver)
		;
	}
};


}} // namespaces
