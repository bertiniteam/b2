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
// Copyright(C) 2017-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Dani Brake
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

	void Observe(AnyEvent const& e) { this->get_override("Observe")(e);}
	
}; // re: ObserverWrapper

void ExportObserver();

}} // namespaces
