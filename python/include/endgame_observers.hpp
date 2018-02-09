//This file is part of Bertini 2.
//
//python/endgame_observers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/endgame_observers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/endgame_observers.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/endgame_observers.hpp:  source file for exposing endgame observers to python.

#pragma once
#include "python_common.hpp"
#include "endgame_export.hpp"

#include <bertini2/endgames/observers.hpp>


namespace bertini{
	namespace python{


void ExportEndgameObservers();


using namespace bertini::endgame;


template<typename ObsT>
struct EndgameObserverVisitor: public def_visitor<EndgameObserverVisitor<ObsT> >
{
	friend class def_visitor_access;

public:

	template<class PyClass>
	void visit(PyClass& cl) const{
	}
};



}} // namespaces
