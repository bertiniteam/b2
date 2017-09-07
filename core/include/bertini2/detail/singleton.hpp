//This file is part of Bertini 2.
//
//bertini2/detail/singleton.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/detail/singleton.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/detail/singleton.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/detail/singleton.hpp 

\brief Provides types and utilities for dealing with global defaults for the detail routines
*/


#pragma once

namespace bertini{

namespace detail{


/**
\brief A virtual singleton base class which can be used in a heirarchy to require a class to be singleton.

\note
This class is virtual, so make sure you mark destructors as virtual.
*/
struct Singleton {


	virtual ~Singleton() = default;
private:
	// mark both the default and copy constructor as private
	Singleton(){}
	Singleton(Singleton const&){}

};


} // namespace detail
} // namespace bertini






