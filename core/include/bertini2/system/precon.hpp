//This file is part of Bertini 2.
//
//bertini2/system/precon.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/system/precon.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/system/precon.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/system/precon.hpp 

\brief Provides preconstructed systems for use around the home and office.
*/


#pragma once


#include "bertini2/system/system.hpp"



namespace bertini{
	namespace system{

		struct Precon{

/**
Creates the Griewank-Osborn system.

From the returned system, you can unpack the variables, etc.
*/
static
System GriewankOsborn();
			
			
/**
 Creates system with crossed paths at t = 1/3 and t = 2/3.
 
 From the returned system, you can unpack the variables, etc.
 */
static
System CrossedPaths();



/**
\brief Constructs a system describing the 2-sphere, in variables x, y, z.

System has one affine variable group: x, y, z;
One function: x^2 + y^2 + z^2 - 1;
*/
static 
System Sphere();


};
		
}} // namespaces

