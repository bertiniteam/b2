//This file is part of Bertini 2.
//
//bertini2/systems/precon.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/systems/precon.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/systems/precon.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/systems/precon.hpp 

\brief Provides preconstructed systems for use around the home and office.
*/


#pragma once


#include "bertini2/systems/system.hpp"



namespace bertini{
	namespace system{

		struct Precon{

/**
Creates the Griewank-Osborn system.

From the returned system, you can unpack the variables, etc.
*/
static
System GriewankOsborn()
{
	using Var = std::shared_ptr<node::Variable>;

	bertini::System griewank_osborn_sys;
	Var x = MakeVariable("x"), t = MakeVariable("t"), y = MakeVariable("y");
	VariableGroup vars{x,y};
	griewank_osborn_sys.AddVariableGroup(vars); 

	griewank_osborn_sys.AddFunction((mpq_rational(29,16))*pow(x,3)-2*x*y);
	griewank_osborn_sys.AddFunction((y - pow(x,2)));

	return griewank_osborn_sys;
}
			
			
			/**
			 Creates system with crossed paths at t = 1/3 and t = 2/3.
			 
			 From the returned system, you can unpack the variables, etc.
			 */
			static
			System CrossedPaths()
			{
				using Var = std::shared_ptr<node::Variable>;
				
				bertini::System crossed_paths_sys;
				Var x = MakeVariable("x"), t = MakeVariable("t"), y = MakeVariable("y");
				auto two = MakeInteger(2);
				auto half = MakeRational("1/2");
				VariableGroup vars{x,y};
				crossed_paths_sys.AddVariableGroup(vars);
				
				crossed_paths_sys.AddFunction(pow(x,3)+ two);
				crossed_paths_sys.AddFunction(pow(y,2) + half);
				
				return crossed_paths_sys;
			}


};
		
}} // namespaces

