//This file is part of Bertini 2.
//
//src/system/precon.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/system/precon.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/system/precon.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame


#include "bertini2/system/precon.hpp"


namespace bertini{
	namespace system {

// has multiplicity 3 root at origin, and 3 infinite roots
System Precon::GriewankOsborn()
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



System Precon::CrossedPaths()
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

 
System Precon::Sphere()
{
	using Var = std::shared_ptr<node::Variable>;
	
	bertini::System sphere_sys;
	Var x = MakeVariable("x"), y = MakeVariable("y"), z = MakeVariable("z");
	auto one = MakeInteger(1);

	VariableGroup vars{x,y,z};
	sphere_sys.AddVariableGroup(vars);
	
	sphere_sys.AddFunction(pow(x,2) + pow(y,2) + pow(z,2) - one);

	return sphere_sys;
}


	} // namespace system
}//bertini