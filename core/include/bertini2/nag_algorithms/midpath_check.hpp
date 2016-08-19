//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/midpath_check.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/midpath_check.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/midpath_check.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_algorithms/midpath_check.hpp 

\brief Provides methods for checking for path crossings.  

This essentially amounts to being certain that no two points are the same.
*/


#pragma once

namespace bertini{
	namespace algorithm{

struct Midpath
{
	struct Data{
		bool Passed()
		{
			return pass_;
		}


	private:
		bool pass_ = true;


	};

	template<typename T>
	static
	Data Check(T const&)
	{
		return Data();
	}
};

	}
}
