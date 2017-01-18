//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/sharpening.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/sharpening.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/sharpening.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame

/**
\file bertini2/nag_algorithms/sharpening.hpp 

\brief Provides the algorithms for sharpening points.  
*/


#pragma once

#include "bertini2/num_traits.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"
#include "bertini2/detail/visitable.hpp"
#include "bertini2/tracking.hpp"
#include "bertini2/nag_algorithms/midpath_check.hpp"
#include "bertini2/io/generators.hpp"


namespace bertini {

	namespace algorithm {

		namespace output {
			struct Sharpen
			{
				
			};
		} //re: namespace output

		template<class SystemType>
		struct Sharpen
		{

		};
		
	} // algo

} // bertini



