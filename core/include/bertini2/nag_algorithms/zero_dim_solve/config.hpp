//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/zero_dim_solve/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/zero_dim_solve/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/zero_dim_solve/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_algorithms/zero_dim_solve/config.hpp 

\brief Provides the config for the zero dimensional algorithm.
*/

#pragma once

#include <string>

namespace bertini{
	namespace algorithm{
		namespace config{



template<typename ComplexT>
struct ZeroDim
{
	unsigned initial_ambient_precision = DoublePrecision();
	unsigned max_num_crossed_path_resolve_attempts = 2; ///< The maximum number of times to attempt to re-solve crossed paths at the endgame boundary.

	ComplexT start_time = ComplexT(1);
	ComplexT endgame_boundary = ComplexT("0.1");
	ComplexT target_time = ComplexT(0);

	std::string path_variable_name = "ZERO_DIM_PATH_VARIABLE";
};


} } } // namespaces