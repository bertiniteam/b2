//This file is part of Bertini 2.0.
//
//correct.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//correct.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with correct.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  correct.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#ifndef correct_hpp
#define correct_hpp

#include "tracking/newton_correct.hpp"

namespace bertini{
	namespace tracking{
		

		template<typename NumType>
		SuccessCode Correct(Vec<NumType> & next_space,
					               System & S,
					               Vec<NumType> const& current_space, // pass by value to get a copy of it
					               NumType const& current_time, 
					               config::PrecisionType PrecType, 
					               NumType tracking_tolerance,
					               NumType path_truncation_threshold,
					               unsigned max_num_newton_iterations,
					               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			return correct::NewtonLoop(next_space,
						               S,
						               current_space, // pass by value to get a copy of it
						               current_time, 
						               PrecType, 
						               tracking_tolerance,
						               path_truncation_threshold,
						               max_num_newton_iterations,
						               AMP_config);
		}
	}
	
}


#endif


