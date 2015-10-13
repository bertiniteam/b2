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

#ifndef BERTINI_CORRECT_HPP
#define BERTINI_CORRECT_HPP

#include "tracking/newton_correct.hpp"

namespace bertini{
	namespace tracking{
			

		template <typename ComplexType, typename RealType>
		SuccessCode Correct(Vec<ComplexType> & next_space,
					               System & S,
					               Vec<ComplexType> const& current_space, 
					               ComplexType const& current_time, 
					               RealType tracking_tolerance,
					               RealType path_truncation_threshold,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations,
					               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			return correct::NewtonLoop(next_space,
						               S,
						               current_space, 
						               current_time, 
						               tracking_tolerance,
						               path_truncation_threshold,
						               min_num_newton_iterations,
						               max_num_newton_iterations,
						               AMP_config);
		}



		template <typename ComplexType, typename RealType>
		SuccessCode Correct(Vec<ComplexType> & next_space,
		                    RealType & norm_delta_z,
			                       RealType & norm_J,
			                       RealType & norm_J_inverse,
			                       RealType & condition_number_estimate,
					               System & S,
					               Vec<ComplexType> const& current_space, 
					               ComplexType const& current_time, 
					               RealType tracking_tolerance,
					               RealType path_truncation_threshold,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations,
					               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			return correct::NewtonLoop(next_space,
			                           norm_delta_z,
			                           norm_J,
			                           norm_J_inverse,
			                           condition_number_estimate,
						               S,
						               current_space, 
						               current_time, 
						               tracking_tolerance,
						               path_truncation_threshold,
						               min_num_newton_iterations,
						               max_num_newton_iterations,
						               AMP_config);
		}


	}
	
}


#endif


