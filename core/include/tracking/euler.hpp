//This file is part of Bertini 2.0.
//
//euler.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//euler.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with euler.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  euler.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#ifndef euler_hpp
#define euler_hpp

#include "tracking/tracking_config.hpp"

#include "system.hpp"
#include <eigen3/Eigen/LU>


namespace bertini{
	namespace tracking{
		namespace predict{

			/**
			Perform an euler-prediction step

			\param current_time The current value of the path variable.
			\param PType The mode for multiple precision.
			\param num_steps_since_last_condition_number_computation Obvious, hopefully.
			*/
			template <typename T>
			SuccessCode EulerStep(Vec<T> & next_space, T & next_time,
			               System const& S,
			               Vec<T> const& current_space, T current_time, 
			               T const& dt,
			               unsigned & num_steps_since_last_condition_number_computation, 
			               unsigned frequency_of_CN_estimation, PrecisionType PType)
			{
				auto dh_dt = -S.TimeDerivative(current_time);
				auto dh_dx = S.Jacobian(current_space, current_time); // this will complain if the system does not depend on time.


				if (PType==PrecisionType::Adaptive)
				{
					if (num_steps_since_last_condition_number_computation > frequency_of_CN_estimation)
					{
						// compute the condition number of the ...
						num_steps_since_last_condition_number_computation = 0; // reset the counter to 0
					}
					else // no need to compute the condition number
						num_steps_since_last_condition_number_computation++;
				}



				return SuccessCode::Success;
			}
		}
	}
}



#endif

