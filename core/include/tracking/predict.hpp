//This file is part of Bertini 2.0.
//
//predict.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//predict.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with predict.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  predict.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#ifndef predict_hpp
#define predict_hpp

#include "tracking/ode_predictors.hpp"

namespace bertini{
	namespace tracking{

		/**
		Wrapper class for calling an ODE predictor.

		\param predictor_choice The enum class selecting the predictor to be used.
		\param next_space The computed prediction.
		\param next_time The next time.
		\param sys The system being solved.
		\param current_space The current space variable vector.
		\param current_time The current time.
		\param dt The size of the time step.
		\param condition_number_estimate The computed estimate of the condition number of the Jacobian.
		\param num_steps_since_last_condition_number_computation.  Updated in this function.
		\param frequency_of_CN_estimation How many steps to take between condition number estimates.
		\param prec_type The operating precision type.  
		\param tracking_tolerance How tightly to track the path.
		\param AMP_config The settings for adaptive multiple precision.

		\tparam T The number type for evaluation.
		*/
		template <typename T>
		SuccessCode Predict(config::Predictor predictor_choice,
		                    Vec<T> & next_space, T & next_time,
				               System & sys,
				               Vec<T> const& current_space, T current_time, 
				               T const& dt,
				               T & condition_number_estimate,
				               unsigned & num_steps_since_last_condition_number_computation, 
				               unsigned frequency_of_CN_estimation, config::PrecisionType prec_type, 
				               T const& tracking_tolerance,
				               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{

			switch (predictor_choice)
			{
				case config::Predictor::Euler:
				{
					return predict::Euler(next_space, next_time,
				               		sys,
				               		current_space, current_time, 
				               		dt,
				               		condition_number_estimate,
				               		num_steps_since_last_condition_number_computation, 
				               		frequency_of_CN_estimation, prec_type, 
				               		tracking_tolerance,
				               		AMP_config);
					break;
				}

				default:
				{
					throw std::runtime_error("unknown predictor choice in Predict");
				}
			}

			return SuccessCode::Failure;
		}

	

	}

	
}



#endif
