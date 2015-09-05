//This file is part of Bertini 2.0.
//
//tracking/step.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracking/step.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/step.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  tracking/step.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015


#ifndef BERTINI_STEP_HPP
#define BERTINI_STEP_HPP

#include "tracking/predict.hpp"
#include "tracking/correct.hpp"

namespace bertini {

	namespace tracking {

		/**
	

		\tparam NumType The type of number being used in the algorithm

		*/
		template<typename NumType>
		SuccessCode Step(config::PredictorChoice predictor_choice,
		                    Vec<T> & next_space, T & next_time,
				               System & sys,
				               Vec<T> const& current_space, T current_time, 
				               T const& dt,
				               T & condition_number_estimate,
				               unsigned & num_steps_since_last_condition_number_computation, 
				               unsigned frequency_of_CN_estimation, PrecisionType prec_type, 
				               T const& tracking_tolerance,
				               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{

			SuccessCode predictor_code = Predict(next_space, next_time,
							               		sys,
							               		current_space, current_time, 
							               		dt,
							               		condition_number_estimate,
							               		num_steps_since_last_condition_number_computation, 
							               		frequency_of_CN_estimation, prec_type, 
							               		tracking_tolerance,
							               		AMP_config);

			if (predictor_code!=SuccessCode::Success)
				return predictor_code;


			SuccessCode corrector_code = Correct(next_space,
								               S,
								               current_space, // pass by value to get a copy of it
								               current_time, 
								               PrecType, 
								               tracking_tolerance,
								               path_truncation_threshold,
								               max_num_newton_iterations,
								               AMP_config);
			if (corrector_code!=SuccessCode::Success)
				return corrector_code;

			return SuccessCode::Success;
		}

	} // re: namespace tracking

} // re: namespace bertini





#endif

