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
	

		\tparam ComplexType The type of number being used in the algorithm

		*/
		template <typename ComplexType, typename RealType>
		SuccessCode Step(config::Predictor predictor_choice,
		                    Vec<ComplexType> & next_space, ComplexType & next_time,
				               System & sys,
				               Vec<ComplexType> const& current_space, ComplexType current_time, 
				               ComplexType const& delta_t,
				               RealType & condition_number_estimate,
				               unsigned & num_steps_since_last_condition_number_computation, 
				               unsigned frequency_of_CN_estimation, config::PrecisionType prec_type, 
				               RealType const& tracking_tolerance,
				               RealType const& path_truncation_threshold,
				               unsigned min_num_newton_iterations,
				               unsigned max_num_newton_iterations,
				               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{

			SuccessCode predictor_code = Predict(next_space,
							               		sys,
							               		current_space, current_time, 
							               		delta_t,
							               		condition_number_estimate,
							               		num_steps_since_last_condition_number_computation, 
							               		frequency_of_CN_estimation, prec_type, 
							               		tracking_tolerance,
							               		AMP_config);

			if (predictor_code!=SuccessCode::Success)
				return predictor_code;

			next_time = current_time + delta_t;

			SuccessCode corrector_code = Correct(next_space,
								               sys,
								               current_space, // pass by value to get a copy of it
								               next_time, 
								               prec_type, 
								               tracking_tolerance,
								               path_truncation_threshold,
								               min_num_newton_iterations,
								               max_num_newton_iterations,
								               AMP_config);


			if (corrector_code!=SuccessCode::Success)
				return corrector_code;

			return SuccessCode::Success;
		}

	} // re: namespace tracking

} // re: namespace bertini





#endif

