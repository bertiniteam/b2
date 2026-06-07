//This file is part of Bertini 2.
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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire


/**
\file step.hpp 

\brief Provides Step wrappers for taking a generic predict-correct step.  

This provides no calls for adaptive precision, those living with the AMPTracker type instead.
*/


#ifndef BERTINI_STEP_HPP
#define BERTINI_STEP_HPP

#include "bertini2/trackers/predict.hpp"
#include "bertini2/trackers/correct.hpp"

namespace bertini {

	namespace tracking {

		/**
	

		\tparam ComplexT The type of number being used in the algorithm

		*/
		template <typename ComplexT, typename RealT>
		SuccessCode Step(Predictor predictor_choice,
		                    Vec<ComplexT> & next_space, ComplexT & next_time,
				               System & sys,
				               Vec<ComplexT> const& current_space, ComplexT current_time, 
				               ComplexT const& delta_t,
				               RealT & condition_number_estimate,
				               unsigned & num_steps_since_last_condition_number_computation, 
				               unsigned frequency_of_CN_estimation, PrecisionType prec_type, 
				               RealT const& tracking_tolerance,
				               RealT const& path_truncation_threshold,
				               unsigned min_num_newton_iterations,
				               unsigned max_num_newton_iterations,
				               AdaptiveMultiplePrecisionConfig const& AMP_config)
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

