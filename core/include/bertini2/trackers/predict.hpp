//This file is part of Bertini 2.
//
//tracking/predict.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracking/predict.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/predict.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file predict.hpp 

\brief Wrapper functions for calling ODE predictors for systems.
*/

#ifndef BERTINI_PREDICT_HPP
#define BERTINI_PREDICT_HPP

#include "bertini2/trackers/ode_predictors.hpp"

namespace bertini{
	namespace tracking{

		/**
		Wrapper class for calling an ODE predictor, using fixed precision

		\param predictor_choice The enum class selecting the predictor to be used.
		\param next_space The computed prediction.
		\param sys The system being solved.
		\param current_space The current space variable vector.
		\param current_time The current time.
		\param delta_t The size of the time step.
		\param condition_number_estimate The computed estimate of the condition number of the Jacobian.
		\param num_steps_since_last_condition_number_computation.  Updated in this function.
		\param frequency_of_CN_estimation How many steps to take between condition number estimates.
		\param prec_type The operating precision type.  
		\param tracking_tolerance How tightly to track the path.
		\param AMP_config The settings for adaptive multiple precision.

		\tparam ComplexT The complex number type for evaluation.
		\tparam TolT The numeric type for tolerances and state.
		*/
		template <typename ComplexT, typename TolT>
		SuccessCode Predict(Predictor predictor_choice,
							Vec<ComplexT> & next_space,
							   System const& sys,
							   Vec<ComplexT> const& current_space, ComplexT const& current_time, 
							   ComplexT const& delta_t,
							   TolT & condition_number_estimate,
							   unsigned & num_steps_since_last_condition_number_computation, 
							   unsigned frequency_of_CN_estimation, 
							   TolT const& tracking_tolerance)
		{
			predict::ExplicitRKPredictor predictor(predictor_choice);

			return predictor.Predict(next_space,
									sys,
									current_space, current_time, 
									delta_t,
									condition_number_estimate,
									num_steps_since_last_condition_number_computation, 
									frequency_of_CN_estimation, 
									tracking_tolerance);
		}







		




		/**
		An overload of Predict which returns (by reference) error estimate, and estimates of the norm of \f$J\f$ and \f$J^{-1}\f$ from AMP2 paper \cite AMP2.

		\see Predict
		*/
		template <typename ComplexT, typename TolT>
		SuccessCode Predict(Predictor predictor_choice,
							Vec<ComplexT> & next_space,
							TolT & size_proportion, /*\f$a\f$ from the AMP2 paper */
							TolT & norm_J,
							TolT & norm_J_inverse,
							System const& sys,
							Vec<ComplexT> const& current_space, ComplexT current_time, 
							ComplexT const& delta_t,
							TolT & condition_number_estimate,
							unsigned & num_steps_since_last_condition_number_computation, 
							unsigned frequency_of_CN_estimation, 
							TolT const& tracking_tolerance,
							AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			predict::ExplicitRKPredictor<ComplexT, TolT> predictor(predictor_choice);

			return predictor.Predict(next_space,
										size_proportion,
										norm_J,
										norm_J_inverse,
										sys,
										current_space, current_time, 
										delta_t,
										condition_number_estimate,
										num_steps_since_last_condition_number_computation, 
										frequency_of_CN_estimation, 
										tracking_tolerance,
										AMP_config);
		}







		/**
		An overload of Predict which returns (by reference) error estimate, and estimates of the norm of \f$J\f$ and \f$J^{-1}\f$, and the size proportion \f$a\f$ from AMP2 paper \cite AMP2.

		\see Predict
		*/
		template <typename ComplexT, typename TolT>
		SuccessCode Predict(Predictor predictor_choice,
							Vec<ComplexT> & next_space,
							TolT & error_estimate,
							TolT & size_proportion, /*\f$a\f$ from the AMP2 paper */
							TolT & norm_J,
							TolT & norm_J_inverse,
							System const& sys,
							Vec<ComplexT> const& current_space, ComplexT current_time, 
							ComplexT const& delta_t,
							TolT & condition_number_estimate,
							unsigned & num_steps_since_last_condition_number_computation, 
							unsigned frequency_of_CN_estimation, 
							TolT const& tracking_tolerance,
							AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			predict::ExplicitRKPredictor<ComplexT, TolT> predictor(predictor_choice);

			return predictor.Predict(next_space,
												error_estimate,
												size_proportion,
												norm_J,
												norm_J_inverse,
												sys,
												current_space, current_time, 
												delta_t,
												condition_number_estimate,
												num_steps_since_last_condition_number_computation, 
												frequency_of_CN_estimation, 
												tracking_tolerance,
												AMP_config);
		}



		/**
		Wrapper class for calling an ODE predictor, using adaptive precision, not returning some meta-data about the step.

		\param predictor_choice The enum class selecting the predictor to be used.
		\param next_space The computed prediction.
		\param sys The system being solved.
		\param current_space The current space variable vector.
		\param current_time The current time.
		\param delta_t The size of the time step.
		\param condition_number_estimate The computed estimate of the condition number of the Jacobian.
		\param num_steps_since_last_condition_number_computation.  Updated in this function.
		\param frequency_of_CN_estimation How many steps to take between condition number estimates.
		\param prec_type The operating precision type.  
		\param tracking_tolerance How tightly to track the path.
		\param AMP_config The settings for adaptive multiple precision.

		\tparam ComplexT The complex number type for evaluation.
		\tparam TolT The complex number type for evaluation.
		*/
		template <typename ComplexT, typename TolT>
		SuccessCode Predict(Predictor predictor_choice,
							Vec<ComplexT> & next_space,
							   System const& sys,
							   Vec<ComplexT> const& current_space, ComplexT current_time, 
							   ComplexT const& delta_t,
							   TolT & condition_number_estimate,
							   unsigned & num_steps_since_last_condition_number_computation, 
							   unsigned frequency_of_CN_estimation, 
							   TolT const& tracking_tolerance,
							   AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			TolT size_proportion, norm_J, norm_J_inverse;

			return Predict(predictor_choice,
			               next_space,
							size_proportion,
							norm_J,
							norm_J_inverse,
							sys,
							current_space, current_time, 
							delta_t,
							condition_number_estimate,
							num_steps_since_last_condition_number_computation, 
							frequency_of_CN_estimation, 
							tracking_tolerance,
							AMP_config);
		}


		



		

	} // re: namespace tracking

	
} // re: namespace bertini



#endif
