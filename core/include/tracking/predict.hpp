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

#ifndef BERTINI_PREDICT_HPP
#define BERTINI_PREDICT_HPP

#include "tracking/ode_predictors.hpp"

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

		\tparam ComplexType The complex number type for evaluation.
		\tparam RealType The complex number type for evaluation.
		*/
		template <typename ComplexType, typename RealType>
		SuccessCode Predict(config::Predictor predictor_choice,
							Vec<ComplexType> & next_space,
							   System & sys,
							   Vec<ComplexType> const& current_space, ComplexType current_time, 
							   ComplexType const& delta_t,
							   RealType & condition_number_estimate,
							   unsigned & num_steps_since_last_condition_number_computation, 
							   unsigned frequency_of_CN_estimation, 
							   RealType const& tracking_tolerance)
		{
			static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");

			switch (predictor_choice)
			{
				case config::Predictor::Euler:
				{
					return predict::Euler(next_space,
									sys,
									current_space, current_time, 
									delta_t,
									condition_number_estimate,
									num_steps_since_last_condition_number_computation, 
									frequency_of_CN_estimation, 
									tracking_tolerance);
					break;
				}

				default:
				{
					throw std::runtime_error("incompatible predictor choice in Predict");
				}
			}

			return SuccessCode::Failure;
		}







		




		/**
		An overload of Predict which returns (by reference) error estimate, and estimates of the norm of \f$J\f$ and \f$J^{-1}\f$ from AMP2 paper \cite AMP2.

		\see Predict
		*/
		template <typename ComplexType, typename RealType>
		SuccessCode Predict(config::Predictor predictor_choice,
							Vec<ComplexType> & next_space,
							RealType & size_proportion, /*\f$a\f$ from the AMP2 paper */
							RealType & norm_J,
							RealType & norm_J_inverse,
							System & sys,
							Vec<ComplexType> const& current_space, ComplexType current_time, 
							ComplexType const& delta_t,
							RealType & condition_number_estimate,
							unsigned & num_steps_since_last_condition_number_computation, 
							unsigned frequency_of_CN_estimation, 
							RealType const& tracking_tolerance,
							config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");

			switch (predictor_choice)
			{
				case config::Predictor::Euler:
				{
					return predict::Euler(next_space,
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
					break;
				}

				default:
				{
					throw std::runtime_error("incompatible predictor choice in Predict with norms returned, but no error estimate");
				}
			}

			return SuccessCode::Failure;
		}







		/**
		An overload of Predict which returns (by reference) error estimate, and estimates of the norm of \f$J\f$ and \f$J^{-1}\f$, and the size proportion \f$a\f$ from AMP2 paper \cite AMP2.

		\see Predict
		*/
		template <typename ComplexType, typename RealType>
		SuccessCode Predict(config::Predictor predictor_choice,
							Vec<ComplexType> & next_space,
							RealType & error_estimate,
							RealType & size_proportion, /*\f$a\f$ from the AMP2 paper */
							RealType & norm_J,
							RealType & norm_J_inverse,
							System & sys,
							Vec<ComplexType> const& current_space, ComplexType current_time, 
							ComplexType const& delta_t,
							RealType & condition_number_estimate,
							unsigned & num_steps_since_last_condition_number_computation, 
							unsigned frequency_of_CN_estimation, 
							RealType const& tracking_tolerance,
							config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");

			switch (predictor_choice)
			{
				case config::Predictor::HeunEuler:
				{
					return predict::HeunEuler(next_space,
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
					break;
				}

				default:
				{
					throw std::runtime_error("incompatible predictor choice in Predict, with return of error estimate");
				}
			}

			return SuccessCode::Failure;
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

		\tparam ComplexType The complex number type for evaluation.
		\tparam RealType The complex number type for evaluation.
		*/
		template <typename ComplexType, typename RealType>
		SuccessCode Predict(config::Predictor predictor_choice,
							Vec<ComplexType> & next_space,
							   System & sys,
							   Vec<ComplexType> const& current_space, ComplexType current_time, 
							   ComplexType const& delta_t,
							   RealType & condition_number_estimate,
							   unsigned & num_steps_since_last_condition_number_computation, 
							   unsigned frequency_of_CN_estimation, 
							   RealType const& tracking_tolerance,
							   config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");

			RealType size_proportion, norm_J, norm_J_inverse;

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


		


		namespace predict{

			/**
			The lowest order of the predictor.  The order of the error estimate is this plus one.
			*/
			inline
			unsigned Order(config::Predictor predictor_choice)
			{
				switch (predictor_choice)
				{
					case (config::Predictor::Euler):
						return 1;
					case (config::Predictor::HeunEuler):
						return 1;
				}
			}
			

			inline bool HasErrorEstimate(config::Predictor predictor_choice)
			{
				switch (predictor_choice)
				{
					case (config::Predictor::Euler):
						return false;
					case (config::Predictor::HeunEuler):
						return true;
				}
			}


		} // re: namespace predict

		

	} // re: namespace tracking

	
} // re: namespace bertini



#endif
