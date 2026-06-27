//This file is part of Bertini 2.0.
//
//heun_euler.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//heun_euler.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with heun_euler.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  predictor.hpp
//
//  copyright 2015
//  James B. Collins
//  West Texas A&M University
//  Department of Mathematics
//  Spring 2016


/**
 \file base_predictor.hpp
 
 \brief Contains a base class for all ODE predictors.
 */

#ifndef BERTINI_BASE_PREDICTORS_HPP
#define BERTINI_BASE_PREDICTORS_HPP

#include "bertini2/trackers/amp_criteria.hpp"
#include "bertini2/trackers/config.hpp"

#include "system.hpp"
#include "mpfr_extensions.hpp"
#include <Eigen/LU>

#include <boost/type_index.hpp>


namespace bertini{
	namespace tracking{
		namespace predict{
			
			
			using Predictor = Predictor;
			
			/**
			 \brief Get the Bertini2 default predictor.
			 
			 Currently set to Euler, though this will change in future versions.
			 */
			inline
			Predictor DefaultPredictor()
			{
				// RKF45 (Runge-Kutta-Fehlberg 4(5)): order-4 prediction with an embedded error
				// estimate.  Was Euler (order 1) -- an ODR clash with the RKF45 definition in
				// explicit_predictors.hpp, and atrocious for performance (tiny steps), and it forced
				// SetSizeProportion onto the no-error-estimate fallback maxCoeff(K)/|delta_t|^p, which
				// blows up at small steps and spuriously inflates AMP precision.  RKF45 has an error
				// estimate, so size_proportion = err_est/|delta_t|^(p+1) stays bounded.
				return Predictor::RKF45;
			}
			
			
			/**
			 The lowest order of the predictor.  The order of the error estimate is this plus one.
			 */
			inline
			unsigned Order(Predictor predictor_choice)
			{
				switch (predictor_choice)
				{
					case (Predictor::Euler):
						return 1;
					case (Predictor::HeunEuler):
						return 1;
					case (Predictor::RK4):
						return 4;
					case (Predictor::RKF45):
						return 4;
					case (Predictor::RKCashKarp45):
						return 4;
					case (Predictor::RKDormandPrince56):
						return 5;
					case (Predictor::RKVerner67):
						return 6;
					default:
					{
						throw std::runtime_error("incompatible predictor choice in Order");
					}
				}
			}
			
			
			inline bool HasErrorEstimate(Predictor predictor_choice)
			{
				switch (predictor_choice)
				{
					case (Predictor::Euler):
						return false;
					case (Predictor::HeunEuler):
						return true;
					case (Predictor::RK4):
						return false;
					case (Predictor::RKF45):
						return true;
					case (Predictor::RKCashKarp45):
						return true;
					case (Predictor::RKDormandPrince56):
						return true;
					case (Predictor::RKVerner67):
						return true;
					default:
					{
						throw std::runtime_error("incompatible predictor choice in HasErrorEstimate");
					}
				}
			}
			
			
			
			
			
			
			
			
			
			/**
			 /class BasePredictor
			 
			 \brief An interface for all predictors.
			 
			 ## Purpose
			 Stores all the information needed to implement the predictor method
				- Butcher Table
				- Number of Stages
				- Order of the method.
			 
			 Also stores information computed during implementation of the method.
			 
			 
			 ## Use
			 Implement the following pure functions:
				- FullStep
				- SetErrorEstimate
				- SetSizeProportion
			 
			 
			 
			 */
			
			class BasePredictor
			{
			public:
				
				
				BasePredictor(Predictor method){};
				
				/**
				 \brief Perform a generic predictor step.
				 
				 \param next_space The computed prediction.
				 \param method An enum class selecting the predictor method to use.
				 \param S The system being solved.
				 \param current_space The current space variable vector.
				 \param current_time The current time.
				 \param delta_t The size of the time step.
				 \param condition_number_estimate The computed estimate of the condition number of the Jacobian.
				 \param num_steps_since_last_condition_number_computation.  Updated in this function.
				 \param frequency_of_CN_estimation How many steps to take between condition number estimates.
				 \param prec_type The operating precision type.
				 \param tracking_tolerance How tightly to track the path.
				 */

				template<typename ComplexT, typename RealT>
				SuccessCode Predict(Vec<ComplexT> & next_space,
										 System const& S,
										 Vec<ComplexT> const& current_space, ComplexT current_time,
										 ComplexT const& delta_t,
										 RealT & condition_number_estimate,
										 unsigned & num_steps_since_last_condition_number_computation,
										 unsigned frequency_of_CN_estimation,
										 RealT const& tracking_tolerance)
				{
					
					return FullStep<ComplexT, RealT>(next_space, S, current_space, current_time, delta_t);
					
					
				}

				
				
				/**
				 \brief Perform a generic predictor step and return size_proportion and condition number information
				 
				 \param next_space The computed prediction.
				 \param method An enum class selecting the predictor method to use.
				 \param size_proportion $a$ in AMP2 paper.
				 \param norm_J The computed estimate of the norm of the Jacobian matrix.
				 \param norm_J_inverse The computed estimate of the norm of the inverse of the Jacobian matrix.
				 \param S The system being solved.
				 \param current_space The current space variable vector.
				 \param current_time The current time.
				 \param delta_t The size of the time step.
				 \param condition_number_estimate The computed estimate of the condition number of the Jacobian.
				 \param num_steps_since_last_condition_number_computation.  Updated in this function.
				 \param frequency_of_CN_estimation How many steps to take between condition number estimates.
				 \param prec_type The operating precision type.
				 \param tracking_tolerance How tightly to track the path.
				 \param AMP_config The settings for adaptive multiple precision.
				 */
				
				template<typename ComplexT, typename RealT>
				SuccessCode Predict(Vec<ComplexT> & next_space,
									RealT & size_proportion,
									RealT & norm_J,
									RealT & norm_J_inverse,
									System const& S,
									Vec<ComplexT> const& current_space, ComplexT current_time,
									ComplexT const& delta_t,
									RealT & condition_number_estimate,
									unsigned & num_steps_since_last_condition_number_computation,
									unsigned frequency_of_CN_estimation,
									RealT const& tracking_tolerance,
											AdaptiveMultiplePrecisionConfig const& AMP_config)
				{
					
					
					auto success_code = Predict<ComplexT, RealT>(next_space, S, current_space, current_time, delta_t,
												condition_number_estimate, num_steps_since_last_condition_number_computation,
												frequency_of_CN_estimation, tracking_tolerance);
					
					if(success_code != SuccessCode::Success)
						return success_code;
					
					// Calculate condition number and updated if needed
					Eigen::PartialPivLU<Mat<ComplexT>>& dhdxref = std::get< Eigen::PartialPivLU<Mat<ComplexT>> >(dh_dx_);
					Mat<ComplexT>& LUref = std::get< Mat<ComplexT> >(LU_);

					Vec<ComplexT> randy = RandomOfUnits<ComplexT>(S.NumVariables());
					Vec<ComplexT> temp_soln = LUref.solve(randy);
					
					norm_J = dhdxref.norm();
					norm_J_inverse = temp_soln.norm();
					
					if (num_steps_since_last_condition_number_computation >= frequency_of_CN_estimation)
					{
						condition_number_estimate = norm_J * norm_J_inverse;
						num_steps_since_last_condition_number_computation = 1; // reset the counter to 1
					}
					else // no need to compute the condition number
						num_steps_since_last_condition_number_computation++;
					
					
					// Set size_proportion
					SetSizeProportion<ComplexT,RealT>(size_proportion, delta_t);
					
					
					
					//AMP Criteria
					if (!amp::CriterionA(norm_J, norm_J_inverse, AMP_config)) // AMP_criterion_A != ok
						return SuccessCode::HigherPrecisionNecessary;
					else if (!amp::CriterionC(norm_J_inverse, current_space, tracking_tolerance, AMP_config)) // AMP_criterion_C != ok
						return SuccessCode::HigherPrecisionNecessary;
					
					
					return success_code;
				}

				
				/**
				 \brief Perform a generic predictor step and return error estimate, size_proportion and condition number information
				 
				 \param next_space The computed prediction.
				 \param method An enum class selecting the predictor method to use.
				 \param error_estimate Estimate of the error from an embedded method.
				 \param size_proportion $a$ in AMP2 paper.
				 \param norm_J The computed estimate of the norm of the Jacobian matrix.
				 \param norm_J_inverse The computed estimate of the norm of the inverse of the Jacobian matrix.
				 \param S The system being solved.
				 \param current_space The current space variable vector.
				 \param current_time The current time.
				 \param delta_t The size of the time step.
				 \param condition_number_estimate The computed estimate of the condition number of the Jacobian.
				 \param num_steps_since_last_condition_number_computation.  Updated in this function.
				 \param frequency_of_CN_estimation How many steps to take between condition number estimates.
				 \param prec_type The operating precision type.
				 \param tracking_tolerance How tightly to track the path.
				 \param AMP_config The settings for adaptive multiple precision.
				 */
				
				template<typename ComplexT, typename RealT>
				SuccessCode Predict(Vec<ComplexT> & next_space,
									RealT & error_estimate,
									RealT & size_proportion,
									RealT & norm_J,
									RealT & norm_J_inverse,
									System const& S,
									Vec<ComplexT> const& current_space, ComplexT current_time,
									ComplexT const& delta_t,
									RealT & condition_number_estimate,
									unsigned & num_steps_since_last_condition_number_computation,
									unsigned frequency_of_CN_estimation,
									RealT const& tracking_tolerance,
									AdaptiveMultiplePrecisionConfig const& AMP_config)
				{
					// If this is a method without an error estimator, then can't calculate size proportion and should throw an error
					
					if(!predict::HasErrorEstimate(predictor_))
					{
						throw std::runtime_error("incompatible predictor choice in ExplicitPredict, no error estimator");
					}
					
					
					
					
					auto success_code = Predict<ComplexT,RealT>(next_space, size_proportion, norm_J, norm_J_inverse,
												S, current_space, current_time, delta_t,
												condition_number_estimate, num_steps_since_last_condition_number_computation,
												frequency_of_CN_estimation, tracking_tolerance, AMP_config);
					
					if(success_code != SuccessCode::Success)
						return success_code;
					
					SetErrorEstimate<ComplexT,RealT>(error_estimate, delta_t);
					
					
					return success_code;
				}

				
				
				
				/**
				 /brief Sets the local variables to correspond to a particular predictor method
				 
				 \param method Enum class that determines the predictor method
				 
				 */
				
				virtual void PredictorMethod(Predictor method) = 0;
				
				Predictor PredictorMethod()
				{
					return predictor_;
				}
				
				
				
				
				/**
				 The lowest order of the predictor.  The order of the error estimate is this plus one.
				 */
				inline
				unsigned Order()
				{
					return p_;
				}
				
				
				inline bool HasErrorEstimate()
				{
					return predict::HasErrorEstimate(predictor_);
				}


				
				
				
			protected:
				///////////////////////////
				//
				//   Protected Members Methods
				//
				////////////////////
				
				/**
				 \brief Performs a full prediction step from current_time to current_time + delta_t
				 
				 \param next_space The computed prediction space
				 \param S The homotopy system
				 \param current_space The current space values
				 \param current_time The current time values
				 \param delta_t The time step
				 
				 \return SuccessCode determining result of the computation
				 */
				
				template<typename ComplexT, typename RealT>
				virtual SuccessCode FullStep(Vec<ComplexT> & next_space,
											 System const& S,
											 Vec<ComplexT> const& current_space, ComplexT current_time,
											 ComplexT const& delta_t) = 0;
				
				
				/**
				 \brief Computes the error estimate of this prediction step.
				 
				 \param error_estimate Computed error estimate
				 \param delta_t The time step
				 
				 \return Success code or the computation
				 
				 */
				
				template<typename ComplexT, typename RealT>
				virtual SuccessCode SetErrorEstimate(RealT & error_estimate, ComplexT const& delta_t) = 0;

				
				
				/**
				 \brief Compute the size proportion variable for AMP computation
				 
				 \param size_proportion Computed size proportion
				 \param delta_t The time step
				 
				 \return Success code of the computation
				 
				 */
				
				template<typename ComplexT, typename RealT>
				virtual SuccessCode SetSizeProportion(RealT & size_proportion, ComplexT const& delta_t) = 0;

				
				
				/**
				 \brief Evaluates the RHS of the Davidenko differential equation at a particular time and space
				 
				 \param S The homotopy system
				 \param space The space variable used to evaluate RHS
				 \param time The time variable used to evaluate RHS
				 \param K Matrix of stage variables
				 \param stage Which stage variable(column of K) should be filled by this computation
				 
				 \return Success code of this computation
				 */
				
				template<typename ComplexT>
				virtual SuccessCode EvalRHS(System const& S,
									Vec<ComplexT> const& space, ComplexT time, Mat<ComplexT> & K, int stage) = 0;


				
				
				
				
				
				
				
				
				///////////////////////////
				//
				//   Protected Data Members
				//
				////////////////////
				
				Predictor predictor_;
				unsigned p_;
				std::tuple< Mat<complex_dbl>, Mat<mpfr> > dh_dx_;
				std::tuple< Eigen::PartialPivLU<Mat<complex_dbl>>, Eigen::PartialPivLU<Mat<mpfr>> > LU_;
				

				
			}; // re: class BasePredictor
			
		} // re: namespace predict
	}// re: namespace tracking
}// re: namespace bertini







#endif

