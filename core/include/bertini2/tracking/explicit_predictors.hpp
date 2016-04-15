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
 \file predictor.hpp
 
 \brief Contains a base class for all ODE predictors.
 */

#ifndef BERTINI_EXPLICIT_PREDICTORS_HPP
#define BERTINI_EXPLICIT_PREDICTORS_HPP

#include "tracking/amp_criteria.hpp"
#include "tracking/tracking_config.hpp"

#include "system.hpp"
#include "mpfr_extensions.hpp"
#include <Eigen/LU>

#include <boost/type_index.hpp>


namespace bertini{
	namespace tracking{
		namespace predict{
			
						
			using Predictor = config::Predictor;
			
			/**
			 \brief Get the Bertini2 default predictor.
			 
			 Currently set to Euler, though this will change in future versions.
			 */
			inline
			config::Predictor DefaultPredictor()
			{
				return Predictor::Euler;
			}
			
			
			/**
			 The lowest order of the predictor.  The order of the error estimate is this plus one.
			 */
			inline
			unsigned Order(config::Predictor predictor_choice)
			{
				switch (predictor_choice)
				{
					case (Predictor::Euler):
						return 1;
					case (Predictor::HeunEuler):
						return 1;
					default:
					{
						throw std::runtime_error("incompatible predictor choice in Order");
					}
				}
			}
			
			
			inline bool HasErrorEstimate(config::Predictor predictor_choice)
			{
				switch (predictor_choice)
				{
					case (Predictor::Euler):
						return false;
					case (Predictor::HeunEuler):
						return true;
					default:
					{
						throw std::runtime_error("incompatible predictor choice in HasErrorEstimate");
					}
				}
			}
			
			
			
			
			
			/**
			 /class ExplicitRKPredictor
			 
			 \brief A static class which stores all the explicit ODE predictor methods.
			 
			 ## Purpose
			 Stores all the information needed to implement the predictor method
				- Butcher Table
				- Number of Stages
				- Order of the method.  
			 
			 Also stores information computed during implementation of the method.
			 
			 
			 ## Use
			 Each predictor method is a static method of this class, therefore all that is required to run the predictor is to call the correct method:
			 
			 \code
			 ExplicitRKPredictors<Complex,Real> euler(config::Predictor::Euler)
			 success_code = euler.Predict( ... )
			 \endcode
			 
			 
			 
			 */
			
			template <typename ComplexType, typename RealType>
			class ExplicitRKPredictor
			{
			public:
				
				
				
				
				/**
				 \brief Constructor for a particular predictor method
				 
				 
				 \param method The predictor method to be implemented.
				 
				 */
				
				
				ExplicitRKPredictor(Predictor method)
				{
					predictor_ = method;
					p_ = Order(method);
					switch(method)
					{
						case Predictor::Euler:
						{
							s_ = 1;
							c_ = Vec<RealType>(s_); c_(0) = 0;
							a_ = Mat<RealType>(s_,s_); a_(0,0) = 0;
							b_ = Vec<RealType>(s_); b_(0) = 1;
							break;
						}
						case Predictor::HeunEuler:
						{
							mpq_rational half = mpq_rational(1,2);
							s_ = 2;
							c_ = Vec<RealType>(s_); c_ << RealType(0), RealType(1);
							a_ = Mat<RealType>(s_, s_); a_ << RealType(0), RealType(0), RealType(1), RealType(0);
							b_ = Vec<RealType>(s_); b_ << RealType(.5), RealType(.5);
							bstar_ = Vec<RealType>(s_); bstar_ << RealType(1), RealType(0);
							break;
						}
						default:
						{
							throw std::runtime_error("incompatible predictor choice in ExplicitPredict");
						}
					}

				}
				
				
				
				
				
				
				
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
				
				SuccessCode Predict(Vec<ComplexType> & next_space, Predictor method,
										 System const& S,
										 Vec<ComplexType> const& current_space, ComplexType current_time,
										 ComplexType const& delta_t,
										 RealType & condition_number_estimate,
										 unsigned & num_steps_since_last_condition_number_computation,
										 unsigned frequency_of_CN_estimation,
										 RealType const& tracking_tolerance)
				{
					
					K_ = Mat<ComplexType>(S.NumTotalFunctions(), s_);
					return FullStep(next_space, S, current_space, current_time, delta_t);
					
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
				
				SuccessCode Predict(Vec<ComplexType> & next_space,
												   Predictor method,
												   RealType & size_proportion,
												   RealType & norm_J,
												   RealType & norm_J_inverse,
												   System const& S,
												   Vec<ComplexType> const& current_space, ComplexType current_time,
												   ComplexType const& delta_t,
												   RealType & condition_number_estimate,
												   unsigned & num_steps_since_last_condition_number_computation,
												   unsigned frequency_of_CN_estimation,
												   RealType const& tracking_tolerance,
												   config::AdaptiveMultiplePrecisionConfig const& AMP_config)
				{
					
					
					auto success_code = Predict(next_space, method, S, current_space, current_time, delta_t,
									condition_number_estimate, num_steps_since_last_condition_number_computation,
														frequency_of_CN_estimation, tracking_tolerance);

					
					// Calculate condition number and updated if needed
					Vec<ComplexType> randy = RandomOfUnits<ComplexType>(S.NumVariables());
					Vec<ComplexType> temp_soln = LU_.solve(randy);
					
					norm_J = dh_dx_.norm();
					norm_J_inverse = temp_soln.norm();
					
					if (num_steps_since_last_condition_number_computation >= frequency_of_CN_estimation)
					{
						condition_number_estimate = norm_J * norm_J_inverse;
						num_steps_since_last_condition_number_computation = 1; // reset the counter to 1
					}
					else // no need to compute the condition number
						num_steps_since_last_condition_number_computation++;

					
					// Set size_proportion
					SetSizeProportion(size_proportion, delta_t);
					
					
					
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
				
				SuccessCode Predict(Vec<ComplexType> & next_space,
												   Predictor method,
												   RealType & error_estimate,
												   RealType & size_proportion,
												   RealType & norm_J,
												   RealType & norm_J_inverse,
												   System const& S,
												   Vec<ComplexType> const& current_space, ComplexType current_time,
												   ComplexType const& delta_t,
												   RealType & condition_number_estimate,
												   unsigned & num_steps_since_last_condition_number_computation,
												   unsigned frequency_of_CN_estimation,
												   RealType const& tracking_tolerance,
												   config::AdaptiveMultiplePrecisionConfig const& AMP_config)
				{
					// If this is a method without an error estimator, then can't calculate size proportion and should throw an error
					
					if(!HasErrorEstimate(method))
					{
						throw std::runtime_error("incompatible predictor choice in ExplicitPredict, no error estimator");
					}

					
					
					
					auto success_code = Predict(next_space, method, size_proportion, norm_J, norm_J_inverse,
									S, current_space, current_time, delta_t,
									condition_number_estimate, num_steps_since_last_condition_number_computation,
														frequency_of_CN_estimation, tracking_tolerance, AMP_config);
					
					SetErrorEstimate(error_estimate, delta_t);
					
					
//					std::cout << K_.col(0).norm() << std::endl;
					
					return success_code;
				}

				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
			private:
				
				
				
				
				
				
				void MethodSetup(Predictor method)
				{
					p_ = Order(method);
					switch(method)
					{
						case Predictor::Euler:
						{
							s_ = 1;
							c_ = Vec<RealType>(s_); c_(0) = 0;
							a_ = Mat<RealType>(s_,s_); a_(0,0) = 0;
							b_ = Vec<RealType>(s_); b_(0) = 1;
							break;
						}
						case Predictor::HeunEuler:
						{
							mpq_rational half = mpq_rational(1,2);
							s_ = 2;
							c_ = Vec<RealType>(s_); c_ << RealType(0), RealType(1);
							a_ = Mat<RealType>(s_, s_); a_ << RealType(0), RealType(0), RealType(1), RealType(0);
							b_ = Vec<RealType>(s_); b_ << RealType(.5), RealType(.5);
							bstar_ = Vec<RealType>(s_); bstar_ << RealType(1), RealType(0);
							break;
						}
						default:
						{
							throw std::runtime_error("incompatible predictor choice in ExplicitPredict");
						}
					}
				}
				
				
				

				
				
				SuccessCode FullStep(Vec<ComplexType> & next_space,
									System const& S,
									Vec<ComplexType> const& current_space, ComplexType current_time,
									 ComplexType const& delta_t)
				{
					Vec<ComplexType> temp = Vec<ComplexType>(S.NumTotalFunctions());
					
					
					if(EvalRHS(S, current_space, current_time, K_, 0) != SuccessCode::Success)
					{
						return SuccessCode::MatrixSolveFailureFirstPartOfPrediction;
					}
					
					for(int ii = 1; ii < s_; ++ii)
					{
						temp.setZero();
						for(int jj = 0; jj < ii; ++jj)
						{
							temp += a_(ii,jj)*K_.col(jj);
							
						}
						
						if(EvalRHS(S, current_space + delta_t*temp, current_time + c_(ii)*delta_t, K_, ii) != SuccessCode::Success)
						{
							return SuccessCode::MatrixSolveFailure;
						}
					}
					
					
					temp.setZero();
					for(int ii = 0; ii < s_; ++ii)
					{
						temp += b_(ii)*K_.col(ii);
					}
										
					next_space = current_space + delta_t*temp;
					
					
					return SuccessCode::Success;
				};

				
				
				
				
				
				
				SuccessCode SetErrorEstimate(RealType & error_estimate, ComplexType const& delta_t)
				{
					auto numFuncs = K_.cols();
					Vec<ComplexType> err = Vec<ComplexType>(numFuncs);
					
					err.setZero();
					for(int ii = 0; ii < s_; ++ii)
					{
						err += (b_(ii) - bstar_(ii))*K_.col(ii);
					}
					
					err *= delta_t;
					
					error_estimate = err.norm();
					
					return SuccessCode::Success;
				};
				
				

				
				
				
				
				SuccessCode SetSizeProportion(RealType & size_proportion, ComplexType const& delta_t)
				{
					if(HasErrorEstimate(predictor_))
					{
						RealType err_est;
						SetErrorEstimate(err_est, delta_t);
						
						using std::pow;
						size_proportion = err_est/(pow(abs(delta_t), p_+1));
						
						return SuccessCode::Success;
					}
					else
					{
						size_proportion = K_.array().abs().maxCoeff();
						return SuccessCode::Success;
					}
				};
				
				
				
				
				

				
				SuccessCode EvalRHS(System const& S,
									 Vec<ComplexType> const& space, ComplexType time, Mat<ComplexType> & K, int stage)
				{
					if(stage == 0)
					{
						dh_dx_ = S.Jacobian(space, time);
						LU_ = dh_dx_.lu();
						
						if (LUPartialPivotDecompositionSuccessful(LU_.matrixLU())!=MatrixSuccessCode::Success)
							return SuccessCode::MatrixSolveFailureFirstPartOfPrediction;
						
						K.col(stage) = LU_.solve(-S.TimeDerivative(space, time));
						
						return SuccessCode::Success;

					}
					else
					{
						auto dh_dx = S.Jacobian(space, time);
						auto LU = dh_dx.lu();
						
						if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
							return SuccessCode::MatrixSolveFailureFirstPartOfPrediction;
						
						K.col(stage) = LU.solve(-S.TimeDerivative(space, time));
						
						return SuccessCode::Success;
					}
				}
				
				
				
				
				
				
				
				
				
				
				// Data Members
				config::Predictor predictor_ = Predictor::None;
				
				
				Mat<ComplexType> dh_dx_;
				Eigen::PartialPivLU<Mat<ComplexType>> LU_;
				Mat<ComplexType> K_;
				
				// Butcher Table (notation from https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)
				int s_; // Number of stages
				Vec<RealType> c_;
				Vec<RealType> b_;
				Vec<RealType> bstar_;
				Mat<RealType> a_;
				int p_;
				
				
				
				
				
				
				// static const variables that store the Butcher table in mpq_rational form
				static const mpq_rational cEulerPtr_[];
				static const Vec<mpq_rational> cEuler_;
				
				
			};
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::cEulerPtr_[] = {mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Vec<mpq_rational> ExplicitRKPredictor<CType, RType>::cEuler_(cEulerPtr_);
			
			
			
			
			
			
			
			
			
			
			

			
			
			
		} // re: predict
	}// re: tracking
}// re: bertini

#endif
