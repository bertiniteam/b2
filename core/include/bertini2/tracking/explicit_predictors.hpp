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

//  base_predictor.hpp
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
			Predictor DefaultPredictor()
			{
				return Predictor::Euler;
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
			 /class ExplicitRKPredictor
			 
			 \brief A class which stores all the explicit single-step multi-stage ODE predictor methods.
			 
			 ## Purpose
			 Stores all the information needed to implement the predictor method
				- Butcher Table
				- Number of Stages
				- Order of the method.  
			 
			 Also stores information computed during implementation of the method.
			 
			 
			 ## Use
			 Each predictor method is stored as a static Butcher table.  To perform a predict step, you must instantiate an object with a particular predictor method and call Predict:
			 
			 \code
			 ExplicitRKPredictors<Complex,Real> euler(config::Predictor::Euler, sys)
			 success_code = euler.Predict( ... )
			 \endcode
			 
			 
			 
			 */
			
			class ExplicitRKPredictor
			{
			public:
				
				
				ExplicitRKPredictor(const System& S)
				{
					numTotalFunctions_ = S.NumTotalFunctions();
					numVariables_ = S.NumVariables();
					std::get< Mat<dbl> >(dh_dx_0_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<mpfr> >(dh_dx_0_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<dbl> >(dh_dx_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<mpfr> >(dh_dx_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Vec<dbl> >(dh_dt_temp_).resize(numTotalFunctions_);
					std::get< Vec<mpfr> >(dh_dt_temp_).resize(numTotalFunctions_);
					PredictorMethod(DefaultPredictor());
				}
				
				/**
				 \brief Constructor for a particular predictor method
				 
				 
				 \param method The predictor method to be implemented.
				 
				 */
				
				ExplicitRKPredictor(Predictor method, const System& S)
				{
					numTotalFunctions_ = S.NumTotalFunctions();
					numVariables_ = S.NumVariables();
					std::get< Mat<dbl> >(dh_dx_0_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<mpfr> >(dh_dx_0_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<dbl> >(dh_dx_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<mpfr> >(dh_dx_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Vec<dbl> >(dh_dt_temp_).resize(numTotalFunctions_);
					std::get< Vec<mpfr> >(dh_dt_temp_).resize(numTotalFunctions_);
					PredictorMethod(method);
				}
				
				
				
				
				/**
				 /brief Sets the local variables to correspond to a particular predictor method
				 
				 \param method Enum class that determines the predictor method
				 
				 */
				
				void PredictorMethod(Predictor method)
				{
					predictor_ = method;
					p_ = predict::Order(method);
					switch(method)
					{
						case Predictor::Euler:
						{
							s_ = 1;
							Mat<double>& arefd = std::get< Mat<double> >(a_);
							Vec<double>& brefd = std::get< Vec<double> >(b_);
							Vec<double>& crefd = std::get< Vec<double> >(c_);
							crefd.resize(s_); crefd(0) = static_cast<double>(cEuler_(0));
							arefd.resize(s_,s_); arefd(0,0) = static_cast<double>(aEuler_(0,0));
							brefd.resize(s_); brefd(0) = static_cast<double>(bEuler_(0));
							Mat<mpfr_float>& arefmp = std::get< Mat<mpfr_float> >(a_);
							Vec<mpfr_float>& brefmp = std::get< Vec<mpfr_float> >(b_);
							Vec<mpfr_float>& crefmp = std::get< Vec<mpfr_float> >(c_);
							crefmp.resize(s_); crefmp(0) = static_cast<mpfr_float>(cEuler_(0));
							arefmp.resize(s_,s_); arefmp(0,0) = static_cast<mpfr_float>(aEuler_(0,0));
							brefmp.resize(s_); brefmp(0) = static_cast<mpfr_float>(bEuler_(0));
							
							break;
						}
						case Predictor::HeunEuler:
						{
							s_ = 2;
							
							FillButcherTable<double>(s_, aHeunEuler_, bHeunEuler_, b_minus_bstarHeunEuler_, cHeunEuler_);
							FillButcherTable<mpfr_float>(s_, aHeunEuler_, bHeunEuler_, b_minus_bstarHeunEuler_, cHeunEuler_);
							
							break;
						}
						case Predictor::RK4:
						{
							s_ = 4;
							
							FillButcherTable<double>(s_, aRK4_, bRK4_, cRK4_);
							FillButcherTable<mpfr_float>(s_, aRK4_, bRK4_, cRK4_);
							
							break;
						}
							
						case Predictor::RKF45:
						{
							s_ = 6;
							
							FillButcherTable<double>(s_, aRKF45_, bRKF45_, b_minus_bstarRKF45_, cRKF45_);
							FillButcherTable<mpfr_float>(s_, aRKF45_, bRKF45_, b_minus_bstarRKF45_, cRKF45_);
							
							break;
						}
							
						case Predictor::RKCashKarp45:
						{
							s_ = 6;
							
							FillButcherTable<double>(s_, aRKCK45_, bRKCK45_, b_minus_bstarRKCK45_, cRKCK45_);
							FillButcherTable<mpfr_float>(s_, aRKCK45_, bRKCK45_, b_minus_bstarRKCK45_, cRKCK45_);
							
							break;
						}
							
						case Predictor::RKDormandPrince56:
						{
							s_ = 8;
							
							FillButcherTable<double>(s_, aRKDP56_, bRKDP56_, b_minus_bstarRKDP56_, cRKDP56_);
							FillButcherTable<mpfr_float>(s_, aRKDP56_, bRKDP56_, b_minus_bstarRKDP56_, cRKDP56_);
							
							break;
						}
							
						case Predictor::RKVerner67:
						{
							s_ = 10;
							
							FillButcherTable<double>(s_, aRKV67_, bRKV67_, b_minus_bstarRKV67_, cRKV67_);
							FillButcherTable<mpfr_float>(s_, aRKV67_, bRKV67_, b_minus_bstarRKV67_, cRKV67_);
							
							break;
						}
							
						default:
						{
							throw std::runtime_error("incompatible predictor choice in ExplicitPredict");
						}
					}
						
					std::get< Mat<dbl> >(K_).resize(numTotalFunctions_, s_);
					std::get< Mat<mpfr> >(K_).resize(numTotalFunctions_, s_);

					
					
				}; // re: PredictorMethod
				
				
				
				
				
				/**
				 \brief Change the system(number of total functions) that the predictor uses.
				 
				 \param S New system
				 
				 */
				void PredictorSystem(const System& S)
				{
					numTotalFunctions_ = S.NumTotalFunctions();
					numVariables_ = S.NumVariables();
					std::get< Mat<dbl> >(K_).resize(numTotalFunctions_, s_);
					std::get< Mat<mpfr> >(K_).resize(numTotalFunctions_, s_);
					std::get< Mat<dbl> >(dh_dx_0_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<mpfr> >(dh_dx_0_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<dbl> >(dh_dx_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<mpfr> >(dh_dx_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Vec<dbl> >(dh_dt_temp_).resize(numTotalFunctions_);
					std::get< Vec<mpfr> >(dh_dt_temp_).resize(numTotalFunctions_);
				}
				
				
				
				
				
				
				/** 
				 /brief Change the precision of the predictor variables and reassign the Butcher table variables.
				 
				 \param new_precision The new precision.
				 
				 */
				void ChangePrecision(unsigned new_precision)
				{
					for(int ii = 0; ii < numTotalFunctions_; ++ii)
					{
						for(int jj = 0; jj < s_; ++jj)
						{
							std::get< Mat<mpfr> >(K_)(ii,jj).precision(new_precision);
						}
					}
					
					for(int ii = 0; ii < numTotalFunctions_; ++ii)
					{
						std::get< Vec<mpfr> >(dh_dt_temp_)(ii).precision(new_precision);
						for(int jj = 0; jj < numVariables_; ++jj)
						{
							std::get< Mat<mpfr> >(dh_dx_0_)(ii,jj).precision(new_precision);
							std::get< Mat<mpfr> >(dh_dx_temp_)(ii,jj).precision(new_precision);
						}
					}

					
					PredictorMethod(predictor_);
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
				
				template<typename ComplexType, typename RealType, typename Derived>
				SuccessCode Predict(Vec<ComplexType> & next_space,
									System const& S,
									const Eigen::MatrixBase<Derived>& current_space, ComplexType current_time,
									ComplexType const& delta_t,
									RealType & condition_number_estimate,
									unsigned & num_steps_since_last_condition_number_computation,
									unsigned frequency_of_CN_estimation,
									RealType const& tracking_tolerance)
				{
					static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");
					static_assert(std::is_same<typename Derived::Scalar, ComplexType>::value, "scalar types must match");

					
					return FullStep<ComplexType, RealType>(next_space, S, current_space, current_time, delta_t);
					
					
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
				
				template<typename ComplexType, typename RealType, typename Derived>
				SuccessCode Predict(Vec<ComplexType> & next_space,
									RealType & size_proportion,
									RealType & norm_J,
									RealType & norm_J_inverse,
									System const& S,
									const Eigen::MatrixBase<Derived>& current_space, ComplexType current_time,
									ComplexType const& delta_t,
									RealType & condition_number_estimate,
									unsigned & num_steps_since_last_condition_number_computation,
									unsigned frequency_of_CN_estimation,
									RealType const& tracking_tolerance,
									config::AdaptiveMultiplePrecisionConfig const& AMP_config)
				{
					static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");
					static_assert(std::is_same<typename Derived::Scalar, ComplexType>::value, "scalar types must match");

					
					auto success_code = Predict<ComplexType, RealType>(next_space, S, current_space, current_time, delta_t,
																	   condition_number_estimate, num_steps_since_last_condition_number_computation,
																	   frequency_of_CN_estimation, tracking_tolerance);
					
					if(success_code != SuccessCode::Success)
						return success_code;
					
					// Calculate condition number and updated if needed
					Eigen::PartialPivLU<Mat<ComplexType>>& LUref = std::get< Eigen::PartialPivLU<Mat<ComplexType>> >(LU_0_);
					Mat<ComplexType>& dhdxref = std::get< Mat<ComplexType> >(dh_dx_0_);
					
					Vec<ComplexType> randy = RandomOfUnits<ComplexType>(S.NumVariables());
					Vec<ComplexType> temp_soln = LUref.solve(randy);
					
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
					SetSizeProportion<ComplexType,RealType>(size_proportion, delta_t);
					
					
					
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
				
				template<typename ComplexType, typename RealType, typename Derived>
				SuccessCode Predict(Vec<ComplexType> & next_space,
									RealType & error_estimate,
									RealType & size_proportion,
									RealType & norm_J,
									RealType & norm_J_inverse,
									System const& S,
									const Eigen::MatrixBase<Derived>& current_space, ComplexType current_time,
									ComplexType const& delta_t,
									RealType & condition_number_estimate,
									unsigned & num_steps_since_last_condition_number_computation,
									unsigned frequency_of_CN_estimation,
									RealType const& tracking_tolerance,
									config::AdaptiveMultiplePrecisionConfig const& AMP_config)
				{
					static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");
					static_assert(std::is_same<typename Derived::Scalar, ComplexType>::value, "scalar types must match");

					// If this is a method without an error estimator, then can't calculate size proportion and should throw an error
					
					if(!predict::HasErrorEstimate(predictor_))
					{
						throw std::runtime_error("incompatible predictor choice in ExplicitPredict, no error estimator");
					}
					
					
					
					
					auto success_code = Predict<ComplexType,RealType>(next_space, size_proportion, norm_J, norm_J_inverse,
																	  S, current_space, current_time, delta_t,
																	  condition_number_estimate, num_steps_since_last_condition_number_computation,
																	  frequency_of_CN_estimation, tracking_tolerance, AMP_config);
					
					if(success_code != SuccessCode::Success)
						return success_code;
					
					SetErrorEstimate<ComplexType,RealType>(error_estimate, delta_t);
					
					
					return success_code;
				}
				
				
				
				
				
				
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
				// Protected Methods
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
				
				template<typename ComplexType, typename RealType, typename Derived>
				SuccessCode FullStep(Vec<typename Derived::Scalar> & next_space,
									System const& S,
									 Eigen::MatrixBase<Derived> const& current_space, ComplexType const& current_time,
									 ComplexType const& delta_t)
				{
					static_assert(std::is_same<typename Derived::Scalar, ComplexType>::value, "scalar types must match");
					
					Mat<ComplexType>& Kref = std::get< Mat<ComplexType> >(K_);
					Mat<RealType>& aref = std::get< Mat<RealType> >(a_);
					Vec<RealType>& bref = std::get< Vec<RealType> >(b_);
					Vec<RealType>& cref = std::get< Vec<RealType> >(c_);
					Kref.fill(ComplexType(0));
					Vec<ComplexType> temp(S.NumTotalFunctions());
					
					if(EvalRHS(S, current_space, current_time, Kref, 0) != SuccessCode::Success)
					{
						return SuccessCode::MatrixSolveFailureFirstPartOfPrediction;
					}
					
					for(int ii = 1; ii < s_; ++ii)
					{
						temp.setZero();
						for(int jj = 0; jj < ii; ++jj)
						{
							temp += aref(ii,jj)*Kref.col(jj);
							
						}
						
						if(EvalRHS(S, current_space + delta_t*temp, current_time + cref(ii)*delta_t, Kref, ii) != SuccessCode::Success)
						{
							return SuccessCode::MatrixSolveFailure;
						}
					}
					
					
					temp.setZero();
					for(int ii = 0; ii < s_; ++ii)
					{
						temp += bref(ii)*Kref.col(ii);
					}
										
					next_space = current_space + delta_t*temp;
					
					return SuccessCode::Success;
				};

				
				
				
				
				/**
				 \brief Computes the error estimate of this prediction step.
				 
				 \param error_estimate Computed error estimate
				 \param delta_t The time step
				 
				 \return Success code or the computation
				 
				 */
				
				template<typename ComplexType, typename RealType>
				SuccessCode SetErrorEstimate(RealType & error_estimate, ComplexType const& delta_t)
				{
					Mat<ComplexType>& Kref = std::get< Mat<ComplexType> >(K_);
					Vec<RealType>& b_minus_bstar_ref = std::get< Vec<RealType> >(b_minus_bstar_);
					
					auto numFuncs = Kref.rows();
					Vec<ComplexType> err(numFuncs);
					
					err.setZero();
					for(int ii = 0; ii < s_; ++ii)
					{
						err += (b_minus_bstar_ref(ii))*Kref.col(ii);
					}
					
					err *= delta_t;
					
					error_estimate = err.norm();
					
					return SuccessCode::Success;
				};
				
				

				
				
				
				/**
				 \brief Compute the size proportion variable for AMP computation
				 
				 \param size_proportion Computed size proportion
				 \param delta_t The time step
				 
				 \return Success code of the computation
				 
				 */
				
				template<typename ComplexType, typename RealType>
				SuccessCode SetSizeProportion(RealType & size_proportion, ComplexType const& delta_t)
				{
					if(predict::HasErrorEstimate(predictor_))
					{
						RealType err_est;
						SetErrorEstimate<ComplexType,RealType>(err_est, delta_t);
						
						using std::pow;
						size_proportion = err_est/(pow(abs(delta_t), p_+1));
						
						return SuccessCode::Success;
					}
					else
					{
						Mat<ComplexType>& Kref = std::get< Mat<ComplexType> >(K_);
						using std::pow;
						size_proportion = Kref.array().abs().maxCoeff()/(pow(abs(delta_t), p_));
						return SuccessCode::Success;
					}
				};
				
				
				
				/**
				 \brief Evaluates the RHS of the Davidenko differential equation at a particular time and space
				 
				 \param S The homotopy system
				 \param space The space variable used to evaluate RHS
				 \param time The time variable used to evaluate RHS
				 \param K Matrix of stage variables
				 \param stage Which stage variable(column of K) should be filled by this computation
				 
				 \return Success code of this computation
				 */
				
				template< typename Derived, typename ComplexType>
				SuccessCode EvalRHS(System const& S,
									const Eigen::MatrixBase<Derived>& space, const ComplexType& time, Mat<ComplexType> & K, unsigned stage)
				{
					static_assert(std::is_same<typename Derived::Scalar, ComplexType>::value, "scalar types must match");

					if(stage == 0)
					{
						Eigen::PartialPivLU<Mat<ComplexType>>& LUref = std::get< Eigen::PartialPivLU<Mat<ComplexType>> >(LU_0_);
						Mat<ComplexType>& dhdxref = std::get< Mat<ComplexType> >(dh_dx_0_);
						S.JacobianInPlace(dhdxref,space, time);
						LUref = dhdxref.lu();
						
						if (LUPartialPivotDecompositionSuccessful(LUref.matrixLU())!=MatrixSuccessCode::Success)
							return SuccessCode::MatrixSolveFailureFirstPartOfPrediction;
						
						Vec<ComplexType>& dhdtref = std::get< Vec<ComplexType> >(dh_dt_temp_);
						S.TimeDerivativeInPlace(dhdtref, space, time);
						K.col(stage) = LUref.solve(-dhdtref);
						
						return SuccessCode::Success;
						
					}
					else
					{
						Mat<ComplexType>& dhdxtempref = std::get< Mat<ComplexType> >(dh_dx_temp_);
						S.JacobianInPlace(dhdxtempref,space, time);
						auto LU = dhdxtempref.lu();
						
						if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
							return SuccessCode::MatrixSolveFailureFirstPartOfPrediction;
						
						Vec<ComplexType>& dhdtref = std::get< Vec<ComplexType> >(dh_dt_temp_);
						S.TimeDerivativeInPlace(dhdtref, space, time);
						K.col(stage) = LU.solve(-dhdtref);
						
						return SuccessCode::Success;
					}
				}
				
				
			
				
				
				
				
				
			private:
				///////////////////////////
				//
				// Private Methods
				//
				////////////////////

				
				
				/**
				 /brief Fills the local embedded butcher table variables a,b,bstar and c with the constant static values stored in the class.
				 
				 \param stages Number of stages.  Used to create correct size on variables
				 \param a Matrix of the Butcher table
				 \param b Weights in the Butcher table
				 \param bstar Weights of embedded method in Butcher table
				 \param c Time offsets in Butcher table
				 
				 */
				
				template<typename RealType>
				void FillButcherTable(int stages, const Mat<mpq_rational>& a,
								 const Mat<mpq_rational> & b,
								 const Mat<mpq_rational> & b_minus_bstar,
								 const Mat<mpq_rational> & c)
				{
					Mat<RealType>& aref = std::get< Mat<RealType> >(a_);
					aref.resize(stages, stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						for(int jj = 0; jj < s_; ++jj)
						{
							aref(ii,jj) = static_cast<RealType>(a(ii,jj));
						}
					}
					
					Vec<RealType>& bref = std::get< Vec<RealType> >(b_);
					bref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						bref(ii) = static_cast<RealType>(b(ii));
					}
					
					Vec<RealType>& b_minus_bstar_ref = std::get< Vec<RealType> >(b_minus_bstar_);
					b_minus_bstar_ref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						b_minus_bstar_ref(ii) = static_cast<RealType>(b_minus_bstar(ii));
					}

					Vec<RealType>& cref = std::get< Vec<RealType> >(c_);
					cref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						cref(ii) = static_cast<RealType>(c(ii));
						
					}

				}
				
				

				
				
				/**
				 /brief Fills the local butcher table variables a,b and c with the constant static values stored in the class.
				 
				 \param stages Number of stages.  Used to create correct size on variables
				 \param a Matrix of the Butcher table
				 \param b Weights in the Butcher table
				 \param c Time offsets in Butcher table
				 
				 */
				
				template<typename RealType>
				void FillButcherTable(int stages, const Mat<mpq_rational>& a,
									  const Mat<mpq_rational> & b,
									  const Mat<mpq_rational> & c)
				{
					Mat<RealType>& aref = std::get< Mat<RealType> >(a_);
					aref.resize(stages, stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						for(int jj = 0; jj < s_; ++jj)
						{
							aref(ii,jj) = static_cast<RealType>(a(ii,jj));
						}
					}
					
					Vec<RealType>& bref = std::get< Vec<RealType> >(b_);
					bref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						bref(ii) = static_cast<RealType>(b(ii));
					}
					
					Vec<RealType>& cref = std::get< Vec<RealType> >(c_);
					cref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						cref(ii) = static_cast<RealType>(c(ii));
						
					}
					
				}

				
				
				
				
				
				///////////////////////////
				//
				// Private Data Members
				//
				////////////////////
				
				unsigned numTotalFunctions_; // Number of total functions for the current system
				unsigned numVariables_;  // Number of variables for the current system
				std::tuple< Mat<dbl>, Mat<mpfr> > K_;  // All the stage variables.  Each column represents a different stage.
				Predictor predictor_;  // Method for prediction
				unsigned p_;  //Order of the prediction method
				std::tuple< Mat<dbl>, Mat<mpfr> > dh_dx_0_;  // Jacobian for the initial stage.  Use for AMP testing
				std::tuple< Mat<dbl>, Mat<mpfr> > dh_dx_temp_;  // Temporary jacobian for all other stages
				std::tuple< Vec<dbl>, Vec<mpfr> > dh_dt_temp_;  // Temporary time derivative used for all stages
				std::tuple< Eigen::PartialPivLU<Mat<dbl>>, Eigen::PartialPivLU<Mat<mpfr>> > LU_0_;  // LU from the intial stage used for AMP testing
				
				
				
				// Butcher Table (notation from https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)
				unsigned s_; // Number of stages
				std::tuple< Mat<double>, Mat<mpfr_float> > a_;
				std::tuple< Vec<double>, Vec<mpfr_float> > b_;
				std::tuple< Vec<double>, Vec<mpfr_float> > b_minus_bstar_;
				std::tuple< Vec<double>, Vec<mpfr_float> > c_;
				
				
				
				
				
				
				
				
				// static const variables that store the Butcher table in mpq_rational form
				// Euler
				static const mpq_rational aEulerPtr_[];
				static const Eigen::Matrix<mpq_rational,1,1> aEuler_;
				static const mpq_rational bEulerPtr_[];
				static const Eigen::Matrix<mpq_rational,1,1> bEuler_;
				static const mpq_rational cEulerPtr_[];
				static const Eigen::Matrix<mpq_rational,1,1> cEuler_;

				// Heun-Euler
				static const mpq_rational aHeunEulerPtr_[];
				static const Eigen::Matrix<mpq_rational,2,2> aHeunEuler_;
				static const mpq_rational bHeunEulerPtr_[];
				static const Eigen::Matrix<mpq_rational,2,1> bHeunEuler_;
				static const mpq_rational b_minus_bstarHeunEulerPtr_[];
				static const Eigen::Matrix<mpq_rational,2,1> b_minus_bstarHeunEuler_;
				static const mpq_rational cHeunEulerPtr_[];
				static const Eigen::Matrix<mpq_rational,2,1> cHeunEuler_;

				// RK4
				static const mpq_rational aRK4Ptr_[];
				static const Eigen::Matrix<mpq_rational,4,4> aRK4_;
				static const mpq_rational bRK4Ptr_[];
				static const Eigen::Matrix<mpq_rational,4,1> bRK4_;
				static const mpq_rational cRK4Ptr_[];
				static const Eigen::Matrix<mpq_rational,4,1> cRK4_;

				// RKF45
				static const mpq_rational aRKF45Ptr_[];
				static const Eigen::Matrix<mpq_rational,6,6> aRKF45_;
				static const mpq_rational bRKF45Ptr_[];
				static const Eigen::Matrix<mpq_rational,6,1> bRKF45_;
				static const mpq_rational b_minus_bstarRKF45Ptr_[];
				static const Eigen::Matrix<mpq_rational,6,1> b_minus_bstarRKF45_;
				static const mpq_rational cRKF45Ptr_[];
				static const Eigen::Matrix<mpq_rational,6,1> cRKF45_;
				
				// RK Cash-Karp45
				static const mpq_rational aRKCK45Ptr_[];
				static const Eigen::Matrix<mpq_rational,6,6> aRKCK45_;
				static const mpq_rational bRKCK45Ptr_[];
				static const Eigen::Matrix<mpq_rational,6,1> bRKCK45_;
				static const mpq_rational b_minus_bstarRKCK45Ptr_[];
				static const Eigen::Matrix<mpq_rational,6,1> b_minus_bstarRKCK45_;
				static const mpq_rational cRKCK45Ptr_[];
				static const Eigen::Matrix<mpq_rational,6,1> cRKCK45_;

				// RK Dormand-Prince 56
				static const mpq_rational aRKDP56Ptr_[];
				static const Eigen::Matrix<mpq_rational,8,8> aRKDP56_;
				static const mpq_rational bRKDP56Ptr_[];
				static const Eigen::Matrix<mpq_rational,8,1> bRKDP56_;
				static const mpq_rational b_minus_bstarRKDP56Ptr_[];
				static const Eigen::Matrix<mpq_rational,8,1> b_minus_bstarRKDP56_;
				static const mpq_rational cRKDP56Ptr_[];
				static const Eigen::Matrix<mpq_rational,8,1> cRKDP56_;

				// RK Verner 67
				static const mpq_rational aRKV67Ptr_[];
				static const Eigen::Matrix<mpq_rational,10,10> aRKV67_;
				static const mpq_rational bRKV67Ptr_[];
				static const Eigen::Matrix<mpq_rational,10,1> bRKV67_;
				static const mpq_rational b_minus_bstarRKV67Ptr_[];
				static const Eigen::Matrix<mpq_rational,10,1> b_minus_bstarRKV67_;
				static const mpq_rational cRKV67Ptr_[];
				static const Eigen::Matrix<mpq_rational,10,1> cRKV67_;


				
			}; //re: ExplicitRKPredictor class
			
			
			

			
			
			
		} // re: predict
	}// re: tracking
}// re: bertini

#endif
