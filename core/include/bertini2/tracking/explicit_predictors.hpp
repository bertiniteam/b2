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
					case (Predictor::RK4):
						return 4;
					case (Predictor::RKF45):
						return 4;
					case (Predictor::RKCashKarp45):
						return 4;
					case (Predictor::RKDormandPrince56):
						return 5;
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
					case (Predictor::RK4):
						return false;
					case (Predictor::RKF45):
						return true;
					case (Predictor::RKCashKarp45):
						return true;
					case (Predictor::RKDormandPrince56):
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
							c_ = Vec<RealType>(s_); c_(0) = static_cast<RealType>(cEuler_(0));
							a_ = Mat<RealType>(s_,s_); a_(0,0) = static_cast<RealType>(aEuler_(0,0));
							b_ = Vec<RealType>(s_); b_(0) = static_cast<RealType>(bEuler_(0));
							break;
						}
						case Predictor::HeunEuler:
						{
							s_ = 2;
							
							FillButcherTable(s_, aHeunEuler_, bHeunEuler_, b_minus_bstarHeunEuler_, cHeunEuler_);
							
							break;
						}
						case Predictor::RK4:
						{
							s_ = 4;
							
							FillButcherTable(s_, aRK4_, bRK4_, cRK4_);
							
							break;
						}

						case Predictor::RKF45:
						{
							s_ = 6;
							
							FillButcherTable(s_, aRKF45_, bRKF45_, b_minus_bstarRKF45_, cRKF45_);
							
							break;
						}

						case Predictor::RKCashKarp45:
						{
							s_ = 6;
							
							FillButcherTable(s_, aRKCK45_, bRKCK45_, b_minus_bstarRKCK45_, cRKCK45_);
							
							break;
						}

						case Predictor::RKDormandPrince56:
						{
							s_ = 8;
							
							FillButcherTable(s_, aRKDP56_, bRKDP56_, b_minus_bstarRKDP56_, cRKDP56_);
							
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

					if(success_code != SuccessCode::Success)
						return success_code;
					
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
					
					if(success_code != SuccessCode::Success)
						return success_code;
					
					SetErrorEstimate(error_estimate, delta_t);
					
					
					
					return success_code;
				}

				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
				
			private:
				
				
				
				
				
				
//				void MethodSetup(Predictor method)
//				{
//					p_ = Order(method);
//					switch(method)
//					{
//						case Predictor::Euler:
//						{
//							s_ = 1;
//							c_ = Vec<RealType>(s_); c_(0) = 0;
//							a_ = Mat<RealType>(s_,s_); a_(0,0) = 0;
//							b_ = Vec<RealType>(s_); b_(0) = 1;
//							break;
//						}
//						case Predictor::HeunEuler:
//						{
//							mpq_rational half = mpq_rational(1,2);
//							s_ = 2;
//							c_ = Vec<RealType>(s_); c_ << RealType(0), RealType(1);
//							a_ = Mat<RealType>(s_, s_); a_ << RealType(0), RealType(0), RealType(1), RealType(0);
//							b_ = Vec<RealType>(s_); b_ << RealType(.5), RealType(.5);
//							bstar_ = Vec<RealType>(s_); bstar_ << RealType(1), RealType(0);
//							break;
//						}
//						default:
//						{
//							throw std::runtime_error("incompatible predictor choice in ExplicitPredict");
//						}
//					}
//				}
				
				
				

				
				
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
					auto numFuncs = K_.rows();
					Vec<ComplexType> err = Vec<ComplexType>(numFuncs);
					
					err.setZero();
					for(int ii = 0; ii < s_; ++ii)
					{
						err += (b_minus_bstar_(ii))*K_.col(ii);
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
				
				
				
				/**
				 /brief Fills the local embedded butcher table variables a,b,bstar and c with the constant static values stored in the class.
				 
				 \param stages Number of stages.  Used to create correct size on variables
				 \param a Matrix of the Butcher table
				 \param b Weights in the Butcher table
				 \param bstar Weights of embedded method in Butcher table
				 \param c Time offsets in Butcher table
				 
				 */
				
				template<typename MatDerived, typename VecDerived>
				void FillButcherTable(int stages, const Eigen::MatrixBase<MatDerived>& a,
								 const Eigen::MatrixBase<VecDerived> & b,
								 const Eigen::MatrixBase<VecDerived> & b_minus_bstar,
								 const Eigen::MatrixBase<VecDerived> & c)
				{
					a_ = Mat<RealType>(stages, stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						for(int jj = 0; jj < s_; ++jj)
						{
							a_(ii,jj) = static_cast<RealType>(a(ii,jj));
						}
					}
					
					b_ = Vec<RealType>(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						b_(ii) = static_cast<RealType>(b(ii));
					}
					
					b_minus_bstar_ = Vec<RealType>(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						b_minus_bstar_(ii) = static_cast<RealType>(b_minus_bstar(ii));
					}

					c_ = Vec<RealType>(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						c_(ii) = static_cast<RealType>(c(ii));
						
					}

				}
				
				

				
				
				/**
				 /brief Fills the local butcher table variables a,b and c with the constant static values stored in the class.
				 
				 \param stages Number of stages.  Used to create correct size on variables
				 \param a Matrix of the Butcher table
				 \param b Weights in the Butcher table
				 \param c Time offsets in Butcher table
				 
				 */
				
				template<typename MatDerived, typename VecDerived>
				void FillButcherTable(int stages, const Eigen::MatrixBase<MatDerived>& a,
									  const Eigen::MatrixBase<VecDerived> & b,
									  const Eigen::MatrixBase<VecDerived> & c)
				{
					a_ = Mat<RealType>(stages, stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						for(int jj = 0; jj < s_; ++jj)
						{
							a_(ii,jj) = static_cast<RealType>(a(ii,jj));
						}
					}
					
					b_ = Vec<RealType>(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						b_(ii) = static_cast<RealType>(b(ii));
					}
					
					c_ = Vec<RealType>(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						c_(ii) = static_cast<RealType>(c(ii));
						
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
				Vec<RealType> b_minus_bstar_;
				Mat<RealType> a_;
				int p_;
				
				
				
				
				
				
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

				// RK Dormand-Prince45
				static const mpq_rational aRKDP56Ptr_[];
				static const Eigen::Matrix<mpq_rational,8,8> aRKDP56_;
				static const mpq_rational bRKDP56Ptr_[];
				static const Eigen::Matrix<mpq_rational,8,1> bRKDP56_;
				static const mpq_rational b_minus_bstarRKDP56Ptr_[];
				static const Eigen::Matrix<mpq_rational,8,1> b_minus_bstarRKDP56_;
				static const mpq_rational cRKDP56Ptr_[];
				static const Eigen::Matrix<mpq_rational,8,1> cRKDP56_;


				
			};
			
			// Ptr must be filled row first.
			// i.e. aPtr[] = {a_11, a_21, a_12, a_22}
			
			/* Euler Butcher Table
			  |
			 0|0
			 -----
			  |1
			*/
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::aEulerPtr_[] = {mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,1,1> ExplicitRKPredictor<CType, RType>::aEuler_(aEulerPtr_);
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::bEulerPtr_[] = {mpq_rational(1,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,1,1> ExplicitRKPredictor<CType, RType>::bEuler_(bEulerPtr_);
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::cEulerPtr_[] = {mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,1,1> ExplicitRKPredictor<CType, RType>::cEuler_(cEulerPtr_);
			

			
			/* Heun-Euler Butcher Table
			 0 |
			 1 | 1
			 -----------
			   |1/2  1/2
			   | 1    0
			 */
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::aHeunEulerPtr_[] = {mpq_rational(0,1), mpq_rational(1,1),
																					mpq_rational(0,1), mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,2,2> ExplicitRKPredictor<CType, RType>::aHeunEuler_(aHeunEulerPtr_);
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::bHeunEulerPtr_[] = {mpq_rational(1,2), mpq_rational(1,2)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,2,1> ExplicitRKPredictor<CType, RType>::bHeunEuler_(bHeunEulerPtr_);
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::b_minus_bstarHeunEulerPtr_[] = {mpq_rational(-1,2), mpq_rational(1,2)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,2,1> ExplicitRKPredictor<CType, RType>::b_minus_bstarHeunEuler_(b_minus_bstarHeunEulerPtr_);
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::cHeunEulerPtr_[] = {mpq_rational(0,1), mpq_rational(1,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,2,1> ExplicitRKPredictor<CType, RType>::cHeunEuler_(cHeunEulerPtr_);

			
			
			/* RK4 Butcher Table
			 0   | 0   0  0  0
			 1/2 | 1/2 0  0  0
			 1/2 | 0  1/2 0  0
			 1   | 0   0  1  0
			 ------------------
				 | 1/6 1/3 1/3 1/6
			 */
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::aRK4Ptr_[] = {mpq_rational(0,1), mpq_rational(1,2),
				mpq_rational(0,1), mpq_rational(0,1),
			mpq_rational(0,1), mpq_rational(0,1), mpq_rational(1,2), mpq_rational(0,1),
			mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(1,1),
			mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,4,4> ExplicitRKPredictor<CType, RType>::aRK4_(aRK4Ptr_);
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::bRK4Ptr_[] = {mpq_rational(1,6), mpq_rational(1,3),
			mpq_rational(1,3),mpq_rational(1,6)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,4,1> ExplicitRKPredictor<CType, RType>::bRK4_(bRK4Ptr_);
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::cRK4Ptr_[] = {mpq_rational(0,1), mpq_rational(1,2),
			mpq_rational(1,2),mpq_rational(1,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,4,1> ExplicitRKPredictor<CType, RType>::cRK4_(cRK4Ptr_);


			
			
			/* RKF45 Butcher Table
			 https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
			 */
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::aRKF45Ptr_[] =
			{mpq_rational(0,1), mpq_rational(1,4),mpq_rational(3,32), mpq_rational(1932,2197),mpq_rational(439,216),mpq_rational(-8,27),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(9,32), mpq_rational(-7200,2197),mpq_rational(-8,1),mpq_rational(2,1),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(7296,2197),mpq_rational(3680,513),mpq_rational(-3544,2565),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(-845,4104),mpq_rational(1859,4104),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(0,1),mpq_rational(-11,40),
				
				mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,6,6> ExplicitRKPredictor<CType, RType>::aRKF45_(aRKF45Ptr_);
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::bRKF45Ptr_[] =
			{mpq_rational(16,135), mpq_rational(0,1), mpq_rational(6656,12825), mpq_rational(28561,56430),mpq_rational(-9,50),mpq_rational(2,55)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor<CType, RType>::bRKF45_(bRKF45Ptr_);

			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::b_minus_bstarRKF45Ptr_[] =
			{mpq_rational(1,360), mpq_rational(0,1), mpq_rational(-128,4275), mpq_rational(-2197,75240),mpq_rational(1,50),mpq_rational(2,55)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor<CType, RType>::b_minus_bstarRKF45_(b_minus_bstarRKF45Ptr_);

			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::cRKF45Ptr_[] =
			{mpq_rational(0,1), mpq_rational(1,4), mpq_rational(3,8), mpq_rational(12,13),mpq_rational(1,1),mpq_rational(1,2)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor<CType, RType>::cRKF45_(cRKF45Ptr_);


			
			
			/* RK Cash-Karp45 Butcher Table
			 https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
			 */
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::aRKCK45Ptr_[] =
			{mpq_rational(0,1), mpq_rational(1,5),mpq_rational(3,40), mpq_rational(3,10),mpq_rational(-11,54),mpq_rational(1631,55296),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(9,40), mpq_rational(-9,10),mpq_rational(5,2),mpq_rational(175,512),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(6,5),mpq_rational(-70,27),mpq_rational(575,13824),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(35,27),mpq_rational(44275,110592),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(0,1),mpq_rational(253,4096),
				
				mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,6,6> ExplicitRKPredictor<CType, RType>::aRKCK45_(aRKCK45Ptr_);
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::bRKCK45Ptr_[] =
			{mpq_rational(37,378), mpq_rational(0,1), mpq_rational(250,621), mpq_rational(125,594),mpq_rational(0,1),mpq_rational(512,1771)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor<CType, RType>::bRKCK45_(bRKCK45Ptr_);
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::b_minus_bstarRKCK45Ptr_[] =
			{mpq_rational(-277,64512), mpq_rational(0,1), mpq_rational(6925,370944), mpq_rational(-6925,202752),mpq_rational(-277,14336),mpq_rational(277,7084)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor<CType, RType>::b_minus_bstarRKCK45_(b_minus_bstarRKCK45Ptr_);
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::cRKCK45Ptr_[] =
			{mpq_rational(0,1), mpq_rational(1,5), mpq_rational(3,10), mpq_rational(3,5),mpq_rational(1,1),mpq_rational(7,8)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor<CType, RType>::cRKCK45_(cRKCK45Ptr_);

			
			
			
			/* RK Dormand-Prince56 Butcher Table
			 Prince and Dormand.  High order embedded Runge-Kutta formulae.  J. Comput. Appl. Math., 7(1):67-75, 1981
			 */
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::aRKDP56Ptr_[] =
			{mpq_rational(0,1), mpq_rational(1,10),mpq_rational(-2,81), mpq_rational(615,1372),mpq_rational(3243,5500),mpq_rational(-26492,37125),mpq_rational(5561,2376),mpq_rational(465467,266112),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(20,81), mpq_rational(-270,343),mpq_rational(-54,55),mpq_rational(72,55),mpq_rational(-35,11),mpq_rational(-2945,1232),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(1053,1372),mpq_rational(50949,71500),mpq_rational(2808,23375),mpq_rational(-24117,31603),mpq_rational(-5610201,14158144),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(4998,17875),mpq_rational(-24206,37125),mpq_rational(899983,200772),mpq_rational(10513573,3212352),
				
				mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(0,1),mpq_rational(338,459),mpq_rational(-5225,1836),mpq_rational(-424325,205632),
				
				mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(3925,4056),mpq_rational(376225,454272),
			
				mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),
			
				mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,8,8> ExplicitRKPredictor<CType, RType>::aRKDP56_(aRKDP56Ptr_);
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::bRKDP56Ptr_[] =
			{mpq_rational(821,10800), mpq_rational(0,1), mpq_rational(19683,71825), mpq_rational(175273,912600),mpq_rational(395,3672),mpq_rational(785,2704),mpq_rational(3,50),mpq_rational(0,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,8,1> ExplicitRKPredictor<CType, RType>::bRKDP56_(bRKDP56Ptr_);
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::b_minus_bstarRKDP56Ptr_[] =
			{mpq_rational(13,2400), mpq_rational(0,1), mpq_rational(-19683,618800), mpq_rational(2401,31200),mpq_rational(-65,816),mpq_rational(15,416),mpq_rational(521,5600),mpq_rational(-1,10)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,8,1> ExplicitRKPredictor<CType, RType>::b_minus_bstarRKDP56_(b_minus_bstarRKDP56Ptr_);
			
			template<typename CType, typename RType>
			const mpq_rational ExplicitRKPredictor<CType, RType>::cRKDP56Ptr_[] =
			{mpq_rational(0,1), mpq_rational(1,10), mpq_rational(2,9), mpq_rational(3,7),mpq_rational(3,5),mpq_rational(4,5),mpq_rational(1,1),mpq_rational(1,1)};
			template<typename CType, typename RType>
			const Eigen::Matrix<mpq_rational,8,1> ExplicitRKPredictor<CType, RType>::cRKDP56_(cRKDP56Ptr_);
			
			
			
			

			
			
			
		} // re: predict
	}// re: tracking
}// re: bertini

#endif
