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

#ifndef Xcode_b2_predictor_utility_hpp
#define Xcode_b2_predictor_utility_hpp

#include "tracking/amp_criteria.hpp"
#include "tracking/tracking_config.hpp"

#include "system.hpp"
#include "mpfr_extensions.hpp"
#include <Eigen/LU>

#include <boost/type_index.hpp>


namespace bertini{
	namespace tracking{
		namespace predict{
			
			
			using mpq_rational = bertini::mpq_rational;
			
			using Predictor = config::Predictor;
			
			
			
			
			
			
			template <typename ComplexType, typename RealType>
			class ExplicitPredictors
			{
			public:
				
				
				
				
				
				
				/*
				Butcher Table for Forward Euler
				c|    a
				0| 0
				-------
				 | 1  b
				 
				 */
				
				
				/*
				 Euler method
				 
				 To call: ExplicitPredictors::Euler( ... )
				 */
				static SuccessCode Euler(Vec<ComplexType> & next_space,
									System const& S,
									Vec<ComplexType> const& current_space, ComplexType current_time,
									ComplexType const& delta_t,
									RealType & condition_number_estimate,
									unsigned & num_steps_since_last_condition_number_computation,
									unsigned frequency_of_CN_estimation,
									RealType const& tracking_tolerance)
				{
					if(predictor_ != Predictor::Euler)
					{
						EulerSetup();
						predictor_ = Predictor::Euler;
					}
					
					K_ = Mat<ComplexType>(S.NumTotalFunctions(), s_);
					return FullStep(next_space, S, current_space, current_time, delta_t);
				}
				
				
				
				
				/* 
				 Euler with norms and size_proportion being calculated
				 
				 To call: ExplicitPredictors<Complex, Real>::Euler( ... )
				 */
				static SuccessCode Euler(Vec<ComplexType> & next_space,
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
					auto success_code = ExplicitPredictors<ComplexType, RealType>::Euler(next_space,
											  S,
											  current_space, current_time,
											  delta_t,
											  condition_number_estimate,
											  num_steps_since_last_condition_number_computation,
											  frequency_of_CN_estimation,
											  tracking_tolerance);
					
					SetNorms(norm_J, norm_J_inverse, S.NumVariables());
					size_proportion = K_.cols(0).array().abs().maxCoeff();
					
					return success_code;
				}
				
				
				
				

				
				
				
				/*
				 Butcher Table for Heun
				 c|     a
				 0| 0 0
				 1| 1 0
				 -------
				  | .5 .5  b
				  | 1  0   bstar
				 
				 */
				
				
				/*
				 Heun method
				 
				 To call: ExplicitPredictors::HeunEuler( ... )
				 */
				static SuccessCode HeunEuler(Vec<ComplexType> & next_space,
										 System const& S,
										 Vec<ComplexType> const& current_space, ComplexType current_time,
										 ComplexType const& delta_t,
										 RealType & condition_number_estimate,
										 unsigned & num_steps_since_last_condition_number_computation,
										 unsigned frequency_of_CN_estimation,
										 RealType const& tracking_tolerance)
				{
					if(predictor_ != Predictor::HeunEuler)
					{
						HeunSetup();
						predictor_ = Predictor::HeunEuler;
					}
					
					K_ = Mat<ComplexType>(S.NumTotalFunctions(), s_);
					return FullStep(next_space, S, current_space, current_time, delta_t);
				}
				
				
				
				
				/*
				 Heun with norms and size_proportion being calculated
				 
				 To call: ExplicitPredictors::HeunEuler( ... )
				 */
				static SuccessCode HeunEuler(Vec<ComplexType> & next_space,
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
					auto success_code = HeunEuler(next_space,
											  S,
											  current_space, current_time,
											  delta_t,
											  condition_number_estimate,
											  num_steps_since_last_condition_number_computation,
											  frequency_of_CN_estimation,
											  tracking_tolerance);
					
					SetNorms(norm_J, norm_J_inverse, S.NumVariables());
					SetSizeProportion(size_proportion, delta_t);
					
					return success_code;
				}
				
				
				/*
				 Heun with error estimate, norms and size_proportion being calculated
				 
				 To call: ExplicitPredictors::HeunEuler( ... )
				 */
				static SuccessCode HeunEuler(Vec<ComplexType> & next_space,
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
					auto success_code = HeunEuler(next_space, size_proportion, norm_J, norm_J_inverse, S,
											 current_space, current_time, delta_t, condition_number_estimate,
											 num_steps_since_last_condition_number_computation,
											 frequency_of_CN_estimation, tracking_tolerance, AMP_config);
					
					SetErrorEstimate(error_estimate, delta_t);
					
					return success_code;
				}


			private:
				
				// This is a static class, no instances should be created.
				ExplicitPredictors() = default;
				
				
				/** 
				 Setup Butcher table and order for Euler's method
				 */
				static void EulerSetup()
				{
					s_ = 1;
					c_ = Vec<RealType>(s_); c_(0) = 0;
					a_ = Mat<RealType>(s_,s_); a_(0,0) = 0;
					b_ = Vec<RealType>(s_); b_(0) = 1;
					p_ = 1;
				}

				
				
				
				/** 
				 Setup Butcher table and order for Heun's method
				 */
				static void HeunSetup()
				{
					s_ = 2;
					c_ = Vec<RealType>(s_); c_ << RealType(0), RealType(1);
					a_ = Mat<RealType>(s_, s_); a_ << RealType(0), RealType(0), RealType(1), RealType(0);
					b_ = Vec<RealType>(s_); b_ << RealType(.5), RealType(.5);
					bstar_ = Vec<RealType>(s_); bstar_ << RealType(1), RealType(0);
					p_ = 1;
				}

				
				
				static SuccessCode FullStep(Vec<ComplexType> & next_space,
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

				
				
				static SuccessCode SetErrorEstimate(RealType & error_estimate, ComplexType const& delta_t)
				{
					int numFuncs = K_.cols();
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
				
				

				
				static SuccessCode SetNorms(RealType & norm_J,RealType & norm_J_inverse, int NumVariables)
				{
					Vec<ComplexType> randy = RandomOfUnits<ComplexType>(NumVariables);
					Vec<ComplexType> temp_soln = LU_.solve(randy);
					
					norm_J = dh_dx_.norm();
					norm_J_inverse = temp_soln.norm();
					
					return SuccessCode::Success;
				};
				
				
				
				static SuccessCode SetSizeProportion(RealType & size_proportion, ComplexType const& delta_t)
				{
					RealType err_est;
					SetErrorEstimate(err_est, delta_t);
					
					size_proportion = err_est/(std::pow(abs(delta_t), p_+1));
					
					return SuccessCode::Success;
				};
				
				
				

				
				static SuccessCode EvalRHS(System const& S,
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
				static config::Predictor predictor_;
				
				
				static Mat<ComplexType> dh_dx_;
				static Eigen::PartialPivLU<Mat<ComplexType>> LU_;
				
				// Butcher Table (notation from https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)
				static int s_; // Number of stages
				static Vec<RealType> c_;
				static Vec<RealType> b_;
				static Vec<RealType> bstar_;
				static Mat<RealType> a_;
				static int p_;
				
				static Mat<ComplexType> K_;
			};
			template<typename Complex, typename Real>
			config::Predictor ExplicitPredictors<Complex,Real>::predictor_ = Predictor::None;
			
			
			template<typename Complex, typename Real>
			int ExplicitPredictors<Complex,Real>::s_;
			template<typename Complex, typename Real>
			Vec<Real> ExplicitPredictors<Complex,Real>::c_;
			template<typename Complex, typename Real>
			Vec<Real> ExplicitPredictors<Complex,Real>::b_;
			template<typename Complex, typename Real>
			Vec<Real> ExplicitPredictors<Complex,Real>::bstar_;
			template<typename Complex, typename Real>
			Mat<Real> ExplicitPredictors<Complex,Real>::a_;
			template<typename Complex, typename Real>
			int ExplicitPredictors<Complex,Real>::p_;
			
			template<typename Complex, typename Real>
			Mat<Complex> ExplicitPredictors<Complex,Real>::K_;
			template<typename Complex, typename Real>
			Mat<Complex> ExplicitPredictors<Complex,Real>::dh_dx_;
			template<typename Complex, typename Real>
			Eigen::PartialPivLU<Mat<Complex>> ExplicitPredictors<Complex,Real>::LU_;
			
			
			
			
			
			
			
			
			
			
			
			
//			template <typename ComplexType, typename RealType>
//			class EulerPred : public Predictor<ComplexType, RealType>
//			{
//				using Predictor<ComplexType, RealType>::s_;
//				using Predictor<ComplexType, RealType>::c_;
//				using Predictor<ComplexType, RealType>::b_;
//				using Predictor<ComplexType, RealType>::bstar_;
//				using Predictor<ComplexType, RealType>::a_;
//				using Predictor<ComplexType, RealType>::K_;
//				
////				using Predictor<ComplexType, RealType>::EvalRHS;
////				using Predictor<ComplexType, RealType>::ErrorEstimate;
////				using Predictor<ComplexType, RealType>::NormAndSize;
//				using Predictor<ComplexType, RealType>::FullStep;
//				
//			public:
//				EulerPred()
//				{
//					s_ = 1;
//					
//				}
//				
//				
//				
//				SuccessCode Predict(Vec<ComplexType> & next_space,
//									System const& S,
//									Vec<ComplexType> const& current_space, ComplexType current_time,
//									ComplexType const& delta_t,
//									RealType & condition_number_estimate,
//									unsigned & num_steps_since_last_condition_number_computation,
//									unsigned frequency_of_CN_estimation,
//									RealType const& tracking_tolerance)
//				{
//					K_ = Mat<ComplexType>(S.NumTotalFunctions(), s_);
//					return FullStep(next_space, S, current_space, current_time, delta_t);
//				}
//				
//				
//			};

			
			
			
		} // re: predict
	}// re: tracking
}// re: bertini

#endif
