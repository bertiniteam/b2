//This file is part of Bertini 2.0.
//
//bertini2/trackers/explicit_predictors.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/trackers/explicit_predictors.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/trackers/explicit_predictors.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  copyright 2015
//  James B. Collins
//  West Texas A&M University
//  Department of Mathematics
//  Spring 2016
//
// silviana amethyst, university of wisconsin-eau claire



/**
 \file explicit_predictors.hpp
 
 \brief Contains a base class for all ODE predictors.
 */

#ifndef BERTINI_EXPLICIT_PREDICTORS_HPP
#define BERTINI_EXPLICIT_PREDICTORS_HPP

#include "bertini2/trackers/amp_criteria.hpp"

#include "bertini2/system/system.hpp"
#include "bertini2/mpfr_extensions.hpp"
#include <Eigen/LU>

#include <boost/type_index.hpp>

#include "bertini2/eigen_extensions.hpp"
#include "bertini2/linalg/lu_solver.hpp"

namespace bertini{
	namespace tracking{
		namespace predict{
			
			
		// Constant,
		// Euler,
		// Heun,
		// RK4,
		// HeunEuler,
		// RKNorsett34,
		// RKF45,
		// RKCashKarp45,
		// RKDormandPrince56,
		// RKVerner67

			/**
			 \brief Get the Bertini2 default predictor.
			 
			 Currently set to RKF45.

			 \return The default predictor method to use.
			 */
			inline
			Predictor DefaultPredictor()
			{
				return Predictor::RKF45;
			}
			
			
			/**
			\brief The order of the predictor.  
			
			The order of the error estimate is this plus one.

			 \return The order of the predictor method.
			 \param predictor_choice The predictor method to query.
			 */
			inline
			unsigned Order(Predictor predictor_choice)
			{
				switch (predictor_choice)
				{
					case (Predictor::Constant):
						return 0;
					case (Predictor::Euler):
						return 1;
					case (Predictor::Heun):
						return 2;
					case (Predictor::HeunEuler):
						return 2;
					case (Predictor::RKNorsett34):
						return 3;
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
			
			/**
			\brief Ask whether a predictor method provides an error estimate.

			\return Yes or no, does it or does it not.
			\param predictor_choice The predictor method to query.
			*/
			inline bool HasErrorEstimate(Predictor predictor_choice)
			{
				switch (predictor_choice)
				{
					case (Predictor::Constant):
						return false;
					case (Predictor::Euler):
						return false;
					case (Predictor::Heun):
						return false;
					case (Predictor::HeunEuler):
						return true;
					case (Predictor::RKNorsett34):
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
			 \class ExplicitRKPredictor
			 
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
			 ExplicitRKPredictors<Complex,Real> euler(Predictor::Euler, sys)
			 success_code = euler.Predict( ... )
			 \endcode
			 */
			class ExplicitRKPredictor
			{
			public:
				
				/**
				\brief Construct a predictor to work on a system.

				\param S the system the predictor will be predicting on.
				*/
				ExplicitRKPredictor(const System& S) : s_(0), current_precision_(DefaultPrecision())
				{
					ChangeSystem(S);
					PredictorMethod(DefaultPredictor());
				}
				
				/**
				 \brief Constructor for a particular predictor method
				 
				 \param method The predictor method to be implemented.
				 \param S the system to be predicting on.
				 */
				ExplicitRKPredictor(Predictor method, const System& S) : s_(0), current_precision_(DefaultPrecision())
				{
					ChangeSystem(S);
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
						case Predictor::Constant:
						{
							s_ = 1;
							Mat<double>& arefd = std::get< Mat<double> >(a_);
							Vec<double>& brefd = std::get< Vec<double> >(b_);
							Vec<double>& crefd = std::get< Vec<double> >(c_);
							crefd.resize(s_); crefd(0) = 0;
							arefd.resize(s_,s_); arefd(0,0) = 0;
							brefd.resize(s_); brefd(0) = 0;
							Mat<real_mp>& arefmp = std::get< Mat<real_mp> >(a_);
							Vec<real_mp>& brefmp = std::get< Vec<real_mp> >(b_);
							Vec<real_mp>& crefmp = std::get< Vec<real_mp> >(c_);
							crefmp.resize(s_); crefmp(0) = 0;
							arefmp.resize(s_,s_); arefmp(0,0) = 0;
							brefmp.resize(s_); brefmp(0) = 0;
							uses_embedded_ = false;
							
							break;
						}
						case Predictor::Euler:
						{
							s_ = 1;
							Mat<double>& arefd = std::get< Mat<double> >(a_);
							Vec<double>& brefd = std::get< Vec<double> >(b_);
							Vec<double>& crefd = std::get< Vec<double> >(c_);
							crefd.resize(s_); crefd(0) = static_cast<double>(cEuler_(0));
							arefd.resize(s_,s_); arefd(0,0) = static_cast<double>(aEuler_(0,0));
							brefd.resize(s_); brefd(0) = static_cast<double>(bEuler_(0));
							Mat<real_mp>& arefmp = std::get< Mat<real_mp> >(a_);
							Vec<real_mp>& brefmp = std::get< Vec<real_mp> >(b_);
							Vec<real_mp>& crefmp = std::get< Vec<real_mp> >(c_);
							crefmp.resize(s_); crefmp(0) = static_cast<real_mp>(cEuler_(0));
							arefmp.resize(s_,s_); arefmp(0,0) = static_cast<real_mp>(aEuler_(0,0));
							brefmp.resize(s_); brefmp(0) = static_cast<real_mp>(bEuler_(0));
							uses_embedded_ = false;
							break;
						}
						case Predictor::HeunEuler:
						{
							s_ = 2;
							
							FillButcherTable<double>(static_cast<int>(s_),aHeunEuler_, bHeunEuler_, b_minus_bstarHeunEuler_, cHeunEuler_);
							FillButcherTable<real_mp>(static_cast<int>(s_),aHeunEuler_, bHeunEuler_, b_minus_bstarHeunEuler_, cHeunEuler_);
							
							break;
						}
						case Predictor::RK4:
						{
							s_ = 4;
							
							FillButcherTable<double>(static_cast<int>(s_),aRK4_, bRK4_, cRK4_);
							FillButcherTable<real_mp>(static_cast<int>(s_),aRK4_, bRK4_, cRK4_);
							
							break;
						}
							
						case Predictor::RKF45:
						{
							s_ = 6;
							
							FillButcherTable<double>(static_cast<int>(s_),aRKF45_, bRKF45_, b_minus_bstarRKF45_, cRKF45_);
							FillButcherTable<real_mp>(static_cast<int>(s_),aRKF45_, bRKF45_, b_minus_bstarRKF45_, cRKF45_);
							
							break;
						}
							
						case Predictor::RKCashKarp45:
						{
							s_ = 6;
							
							FillButcherTable<double>(static_cast<int>(s_),aRKCK45_, bRKCK45_, b_minus_bstarRKCK45_, cRKCK45_);
							FillButcherTable<real_mp>(static_cast<int>(s_),aRKCK45_, bRKCK45_, b_minus_bstarRKCK45_, cRKCK45_);
							
							break;
						}
							
						case Predictor::RKDormandPrince56:
						{
							s_ = 8;
							
							FillButcherTable<double>(static_cast<int>(s_),aRKDP56_, bRKDP56_, b_minus_bstarRKDP56_, cRKDP56_);
							FillButcherTable<real_mp>(static_cast<int>(s_),aRKDP56_, bRKDP56_, b_minus_bstarRKDP56_, cRKDP56_);
							
							break;
						}
							
						case Predictor::RKVerner67:
						{
							s_ = 10;
							
							FillButcherTable<double>(static_cast<int>(s_),aRKV67_, bRKV67_, b_minus_bstarRKV67_, cRKV67_);
							FillButcherTable<real_mp>(static_cast<int>(s_),aRKV67_, bRKV67_, b_minus_bstarRKV67_, cRKV67_);
							
							break;
						}
							
						default:
						{
							throw std::runtime_error("incompatible predictor choice in ExplicitPredict");
						}
					}
					ResizeK();
				}; // re: PredictorMethod
				
				
				
				
				
				/**
				 \brief Change the system (number of total functions) that the predictor uses.
				 
				 \param S New system to switch to.
				 */
				void ChangeSystem(const System& S)
				{
					numTotalFunctions_ = static_cast<unsigned>(S.NumTotalFunctions());
					numVariables_ = static_cast<unsigned>(S.NumVariables());
					// you cannot set K_ here, because s_ may not have been set
					std::get< Mat<complex_dbl> >(dh_dx_0_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<complex_mp> >(dh_dx_0_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<complex_dbl> >(dh_dx_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<complex_mp> >(dh_dx_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Vec<complex_dbl> >(dh_dt_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_mp> >(dh_dt_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_dbl> >(step_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_mp> >(step_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_dbl> >(rand_temp_) = RandomOfUnits<complex_dbl>(numVariables_);
					std::get< Vec<complex_mp> >(rand_temp_) = RandomOfUnits<complex_mp>(numVariables_);
					std::get< Vec<complex_dbl> >(solve_temp_).resize(numVariables_);
					std::get< Vec<complex_mp> >(solve_temp_).resize(numVariables_);

					std::get< Vec<complex_dbl> >(stage_pt_temp_).resize(numVariables_);
					std::get< Vec<complex_mp> >(stage_pt_temp_).resize(numVariables_);
					std::get< Vec<complex_dbl> >(err_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_mp> >(err_temp_).resize(numTotalFunctions_);

					ResizeK();
				}
				
				
				/// \brief Resize the internal Runge-Kutta stage matrix K to match the current system and stage count.
				void ResizeK()
				{
					std::get< Mat<complex_dbl> >(K_).resize(numTotalFunctions_, s_);
					std::get< Mat<complex_mp> >(K_).resize(numTotalFunctions_, s_);

					std::get< linalg::PartialPivLU<complex_dbl> >(LU_).ChangeSize(numVariables_);
					std::get< linalg::PartialPivLU<complex_mp> >(LU_).ChangeSize(numVariables_);
				}


				/**
				 \brief Adopt an externally-owned condition-number probe direction (both precisions).

				 The tracker calls this once per path so the predictor estimates ||J^{-1}|| against the
				 SAME random direction as the corrector (ADR-0024: one tracker-owned, per-path probe).
				 Standalone use (e.g. unit tests with no tracker) keeps the per-system probe drawn in
				 ChangeSystem.
				 */
				void SetConditionProbe(std::tuple< Vec<complex_dbl>, Vec<complex_mp> > const& probe)
				{
					rand_temp_ = probe;
				}
				
				
				/**
				\brief get the current precision of the predictor
				*/
				unsigned precision() const
				{
					return current_precision_;
				}

				/** 
				 /brief Change the precision of the predictor variables and reassign the Butcher table variables.
				 
				 \param new_precision The new precision.
				 
				 */
				void ChangePrecision(unsigned new_precision)
				{
					Precision(std::get< Mat<complex_mp> >(K_),new_precision);

					Precision(std::get< Vec<complex_mp> >(dh_dt_temp_),new_precision);
					Precision(std::get< Mat<complex_mp> >(dh_dx_0_),new_precision);
					Precision(std::get< Mat<complex_mp> >(dh_dx_temp_),new_precision);
					Precision(std::get< Vec<complex_mp> >(step_temp_),new_precision);
					Precision(std::get< Vec<complex_mp> >(rand_temp_),new_precision);
					Precision(std::get< Vec<complex_mp> >(solve_temp_),new_precision);
					Precision(std::get< Vec<complex_mp> >(stage_pt_temp_),new_precision);
					Precision(std::get< Vec<complex_mp> >(err_temp_),new_precision);
					std::get< linalg::PartialPivLU<complex_mp> >(LU_).ChangePrecision(new_precision);

					Precision(std::get< Mat<real_mp> >(a_),new_precision);
					Precision(std::get< Vec<real_mp> >(b_),new_precision);
					Precision(std::get< Vec<real_mp> >(b_minus_bstar_),new_precision);
					Precision(std::get< Vec<real_mp> >(c_),new_precision);

					PredictorMethod(predictor_);

					current_precision_ = new_precision;

					PrecisionSanityCheck();
				}
				
				/// \brief Assert (in debug builds) that the predictor's state is all at the expected precision.
				void PrecisionSanityCheck() const
				{
#ifndef NDEBUG
					// ThreadPrecision: correct when running on a std::thread worker,
					// where precision is set via SetThreadPrecision (thread-local only).
					assert(current_precision_==ThreadPrecision());

					Vec<complex_mp>& dhdttemp = std::get< Vec<complex_mp> >(dh_dt_temp_);
					Mat<complex_mp>& dhdx0 = std::get< Mat<complex_mp> >(dh_dx_0_); 
					Mat<complex_mp>& dhdxtemp = std::get< Mat<complex_mp> >(dh_dx_temp_); 

					Mat<real_mp>& a = std::get< Mat<real_mp> >(a_); 
					Vec<real_mp>& b = std::get< Vec<real_mp> >(b_); 
					Vec<real_mp>& bstar = std::get< Vec<real_mp> >(b_minus_bstar_); 
					Vec<real_mp>& c = std::get< Vec<real_mp> >(c_);



					assert(Precision(dhdttemp)==current_precision_);
					assert(Precision(dhdx0)==current_precision_);
					assert(Precision(dhdxtemp)==current_precision_);

					assert(Precision(a)==current_precision_);
					assert(Precision(b)==current_precision_);
					if (uses_embedded_)
						assert(Precision(bstar)==current_precision_);
					assert(Precision(c)==current_precision_);
#endif
				}
				
				
				/**
				 \brief Perform a generic predictor step.
				 
				 \param next_space The computed prediction.
				 \param meta Step metadata, populated during the step with the Jacobian norms and condition number estimate.
				 \param S The system being solved.
				 \param current_space The current space variable vector.
				 \param current_time The current time.
				 \param delta_t The size of the time step.
				 \param num_steps_since_last_condition_number_computation Updated in this function.
				 \param frequency_of_CN_estimation How many steps to take between condition number estimates.
				 \param tracking_tolerance How tightly to track the path.
				 \param AMP_config Optional adaptive-multiple-precision configuration; when null, fixed-precision behaviour is used.

				 \return SuccessCode indicating how the prediction went.
				 */
				
				template<typename ComplexT>
				SuccessCode Predict(Vec<ComplexT> & next_space,
									StepMetadata & meta,
									System const& S,
									const Vec<ComplexT>& current_space, ComplexT current_time,
									ComplexT const& delta_t,
									unsigned & num_steps_since_last_condition_number_computation,
									unsigned frequency_of_CN_estimation,
									NumErrorT const& tracking_tolerance,
									AdaptiveMultiplePrecisionConfig const* AMP_config = nullptr)
				{
					auto step_success = FullStep(next_space, S, current_space, current_time, delta_t);

					// Condition estimate (norm_J, norm_J_inverse, condition_number_estimate) is always
					// computed -- it is reported and used by the fixed-precision path too.
					SetNormsCond<ComplexT>(meta.norm_J, meta.norm_J_inverse, meta.condition_number_estimate,
					                       num_steps_since_last_condition_number_computation, frequency_of_CN_estimation);

					// Fixed-precision (nullptr): just the prediction + condition estimate, no AMP criteria.
					if (AMP_config == nullptr)
						return step_success;

					SetSizeProportion<ComplexT>(meta.size_proportion, delta_t);
					if (predict::HasErrorEstimate(predictor_))
						SetErrorEstimate<ComplexT>(meta.error_estimate, delta_t);

					if (step_success != SuccessCode::Success)
						return step_success;

					if (!amp::CriterionA<ComplexT>(meta.norm_J, meta.norm_J_inverse, *AMP_config))
						return SuccessCode::HigherPrecisionNecessary;
					if (!amp::CriterionC<ComplexT>(meta.norm_J_inverse, current_space, tracking_tolerance, *AMP_config))
						return SuccessCode::HigherPrecisionNecessary;

					return step_success;
				}
				
				
				
				
				
				
				
				
				/**
				\brief Get the currently used prediction method.

				\return The current method.
				*/
				Predictor PredictorMethod()
				{
					return predictor_;
				}
				
				
				
				
				/**
				\brief Get the order of the currently used prediction method.

				This is the lowest order of the predictor.  The order of the error estimate is this plus one.

				\return The aforementioned order.
				 */
				inline
				unsigned Order()
				{
					return p_;
				}
				
				/**
				\brief Get whether the current prediction method provides an error estimate.

				\return Yes or no.
				*/
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
				
				template<typename ComplexT>
				SuccessCode FullStep(Vec<ComplexT> & next_space,
									System const& S,
									 Vec<ComplexT> const& current_space, ComplexT const& current_time,
									 ComplexT const& delta_t)
				{
					
					// If using constant predictor
					if(s_ == 0)
					{
						next_space = current_space;
						return SuccessCode::Success;
					}
					
					using RealT = typename Eigen::NumTraits<ComplexT>::Real;

					Mat<ComplexT>& Kref = std::get< Mat<ComplexT> >(K_);
					Mat<RealT>& aref = std::get< Mat<RealT> >(a_);
					Vec<RealT>& bref = std::get< Vec<RealT> >(b_);
					Vec<RealT>& cref = std::get< Vec<RealT> >(c_);
					Kref.fill(ComplexT(0));
					Vec<ComplexT>& temp = std::get< Vec<ComplexT> >(step_temp_);

					if(EvalRHS(S, current_space, current_time, Kref, 0) != SuccessCode::Success)
					{
						return SuccessCode::MatrixSolveFailureFirstPartOfPrediction;
					}
					
					Vec<ComplexT>& stage_pt = std::get< Vec<ComplexT> >(stage_pt_temp_);
					for(unsigned ii = 1; ii < s_; ++ii)
					{
						temp.setZero(); // see https://github.com/bertiniteam/b2/issues/198
						for(unsigned jj = 0; jj < ii; ++jj)
							temp += aref(ii,jj)*Kref.col(jj);

						// Evaluate into the preallocated stage-point scratch rather than passing the
						// expression (current_space + delta_t*temp) to EvalRHS's const Vec& param, which
						// would materialize a fresh temporary Vec every stage.
						stage_pt.noalias() = current_space + delta_t*temp;
						if(EvalRHS<ComplexT>(S, stage_pt, current_time + cref(ii)*delta_t, Kref, ii) != SuccessCode::Success)
							return SuccessCode::MatrixSolveFailure;
					}
					
					
					temp.setZero();
					for(unsigned ii = 0; ii < s_; ++ii)
						temp += bref(ii)*Kref.col(ii);

					// next_space is a distinct buffer from current_space/temp, so noalias avoids a
					// materialized temporary for the axpy.
					next_space.noalias() = current_space + delta_t*temp;
					
					return SuccessCode::Success;
				};

				
				/// \brief Compute and (when due) refresh the Jacobian norms and condition-number estimate.
				/// \tparam ComplexT The complex number type to compute at.
				/// \param[out] norm_J Set to ||J||.
				/// \param[out] norm_J_inverse Set to the estimate of ||J^{-1}||.
				/// \param[out] condition_number_estimate Set to the product of the two norms.
				/// \param num_steps_since_last_condition_number_computation Steps elapsed since the last estimate.
				/// \param frequency_of_CN_estimation Recompute the estimate once this many steps have passed.
				template<typename ComplexT>
				void SetNormsCond(NumErrorT & norm_J, NumErrorT & norm_J_inverse, NumErrorT & condition_number_estimate, unsigned num_steps_since_last_condition_number_computation, unsigned frequency_of_CN_estimation)
				{
					// Calculate condition number and update if needed
					linalg::PartialPivLU<ComplexT>& LUref = std::get< linalg::PartialPivLU<ComplexT> >(LU_);
					Mat<ComplexT>& dhdxref = std::get< Mat<ComplexT> >(dh_dx_0_);

					Vec<ComplexT> const& randy = std::get< Vec<ComplexT> >(rand_temp_);
					Vec<ComplexT>& solve_ref = std::get< Vec<ComplexT> >(solve_temp_);
					LUref.Solve(randy, solve_ref);

					norm_J = NumErrorT(dhdxref.norm());
					norm_J_inverse = NumErrorT(solve_ref.norm());


					if (num_steps_since_last_condition_number_computation >= frequency_of_CN_estimation)
					{
						condition_number_estimate = NumErrorT(norm_J * norm_J_inverse);
						num_steps_since_last_condition_number_computation = 1; // reset the counter to 1
					}
					else // no need to compute the condition number
						num_steps_since_last_condition_number_computation++;
				}
				
				
				/**
				 \brief Computes the error estimate of this prediction step.
				 
				 \param error_estimate Computed error estimate
				 \param delta_t The time step
				 
				 \return Success code or the computation
				 
				 */
				
				template<typename ComplexT>
				SuccessCode SetErrorEstimate(NumErrorT & error_estimate, ComplexT const& delta_t)
				{
					using RealT = typename Eigen::NumTraits<ComplexT>::Real;

					Mat<ComplexT>& Kref = std::get< Mat<ComplexT> >(K_);
					Vec<RealT>& b_minus_bstar_ref = std::get< Vec<RealT> >(b_minus_bstar_);
					
					Vec<ComplexT>& err = std::get< Vec<ComplexT> >(err_temp_);  // reused scratch, not a fresh per-step alloc

					err.setZero();
					for(unsigned ii = 0; ii < s_; ++ii)
					{
						err += (b_minus_bstar_ref(ii))*Kref.col(ii);
					}

					err *= delta_t;
					
					error_estimate = NumErrorT(err.norm());
					
					return SuccessCode::Success;
				};
				
				

				
				
				
				/**
				 \brief Compute the size proportion variable for AMP computation
				 
				 \param size_proportion Computed size proportion
				 \param delta_t The time step
				 
				 \return Success code of the computation
				 
				 */
				
				template<typename ComplexT>
				SuccessCode SetSizeProportion(NumErrorT & size_proportion, ComplexT const& delta_t)
				{
					if(predict::HasErrorEstimate(predictor_))
					{
						NumErrorT err_est;
						SetErrorEstimate(err_est, delta_t);
						
						using std::pow;
						size_proportion = err_est/NumErrorT(pow(abs(delta_t), p_+1));
						
						return SuccessCode::Success;
					}
					else
					{
						// No embedded error estimate (e.g. Euler, RK4).  AMP2 (bhswAMP2, Eqs 9-10)
						// derives the size proportion $a$ from the prediction step itself, written as
						// the initial Newton residual: ||d|| = a|s|.  The predicted step in z is
						// delta_t*(sum_i b_i K_i), so ||d|| = |delta_t|*||K|| and therefore
						//     a = ||d|| / |s| = ||K|| ~ maxCoeff(|K|),
						// with NO further division by |delta_t|.  The previous
						// maxCoeff(K)/|delta_t|^p over-divided by the step: it inflated $a$ like
						// 1/|delta_t|^p as the step shrank, spuriously escalating AMP precision
						// (DigitsB) on small steps.  $a$ is meant to be an O(1), step-independent
						// proportionality constant -- see SetSizeProportion's error-estimate branch,
						// which divides err_est by |delta_t|^(p+1) for exactly the same reason
						// (AMP3 / bhsODEAMP Eq. 6).
						Mat<ComplexT>& Kref = std::get< Mat<ComplexT> >(K_);
						size_proportion = NumErrorT(Kref.array().abs().maxCoeff());
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
				
				template<typename ComplexT>
				SuccessCode EvalRHS(System const& S,
									const Vec<ComplexT>& space, const ComplexT& time, Mat<ComplexT> & K, unsigned stage)
				{


					if (std::is_same<ComplexT, complex_mp>::value)
						PrecisionSanityCheck();

					if(stage == 0)
					{
						linalg::PartialPivLU<ComplexT>& LUref = std::get< linalg::PartialPivLU<ComplexT> >(LU_);
						Mat<ComplexT>& dhdxref = std::get< Mat<ComplexT> >(dh_dx_0_);

						if (!std::is_same<ComplexT,complex_dbl>::value)
						{
							assert(ThreadPrecision()==current_precision_);

							assert(Precision(space)==current_precision_);
							assert(Precision(time)==current_precision_);
							assert(Precision(dhdxref)==current_precision_);
							assert(Precision(K)==current_precision_);
						}
						S.SetAndReset<ComplexT>(space, time);
						S.JacobianInPlace(dhdxref);
						// Factor a copy of dh/dx (dh_dx_0_ is read again in SetNormsCond); health check folded in.
						auto lu_code = LUref.Factor(dhdxref);
						if (!std::is_same<ComplexT,complex_dbl>::value)
						{
							assert(Precision(dhdxref)==current_precision_);
							assert(Precision(LUref.Factors())==current_precision_);
						}

						if (lu_code!=MatrixSuccessCode::Success)
							return SuccessCode::MatrixSolveFailureFirstPartOfPrediction;

						Vec<ComplexT>& dhdtref = std::get< Vec<ComplexT> >(dh_dt_temp_);
						S.TimeDerivativeInPlace(dhdtref);
						dhdtref = -dhdtref;                     // in place; Solve needs materialized rhs
						LUref.Solve(dhdtref, K.col(stage));
						
						return SuccessCode::Success;
						
					}
					else
					{
						S.SetAndReset<ComplexT>(space, time);

						Mat<ComplexT>& dhdxtempref = std::get< Mat<ComplexT> >(dh_dx_temp_);
						S.JacobianInPlace(dhdxtempref);
						linalg::PartialPivLU<ComplexT>& LU_temp = std::get< linalg::PartialPivLU<ComplexT> >(LU_);
						// dh_dx_temp_ is pure scratch here -> factor destructively, skipping the copy.
						if (LU_temp.FactorDestructive(dhdxtempref)!=MatrixSuccessCode::Success)
							return SuccessCode::MatrixSolveFailure;

						Vec<ComplexT>& dhdtref = std::get< Vec<ComplexT> >(dh_dt_temp_);
						S.TimeDerivativeInPlace(dhdtref);
						dhdtref = -dhdtref;                     // in place; Solve needs materialized rhs
						LU_temp.Solve(dhdtref, K.col(stage));
						
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
				
				template<typename RealT>
				void FillButcherTable(int stages, const Mat<mpq_rational>& a,
								 const Mat<mpq_rational> & b,
								 const Mat<mpq_rational> & b_minus_bstar,
								 const Mat<mpq_rational> & c)
				{
					Mat<RealT>& aref = std::get< Mat<RealT> >(a_);
					aref.resize(stages, stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						for(unsigned jj = 0; jj < s_; ++jj)
						{
							aref(ii,jj) = static_cast<RealT>(a(ii,jj));
						}
					}
					
					Vec<RealT>& bref = std::get< Vec<RealT> >(b_);
					bref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						bref(ii) = static_cast<RealT>(b(ii));
					}
					
					Vec<RealT>& b_minus_bstar_ref = std::get< Vec<RealT> >(b_minus_bstar_);
					b_minus_bstar_ref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						b_minus_bstar_ref(ii) = static_cast<RealT>(b_minus_bstar(ii));
					}

					Vec<RealT>& cref = std::get< Vec<RealT> >(c_);
					cref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						cref(ii) = static_cast<RealT>(c(ii));
						
					}
					uses_embedded_ = true;
				}
				
				

				
				
				/**
				 /brief Fills the local butcher table variables a,b and c with the constant static values stored in the class.
				 
				 \param stages Number of stages.  Used to create correct size on variables
				 \param a Matrix of the Butcher table
				 \param b Weights in the Butcher table
				 \param c Time offsets in Butcher table
				 
				 */
				
				template<typename RealT>
				void FillButcherTable(int stages, const Mat<mpq_rational>& a,
									  const Mat<mpq_rational> & b,
									  const Mat<mpq_rational> & c)
				{
					Mat<RealT>& aref = std::get< Mat<RealT> >(a_);
					aref.resize(stages, stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						for(unsigned jj = 0; jj < s_; ++jj)
						{
							aref(ii,jj) = static_cast<RealT>(a(ii,jj));
						}
					}
					
					Vec<RealT>& bref = std::get< Vec<RealT> >(b_);
					bref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						bref(ii) = static_cast<RealT>(b(ii));
					}
					
					Vec<RealT>& cref = std::get< Vec<RealT> >(c_);
					cref.resize(stages);
					for(int ii = 0; ii < stages; ++ii)
					{
						cref(ii) = static_cast<RealT>(c(ii));
						
					}
					uses_embedded_ = false;
				}

				
				
				
				
				
				///////////////////////////
				//
				// Private Data Members
				//
				////////////////////
				
				unsigned numTotalFunctions_; // Number of total functions for the current system
				unsigned numVariables_;  // Number of variables for the current system
				mutable std::tuple< Mat<complex_dbl>, Mat<complex_mp> > K_;  // All the stage variables.  Each column represents a different stage.
				Predictor predictor_;  // Method for prediction
				unsigned p_;  //Order of the prediction method
				mutable std::tuple< Mat<complex_dbl>, Mat<complex_mp> > dh_dx_0_;  // Jacobian for the initial stage.  Use for AMP testing
				mutable std::tuple< Mat<complex_dbl>, Mat<complex_mp> > dh_dx_temp_;  // Temporary jacobian for all other stages
				mutable std::tuple< Vec<complex_dbl>, Vec<complex_mp> > dh_dt_temp_;  // Temporary time derivative used for all stages
				// std::tuple< Eigen::PartialPivLU<Mat<complex_dbl>>, Eigen::PartialPivLU<Mat<complex_mp>> > LU_0_;  // LU from the intial stage used for AMP testing

				mutable std::tuple< linalg::PartialPivLU<complex_dbl>, linalg::PartialPivLU<complex_mp> > LU_;

				mutable std::tuple< Vec<complex_dbl>, Vec<complex_mp> > step_temp_;  // reused scratch for FullStep stage accumulation
				mutable std::tuple< Vec<complex_dbl>, Vec<complex_mp> > rand_temp_;  // reused scratch: random RHS for norm_J_inverse
				mutable std::tuple< Vec<complex_dbl>, Vec<complex_mp> > solve_temp_; // reused scratch: LU solve result
				mutable std::tuple< Vec<complex_dbl>, Vec<complex_mp> > stage_pt_temp_; // reused scratch: RK stage point (current_space + delta_t*temp)
				mutable std::tuple< Vec<complex_dbl>, Vec<complex_mp> > err_temp_;   // reused scratch: embedded-RK error estimate vector
				
				
				// Butcher Table (notation from https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods )
				mutable unsigned s_; // Number of stages
				mutable std::tuple< Mat<double>, Mat<real_mp> > a_;
				mutable std::tuple< Vec<double>, Vec<real_mp> > b_;
				mutable std::tuple< Vec<double>, Vec<real_mp> > b_minus_bstar_;
				mutable std::tuple< Vec<double>, Vec<real_mp> > c_;
				
				mutable bool uses_embedded_;
				mutable unsigned current_precision_;
				
				
				
				
				
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
			
			
			

			
			
			
		/// \cond INTERNAL
		// Explicit instantiation declarations — suppress re-instantiation in every
		// including TU.  The definitions live in core/src/tracking/explicit_predictors.cpp.
		// Concrete types: complex_dbl = std::complex<double>, complex_mp (multiprecision).
		// NumErrorT = double (from bertini2/common/config.hpp).

		extern template SuccessCode ExplicitRKPredictor::Predict<complex_dbl>(
		    Vec<complex_dbl>&, StepMetadata&, System const&, Vec<complex_dbl> const&, complex_dbl, complex_dbl const&,
		    unsigned&, unsigned, NumErrorT const&, AdaptiveMultiplePrecisionConfig const*);
		extern template SuccessCode ExplicitRKPredictor::Predict<complex_mp>(
		    Vec<complex_mp>&, StepMetadata&, System const&, Vec<complex_mp> const&, complex_mp, complex_mp const&,
		    unsigned&, unsigned, NumErrorT const&, AdaptiveMultiplePrecisionConfig const*);

		extern template SuccessCode ExplicitRKPredictor::FullStep<complex_dbl>(
		    Vec<complex_dbl>&, System const&, Vec<complex_dbl> const&, complex_dbl const&, complex_dbl const&);
		extern template SuccessCode ExplicitRKPredictor::FullStep<complex_mp>(
		    Vec<complex_mp>&, System const&, Vec<complex_mp> const&, complex_mp const&, complex_mp const&);

		extern template void ExplicitRKPredictor::SetNormsCond<complex_dbl>(
		    double&, double&, double&, unsigned, unsigned);
		extern template void ExplicitRKPredictor::SetNormsCond<complex_mp>(
		    double&, double&, double&, unsigned, unsigned);

		extern template SuccessCode ExplicitRKPredictor::SetErrorEstimate<complex_dbl>(double&, complex_dbl const&);
		extern template SuccessCode ExplicitRKPredictor::SetErrorEstimate<complex_mp>(double&, complex_mp const&);

		extern template SuccessCode ExplicitRKPredictor::SetSizeProportion<complex_dbl>(double&, complex_dbl const&);
		extern template SuccessCode ExplicitRKPredictor::SetSizeProportion<complex_mp>(double&, complex_mp const&);

		extern template SuccessCode ExplicitRKPredictor::EvalRHS<complex_dbl>(
		    System const&, Vec<complex_dbl> const&, complex_dbl const&, Mat<complex_dbl>&, unsigned);
		extern template SuccessCode ExplicitRKPredictor::EvalRHS<complex_mp>(
		    System const&, Vec<complex_mp> const&, complex_mp const&, Mat<complex_mp>&, unsigned);

		extern template void ExplicitRKPredictor::FillButcherTable<double>(
		    int, Mat<mpq_rational> const&, Mat<mpq_rational> const&,
		    Mat<mpq_rational> const&, Mat<mpq_rational> const&);
		extern template void ExplicitRKPredictor::FillButcherTable<real_mp>(
		    int, Mat<mpq_rational> const&, Mat<mpq_rational> const&,
		    Mat<mpq_rational> const&, Mat<mpq_rational> const&);

		extern template void ExplicitRKPredictor::FillButcherTable<double>(
		    int, Mat<mpq_rational> const&, Mat<mpq_rational> const&, Mat<mpq_rational> const&);
		extern template void ExplicitRKPredictor::FillButcherTable<real_mp>(
		    int, Mat<mpq_rational> const&, Mat<mpq_rational> const&, Mat<mpq_rational> const&);
		/// \endcond

		} // re: predict
	}// re: tracking
}// re: bertini

#endif
