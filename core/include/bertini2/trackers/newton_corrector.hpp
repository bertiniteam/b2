//
//  newton_corrector.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 4/27/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef BERTINI_NEWTON_CORRECTOR_HPP
#define BERTINI_NEWTON_CORRECTOR_HPP

#include "bertini2/trackers/amp_criteria.hpp"
#include "bertini2/trackers/config.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/linalg/lu_solver.hpp"


namespace bertini{
	namespace tracking{
		namespace correct{


			/**
			 /class NewtonCorrector
			 
			 \brief A command class which performs a newton correction step.
			 
			 ## Purpose
			 
			 Stores information computed during implementation of the method.
			 
			 
			 ## Use
			 To perform a correction step, you must instantiate an object and call Correct:
			 
			 \code
			 NewtonCorrect<Complex,Real> newton(sys)
			 success_code = newton.Correct( ... )
			 \endcode
			 
			 
			 
			 */

			class NewtonCorrector
			{
			public:
				
				
				/// \brief Construct a Newton corrector for a system, at the current default precision.
				NewtonCorrector(const System& S) : current_precision_(DefaultPrecision())
				{
					ChangeSystem(S);
				}

				
				
				
				
				
				
				/**
				 \brief Sets the Newton configuration settings.
				 
				 \param newton_settings A Newton struct holding configuration settings.
				 
				 */
				
				void Settings(const NewtonConfig& newton_settings)
				{
					newton_config_ = newton_settings;
				}
				
				
	
				
				
				/**
				 \brief Change the precision of the predictor variables and reassign the Butcher table variables.
				 
				 \param new_precision The new precision.
				 
				 */
				void ChangePrecision(unsigned new_precision)
				{
					Precision(std::get< Vec<complex_mp> >(f_temp_), new_precision);
					Precision(std::get< Vec<complex_mp> >(step_temp_), new_precision);
					Precision(std::get< Mat<complex_mp> >(J_temp_), new_precision);
					Precision(std::get< Vec<complex_mp> >(rand_temp_), new_precision);
					Precision(std::get< Vec<complex_mp> >(solve_temp_), new_precision);
					std::get< linalg::PartialPivLU<complex_mp> >(LU_).ChangePrecision(new_precision);

					current_precision_ = new_precision;
				}


				/// \brief Get the corrector's current working precision.
				unsigned precision() const
				{
					return current_precision_;
				}
				
				/**
				 \brief Change the system(number of total functions) that the predictor uses.
				 
				 \param S New system
				 
				 */
				void ChangeSystem(const System& S)
				{
					numTotalFunctions_ = static_cast<unsigned>(S.NumTotalFunctions());
					numVariables_ = static_cast<unsigned>(S.NumVariables());
					std::get< Mat<complex_dbl> >(J_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Mat<complex_mp> >(J_temp_).resize(numTotalFunctions_, numVariables_);
					std::get< Vec<complex_dbl> >(f_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_mp> >(f_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_dbl> >(step_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_mp> >(step_temp_).resize(numTotalFunctions_);
					std::get< Vec<complex_dbl> >(solve_temp_).resize(numVariables_);
					std::get< Vec<complex_mp> >(solve_temp_).resize(numVariables_);
					std::get< linalg::PartialPivLU<complex_dbl> >(LU_).ChangeSize(numVariables_);
					std::get< linalg::PartialPivLU<complex_mp> >(LU_).ChangeSize(numVariables_);
					RefreshRandomDirection();
				}

				/**
				 \brief (Re)draw the random probe direction used to estimate ||J^{-1}|| (the condition
				 number) from this thread's RNG engine.

				 The tracker calls this once at the start of a path track (see Tracker::TrackPath's
				 caller / the per-path reseed point), so the direction is held fixed for the ENTIRE
				 track of that point -- every Newton step and the endgame's sample-circle sub-tracks
				 -- and, given a per-path RNG reseed, is reproducible regardless of how paths were
				 distributed across workers.  It is NOT redrawn per Newton step or per TrackPath, which
				 would both perturb condition estimates and (in the endgame) churn the RNG mid-track.
				 */
				void RefreshRandomDirection()
				{
					std::get< Vec<complex_dbl> >(rand_temp_) = RandomOfUnits<complex_dbl>(numVariables_);
					std::get< Vec<complex_mp> >(rand_temp_) = RandomOfUnits<complex_mp>(numVariables_);
				}

				/**
				 \brief The condition-number probe direction (both precisions), so the tracker can share
				 the SAME per-path probe with the predictor -- making the predictor's and corrector's
				 ||J^{-1}|| estimates use one consistent direction.  \see RefreshRandomDirection.
				 */
				std::tuple< Vec<complex_dbl>, Vec<complex_mp> > const& ConditionProbe() const
				{
					return rand_temp_;
				}

				
				
				


				/**
				 \brief Run Newton's method, optionally with adaptive multiple precision.

				 One method replaces the former fixed/AMP/out-param overloads.  When AMP_config is null
				 it is a plain fixed-precision Newton loop (convergence only).  When AMP_config is
				 non-null it additionally fills the step metadata (norm_delta_z, norm_J, norm_J_inverse,
				 condition_number_estimate) and enforces AMP criteria B and C, returning
				 HigherPrecisionNecessary on violation.

				 \param[out] next_space The computed next space point.
				 \param[out] meta Step metadata (corrector fields; only filled in the AMP path).
				 \param S The system we are tracking on.
				 \param current_space The base point for newton correcting.
				 \param current_time The current time value.
				 \param tracking_tolerance Iterate until the step is shorter than this.
				 \param min_num_newton_iterations Take at least this many steps (>= 1).
				 \param max_num_newton_iterations The maximum number of iterations.
				 \param AMP_config Adaptive-precision settings, or nullptr for fixed precision.
				 */
				template <typename ComplexT>
				SuccessCode Correct(Vec<ComplexT> & next_space,
									   StepMetadata & meta,
									   System const& S,
									   Vec<ComplexT> const& current_space,
									   ComplexT const& current_time,
									   NumErrorT const& tracking_tolerance,
									   unsigned min_num_newton_iterations,
									   unsigned max_num_newton_iterations,
									   AdaptiveMultiplePrecisionConfig const* AMP_config = nullptr)
				{
					#ifndef BERTINI_DISABLE_ASSERTS
					assert(max_num_newton_iterations >= min_num_newton_iterations && "max number newton iterations must be at least the min.");
					#endif

					Vec<ComplexT>& step_ref = std::get< Vec<ComplexT> >(step_temp_);

					next_space = current_space;
					for (unsigned ii = 0; ii < max_num_newton_iterations; ++ii)
					{
						//Update the newton iterate by one iteration
						auto success_code = EvalIterationStep(step_ref, S, next_space, current_time);
						if(success_code != SuccessCode::Success)
							return success_code;

						next_space -= step_ref;  // step_ref = +J^{-1}f = -(Newton step); see EvalIterationStep

						// Fixed precision: cheap convergence-only loop, no norms / probe / criteria.
						if (AMP_config == nullptr)
						{
							if ( (step_ref.template lpNorm<Eigen::Infinity>() < tracking_tolerance) && (ii >= (min_num_newton_iterations-1)) )
								return SuccessCode::Success;
							continue;
						}

						// Adaptive precision: fill metadata + enforce AMP criteria B and C.
						Mat<ComplexT>& J_temp_ref = std::get< Mat<ComplexT> >(J_temp_);
						linalg::PartialPivLU<ComplexT>& LU_ref = std::get< linalg::PartialPivLU<ComplexT> >(LU_);

						meta.norm_delta_z = NumErrorT(step_ref.template lpNorm<Eigen::Infinity>());
						meta.norm_J = NumErrorT(J_temp_ref.norm());
						{
							// Reuse the FIXED probe vector generated once at setup (do NOT regenerate
							// per call): a fresh random probe direction every Newton step occasionally
							// produced an inflated ||J^{-1}|| estimate -> spurious HigherPrecisionNecessary
							// -> precision escalation/grind.  A single fixed direction also makes the
							// condition estimates comparable across steps and keeps tracking deterministic.
							Vec<ComplexT>& rand_ref = std::get< Vec<ComplexT> >(rand_temp_);
							Vec<ComplexT>& solve_ref = std::get< Vec<ComplexT> >(solve_temp_);
							LU_ref.Solve(rand_ref, solve_ref);   // reuse the factorization from EvalIterationStep
							meta.norm_J_inverse = NumErrorT(solve_ref.norm());
						}
						meta.condition_number_estimate = NumErrorT(meta.norm_J*meta.norm_J_inverse);

						if ( (meta.norm_delta_z < tracking_tolerance) && (ii >= (min_num_newton_iterations-1)) )
							return SuccessCode::Success;

						if (!amp::CriterionB<ComplexT>(meta.norm_J, meta.norm_J_inverse, max_num_newton_iterations - ii, tracking_tolerance, meta.norm_delta_z, *AMP_config))
							return SuccessCode::HigherPrecisionNecessary;

						if (!amp::CriterionC<ComplexT>(meta.norm_J_inverse, next_space, tracking_tolerance, *AMP_config))
							return SuccessCode::HigherPrecisionNecessary;
					}

					return SuccessCode::FailedToConverge;
				}


			private:

				///////////////////////////
				//
				// Private Data Methods
				//
				////////////////////
				
				
				/**
				 \brief This function computes the newton step for a system given information about the previous iteration
				 
				 \param newton_step The computed step for Newton's method
				 \param S The system used in the computations
				 \param current_space The space from the previous Newton iteration
				 \param current_time The time from the previous Newton iteration
				 
				 */
				
				template<typename ComplexT, typename Derived>
				SuccessCode EvalIterationStep(Vec<ComplexT> & newton_step,
											  const System& S,
											  const Eigen::MatrixBase<Derived>& current_space, const ComplexT& current_time)
				{
					Vec<ComplexT>& f_temp_ref = std::get< Vec<ComplexT> >(f_temp_);
					Mat<ComplexT>& J_temp_ref = std::get< Mat<ComplexT> >(J_temp_);
					
					linalg::PartialPivLU<ComplexT>& LU_ref = std::get< linalg::PartialPivLU<ComplexT> >(LU_);

					S.SetAndReset<ComplexT>(current_space, current_time);
					S.EvalInPlace(f_temp_ref);
					S.JacobianInPlace(J_temp_ref);
					// Factor a copy of J (J_temp_ref is read again afterward for meta.norm_J).  The health
					// check is folded into Factor().
					if (LU_ref.Factor(J_temp_ref)!=MatrixSuccessCode::Success)
						return SuccessCode::MatrixSolveFailure;

					// Solve J*newton_step = f (NOT -f): this lets us skip materializing the negated RHS
					// temporary.  newton_step therefore holds +J^{-1}f = -(true Newton step), so callers
					// SUBTRACT it (next_space -= newton_step).  The convergence test uses lpNorm, which is
					// sign-insensitive, so it is unaffected.
					LU_ref.Solve(f_temp_ref, newton_step);

					return SuccessCode::Success;
					
				}
				

				
				///////////////////////////
				//
				// Private Data Members
				//
				////////////////////
				
				unsigned numTotalFunctions_; // Number of total functions for the current system
				unsigned numVariables_;  // Number of variables for the current system
				
				std::tuple< Vec<complex_dbl>, Vec<complex_mp> > f_temp_;
				std::tuple< Vec<complex_dbl>, Vec<complex_mp> > step_temp_;
				std::tuple< Mat<complex_dbl>, Mat<complex_mp> > J_temp_;
				std::tuple< Vec<complex_dbl>, Vec<complex_mp> > rand_temp_;  // reused scratch: random RHS for norm_J_inverse
				std::tuple< Vec<complex_dbl>, Vec<complex_mp> > solve_temp_; // reused scratch: LU solve result

				std::tuple< linalg::PartialPivLU<complex_dbl>, linalg::PartialPivLU<complex_mp> > LU_;
				
				unsigned current_precision_;

				NewtonConfig newton_config_; // Hold the settings of the Newton iteration

				
			}; //re: class NewtonCorrector
			
		} //re: namespace correct
	}// re: namespace tracking
}// re: namespace bertini

#endif
