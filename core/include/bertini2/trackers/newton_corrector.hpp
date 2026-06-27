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
				 /brief Change the precision of the predictor variables and reassign the Butcher table variables.
				 
				 \param new_precision The new precision.
				 
				 */
				void ChangePrecision(unsigned new_precision)
				{
					Precision(std::get< Vec<complex_mp> >(f_temp_), new_precision);
					Precision(std::get< Vec<complex_mp> >(step_temp_), new_precision);
					Precision(std::get< Mat<complex_mp> >(J_temp_), new_precision);
					Precision(std::get< Vec<complex_mp> >(rand_temp_), new_precision);
					Precision(std::get< Vec<complex_mp> >(solve_temp_), new_precision);

					current_precision_ = new_precision;
				}


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
				 \brief Run Newton's method in fixed precision.
				 
				 Run Newton's method until it converges (\f$\Delta z\f$ < tol), or the next point's norm exceeds the path truncation threshold.
				 
				 \return The SuccessCode indicating what happened.
				 
				 \tparam ComplexT The complex type for arithmetic
				 
				 \param[out] next_space The computed next space point.
				 \param S The system we are tracking on.
				 \param current_space The base point for newton correcting.
				 \param current_time The current time value.  Note it is complex.
				 \param tracking_tolerance The upper threshold for step size.  Must iterate correcting until the corrector step is less than this threshold in length.
				 \param path_truncation_threshold Correcting stops the the norm of the current solution exceeds this number.
				 \param min_num_newton_iterations The corrector must take at least this many steps.  This should be at least 1.
				 \param max_num_newton_iterations The maximum number of iterations to run Newton's method for.
				 
				 */
				
				template <typename ComplexT>
				SuccessCode Correct(Vec<ComplexT> & next_space,
									   System const& S,
									   Vec<ComplexT> const& current_space, // pass by value to get a copy of it
									   ComplexT const& current_time,
									   NumErrorT const& tracking_tolerance,
									   unsigned min_num_newton_iterations,
									   unsigned max_num_newton_iterations)
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
						
						next_space += step_ref;
						
						if ( (step_ref.template lpNorm<Eigen::Infinity>() < tracking_tolerance) && (ii >= (min_num_newton_iterations-1)) )
							return SuccessCode::Success;
					}
					
					return SuccessCode::FailedToConverge;

					

				}

				
				
				
				
				
				/**
				 \brief Run Newton's method in multiple precision.
				 
				 Run Newton's method until it converges (\f$\Delta z\f$ < tol), an AMP criterion (B or C) is violated, or the next point's norm exceeds the path truncation threshold.
				 
				 \tparam ComplexT The complex type for arithmetic
				 
				 \param[out] next_space The computed next space point.
				 \param S The system we are tracking on.
				 \param current_space The base point for newton correcting.
				 \param current_time The current time value.  Note it is complex.
				 \param tracking_tolerance The upper threshold for step size.  Must iterate correcting until the corrector step is less than this threshold in length.
				 \param path_truncation_threshold Correcting stops the the norm of the current solution exceeds this number.
				 \param min_num_newton_iterations The corrector must take at least this many steps.  This should be at least 1.
				 \param max_num_newton_iterations The maximum number of iterations to run Newton's method for.
				 \param AMP_config Adaptive multiple precision settings.  Using this argument is how Bertini2 knows you want to use adaptive precision.
				 */
				template <typename ComplexT>
				SuccessCode Correct(Vec<ComplexT> & next_space,
									   System const& S,
									   Vec<ComplexT> const& current_space, // pass by value to get a copy of it
									   ComplexT const& current_time,
									   NumErrorT const& tracking_tolerance,
									   unsigned min_num_newton_iterations,
									   unsigned max_num_newton_iterations,
									   AdaptiveMultiplePrecisionConfig const& AMP_config)
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
						
						next_space += step_ref;
						
						Mat<ComplexT>& J_temp_ref = std::get< Mat<ComplexT> >(J_temp_);
						Eigen::PartialPivLU< Mat<ComplexT> >& LU_ref = std::get< Eigen::PartialPivLU< Mat<ComplexT> > >(LU_);
						
						if ( (step_ref.template lpNorm<Eigen::Infinity>() < tracking_tolerance) && (ii >= (min_num_newton_iterations-1)) )
							return SuccessCode::Success;
						
						Vec<ComplexT> const& rand_ref = std::get< Vec<ComplexT> >(rand_temp_);
						Vec<ComplexT>& solve_ref = std::get< Vec<ComplexT> >(solve_temp_);
						solve_ref = LU_ref.solve(rand_ref);
						NumErrorT norm_J_inverse(solve_ref.norm());

						if (!amp::CriterionB<ComplexT>(NumErrorT(J_temp_ref.norm()), norm_J_inverse, max_num_newton_iterations - ii, tracking_tolerance, NumErrorT(step_ref.template lpNorm<Eigen::Infinity>()), AMP_config))
							return SuccessCode::HigherPrecisionNecessary;

						if (!amp::CriterionC<ComplexT>(norm_J_inverse, next_space, tracking_tolerance, AMP_config))
							return SuccessCode::HigherPrecisionNecessary;
					}

					return SuccessCode::FailedToConverge;
				}






				/**
				 \brief Run Newton's method in multiple precision.
				 
				 Run Newton's method until it converges (\f$\Delta z\f$ < tol), an AMP criterion (B or C) is violated, or the next point's norm exceeds the path truncation threshold.
				 
				 \tparam ComplexT The complex type for arithmetic
				 \tparam RealT The underlying real number type, used for comparitors.
				 
				 \param[out] next_space The computed next space point.
				 \param[out] norm_delta_z The norm of the last step size.
				 \param[out] norm_J The matrix norm of the Jacobian matrix.
				 \param[out] norm_J_inverse The matrix norm of the inverse of the Jacobian matrix.
				 \param[out] condition_number_estimate An estimate on the condition number.
				 \param S The system we are tracking on.
				 \param current_space The base point for newton correcting.
				 \param current_time The current time value.  Note it is complex.
				 \param tracking_tolerance The upper threshold for step size.  Must iterate correcting until the corrector step is less than this threshold in length.
				 \param path_truncation_threshold Correcting stops the the norm of the current solution exceeds this number.
				 \param min_num_newton_iterations The corrector must take at least this many steps.  This should be at least 1.
				 \param max_num_newton_iterations The maximum number of iterations to run Newton's method for.
				 \param AMP_config Adaptive multiple precision settings.  Using this argument is how Bertini2 knows you want to use adaptive precision.
				 */
				template <typename ComplexT>
				SuccessCode Correct(Vec<ComplexT> & next_space,
									   NumErrorT & norm_delta_z,
									   NumErrorT & norm_J,
									   NumErrorT & norm_J_inverse,
									   NumErrorT & condition_number_estimate,
									   System const& S,
									   Vec<ComplexT> const& current_space, // pass by value to get a copy of it
									   ComplexT const& current_time,
									   NumErrorT const& tracking_tolerance,
									   unsigned min_num_newton_iterations,
									   unsigned max_num_newton_iterations,
									   AdaptiveMultiplePrecisionConfig const& AMP_config)
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
						
						next_space += step_ref;
						
						Mat<ComplexT>& J_temp_ref = std::get< Mat<ComplexT> >(J_temp_);
						Eigen::PartialPivLU< Mat<ComplexT> >& LU_ref = std::get< Eigen::PartialPivLU< Mat<ComplexT> > >(LU_);
						
						
						norm_delta_z = NumErrorT(step_ref.template lpNorm<Eigen::Infinity>());
						norm_J = NumErrorT(J_temp_ref.norm());
						{
							// Reuse the FIXED probe vector generated once at setup (do NOT regenerate
							// per call): a fresh random probe direction every Newton step occasionally
							// produced an inflated ||J^{-1}|| estimate -> spurious HigherPrecisionNecessary
							// -> precision escalation/grind (confirmed by A/B: regenerating reintroduces
							// the spikes, fixed does not).  A single fixed direction also makes the
							// condition estimates comparable across steps, is cheaper, and keeps tracking
							// deterministic / parallel-bit-identical.
							Vec<ComplexT>& rand_ref = std::get< Vec<ComplexT> >(rand_temp_);
							Vec<ComplexT>& solve_ref = std::get< Vec<ComplexT> >(solve_temp_);
							solve_ref = LU_ref.solve(rand_ref);
							norm_J_inverse = NumErrorT(solve_ref.norm());
						}
						condition_number_estimate = NumErrorT(norm_J*norm_J_inverse);
												
						if ( (norm_delta_z < tracking_tolerance) && (ii >= (min_num_newton_iterations-1)) )
							return SuccessCode::Success;
						
						if (!amp::CriterionB<ComplexT>(norm_J, norm_J_inverse, max_num_newton_iterations - ii, tracking_tolerance, norm_delta_z, AMP_config))
							return SuccessCode::HigherPrecisionNecessary;

						if (!amp::CriterionC<ComplexT>(norm_J_inverse, next_space, tracking_tolerance, AMP_config))
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
					
					Eigen::PartialPivLU< Mat<ComplexT> >& LU_ref = std::get< Eigen::PartialPivLU< Mat<ComplexT> > >(LU_);

					S.SetAndReset<ComplexT>(current_space, current_time);
					S.EvalInPlace(f_temp_ref);
					S.JacobianInPlace(J_temp_ref);
					LU_ref.compute(J_temp_ref);
					
					if (LUPartialPivotDecompositionSuccessful(LU_ref.matrixLU())!=MatrixSuccessCode::Success)
						return SuccessCode::MatrixSolveFailure;
					
					newton_step = LU_ref.solve(-f_temp_ref);
					
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

				std::tuple< Eigen::PartialPivLU<Mat<complex_dbl>>, Eigen::PartialPivLU<Mat<complex_mp>> > LU_;
				
				unsigned current_precision_;

				NewtonConfig newton_config_; // Hold the settings of the Newton iteration

				
			}; //re: class NewtonCorrector
			
		} //re: namespace correct
	}// re: namespace tracking
}// re: namespace bertini

#endif
