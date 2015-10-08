//This file is part of Bertini 2.0.
//
//tracker_loop.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracker_loop.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracker_loop.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  tracker_loop.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015


#ifndef BERTINI_TRACKER_LOOP_HPP
#define BERTINI_TRACKER_LOOP_HPP

#include "tracking/step.hpp"

namespace bertini{

	namespace tracking{


		


		enum class Branch
		{
			A, B
		};

		
		

		template <typename RealType>
		RealType epsilon(unsigned precision)
		{
			return pow( RealType(10), -precision+3);
		}







		template <typename ComplexType>
		void ConvergenceError(unsigned & new_precision, ComplexType & new_stepsize, 
							  unsigned old_precision, ComplexType const& old_stepsize,
							  RealType const& step_adjustment_factor)
		{
			
			new_precision = old_precision;
			new_stepsize = step_adjustment_factor*old_stepsize;

			while (norm(new_stepsize) < epsilon<ComplexType>(new_precision))
				if new_precision==16
					new_precision=30;
				else
					new_precision+=PrecisionIncrement;
		}




		/**
		 \brief Compute the cost function from \cite{amp2}, \f$C(P)\f$.  As currently implemented, this is 
		 \f$ 10.35 + 0.13 P \f$ where P is the precision.

		 \todo Recompute this cost function for boost::multiprecision::mpfr_float
		*/
		double Cost(unsigned precision)
		{
			if (precision==16)
				return 1;
			else
				return 10.35 + 0.13 * precision;
		}






		template <typename RealType>
		struct PrecStep
		{
			RealType stepsize;
			unsigned precision;

			PrecStep(RealType const& ss, unsigned prec) : stepsize(ss), precision(prec)
			{}

			PrecStep(unsigned prec, RealType const& ss) : stepsize(ss), precision(prec)
			{}

		};


		template <RealType> 
		RealType StepsizeSatisfyingNine(unsigned current_precision,
										RealType const& tau, // -log10 of track tolerance
		                                RealType const& D, // from condition B, C
		                                RealType const& a, // size_proportion
		                                unsigned predictor_order = 1,
		                                unsigned num_newton_iterations,
		                                RealType const& d, // norm(delta_z, delta_t)
		                                AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			RealType stepsize;

			auto rhs = num_newton_iterations/(predictor_order+1)*( AMP_config.safety_digits_1 + D + (tau + log10(a))/num_newton_iterations) - current_precision; 
			stepsize = pow(10,rhs);

			return stepsize;
		}








 
		void MinimizeCost(unsigned & new_precision, mpfr_float & new_stepsize, 
						  unsigned min_precision, mpfr_float const& old_stepsize,
						  unsigned max_precision,
						  AdaptiveMultiplePrecisionConfig const& AMP_config_)
		{
			std::set<PrecStep<mpfr_float> > computed_candidates;

			unsigned lowest_mp_precision_to_test = min_precision;

			if (min_precision==16)
			{
				unsigned test_precision = 16;
				
				auto test_stepsize = StepsizeSatisfyingNine(test_precision,
				                                            unsigned current_precision,
															tau, // -log10 of track tolerance
							                                D, // from condition B
							                                a, // size_proportion
							                                predictor_order,
							                                num_newton_iterations,
							                                d, // norm(delta_z, delta_t)
							                                AMP_config);
				
				if (min(abs(test_stepsize),abs(old_stepsize)) > epsilon<RealType>(test_precision))
					computed_candidates.insert(PrecStep<RealType>(test_precision, test_stepsize));

				lowest_mp_precision_to_test = LowestMultiplePrecision;
			}

			
			for (unsigned test_precision = lowest_mp_precision_to_test; test_precision <= max_precision; test_precision+=PrecisionIncrement)
			{

				auto test_stepsize = StepsizeSatisfyingNine(test_precision,
				                                            unsigned current_precision,
															tau, // -log10 of track tolerance
							                                D, // from condition B
							                                a, // size_proportion
							                                predictor_order,
							                                num_newton_iterations,
							                                d, // norm(delta_z, delta_t)
							                                AMP_config);

				if (min(abs(test_stepsize),abs(old_stepsize)) > epsilon<RealType>(test_precision))
					computed_candidates.insert(PrecStep<mpfr_float>(test_precision, test_stepsize));
			}

			if (computed_candidates.empty())
			{
				new_precision = max_precision + 1;
				new_stepsize = old_stepsize;
			}
			else
			{
				// find the minimizer of the Cost function
				double min_cost(1e300);
				for (auto candidate : computed_candidates)
				{	
					double curr_cost = Cost(candidate.precision) / abs(candidate.stepsize)
					if (curr_cost < min_cost)
					{
						min_cost = curr_cost;
						new_precision = candidate.precision;
						new_stepsize = candidate.stepsize;
					}
				}
				
			}
		}






		template <typename RealType, typename ComplexType>
		Branch SafetyError(unsigned new_precision, mpfr_float & new_stepsize, 
		                   Vec<ComplexType> & new_delta_z, ComplexType & new_delta_t,

						 unsigned current_precision, mpfr_float const& current_stepsize, 
						 Vec<ComplexType> const& delta_z, ComplexType const& delta_t,

						 RealType const& norm_J, RealType const& norm_J_inverse,
						 RealType const& size_proportion,
						 unsigned predictor_order
						 )
		{
			new_precision = current_precision;
			new_stepsize = current_stepsize;
			new_delta_t = delta_t;
			new_delta_z = delta_z;

			while(!amp::CriterionC())
			{

			}

			if ()
				return Branch::A;
			else
				return Branch::B;
		}


		/**
		Computes the minimum and maximum precisions necessary for computation to continue.
		*/
		void AMP3_update(unsigned & new_precision, mpfr_float & new_stepsize,
		                 unsigned old_precision, mpfr_float const& old_stepsize,
		                 unsigned digits_tau,
		                 unsigned digits_final = 0,
		                 unsigned P0,
		                 unsigned digits_C,
		                 mpfr const& current_time,
		                 mpfr_float const& eta_min,
		                 mpfr_float const& eta_max,
		                 unsigned max_num_newton_iterations,
		                 AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			// compute the number of digits necessary to trigger a lowering.  
			unsigned precision_decrease_threshold;
			if (old_precision==DoublePrecision)
				precision_decrease_threshold = 0; // if in double, cannot lower precision any further
			else if (old_precision==LowestMultiplePrecision)
				precision_decrease_threshold = LowestMultiplePrecision - DoublePrecision; // if in lowest multiple, have to clear the difference down to double
			else 
				precision_decrease_threshold = PrecisionIncrement; // the threshold is the step between precisions in multiple precision


			unsigned digits_B = ceil(P0 - (predictor_order_+1)* eta_min/max_num_newton_iterations);

			if (digits_B < old_precision)
				digits_B = min(digits_B+precision_decrease_threshold, old_precision);

			if (digits_C < old_precision)
				digits_C = min(digits_C+precision_decrease_threshold, old_precision);

			unsigned digits_stepsize = max(MinDigitsFromEta(eta_min,current_time), MinDigitsFromEta(eta_max,current_time));

			if (digits_stepsize < old_precision)
				digits_stepsize = min(digits_stepsize+precision_decrease_threshold, old_precision);

			min_digits_needed = max(digits_tau, digits_final, digits_B, digits_C, digits_stepsize);
			max_digits_needed = max(min_digits_needed, ceil(P0 - (predictor_order_+1)* eta_max/max_num_newton_iterations));
			

			MinimizeCost(
			             new_precision, new_stepsize, // the computed values
					  	old_precision, old_stepsize,
					  min_digits_needed, max_digits_needed,
					  AMP_config);
		}
		

		unsigned MinDigitsFromEta(mpfr_float eta, mpfr current_time)
		{
			unsigned min_digits = ceil(eta) + ceil(log10(current_time)) + 2;
		}



		/** 
		\brief Functor-like class for tracking paths on a system

		create an instance of this class, feeding it the system to be tracked on, and some configuration.  Then, this this tracker to track paths of the system.
		*/
		class Tracker
		{
		public:



			Tracker(System sys, 
			        config::Predictor new_predictor_choice = config::Predictor::Euler,
			        ) : tracked_system(sys)
			{	
				Predictor(config::Predictor new_predictor_choice);
				
				// initialize to the frequency so guaranteed to compute it the first try 
				num_steps_since_last_condition_number_computation = frequency_of_CN_estimation;
			}


			/**
			\brief Get the tracker set up for tracking.

			Pass the tracker the configuration for tracking, to get it set up.
			*/
			void Setup(System const& sys,
						config::Predictor predictor_choice,
						config::PrecisionType prec_type = config::PrecisionType::Adaptive,
						mpfr_float const& tracking_tolerance,
						mpfr_float const& path_truncation_threshold,
						config::Stepping const& stepping,
						config::Newton const& newton,
						config::Security const& security,
						config::AdaptiveMultiplePrecisionConfig const& AMP 
						)
			{
				tracked_system_ = sys;
				Predictor(predictor_choice);

				tracking_tolerance_ = tracking_tolerance;
				path_truncation_threshold_ = path_truncation_threshold;

				stepping_config_ = stepping;
				newton_config_ = newton;
				security_config_ = security;
				AMP_config_ = AMP_config;
			}







			SuccessCode TrackPathAMP(Vec<mpfr> & solution_at_endtime,
									mpfr const& start_time, mpfr const& endtime,
									Vec<mpfr> const& start_point
									)
			{
				// set up the master current time and the current step size
				current_time_.precision(start_time.precision())
				current_time_ = start_time;

				current_stepsize_.precision(stepping_.initial_step_size.precision())
				current_stepsize_ = stepping_.initial_step_size;

				current_precision_ = start_point.precision(); // get the current precision.
				num_consecutive_successful_steps = 0;
				num_successful_steps_taken = 0;
				num_failed_steps = 0;


				// populate the current space value with the start point, in appropriate precision
				if (current_precision_==16)
					MultipleToDouble( start_point);
				else
					MultipleToMultiple( start_point);


				SuccessCode initial_refinement = Refine();
				if (initial_refinement!=SuccessCode::Success)
				{
					do {
						if (current_precision_ > AMP_config_.maximum_precision)
							return SuccessCode::MaxPrecisionReached;

						if (current_precision_==16)
							initial_refinement = ChangePrecision(LowestMultiplePrecision);
						else
							initial_refinement = ChangePrecision(current_precision_+PrecisionIncrement);
					}
					while (initial_refinement!=SuccessCode::Success);
				}



				// as precondition to this while loop, the correct container, either dbl or mpfr, must have the correct data.
				while (current_time_)!=endtime)
				{
					if (num_successful_steps_taken >= stepping_.max_num_steps)
						return SuccessCode::MaxNumStepsTaken;
					if (current_stepsize_ < stepping_.min_step_size)
						return SuccessCode::MinStepSizeReached;
					



					// compute the next delta_t
					if (abs(endtime-current_time_) < abs(current_stepsize_))
						delta_t_ = endtime-current_time_;
					else
						delta_t_ = current_stepsize_ * (endtime - current_time_)/abs(endtime - current_time_);



					if (current_precision_==16)
						step_success_code = TrackerIteration<dbl>();
					else
						step_success_code = TrackerIteration<mpfr>();



					if (step_success_code==SuccessCode::Success)
					{
						num_steps_taken_++; 
						num_consecutive_successful_steps_++;
						current_time_ += delta_t_;
					}
					else
					{
						num_consecutive_successful_steps=0;
						num_failed_steps++;
					}


					if (num_consecutive_successful_steps_ >= stepping_.consecutive_successful_steps_before_precision_decrease)
					{
						if (current_precision_==LowestMultiplePrecision)
							MultipleToDouble();
						else if (current_precision_ > LowestMultiplePrecision)
							MultipleToMultiple(current_precision_ - PrecisionIncrement);
					}


					if (num_consecutive_successful_steps_ >= stepping_.consecutive_successful_steps_before_stepsize_increase)
					{
						if (current_stepsize_ < stepping_.max_step_size)
							current_stepsize_ = min(stepping_.max_step_size, current_stepsize_*stepping_.step_size_success_factor);
					}


					if (next_precision!=current_precision_)
					{
						SuccessCode precision_change_code = ChangePrecision(next_precision);
						if (precision_change_code!=SuccessCode::Success)
							return precision_change_code;
					}

				}// re: while


				return SuccessCode::Success;
			}








			/**
			Change tracker to use a predictor

			\param new_predictor_choice The new predictor to be used.

			\see config::Predictor
			*/
			void Predictor(config::Predictor new_predictor_choice)
			{
				predictor_choice_ = new_predictor_choice;
				predictor_order_ = PredictorOrder(predictor_choice);
			}


			/**
			\brief Query the currently set predictor
			*/
			config::Predictor Predictor()
			{
				return predictor_choice_;
			}


		private:


			// Vec<ComplexType> & next_space, ComplexType & next_time, 
			// 							ComplexType & next_delta_t, unsigned & next_precision,
			// 							RealType & condition_number_estimate,
			// 							// end the output variables
			// 							Vec<ComplexType> const& current_space, ComplexType const& current_time, 
			// 							ComplexType const& delta_t
										

			template <typename ComplexType>
			SuccessCode TrackerIteration()
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
				#endif

				typedef NumTraits<ComplexType>::Real RealType;

				Vec<ComplexType>& predicted_space = std::get<Vec<ComplexType> >(temporary_space_); // this will be populated in the Predict step
				Vec<ComplexType>& current_space = std::get<Vec<ComplexType> >(current_space_);


				ComplexType current_time = ComplexType(current_time_);
				ComplexType delta_t = ComplexType(delta_t_);

				SuccessCode predictor_code = Predict(predicted_space, current_space, current_time, delta_t);

				if (predictor_code==SuccessCode::MatrixSolveFailure)
				{
					ConvergenceError();
					return predictor_code;
				}	

				if (predictor_code==SuccessCode::HigherPrecisionNecessary)
				{	
					SafetyError();
					return predictor_code;
				}

				Vec<ComplexType>& tentative_next_space = std::get<Vec<ComplexType> >(tentative_space_); // this will be populated in the Correct step

				ComplexType tentative_next_time = current_time + delta_t;

				SuccessCode corrector_code = Correct(tentative_next_space,
				                                     predicted_space
				                                     tentative_next_time);

				if (corrector_code==SuccessCode::MatrixSolveFailure || corrector_code==SuccessCode::FailedToConverge)
				{
						ConvergenceError();
						return corrector_code;
				}
				else if (corrector_code == SuccessCode::HigherPrecisionNecessary)
				{
						SafetyError();
						return corrector_code;
				}

				return SuccessCode::Success;
			}



			/**
			\brief Wrapper function for calling the correct predictor.
			*/
			template<typename ComplexType>
			SuccessCode Predict(Vec<ComplexType> & predicted_space, 
			                    Vec<ComplexType> const& current_space, 
			                    ComplexType const& current_time, ComplexType const& delta_t)
			{
				typedef NumTraits<ComplexType>::Real RealType;

				if (tracking::HasErrorEstimate(predictor_choice_))
					return = Predict(predictor_choice_,
									predicted_space,
									std::get<RealType>(error_estimate),
									std::get<RealType>(size_proportion),
									std::get<RealType>(norm_J),
									std::get<RealType>(norm_J_inverse_),
									tracked_system_,
									current_space, current_time, 
									delta_t,
									std::get<RealType>(condition_number_estimate_),
									num_steps_since_last_condition_number_computation, 
									frequency_of_CN_estimation, 
									RealType(tracking_tolerance_),
									AMP_config);
				else
					return = Predict(predictor_choice_,
									predicted_space,
									std::get<RealType>(size_proportion),
									std::get<RealType>(norm_J),
									std::get<RealType>(norm_J_inverse_),
									tracked_system_,
									current_space, current_time, 
									delta_t,
									std::get<RealType>(condition_number_estimate_),
									num_steps_since_last_condition_number_computation, 
									frequency_of_CN_estimation, PrecisionType::Adaptive, 
									RealType(tracking_tolerance_),
									AMP_config);
			}




			/**
			Wrapper function for calling Correct and getting the error estimates etc directly into the tracker object.
			*/
			template<typename ComplexType>
			SuccessCode Correct(Vec<ComplexType> & corrected_space, 
			                    Vec<ComplexType> const& current_space, 
			                    ComplexType const& current_time)
			{
				return Correct(corrected_space,
						   tracked_system_,
						   current_space,
						   current_time, 
						   PrecisionType::Adaptive, 
						   RealType(tracking_tolerance_),
						   RealType(path_truncation_threshold_),
						   newton_.min_num_newton_iterations_,
						   newton_.max_num_newton_iterations_,
						   AMP_config_);
			}




			/**
			\brief Increase stepsize or decrease precision, because of consecutive successful steps.
			*/
			void StepSuccess()
			{
				if (HasErrorEstimate(predictor_))
					StepSuccessWithErrorEstimate();
				else
					StepSuccessNoErrorEstimate();

			}

			/**
			Computes new stepsize and precision
			*/
			void StepSuccessWithErrorEstimate(double)
			{
				mpfr_float minimum_stepsize = current_stepsize_ * stepping_.step_size_fail_factor;


				unsigned precision_from_criterion_B = amp::CriterionBRHS(norm_J, norm_J_inverse, num_newton_iterations_remaining, tracking_tolerance,  norm_of_latest_newton_residual, AMP_config_);

				unsigned precision_from_criterion_C = round( amp::CriterionCRHS(norm_J_inverse, z, tracking_tolerance, AMP_config));
				precision_from_criterion_C = max(precision_from_criterion_C, current_precision_);


				mpfr_float maximum_stepsize = current_stepsize_;

				if (num_consecutive_successful_steps_ < stepping_.consecutive_successful_steps_before_stepsize_increase)
					// this space left intentionally blank
				else if (num_consecutive_successful_steps_ == stepping_.consecutive_successful_steps_before_stepsize_increase)
					maximum_stepsize = min(current_stepsize_ * stepping_.step_size_success_factor, stepping_.max_step_size);
				else
					num_consecutive_successful_steps_=0;

				mpfr_float eta_min = -log10(minimum_stepsize);
				mpfr_float eta_max = -log10(maximum_stepsize);


				// now call AMP3_update()
			}



			/**
			Computes new stepsize and precision
			*/
			void StepSuccessWithErrorEstimate(mpfr_float)
			{
				mpfr_float minimum_stepsize = current_stepsize_ * stepping_.step_size_fail_factor;


				unsigned precision_from_criterion_B = amp::CriterionBRHS(norm_J, norm_J_inverse, num_newton_iterations_remaining, tracking_tolerance,  norm_of_latest_newton_residual, AMP_config_);

				unsigned precision_from_criterion_C = round( amp::CriterionCRHS(norm_J_inverse, z, tracking_tolerance, AMP_config)); // hopefully this number is smaller than the current precision, allowing us to reduce precision


				mpfr_float maximum_stepsize = current_stepsize_;

				if (num_consecutive_successful_steps_ < stepping_.consecutive_successful_steps_before_stepsize_increase)
				{	
					precision_from_criterion_C = max(precision_from_criterion_C, current_precision_);
				}
				else if (num_consecutive_successful_steps_ == stepping_.consecutive_successful_steps_before_stepsize_increase)
				{
					// allow stepsize increase, up to maximum step size
					maximum_stepsize = min(current_stepsize_ * stepping_.step_size_success_factor, stepping_.max_step_size);
				}	
				else if (num_consecutive_successful_steps_ < stepping_.consecutive_successful_steps_before_precision_decrease)
				{
					precision_from_criterion_C = max(precision_from_criterion_C, current_precision_);
				}
				else if (num_consecutive_successful_steps_ == stepping_.consecutive_successful_steps_before_precision_decrease)
				{
					// allow stepsize increase, up to maximum step size
					maximum_stepsize = min(current_stepsize_ * stepping_.step_size_success_factor, stepping_.max_step_size);
					if (num_precision_decreases >= AMP_config_.max_num_precision_decreases)
						precision_from_criterion_C = max(precision_from_criterion_C, current_precision_);
				}
				else
				{
					// reset number of consecutive steps to 0.
					num_consecutive_successful_steps_=0;
					precision_from_criterion_C = max(precision_from_criterion_C, current_precision_);
				}


				max_num_precision_decreases

				mpfr_float eta_min = -log10(minimum_stepsize);
				mpfr_float eta_max = -log10(maximum_stepsize);


				AMP3_update()

				// now call AMP3_update()
			}


			


			/////////////////
			//
			//  Functions for refining a current space value
			//
			///////////////////////


			/**
			\brief Run Newton's method from the currently stored time and space value, in current precision.

			This overwrites the value of the current space, if successful, with the refined value.
			*/
			SuccessCode Refine()
			{
				SuccessCode code;
				if (current_precision==16)
				{
					code = Refine(std::get<Vec<dbl> >(temporary_space_),std::get<Vec<dbl> >(current_space_), std::get<dbl>(current_time_),double(tracking_tolerance_));
					if (code == SuccessCode::Success)
						std::get<Vec<dbl> >(current_space_) = std::get<Vec<dbl> >(temporary_space_);
				}
				else
				{
					code = Refine(std::get<Vec<mpfr> >(temporary_space_),std::get<Vec<mpfr> >(current_space_), std::get<mpfr>(current_time_),tracking_tolerance_);
					if (code == SuccessCode::Success)
						std::get<Vec<mpfr> >(current_space_) = std::get<Vec<mpfr> >(temporary_space_);
				}
				return code;
			}







			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode Refine(Vec<ComplexType> & new_space
			                    Vec<ComplexType> const& start_point, ComplexType const& current_time,
			                    RealType const& tolerance)
			{
				return Correct(next_space,
							   tracked_system_,
							   start_point,
							   current_time, 
							   PrecisionType::Adaptive, 
							   tolerance,
							   path_truncation_threshold_,
							   newton_.min_num_newton_iterations,
							   newton_.max_num_newton_iterations,
							   AMP_config_);
			}















			/////////////////
			//
			//  Functions for converting between precision types
			//
			///////////////////////


			/**
			Change precision of tracker to next_precision.  Converts the internal temporaries, and adjusts precision of system. Then refines if necessary.

			If the new precision is higher than current precision, a refine step will be called, which runs Newton's method.  This may fail, leaving the tracker in a state with higher precision internals, but garbage digits after the previously known digits.
			*/
			SuccessCode ChangePrecision(unsigned new_precision)
			{
				if (new_precision==current_precision_) // no op
					return;

				num_steps_since_last_condition_number_computation = frequency_of_CN_estimation;

				bool upsampling_needed = new_precision > current_precision_;

				if (new_precision==16 && current_precision_>16)
				{
					// convert from multiple precision to double precision
					MultipleToDouble();
				}
				else if(new_precision > 16 && current_precision_ == 16)
				{
					// convert from double to multiple precision
					DoubleToMultiple();
				}
				else if(new_precision < current_precision_)
				{
					MultipleToMultiple();
				}
				else if (new_precision > current_precision_)
				{
					MultipleToMultiple();
				}


				#ifndef BERTINI_DISABLE_ASSERTS
				assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
				#endif


				if (upsampling_needed)
					return Refine();
				else
					return SuccessCode::Success;
			}





			


			/**
			\brief Converts from multiple to double

			Copies a multiple-precision into the double storage vector, and changes precision of the time and delta_t.
			*/
			void MultipleToDouble(Vec<mpfr> const& source_point, mpfr const& source_time, mpfr const& source_delta_t)
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == system.NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				assert(source_point.size() == std::get<Vec<dbl> >(current_space_).size() && "source point and current space point must be same size when changing precision.  Precondition not met.");
				#endif

				current_precision_ = 16;
				system.precision(16);

				if (std::get<Vec<dbl> >(current_space_).size()!=source_point.size())
					std::get<Vec<dbl> >(current_space_).size(source_point.size());

				for (unsigned ii=0; ii<source_point.size())
					std::get<Vec<dbl> >(current_space_)(ii) = dbl(source_point(ii));
			}

			/**
			\brief Converts from multiple to double

			Changes the precision of the internal temporaries to double precision
			*/
			void MultipleToDouble()
			{
				MultipleToDouble(std::get<Vec<mpfr> >(current_space_), std::get<mpfr>(current_time_), std::get<mpfr>(delta_t_));
			}





			/**
			\brief Converts from double to multiple

			Copies a double-precision into the multiple-precision storage vector.  

			You should call Newton after this to populate the new digits with non-garbage data.
			*/
			void DoubleToMultiple(unsigned new_precision, Vec<dbl> const& source_point, dbl const& source_time, dbl const& source_delta_t)
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == system.NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				assert(new_precision > 16 && "must convert to precision higher than 16 when converting to multiple precision");
				#endif

				current_precision_ = new_precision;
				mpfr_float.precision(new_precision);
				system.precision(new_precision);

				if (std::get<Vec<mpfr> >(current_space_).size()!=source_point.size())
					std::get<Vec<mpfr> >(current_space_).size(source_point.size());

				for (unsigned ii=0; ii<source_point.size())
					std::get<Vec<mpfr> >(current_space_)(ii) = mpfr(source_point(ii));

				AdjustTemporariesPrecision(new_precision);

				#ifndef BERTINI_DISABLE_ASSERTS
				assert(std::get<mpfr>(current_space_)(0).precision() == current_precision_ && "precision of time in mpfr doesn't match tracker");
				#endif
			}



			/**
			\brief Converts from double to multiple

			Copies the double-precision temporaries into the multiple-precision temporaries.  You should call Newton after this to populate the new digits with non-garbage data.
			*/
			void DoubleToMultiple(unsigned new_precision)
			{
				DoubleToMultiple( new_precision, std::get<Vec<dbl> >(current_space_), std::get<dbl>(current_time_), std::get<dbl>(delta_t_)
			}





			/**
			\brief Converts from multiple to different precision multiple precision

			Copies a multiple-precision into the multiple-precision storage vector, and changes precision of the time and delta_t.
			Also resets counter so have to re-compute the condition number on next step attempt.
			*/
			void MultipleToMultiple(unsigned new_precision, Vec<mpfr> const& source_point, mpfr const& source_time, mpfr const& source_delta_t)
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == system.NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				assert(new_precision > 16 && "must convert to precision higher than 16 when converting to multiple precision");
				#endif

				current_precision_ = new_precision;
				mpfr_float.precision(new_precision);
				system.precision(new_precision);

				if (std::get<Vec<mpfr> >(current_space_).size()!=source_point.size())
					std::get<Vec<mpfr> >(current_space_).size(source_point.size());

				for (unsigned ii=0; ii<source_point.size())
					std::get<Vec<mpfr> >(current_space_)(ii) = mpfr(source_point(ii));

				AdjustTemporariesPrecision(new_precision);

				#ifndef BERTINI_DISABLE_ASSERTS
				assert(std::get<mpfr>(current_space_)(0).precision() == current_precision_ && "precision of time in mpfr doesn't match tracker");
				#endif
			}


			/**
			\brief Converts from multiple to different precision multiple precision

			Changes the precision of the internal temporaries to desired precision
			*/
			void MultipleToMultiple(unsigned new_precision)
			{
				MultipleToMultiple( new_precision, std::get<Vec<mpfr> >(current_space_), std::get<mpfr>(current_time_), std::get<mpfr>(delta_t_));
			}



			void AdjustTemporariesPrecision(unsigned new_precision)
			{
				//  the current_space value is adjusted in the appropriate ChangePrecision function
				std::get<Vec<mpfr> >(tentative_space_) = Vec<mpfr>(tracked_system_.NumVariables());
				std::get<Vec<mpfr> >(temporary_space_) = Vec<mpfr>(tracked_system_.NumVariables());

				std::get<mpfr_float>(condition_number_estimate_).precision(new_precision);
				std::get<mpfr_float>(error_estimate_).precision(new_precision);
				std::get<mpfr_float>(norm_J_).precision(new_precision);
				std::get<mpfr_float>(norm_J_inverse_).precision(new_precision);
			}



			

			bool PrecisionSanityCheck() const
			{	
				if (current_precision_==16)
					return true;
				else
				{
					auto checker = [this](bool val){return val==current_precision_};
				
					return checker(tracked_system_.precision()) &&
							checker(std::get<Vec<mpfr> >(current_space_(0))) && 
							checker(std::get<mpfr_float>(condition_number_estimate_));
				}
				
			}


			/////////////////////////////////////////////
			//////////////////////////////////////
			/////////////////////////////
			////////////////////  data members stored in this class
			////////////
			//////
			//


			System tracked_system_;



			////////////
			// state variables
			/////////////


			unsigned current_precision_;

			// tracking the numbers of things
			unsigned num_total_steps_taken_;
			unsigned num_successful_steps_taken_;
			unsigned num_consecutive_successful_steps_;
			unsigned num_failed_steps_;

			// configuration for tracking
			unsigned predictor_order_;

			
			unsigned frequency_of_CN_estimation_;
			unsigned num_steps_since_last_condition_number_computation_;

			mpfr_float tracking_tolerance_;
			mpfr_float path_truncation_threshold_;

			// permanent temporaries
			unsigned next_precision_;
			SuccessCode step_success_code_;

			mpfr current_time_;
			mpfr delta_t_;
			mpfr_float current_stepsize_;

			std::tuple< Vec<dbl>, Vec<mpfr> > current_space_;
			std::tuple< Vec<dbl>, Vec<mpfr> > tentative_space_;
			std::tuple< Vec<dbl>, Vec<mpfr> > temporary_space_;

			std::tuple< double, mpfr_float > condition_number_estimate_;			
			std::tuple< double, mpfr_float > error_estimate_;
			std::tuple< double, mpfr_float > norm_J_;
			std::tuple< double, mpfr_float > norm_J_inverse_;


			config::Stepping stepping_config_;
			config::Newton newton_config_;
			config::Security security_config_;
			config::AdaptiveMultiplePrecisionConfig AMP_config_;
		}; // re: class Tracker


	} // namespace tracking
} // namespace bertini





#endif

