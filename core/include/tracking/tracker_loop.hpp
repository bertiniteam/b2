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


		enum 
		{
			PrecisionIncrement = 10
		}


		enum class Branch
		{
			A, B
		};

		template <typename T>
		using PrecStep = std::pair<unsigned, T>
		

		template <typename RealType>
		RealType epsilon(unsigned precision)
		{
			return pow(10, -precision+3);
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
					new_precision+=30;
				else
					new_precision+=PrecisionIncrement;
		}


		template <typename ComplexType> 
		void MinimizeCost(unsigned & new_precision, ComplexType & new_stepsize, 
						  unsigned old_precision, ComplexType const& old_stepsize,
						  unsigned max_precision)
		{
			std::set G;

			unsigned base_precision = old_precision;


			if (old_precision==16)
			{
				unsigned candidate_precision = 16;
				
				compute |candidate_stepsize| that satisfies Eq 9 as an equation
				
				if (min(abs(candidate_stepsize),abs(old_stepsize)) > epsilon<RealType>(candidate_precision))
					G.insert(PrecStep<RealType>(candidate_precision, candidate_stepsize));

				base_precision = 30;
			}

			for (unsigned candidate_precision = base_precision; base_precision <= max_precision; base_precision+=PrecisionIncrement)
			{

				compute |candidate_stepsize| that satisfies Eq 9 as an equation

				if (min(abs(candidate_stepsize),abs(old_stepsize)) > epsilon<RealType>(candidate_precision))
					G.insert(PrecStep<RealType>(candidate_precision, candidate_stepsize));
			}

			if (G.empty())
			{
				new_precision = max_precision + 1;
				new_stepsize = old_precision;
			}
			else
			{
				// find the minimizer
			}
		}






		template <typename RealType, typename ComplexType>
		Branch SafetyError(unsigned new_precision, RealType const& new_stepsize, Vec<ComplexType> const& new_delta_z, ComplexType const& new_delta_t,
						 unsigned starting_precision, RealType const& starting_stepsize, Vec<ComplexType> const& delta_z, ComplexType const& delta_t)
		{
			new_precision = starting_precision;
			new_stepsize = starting_stepsize;
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




		inline 
		SuccessCode TrackPath(Vec<mpfr> & solution_at_endtime,
								mpfr const& start_time, mpfr const& endtime,
								Vec<mpfr> const& start_point,
								System & sys,
								config::Predictor predictor_choice,
								config::PrecisionType prec_type = config::PrecisionType::Adaptive,
								mpfr_float const& tracking_tolerance,
								mpfr_float const& path_truncation_threshold,
								config::Stepping const& stepping,
								config::Newton const& newton,
								config::Security const& security,
								config::AdaptiveMultiplePrecisionConfig const& AMP_config 
								)
		{

			unsigned current_precision = start_point.precision(); // get the current precision.


			unsigned num_consecutive_successful_steps = 0;
			unsigned num_steps_taken = 0;
			


			unsigned predictor_order = PredictorOrder(predictor_choice);




			RealType condition_number_estimate;
			unsigned num_steps_since_last_condition_number_computation = frequency_of_CN_estimation;
			// initialize to the frequency so guaranteed to compute it the first try 

			// create temporaries for use in either tracking mode, double or multiple.
			Vec<mpfr> x_mpfr(start_point.size());
			Vec<dbl> x_dbl(start_point.size());

			mfpr t_mpfr;
			dbl t_dbl;

			mpfr delta_t_mpfr;
			dbl delta_t_dbl

			if (current_precision==16)
			{
				t_dbl = dbl(start_time);
				delta_t_dbl = dbl(stepping.max_step_size);

				for (unsigned ii=0; ii<start_point.size())
					x_dbl(ii) = dbl(start_point(ii));
			}
			else
			{
				t_mpfr = mpfr(start_time);
				delta_t_mpfr = mpfr(stepping.max_step_size);

				for (unsigned ii=0; ii<start_point.size())
					x_mpfr(ii) = mpfr(start_point(ii));

				mpfr_float.default_precision(current_precision);
				sys.precision(current_precision);
			}



			while (num_steps_taken < stepping.max_num_steps)
			{
				SuccessCode step_success_code;
				if (current_precision==16)
					step_success_code = TrackerIteration<mpfr>();
				else
					step_success_code = TrackerIteration<dbl>();



				if (step_success_code==SuccessCode::Success)
				{
					num_steps_taken++; 
					num_consecutive_successful_steps++;
					current_time += delta_t;
				}
				else
				{
					num_consecutive_successful_steps=0;

					
				}


				if (next_precision==16 && current_precision>16)
				{
					// convert from multiple precision to double precision
				}
				else if(next_precision > 16 && current_precision == 16)
				{
					// convert from double to multiple precision

					upsampling_needed = true;					
				}


				if (next_precision - current_precision > 0)
				{
					// sharpen here!!!
				}

			}// re: while


			return SuccessCode::Success;
		}







		template <typename RealType, typename ComplexType>
		SuccessCode TrackerIteration(next_space, next_time, 
									next_delta_t, next_precision,

									current_space, current_time, 
									delta_t, 
									//precision implied by the number type
									predictor_choice, 
									condition_number_estimate,
									num_steps_since_last_condition_number_computation,
									frequency_of_CN_estimation, prec_type,
									tracking_tolerance, AMP_config
									)
		{

			Vec<ComplexType> tentative_next_space;

			SuccessCode predictor_code = Predict(predictor_choice,
											tentative_next_space,
											sys,
											current_space, current_time, 
											delta_t,
											condition_number_estimate,
											num_steps_since_last_condition_number_computation, 
											frequency_of_CN_estimation, prec_type, 
											tracking_tolerance,
											AMP_config);

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

			ComplexType tentative_next_time = current_time + delta_t;

			SuccessCode corrector_code = Correct(next_space,
										   sys,
										   current_space, // pass by value to get a copy of it
										   next_time, 
										   prec_type, 
										   tracking_tolerance,
										   path_truncation_threshold,
										   min_num_newton_iterations,
										   max_num_newton_iterations,
										   AMP_config);

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


	} // namespace tracking
} // namespace bertini





#endif

