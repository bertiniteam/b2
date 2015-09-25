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

			template <typename T> using PrecStep = std::pair<unsigned, T>

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



		/** 
		\brief Functor-like class for tracking paths on a system

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
				current_precision = start_point.precision(); // get the current precision.

				num_consecutive_successful_steps = 0;
				num_successful_steps_taken = 0;
				num_failed_steps = 0;


				x_mpfr = start_point;
				x_dbl = Vec<dbl>(start_point.size());


				if (current_precision==16)
					MultipleToDouble( start_point, start_time, stepping.max_step_size);


				while (num_steps_taken < stepping.max_num_steps)
				{
					// as precondition to this while loop, the correct container, either dbl or mpfr, must have the correct data.
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
						num_failed_steps++;
					}


					if (next_precision!=current_precision)
					{
						SuccessCode precision_change_code = ChangePrecision(next_precision);
						if (precision_change_code!=SuccessCode::Success)
							return precision_change_code;
					}

				}// re: while


				return SuccessCode::Success;
			}

		private:


			/**
			\brief Converts from double to multiple

			Copies a double-precision into the multiple-precision storage vector.  You should call Newton after this to populate the new digits with non-garbage data.
			*/
			void DoubleToMultiple(unsigned new_precision, Vec<dbl> const& source_point, dbl const& source_time, dbl const& source_delta_t)
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == system.NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				assert(new_precision > 16 && "must convert to precision higher than 16 when converting to multiple precision");
				#endif

				current_precision = new_precision;
				mpfr_float.precision(new_precision);
				system.precision(new_precision);

				t_mpfr = mpfr(source_time);
				delta_t_mpfr = mpfr(source_delta_t);

				for (unsigned ii=0; ii<source_point.size())
					x_mpfr(ii) = mpfr(source_point(ii));

				#ifndef BERTINI_DISABLE_ASSERTS
				assert(t_mpfr.precision() == current_precision)
				#endif
			}


			/**
			\brief Converts from multiple to double

			Copies a multiple-precision into the double storage vector, and changes precision of the time and delta_t.
			*/
			void MultipleToDouble(Vec<mpfr> const& source_point, mpfr const& source_time, mpfr const& source_delta_t)
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == system.NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				#endif

				current_precision = new_precision;
				system.precision(16);

				t_dbl = dbl(source_time);
				delta_t_dbl = dbl(source_delta_t);

				for (unsigned ii=0; ii<source_point.size())
					x_dbl(ii) = dbl(source_point(ii));
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

				current_precision = new_precision;
				mpfr_float.precision(new_precision);
				system.precision(new_precision);

				t_mpfr = mpfr(source_time);
				delta_t_mpfr = mpfr(source_delta_t);

				for (unsigned ii=0; ii<source_point.size())
					x_mpfr(ii) = mpfr(source_point(ii));


				#ifndef BERTINI_DISABLE_ASSERTS
				assert(t_mpfr.precision() == current_precision)
				#endif
			}






			/**
			Change precision of tracker to next_precision.  Converts the internal temporaries, and adjusts precision of system. 
			*/
			SuccessCode ChangePrecision(unsigned next_precision)
			{
				if (next_precision==current_precision)
					return SuccessCode::Success;

				num_steps_since_last_condition_number_computation = frequency_of_CN_estimation;

				if (next_precision==16 && current_precision>16)
				{
					// convert from multiple precision to double precision
					MultipleToDouble( x_mpfr, t_mpfr, delta_t_mpfr);
					upsampling_needed = false;
				}
				else if(next_precision > 16 && current_precision == 16)
				{
					// convert from double to multiple precision
					DoubleToMultiple( x_dbl, t_dbl, delta_t_dbl);
					upsampling_needed = true;					
				}
				else if(next_precision < current_precision)
				{
					MultipleToMultiple(x_mpfr, t_mpfr, delta_t_mpfr);
					upsampling_needed = false;
				}
				else if (next_precision > current_precision)
				{
					MultipleToMultiple(x_mpfr, t_mpfr, delta_t_mpfr);
					upsampling_needed = true;
				}

				current_precision = next_precision;

				if (upsampling_needed)
					return = Sharpen(x_mpfr, t_mpfr, some_tolerance_based_on_precision);
				else
					return SuccessCode::Success;
			}


			void Predictor(config::Predictor new_predictor_choice)
			{
				predictor_choice = new_predictor_choice;
				predictor_order = PredictorOrder(predictor_choice);
			}





			void Sharpen()
			{

			}

			/////////////////////////////////////////////
			//////////////////////////////////////
			/////////////////////////////
			////////////////////  data members stored in this class
			////////////
			//////
			//


			System tracked_system;



			////////////
			// state variables
			/////////////


			unsigned current_precision;

			// tracking the numbers of things
			unsigned num_total_steps_taken;
			unsigned num_successful_steps_taken;
			unsigned num_consecutive_successful_steps;
			unsigned num_failed_steps;



			// configuration for tracking
			unsigned predictor_order;


			double condition_number_estimate_dbl;
			mpfr_float condition_number_estimate_mpfr;
			
			unsigned frequency_of_CN_estimation;
			unsigned num_steps_since_last_condition_number_computation;


			// permanent temporaries
			dbl t_dbl;
			mfpr t_mpfr;

			dbl delta_t_dbl;
			mpfr delta_t_mpfr;
			
			Vec<dbl> x_dbl;
			Vec<mpfr> x_mpfr;
		};



		







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

