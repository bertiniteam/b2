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

/**
\file tracker.hpp

\brief Contains the base Tracker, and AMPTracker types.
*/

#ifndef BERTINI_TRACKER_LOOP_HPP
#define BERTINI_TRACKER_LOOP_HPP

#include <boost/multiprecision/number.hpp>

#include <algorithm>
#include "tracking/step.hpp"
#include "limbo.hpp"
#include "logging.hpp"


namespace bertini{

	namespace tracking{

		using namespace bertini::tracking::config;
		using std::max;
		using std::min;
		using std::pow;

		using bertini::max;

		// one of these needs the current time...  i think it's mintimeforcurrentprecision...


		mpfr_float MinTimeForCurrentPrecision(unsigned precision)
		{
			mpfr_float temp = pow( mpfr_float("10"), -long(precision)+3);
			if (precision==DoublePrecision)
				return max( temp,mpfr_float("1e-150"));
			else
				return pow( mpfr_float("10"), -long(precision)+3);
		}


		mpfr_float MinStepSizeForCurrentPrecision(unsigned precision)
		{
			mpfr_float temp = pow( mpfr_float("10"), -long(precision)+3);
			if (precision==DoublePrecision)
				return max( temp,mpfr_float("1e-150"));
			else
				return pow( mpfr_float("10"), -long(precision)+3);
		}


		unsigned MinDigitsFromEta(mpfr_float const& eta, mpfr const& current_time)
		{
			return unsigned(ceil(eta)+ ceil(log10(abs(current_time)))+2);
		}


		/**
		\brief The convergence_error function from \cite{AMP2}.  

		Decrease stepsize, then increase precision until the new stepsize is greater than a precision-dependent minimum stepsize.  
	
		This function is used in the AMPTracker loop, when Newton's method fails to converge.

		\see MinStepSizeForCurrentPrecision
		*/
		void ConvergenceError(unsigned & new_precision, mpfr_float & new_stepsize, 
							  unsigned old_precision, mpfr_float const& old_stepsize,
							  mpfr_float const& step_adjustment_factor)
		{
			BOOST_LOG_TRIVIAL(severity_level::trace) << "ConvergenceError";

			new_precision = old_precision;
			new_stepsize = step_adjustment_factor*old_stepsize;

			while (new_stepsize < MinStepSizeForCurrentPrecision(new_precision))
			{
				BOOST_LOG_TRIVIAL(severity_level::trace) << new_precision << " " << MinStepSizeForCurrentPrecision(new_precision);
				if (new_precision==DoublePrecision)
					new_precision=LowestMultiplePrecision;
				else
					new_precision+=PrecisionIncrement;
			}
		}




		/**
		 \brief Compute the cost function from \cite{amp2}, \f$C(P)\f$.  As currently implemented, this is 
		 \f$ 10.35 + 0.13 P \f$ where P is the precision.
		
		 This function tells you the relative cost of arithmetic at a given precision. 

		 \todo Recompute this cost function for boost::multiprecision::mpfr_float
		*/
		mpfr_float Cost(unsigned precision)
		{
			if (precision==DoublePrecision)
				return 1;
			else
				return mpfr_float("10.35") + mpfr_float("0.13") * precision;
		}




		/**
		 \Compute a stepsize satisfying AMP Criterion B with a given precision
		*/
		mpfr_float StepsizeSatisfyingCriterionB(unsigned precision,
										unsigned criterion_B_rhs,
										unsigned num_newton_iterations,
										unsigned predictor_order = 0)
		{
			return pow(mpfr_float(10),-(long(criterion_B_rhs) - long(precision))*num_newton_iterations/(predictor_order+1));
		}




		/**
		 \brief Compute precision and stepsize minimizing the Cost() of tracking.
		 
		 For a given range of precisions, an old stepsize, and a maximum stepsize, the Cost of tracking is computed, and a minimizer found.  
		
		  This function is used in the AMPTracker tracking loop, both in case of successful steps and in Criterion errors.


		 \param[out] new_precision The minimizing precision.
		 \param[out] new_stepsize The minimizing stepsize.
		 \param[in] min_precision The minimum considered precision.
		 \param[in] old_stepsize The previously used stepsize.
		 \param[in] max_precision The maximum considered precision.
		 \param[in] max_stepsize The maximum permitted stepsize.  
		 \param[in] criterion_B_rhs The right hand side of CriterionB from \cite{AMP1,AMP2}
		 \param num_newton_iterations The number of allowed Newton corrector iterations.
		 \param AMP_config The configuration of AMP settings for tracking.
		 \param predictor_order The order of the predictor being used.  This is the order itself, not the order of the error estimate.

		 \see Cost
		*/
 		template <typename RealType>
		void MinimizeCost(unsigned & new_precision, mpfr_float & new_stepsize, 
						  unsigned min_precision, mpfr_float const& old_stepsize,
						  unsigned max_precision, mpfr_float const& max_stepsize,
						  unsigned criterion_B_rhs,
						  unsigned num_newton_iterations,
						  AdaptiveMultiplePrecisionConfig const& AMP_config,
						  unsigned predictor_order = 0)
		{
			BOOST_LOG_TRIVIAL(severity_level::trace) << "MinimizeCost";
			BOOST_LOG_TRIVIAL(severity_level::trace) << "max_stepsize " << max_stepsize ;
			BOOST_LOG_TRIVIAL(severity_level::trace) << " min_precision " << min_precision << " N " << num_newton_iterations << " order " << predictor_order;

			mpfr_float min_cost = Eigen::NumTraits<mpfr_float>::highest();

			new_precision = MaxPrecisionAllowed+1; // initialize to an impossible value.
			new_stepsize = old_stepsize; // initialize to original step size.


			unsigned lowest_mp_precision_to_test = min_precision;

			if (min_precision<=DoublePrecision)
			{
				unsigned candidate_precision = DoublePrecision;
				
				mpfr_float candidate_stepsize = min(
				                                    StepsizeSatisfyingCriterionB(candidate_precision,
																		criterion_B_rhs,
																		num_newton_iterations,
																		predictor_order),
				                                    max_stepsize);
				
				BOOST_LOG_TRIVIAL(severity_level::trace) << "for candidate precision " << candidate_precision << ", stepsize is " << candidate_stepsize;

				mpfr_float current_cost = Cost(candidate_precision) / abs(candidate_stepsize);

				if (current_cost < min_cost)
				{
					min_cost = current_cost;
					new_stepsize = candidate_stepsize;
					new_precision = candidate_precision;
				}

				lowest_mp_precision_to_test = LowestMultiplePrecision;
			}

			if (lowest_mp_precision_to_test < LowestMultiplePrecision)
				lowest_mp_precision_to_test = LowestMultiplePrecision;
			else
			{
				lowest_mp_precision_to_test = (lowest_mp_precision_to_test/PrecisionIncrement) * PrecisionIncrement;
			}

			if (max_precision < LowestMultiplePrecision)
				max_precision = LowestMultiplePrecision;
			else
			{
				max_precision = (max_precision/PrecisionIncrement) * PrecisionIncrement;
			}


			for (unsigned candidate_precision = lowest_mp_precision_to_test; candidate_precision <= max_precision; candidate_precision+=PrecisionIncrement)
			{

				mpfr_float candidate_stepsize = min(
				                                    StepsizeSatisfyingCriterionB(candidate_precision,
																		criterion_B_rhs,
																		num_newton_iterations,
																		predictor_order),
				                                    max_stepsize);
				BOOST_LOG_TRIVIAL(severity_level::trace) << "for candidate precision " << candidate_precision << ", stepsize is " << candidate_stepsize;				
				BOOST_LOG_TRIVIAL(severity_level::trace) << "computed stepsize is " << StepsizeSatisfyingCriterionB(candidate_precision,
																		criterion_B_rhs,
																		num_newton_iterations,
																		predictor_order);

				mpfr_float current_cost = Cost(candidate_precision) / abs(candidate_stepsize);

				if (current_cost < min_cost)
				{
					min_cost = current_cost;
					new_stepsize = candidate_stepsize;
					new_precision = candidate_precision;
				}
			}

			BOOST_LOG_TRIVIAL(severity_level::trace) << "final minimizer of cost: " << new_precision << " " << new_stepsize;
		}





		/**
		\brief Computes the next precision and stepsize for tracking to continue with minimal cost.
		
		Computes the minimum and maximum precisions necessary for computation to continue.
		
		\tparam RealType The real number type being used.
		\param[out] new_precision The computed next precision to use.
		\param[out] new_stepsize The computed next stepsize to use.
		\param current_precision The previously used precision being used.
		\param current_stepsize The current stepsize being used for tracking.
		\param max_stepsize The maximum allowable stepsize.
		\param digits_tracking_tolerance The number of digits required by the tracking tolerance.
		\param criterion_B_rhs The right hand side of CriterionB from \cite{AMP1,AMP2}.
		\param digits_C The number of digits required by CriterionC, same papers.
		\param current_time The current time in tracking.
		\param norm_J An estimate on the norm of the Jacobian.
		\param norm_J_inverse An estimate on the norm of the inverse of the Jacobian.
		\param size_proportion The relative size of the previous step.
		\param eta_min Minimum log10 of stepsize.
		\param eta_max Maximum log10 of stepsize.
		\param max_num_newton_iterations The maximum number of newton iterations to take.
		\param AMP_config The settings used for AMP.
		\param digits_final The number of digits required for tracking during the endgame, perhaps.
		\param predictor_order The order of the predictor being used.  This is the order itself, not the order of the error estimate.
		*/
		template <typename RealType>
		void AMP3_update(unsigned & new_precision, mpfr_float & new_stepsize,
						 unsigned current_precision, mpfr_float const& current_stepsize,
						 mpfr_float const& max_stepsize,
						 unsigned digits_tracking_tolerance,
						 unsigned criterion_B_rhs,
						 unsigned digits_C,
						 mpfr const& current_time,
						 RealType const& norm_J,
						 RealType const& norm_J_inverse,
						 RealType const& size_proportion,
						 mpfr_float const& eta_min,
						 mpfr_float const& eta_max,
						 unsigned max_num_newton_iterations,
						 AdaptiveMultiplePrecisionConfig const& AMP_config,
						 unsigned digits_final = 0,
						 unsigned predictor_order = 0)
		{
			BOOST_LOG_TRIVIAL(severity_level::trace) << "AMP3_update";
			BOOST_LOG_TRIVIAL(severity_level::trace) << current_precision << " " << current_stepsize;

			BOOST_LOG_TRIVIAL(severity_level::trace) << " criterion_B_rhs " << criterion_B_rhs << " N " << max_num_newton_iterations << " digits_final " << digits_final << " order " << predictor_order;

			// compute the number of digits necessary to trigger a lowering.  
			unsigned precision_decrease_threshold;
			if (current_precision<=LowestMultiplePrecision)
				precision_decrease_threshold = current_precision - DoublePrecision; // if in lowest multiple, have to clear the difference down to double
			else 
				precision_decrease_threshold = PrecisionIncrement; // the threshold is the step between precisions in multiple precision

			BOOST_LOG_TRIVIAL(severity_level::trace) << "precision_decrease_threshold " << precision_decrease_threshold;

			unsigned digits_B(ceil(criterion_B_rhs - (predictor_order+1)* eta_min/max_num_newton_iterations));

			unsigned digits_stepsize = max(MinDigitsFromEta(eta_min,current_time), MinDigitsFromEta(eta_max,current_time));

			BOOST_LOG_TRIVIAL(severity_level::trace) << digits_tracking_tolerance << " " << digits_final << " " << digits_B << " " << digits_C << " " << digits_stepsize;


			// if (digits_B < current_precision)
			// 	digits_B = min(digits_B+precision_decrease_threshold, current_precision);

			// if (digits_C < current_precision)
			// 	digits_C = min(digits_C+precision_decrease_threshold, current_precision);

			// if (digits_stepsize < current_precision)
			// 	digits_stepsize = min(digits_stepsize+precision_decrease_threshold, current_precision);

			unsigned min_digits_needed = max(digits_tracking_tolerance, digits_final, digits_B, digits_C, digits_stepsize);

			unsigned max_digits_needed = max(min_digits_needed, unsigned(ceil(criterion_B_rhs - (predictor_order+1)* eta_max/max_num_newton_iterations)));


			if ( (min_digits_needed > current_precision - precision_decrease_threshold) 
					 && 
				  (min_digits_needed < current_precision)
				)
				min_digits_needed = current_precision;

			if ( (max_digits_needed > current_precision - precision_decrease_threshold) 
					 && 
				  (max_digits_needed < current_precision)
				)
				max_digits_needed = current_precision;

			BOOST_LOG_TRIVIAL(severity_level::trace) << min_digits_needed << " " << max_digits_needed;
			MinimizeCost<RealType>(new_precision, new_stepsize, 
						min_digits_needed, current_precision,
						max_digits_needed, max_stepsize,
						criterion_B_rhs,
						max_num_newton_iterations,
						AMP_config,
						predictor_order);


			
		}






		/**
		\brief Computes the next precision and stepsize for tracking to continue with minimal cost.
		
		Computes the minimum and maximum precisions necessary for computation to continue.
		
		\tparam RealType The real number type being used.
		\param[out] new_precision The computed next precision to use.
		\param[out] new_stepsize The computed next stepsize to use.
		\param current_precision The previously used precision being used.
		\param current_stepsize The current stepsize being used for tracking.
		\param current_time The current time in tracking.
		\param norm_J An estimate on the norm of the Jacobian.
		\param norm_J_inverse An estimate on the norm of the inverse of the Jacobian.
		\param size_proportion The relative size of the previous step.
		\param tracking_tolerance The current tracking tolerance.
		\param max_num_newton_iterations The maximum number of newton iterations to take.
		\param step_size_fail_factor The ratio by which to decrease the step size in case of stepping failure.
		\param step_size_success_factor The ratio by which to increase the step size in case of continued success.
		\param AMP_config The settings used for AMP.
		\param digits_final The number of digits required for tracking during the endgame, perhaps.
		\param predictor_order The order of the predictor being used.  This is the order itself, not the order of the error estimate.
		*/
		template <typename RealType>
		void SafetyError(unsigned new_precision, mpfr_float & new_stepsize, 
						 unsigned current_precision, mpfr_float const& current_stepsize, 
						 mpfr const& current_time,
						 RealType const& norm_J, RealType const& norm_J_inverse,
						 RealType const& size_proportion,
						 unsigned num_newton_iterations,
						 mpfr_float const& tracking_tolerance,
						 mpfr_float const& step_size_fail_factor,
						 mpfr_float const& step_size_success_factor,
						 RealType const& norm_of_current_solution,
						 AdaptiveMultiplePrecisionConfig const& AMP_config,
						 unsigned digits_final = 0,
						 unsigned predictor_order= 0)
		{
			mpfr_float minimum_stepsize = MinTimeForCurrentPrecision(current_precision);
			mpfr_float maximum_stepsize = current_stepsize * step_size_fail_factor;

			if (minimum_stepsize > maximum_stepsize)
			{
				// stepsizes are incompatible, must increase precision
				new_precision = LowestMultiplePrecision;
				new_stepsize = current_stepsize * (1+step_size_fail_factor)/2;
			}
			else
			{
				mpfr_float eta_min = -log10(minimum_stepsize);
				mpfr_float eta_max = -log10(maximum_stepsize);

				unsigned criterion_B_rhs(ceil(amp::CriterionBRHS(norm_J, norm_J_inverse, num_newton_iterations, RealType(tracking_tolerance), size_proportion, AMP_config)));

				unsigned digits_C(ceil( amp::CriterionCRHS(norm_J_inverse, norm_of_current_solution, RealType(tracking_tolerance), AMP_config))); 

				unsigned digits_tracking_tolerance(ceil(-log10(tracking_tolerance)));
				

				AMP3_update(new_precision, new_stepsize,
						 current_precision, current_stepsize,
						 maximum_stepsize,
						 digits_tracking_tolerance,
						 criterion_B_rhs,
						 digits_C,
						 current_time,
						 norm_J,
						 norm_J_inverse,
						 size_proportion,
						 eta_min,
						 eta_max,
						 num_newton_iterations,
						 AMP_config,
						  digits_final,
						  predictor_order);
			}
		}


		
		
		/**
		\brief The minimum number of digits based on a log10 of a stepsize.

		\param eta log10 of a stepsize
		\param current_time The current time
		*/
		unsigned MinDigitsFromEta(mpfr_float eta, mpfr current_time)
		{
			return unsigned(ceil(eta) + ceil(log10(abs(current_time))) + 2);
		}



		/**
		\class Tracker

		\brief Base tracker class for trackers offered in Bertini2.

		\see AMPTracker
		*/
		class Tracker
		{

		};

		/** 
		\class AMPTracker

		\brief Functor-like class for tracking paths on a system
		
		
		## Explanation
		
		The bertini::AMPTracker class enables tracking using Adaptive Multiple Precision on an arbitrary square homotopy.  

		The intended usage is to:

		1. Create a system, and instantiate some settings.
		2. Create an AMPTracker, associating it to the system you are going to solve or track on.
		3. Run AMPTracker::Setup, getting the settings in line for tracking.
		4. Repeatedly, or as needed, call the AMPTracker::TrackPath function, feeding it a start point, and start and end times.  The initial precision is that of the start point.  

		Working precision and stepsize are adjusted automatically to get around nearby singularities to the path.  If the endpoint is singular, this may very well fail, as prediction and correction get more and more difficult with proximity to singularities.  

		The TrackPath method is intended to allow the user to track to nonsingular endpoints, or to an endgame boundary, from which an appropriate endgame will be called.
		
		## Some notes

		The AMPTracker has internal state.  That is, it stores the current state of the path being tracked as data members.  After tracking has concluded, these statistics may be extracted.  If you need additional accessors for these data, contact the software authors.

		The class additionally uses the Boost.Log library for logging.  At time of this writing, the trivial logger is being used.  As development continues we will move toward using a more sophisticated logging model.  Suggestions are welcome.

		## Example Usage
		
		Below we demonstrate a basic usage of the AMPTracker class to track a single path.  

		The pattern is as described above: create an instance of the class, feeding it the system to be tracked, and some configuration.  Then, use the tracker to track paths of the system.

		\code{.cpp}
		mpfr_float::default_precision(30); // set initial precision.  This is not strictly necessary.

		using namespace bertini::tracking;

		// 1. Create the system
		Var x = std::make_shared<Variable>("x");
		Var y = std::make_shared<Variable>("y");
		Var t = std::make_shared<Variable>("t");

		System sys;

		VariableGroup v{x,y};

		sys.AddFunction(pow(x,2) + (1-t)*x - 1);
		sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
		sys.AddPathVariable(t);
		sys.AddVariableGroup(v);

		auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
		
		//  2. Create the Tracker object, associating the system to it.
		bertini::tracking::AMPTracker tracker(sys);

		config::Stepping stepping_preferences;
		config::Newton newton_preferences;

		// 3. Get the settings into the tracker 
		tracker.Setup(config::Predictor::Euler,
		              	mpfr_float("1e-5"),
						mpfr_float("1e5"),
						stepping_preferences,
						newton_preferences,
						AMP);
	
		//  4. Create a start and end time.  These are complex numbers.
		mpfr t_start("1.0");
		mpfr t_end("0");
		
		//  5. Create a start point, and container for the end point.
		Vec<mpfr> start_point(2);
		start_point << mpfr("1"), mpfr("1.414");  // set the value of the start point.  This is Eigen syntax.

		Vec<mpfr> end_point;

		// 6. actually do the tracking
		SuccessCode tracking_success = tracker.TrackPath(end_point,
		                  t_start, t_end, start_point);
		
		// 7. and then onto whatever processing you are doing to the computed point.
		\endcode
		
		If this documentation is insufficient, please contact the authors with suggestions, or get involved!  Pull requests welcomed.
		*/
		class AMPTracker : Tracker
		{
		public:



			AMPTracker(System sys) : tracked_system_(sys)
			{	
				BOOST_LOG_TRIVIAL(severity_level::trace) << "creating tracker from system " << sys;
			}


			/**
			\brief Get the tracker set up for tracking.

			Pass the tracker the configuration for tracking, to get it set up.
			*/
			void Setup(config::Predictor new_predictor_choice,
			           mpfr_float const& tracking_tolerance,
						mpfr_float const& path_truncation_threshold,
						config::Stepping const& stepping,
						config::Newton const& newton,
						config::AdaptiveMultiplePrecisionConfig const& AMP_config
						)
			{
				Predictor(new_predictor_choice);
				
				tracking_tolerance_ = tracking_tolerance;
				digits_tracking_tolerance_ = unsigned(ceil(-log10(tracking_tolerance)));

				path_truncation_threshold_ = path_truncation_threshold;

				stepping_config_ = stepping;
				newton_config_ = newton;
				AMP_config_ = AMP_config;
			}






			/**
			\brief Track a start point through time, from a start time to a target time.

			\param[out] solution_at_endtime The value of the solution at the end time.
			\param start_time The time at which to start tracking.
			\param endtime The time to track to.
			\param start_point The intial space values for tracking.
			\return A success code indicating whether tracking was successful.  Will be SuccessCode::Success if was successful, and something else otherwise.
			
			The is the fundamental method for the tracker.  First, you create and set up the tracker, telling it what system you will solve, and the settings to use.  Then, you actually do the tracking.
			*/
			SuccessCode TrackPath(Vec<mpfr> & solution_at_endtime,
									mpfr const& start_time, mpfr const& endtime,
									Vec<mpfr> const& start_point
									)
			{	

				BOOST_LOG_TRIVIAL(severity_level::debug) << "tracking path from t=" << start_time << " to t=" << endtime << " from x = " << start_point;

				if (start_point.size()!=tracked_system_.NumVariables())
				{
					throw std::runtime_error("start point size must match the number of variables in the system to be tracked");
				}

				
				TrackerLoopInitialization(start_time, start_point);
				

				SuccessCode initial_refinement_code = InitialRefinement();
				if (initial_refinement_code!=SuccessCode::Success)
					return initial_refinement_code;


				BOOST_LOG_TRIVIAL(severity_level::trace) << "starting while loop";
				// as precondition to this while loop, the correct container, either dbl or mpfr, must have the correct data.
				while (current_time_!=endtime)
				{
					if (num_successful_steps_taken_ >= stepping_config_.max_num_steps)
						return SuccessCode::MaxNumStepsTaken;
					if (current_stepsize_ < stepping_config_.min_step_size)
						return SuccessCode::MinStepSizeReached;
					if (current_precision_ > AMP_config_.maximum_precision)
						return SuccessCode::MaxPrecisionReached;



					// compute the next delta_t
					if (abs(endtime-current_time_) < abs(current_stepsize_))
						delta_t_ = endtime-current_time_;
					else
						delta_t_ = current_stepsize_ * (endtime - current_time_)/abs(endtime - current_time_);



					if (current_precision_==DoublePrecision)
						step_success_code_ = TrackerIteration<dbl, double>();
					else
						step_success_code_ = TrackerIteration<mpfr, mpfr_float>();



					if (step_success_code_==SuccessCode::Success)
					{	
						BOOST_LOG_TRIVIAL(severity_level::trace) << "tracker iteration successful";

						num_successful_steps_taken_++; 
						num_consecutive_successful_steps_++;
						current_time_ += delta_t_;
					}
					else
					{
						BOOST_LOG_TRIVIAL(severity_level::trace) << "tracker iteration unsuccessful";

						num_consecutive_successful_steps_=0;
						num_failed_steps_taken_++;
					}


					current_stepsize_ = next_stepsize_;
					if (next_precision_!=current_precision_)
					{
						SuccessCode precision_change_code = ChangePrecision(next_precision_);
						if (precision_change_code!=SuccessCode::Success)
							return precision_change_code;
					}
				}// re: while


				CopyFinalSolution(solution_at_endtime);

				return SuccessCode::Success;
			}








			/**
			\brief Change tracker to use a predictor

			\param new_predictor_choice The new predictor to be used.

			\see config::Predictor
			*/
			void Predictor(config::Predictor new_predictor_choice)
			{
				predictor_choice_ = new_predictor_choice;
				predictor_order_ = predict::Order(predictor_choice_);
			}


			/**
			\brief Query the currently set predictor
			*/
			config::Predictor Predictor()
			{
				return predictor_choice_;
			}


			/**
			\brief See how many steps have been taken.

			\return The total number of steps taken, including successes and fails.
			*/
			unsigned NumTotalStepsTaken () const
			{
				return num_failed_steps_taken_ + num_successful_steps_taken_;
			}


		private:

			/**
			\brief Set up the internals of the tracker for a fresh start.  

			Copies the start time, current stepsize, and start point.  Adjusts the current precision to match the precision of the start point.  Zeros counters.

			\param start_time The time at which to start tracking.
			\param start_point The space values from which to start tracking.
			*/
			void TrackerLoopInitialization(mpfr const& start_time,
										   Vec<mpfr> const& start_point)
			{
				// set up the master current time and the current step size
				current_time_.precision(start_time.precision());
				current_time_ = start_time;

				current_stepsize_.precision(stepping_config_.initial_step_size.precision());
				current_stepsize_ = stepping_config_.initial_step_size;

				// populate the current space value with the start point, in appropriate precision
				if (start_point(0).precision()==DoublePrecision)
					MultipleToDouble( start_point);
				else
					MultipleToMultiple(start_point(0).precision(), start_point);
				
				ChangePrecision(start_point(0).precision());

				// reset a bunch of counters to 0.
				num_consecutive_successful_steps_ = 0;
				num_successful_steps_taken_ = 0;
				num_failed_steps_taken_ = 0;
				num_precision_decreases_ = 0;

				// initialize to the frequency so guaranteed to compute it the first try 
				num_steps_since_last_condition_number_computation_ = frequency_of_CN_estimation_;
			}


			/** 
			\brief Run an initial refinement of the start point, to ensure in high enough precision to start.

			\return Whether initial refinement was successful.  
			*/
			SuccessCode InitialRefinement()
			{
				SuccessCode initial_refinement = Refine();
				if (initial_refinement!=SuccessCode::Success)
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << "initial refinement at precision " << current_precision_ << " failed";
					do {
						if (current_precision_ > AMP_config_.maximum_precision)
						{
							BOOST_LOG_TRIVIAL(severity_level::trace) << " singular start point, tracking failed.";
							return SuccessCode::SingularStartPoint;
						}

						if (current_precision_==DoublePrecision)
							initial_refinement = ChangePrecision(LowestMultiplePrecision);
						else
							initial_refinement = ChangePrecision(current_precision_+PrecisionIncrement);
					}
					while (initial_refinement!=SuccessCode::Success);
				}
				BOOST_LOG_TRIVIAL(severity_level::trace) << "initial refinement successful at precision " << current_precision_;

				return SuccessCode::Success;
			}		


			/**
			\brief Copy from the internally stored current solution into a final solution.

			\param[out] solution_at_endtime The solution at the end time
			*/
			void CopyFinalSolution(Vec<mpfr> & solution_at_endtime) const
			{
				// the current precision is the precision of the output solution point.
				if (current_precision_==DoublePrecision)
				{
					unsigned num_vars = tracked_system_.NumVariables();
					solution_at_endtime.resize(num_vars);
					for (unsigned ii=0; ii<num_vars; ii++)
						solution_at_endtime(ii) = mpfr(std::get<Vec<dbl> >(current_space_)(ii));
				}
				else
					solution_at_endtime = std::get<Vec<mpfr> >(current_space_);
			}



			/**
			\brief Run an iteration of AMP tracking.

			\return SuccessCode indicating whether the iteration was successful.
			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode TrackerIteration()
			{	
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				BOOST_LOG_TRIVIAL(severity_level::trace) << "TrackerIteration " << NumTotalStepsTaken() << "\ncurrent_precision: " << current_precision_;


				#ifndef BERTINI_DISABLE_ASSERTS
				assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
				#endif

				Vec<ComplexType>& predicted_space = std::get<Vec<ComplexType> >(temporary_space_); // this will be populated in the Predict step
				Vec<ComplexType>& current_space = std::get<Vec<ComplexType> >(current_space_); // the thing we ultimately wish to update
				ComplexType current_time = ComplexType(current_time_);
				ComplexType delta_t = ComplexType(delta_t_);

				BOOST_LOG_TRIVIAL(severity_level::trace) << "t = " << current_time;
				BOOST_LOG_TRIVIAL(severity_level::trace) << "delta_t = " << delta_t;
				BOOST_LOG_TRIVIAL(severity_level::trace) << "x = " << current_space;

				SuccessCode predictor_code = Predict<ComplexType, RealType>(predicted_space, current_space, current_time, delta_t);

				if (predictor_code==SuccessCode::MatrixSolveFailure)
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << "Preditor, matrix solve failure";

					ConvergenceError();
					return predictor_code;
				}	

				if (predictor_code==SuccessCode::HigherPrecisionNecessary)
				{	
					BOOST_LOG_TRIVIAL(severity_level::trace) << "Preditor, higher precision necessary";

					SafetyError<ComplexType, RealType>();
					return predictor_code;
				}

				BOOST_LOG_TRIVIAL(severity_level::trace) << "predicted_space = " << predicted_space;

				Vec<ComplexType>& tentative_next_space = std::get<Vec<ComplexType> >(tentative_space_); // this will be populated in the Correct step

				ComplexType tentative_next_time = current_time + delta_t;

				SuccessCode corrector_code = Correct<ComplexType, RealType>(tentative_next_space,
													 predicted_space,
													 tentative_next_time);

				BOOST_LOG_TRIVIAL(severity_level::trace) << "condition number: " << std::get<mpfr_float>(condition_number_estimate_);

				if (corrector_code==SuccessCode::MatrixSolveFailure || corrector_code==SuccessCode::FailedToConverge)
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << "corrector, matrix solve failure or failure to converge";

					ConvergenceError();
					return corrector_code;
				}
				else if (corrector_code == SuccessCode::HigherPrecisionNecessary)
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << "corrector, higher precision necessary";

					SafetyError<ComplexType, RealType>();
					return corrector_code;
				}

				// copy the tentative vector into the current space vector;
				current_space = tentative_next_space;
				StepSuccess<ComplexType,RealType>();
				return SuccessCode::Success;
			}





			/**
			\brief Adjust step size and precision due to AMP Criterion violation
			
			This function adjusts internals to the tracker object.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template<typename ComplexType, typename RealType>
			void SafetyError()
			{	
				RealType& norm_J = std::get<RealType>(norm_J_);
				RealType& norm_J_inverse = std::get<RealType>(norm_J_inverse_);
				RealType& size_proportion = std::get<RealType>(size_proportion_);
				RealType norm_z = std::get<Vec<ComplexType> >(current_space_).norm();

				bertini::tracking::SafetyError(next_precision_, next_stepsize_, 
						 current_precision_, current_stepsize_, 
						 current_time_,
						 norm_J, norm_J_inverse,
						 size_proportion,
						 newton_config_.max_num_newton_iterations,
						 	RealType(tracking_tolerance_),
						 	RealType(stepping_config_.step_size_fail_factor),
						 	RealType(stepping_config_.step_size_success_factor),
						 	norm_z,
						 	AMP_config_,
						 	digits_final_,
						  predictor_order_);
			}



			/**
			 Adjust precision and stepsize due to Newton's method failure to converge in the maximum number of steps.
			*/
			void ConvergenceError()
			{
				::bertini::tracking::ConvergenceError(next_precision_, next_stepsize_, 
							  current_precision_, current_stepsize_,
							  stepping_config_.step_size_fail_factor);
			}


			/**
			\brief Wrapper function for calling the correct predictor.
			
			This function computes the next predicted space value, and sets some internals based on the prediction, such as the norm of the Jacobian.

			The real type and complex type must be commensurate.

			\param[out] predicted_space The result of the prediction
			\param current_space The current space point.
			\param current_time The current time value.
			\param delta_t The time differential for this step.  Allowed to be complex.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template<typename ComplexType, typename RealType>
			SuccessCode Predict(Vec<ComplexType> & predicted_space, 
								Vec<ComplexType> const& current_space, 
								ComplexType const& current_time, ComplexType const& delta_t)
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				RealType& norm_J = std::get<RealType>(norm_J_);
				RealType& norm_J_inverse = std::get<RealType>(norm_J_inverse_);
				RealType& size_proportion = std::get<RealType>(size_proportion_);
				RealType& error_estimate = std::get<RealType>(error_estimate_);
				RealType& condition_number_estimate = std::get<RealType>(condition_number_estimate_);

				if (predict::HasErrorEstimate(predictor_choice_))
					return ::bertini::tracking::Predict(predictor_choice_,
									predicted_space,
									error_estimate,
									size_proportion,
									norm_J,
									norm_J_inverse,
									tracked_system_,
									current_space, current_time, 
									delta_t,
									condition_number_estimate,
									num_steps_since_last_condition_number_computation_, 
									frequency_of_CN_estimation_, 
									RealType(tracking_tolerance_),
									AMP_config_);
				else
					return ::bertini::tracking::Predict(predictor_choice_,
									predicted_space,
									size_proportion,
									norm_J,
									norm_J_inverse,
									tracked_system_,
									current_space, current_time, 
									delta_t,
									condition_number_estimate,
									num_steps_since_last_condition_number_computation_, 
									frequency_of_CN_estimation_, 
									RealType(tracking_tolerance_),
									AMP_config_);
			}



			/////////////////
			//
			//  Functions for refining a current space value
			//
			///////////////////////


			/**
			\brief Run Newton's method from the currently stored time and space value, in current precision.

			This overwrites the value of the current space, if successful, with the refined value.

			\return Whether the refinement was successful.
			*/
			SuccessCode Refine()
			{
				SuccessCode code;
				if (current_precision_==DoublePrecision)
				{
					code = Refine(std::get<Vec<dbl> >(temporary_space_),std::get<Vec<dbl> >(current_space_), dbl(current_time_),double(tracking_tolerance_));
					if (code == SuccessCode::Success)
						std::get<Vec<dbl> >(current_space_) = std::get<Vec<dbl> >(temporary_space_);
				}
				else
				{
					code = Refine(std::get<Vec<mpfr> >(temporary_space_),std::get<Vec<mpfr> >(current_space_), current_time_, tracking_tolerance_);
					if (code == SuccessCode::Success)
						std::get<Vec<mpfr> >(current_space_) = std::get<Vec<mpfr> >(temporary_space_);
				}
				return code;
			}







			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			\param[out] new_space The result of running the refinement.
			\param start_point The base point for running Newton's method.
			\param current_time The current time value.
			\param tolerance The tolerance for convergence.  This is a tolerance on \f$\Delta x\f$, not on function residuals.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode Refine(Vec<ComplexType> & new_space,
								Vec<ComplexType> const& start_point, ComplexType const& current_time,
								RealType const& tolerance)
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				return bertini::tracking::Correct(new_space,
							   tracked_system_,
							   start_point,
							   current_time, 
							   tolerance,
							   RealType(path_truncation_threshold_),
							   newton_config_.min_num_newton_iterations,
							   newton_config_.max_num_newton_iterations,
							   AMP_config_);
			}



			/**
			\brief Run Newton's method.

			Wrapper function for calling Correct and getting the error estimates etc directly into the tracker object.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			\param corrected_space[out] The spatial result of the correction loop.
			\param current_space The start point in space for running the corrector loop.
			\param current_time The current time value.

			\return A SuccessCode indicating whether the loop was successful in converging in the max number of allowable newton steps, to the current path tolerance.
			*/
			template<typename ComplexType, typename RealType>
			SuccessCode Correct(Vec<ComplexType> & corrected_space, 
								Vec<ComplexType> const& current_space, 
								ComplexType const& current_time)
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");


				RealType& norm_J = std::get<RealType>(norm_J_);
				RealType& norm_J_inverse = std::get<RealType>(norm_J_inverse_);
				RealType& norm_delta_z = std::get<RealType>(norm_delta_z_);
				RealType& condition_number_estimate = std::get<RealType>(condition_number_estimate_);


				return bertini::tracking::Correct(corrected_space,
												norm_delta_z,
												norm_J,
												norm_J_inverse,
												condition_number_estimate,
												tracked_system_,
												current_space,
												current_time, 
												RealType(tracking_tolerance_),
												RealType(path_truncation_threshold_),
												newton_config_.min_num_newton_iterations,
												newton_config_.max_num_newton_iterations,
												AMP_config_);
			}




			/**
			\brief Increase stepsize or decrease precision, because of consecutive successful steps.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			If the most recent step was successful, maybe adjust down precision and up stepsize.  

			The number of consecutive successful steps is recorded as state in this class, and if this number exceeds a user-determine threshold, the precision or stepsize are allowed to favorably change.  If not, then precision can only go up or remain the same.  Stepsize can only decrease.  These changes depend on the AMP criteria and current tracking tolerance.
			*/
			template <typename ComplexType, typename RealType>
			void StepSuccess()
			{
				BOOST_LOG_TRIVIAL(severity_level::trace) << "StepSuccess";
				BOOST_LOG_TRIVIAL(severity_level::trace) << "current precision: " << current_precision_;
				BOOST_LOG_TRIVIAL(severity_level::trace) << "current stepsize: " << current_stepsize_;

				RealType& norm_J = std::get<RealType>(norm_J_);
				RealType& norm_J_inverse = std::get<RealType>(norm_J_inverse_);
				RealType& size_proportion = std::get<RealType>(size_proportion_);
				RealType norm_of_current_solution = std::get<Vec<ComplexType> > (current_space_).norm();

				mpfr_float minimum_stepsize = current_stepsize_ * stepping_config_.step_size_fail_factor;


				unsigned criterion_B_rhs(round(amp::CriterionBRHS(norm_J, norm_J_inverse, newton_config_.max_num_newton_iterations, RealType(tracking_tolerance_),  size_proportion, AMP_config_)));

				unsigned digits_C(round( amp::CriterionCRHS(norm_J_inverse, norm_of_current_solution, RealType(tracking_tolerance_), AMP_config_))); 
				// hopefully this number is smaller than the current precision, allowing us to reduce precision

				BOOST_LOG_TRIVIAL(severity_level::trace) << "criterion_B_rhs: " << criterion_B_rhs << " digits_C: " << digits_C;

				mpfr_float maximum_stepsize = current_stepsize_;

				if (num_consecutive_successful_steps_ < stepping_config_.consecutive_successful_steps_before_stepsize_increase)
				{	
					digits_C = max(digits_C, current_precision_);
				}
				else if (num_consecutive_successful_steps_ == stepping_config_.consecutive_successful_steps_before_stepsize_increase)
				{
					// disallow precision changing 
					digits_C = max(digits_C, current_precision_);

					// allow stepsize increase, up to maximum step size
					maximum_stepsize = min(current_stepsize_ * stepping_config_.step_size_success_factor, stepping_config_.max_step_size);
				}	
				else if (current_precision_ > DoublePrecision && num_consecutive_successful_steps_>stepping_config_.consecutive_successful_steps_before_stepsize_increase) // exclude the next two cases from being executed in double precision, because cannot decrease from double.
				{
					if (num_consecutive_successful_steps_ < stepping_config_.consecutive_successful_steps_before_precision_decrease)
					{
						digits_C = max(digits_C, current_precision_);
					}
					else if (num_consecutive_successful_steps_ == stepping_config_.consecutive_successful_steps_before_precision_decrease)
					{ // this is the only case where decreasing precision is allowed
						// allow stepsize increase, up to maximum step size
						maximum_stepsize = min(current_stepsize_ * stepping_config_.step_size_success_factor, stepping_config_.max_step_size);
						if (num_precision_decreases_ >= AMP_config_.max_num_precision_decreases)
							digits_C = max(digits_C, current_precision_);
					}
				}
				else
				{
					// reset number of consecutive steps to 0.
					num_consecutive_successful_steps_=0;
					digits_C = max(digits_C, current_precision_);
				}
				BOOST_LOG_TRIVIAL(severity_level::trace) << "maximum_stepsize " << maximum_stepsize << " minimum_stepsize " << minimum_stepsize;

				mpfr_float eta_min = -log10(minimum_stepsize);
				mpfr_float eta_max = -log10(maximum_stepsize);

				AMP3_update(next_precision_, next_stepsize_,
						 current_precision_, current_stepsize_,
						 maximum_stepsize,
						 digits_tracking_tolerance_,
						 criterion_B_rhs,
						 digits_C,
						 current_time_,
						 norm_J,
						 norm_J_inverse,
						 size_proportion,
						 eta_min,
						 eta_max,
						 newton_config_.max_num_newton_iterations,
						 AMP_config_,
						 digits_final_,
						 predictor_order_);

				BOOST_LOG_TRIVIAL(severity_level::trace) << "new precision: " << next_precision_;
				BOOST_LOG_TRIVIAL(severity_level::trace) << "new stepsize: " << next_stepsize_;
			}









			















			/////////////////
			//
			//  Functions for converting between precision types
			//
			///////////////////////


			/**
			Change precision of tracker to next_precision.  Converts the internal temporaries, and adjusts precision of system. Then refines if necessary.

			If the new precision is higher than current precision, a refine step will be called, which runs Newton's method.  This may fail, leaving the tracker in a state with higher precision internals, but garbage digits after the previously known digits.

			\param new_precision The precision to change to.
			\return SuccessCode indicating whether the change was successful.  If the precision increases, and the refinement loop fails, this could be not Success.  Changing down is guaranteed to succeed.
			*/
			SuccessCode ChangePrecision(unsigned new_precision)
			{
				if (new_precision==current_precision_) // no op
					return SuccessCode::Success;

				BOOST_LOG_TRIVIAL(severity_level::debug) << "changing precision to " << new_precision;

				num_steps_since_last_condition_number_computation_ = frequency_of_CN_estimation_;

				bool upsampling_needed = new_precision > current_precision_;

				if (new_precision==DoublePrecision && current_precision_>DoublePrecision)
				{
					// convert from multiple precision to double precision
					MultipleToDouble();
				}
				else if(new_precision > DoublePrecision && current_precision_ == DoublePrecision)
				{
					// convert from double to multiple precision
					DoubleToMultiple(new_precision);
				}
				else
				{
					MultipleToMultiple(new_precision);
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

			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void MultipleToDouble(Vec<mpfr> const& source_point)
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == tracked_system_.NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				#endif

				current_precision_ = DoublePrecision;
				mpfr_float::default_precision(DoublePrecision);

				if (std::get<Vec<dbl> >(current_space_).size()!=source_point.size())
					std::get<Vec<dbl> >(current_space_).resize(source_point.size());

				for (unsigned ii=0; ii<source_point.size(); ii++)
					std::get<Vec<dbl> >(current_space_)(ii) = dbl(source_point(ii));
			}

			/**
			\brief Converts from multiple to double

			Changes the precision of the internal temporaries to double precision
			*/
			void MultipleToDouble()
			{
				MultipleToDouble(std::get<Vec<mpfr> >(current_space_));
			}

			//, std::get<mpfr>(current_time_), std::get<mpfr>(delta_t_)



			/**
			\brief Converts from double to multiple

			Copies a double-precision into the multiple-precision storage vector.  

			You should call Refine after this to populate the new digits with non-garbage data.
	
			\param new_precision The new precision.
			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void DoubleToMultiple(unsigned new_precision, Vec<dbl> const& source_point)
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == tracked_system_.NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				assert(new_precision > DoublePrecision && "must convert to precision higher than DoublePrecision when converting to multiple precision");
				#endif

				current_precision_ = new_precision;
				mpfr_float::default_precision(new_precision);
				tracked_system_.precision(new_precision);

				if (std::get<Vec<mpfr> >(current_space_).size()!=source_point.size())
					std::get<Vec<mpfr> >(current_space_).resize(source_point.size());

				for (unsigned ii=0; ii<source_point.size(); ii++)
					std::get<Vec<mpfr> >(current_space_)(ii) = mpfr(source_point(ii));

				AdjustTemporariesPrecision(new_precision);

				#ifndef BERTINI_DISABLE_ASSERTS
				assert(std::get<Vec<mpfr> >(current_space_)(0).precision() == current_precision_ && "precision of time in mpfr doesn't match tracker");
				#endif
			}



			/**
			\brief Converts from double to multiple

			Copies the double-precision temporaries into the multiple-precision temporaries.  You should call Newton after this to populate the new digits with non-garbage data.

			\param new_precision The new precision.
			*/
			void DoubleToMultiple(unsigned new_precision)
			{
				DoubleToMultiple( new_precision, std::get<Vec<dbl> >(current_space_));
			}





			/**
			\brief Converts from multiple to different precision multiple precision

			Copies a multiple-precision into the multiple-precision storage vector, and changes precision of the time and delta_t.
			Also resets counter so have to re-compute the condition number on next step attempt.

			\param new_precision The new precision.
			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void MultipleToMultiple(unsigned new_precision, Vec<mpfr> const& source_point)
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == tracked_system_.NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				assert(new_precision > DoublePrecision && "must convert to precision higher than DoublePrecision when converting to multiple precision");
				#endif

				current_precision_ = new_precision;
				mpfr_float::default_precision(new_precision);
				tracked_system_.precision(new_precision);

				if (std::get<Vec<mpfr> >(current_space_).size()!=source_point.size())
					std::get<Vec<mpfr> >(current_space_).resize(source_point.size());

				for (unsigned ii=0; ii<source_point.size(); ii++)
					std::get<Vec<mpfr> >(current_space_)(ii) = mpfr(source_point(ii));

				AdjustTemporariesPrecision(new_precision);

				#ifndef BERTINI_DISABLE_ASSERTS
				assert(std::get<Vec<mpfr> >(current_space_)(0).precision() == current_precision_ && "precision of time in mpfr doesn't match tracker");
				#endif
			}


			/**
			\brief Converts from multiple to different precision multiple precision

			Changes the precision of the internal temporaries to desired precision

			\param new_precision The new precision.
			*/
			void MultipleToMultiple(unsigned new_precision)
			{
				MultipleToMultiple( new_precision, std::get<Vec<mpfr> >(current_space_));
			}


			/**
			\brief Change precision of all temporary internal state variables.

			This excludes those which canot be re-written without copying -- the current space point most notably.

			\brief new_precision The new precision to adjust to.
			*/
			void AdjustTemporariesPrecision(unsigned new_precision)
			{
				unsigned num_vars = tracked_system_.NumVariables();

				//  the current_space value is adjusted in the appropriate ChangePrecision function
				std::get<Vec<mpfr> >(tentative_space_).resize(num_vars);
				for (unsigned ii = 0; ii < num_vars; ++ii)
					std::get<Vec<mpfr> >(tentative_space_)(ii).precision(new_precision);

				std::get<Vec<mpfr> >(temporary_space_).resize(num_vars);
				for (unsigned ii = 0; ii < num_vars; ++ii)
					std::get<Vec<mpfr> >(temporary_space_)(ii).precision(new_precision);

				std::get<mpfr_float>(condition_number_estimate_).precision(new_precision);
				std::get<mpfr_float>(error_estimate_).precision(new_precision);
				std::get<mpfr_float>(norm_J_).precision(new_precision);
				std::get<mpfr_float>(norm_J_inverse_).precision(new_precision);
				std::get<mpfr_float>(norm_delta_z_).precision(new_precision);
				std::get<mpfr_float>(size_proportion_).precision(new_precision);
			}



			
			/**
			\brief Ensure that all internal state is in uniform precision.

			\return True if all internal state variables are in the same precision as current_precision_, false otherwise.
			*/
			bool PrecisionSanityCheck() const
			{	
				if (current_precision_==DoublePrecision)
				{
					return true;
				}
				else
				{
					assert(mpfr_float::default_precision()==current_precision_ && "current precision differs from the default precision");

					auto checker = [this](bool val) -> bool { return (val==current_precision_); };

					return tracked_system_.precision() == current_precision_
							&&
							std::get<Vec<mpfr> >(current_space_)(0).precision() == current_precision_ &&
							std::get<Vec<mpfr> >(tentative_space_)(0).precision() == current_precision_ &&
							std::get<Vec<mpfr> >(temporary_space_)(0).precision() == current_precision_ &&
							std::get<mpfr_float>(norm_delta_z_).precision() == current_precision_ &&
							std::get<mpfr_float>(condition_number_estimate_).precision() == current_precision_ &&
							std::get<mpfr_float>(error_estimate_).precision() == current_precision_ &&
							std::get<mpfr_float>(norm_J_).precision() == current_precision_ &&
							std::get<mpfr_float>(norm_J_inverse_).precision() == current_precision_ &&
							std::get<mpfr_float>(size_proportion_).precision() == current_precision_
							        ;
				}
				
			}




			/////////////////////////////////////////////
			//////////////////////////////////////
			/////////////////////////////
			////////////////////  data members stored in this class
			////////////
			//////
			//


			System tracked_system_; ///< The system being tracked.



			////////////
			// state variables
			/////////////


			unsigned current_precision_; ///< The current precision of the tracker, the system, and all temporaries.

			// tracking the numbers of things
			unsigned num_total_steps_taken_; ///< The number of steps taken, including failures and successes.
			unsigned num_successful_steps_taken_;  ///< The number of successful steps taken so far.
			unsigned num_consecutive_successful_steps_; ///< The number of CONSECUTIVE successful steps taken in a row.
			unsigned num_failed_steps_taken_; ///< The total number of failed steps taken.
			unsigned num_precision_decreases_; ///< The number of times precision has decreased this track.


			// configuration for tracking
			config::Predictor predictor_choice_; ///< The predictor to use while tracking.
			unsigned predictor_order_; ///< The order of the predictor -- one less than the error estimate order.

			
			unsigned frequency_of_CN_estimation_; ///< How frequently the condition number should be re-estimated.
			unsigned num_steps_since_last_condition_number_computation_; ///< How many steps have passed since the most recent condition number estimate.


			unsigned digits_final_ = 0; ///< The number of digits to track to, due to being in endgame zone.
			unsigned digits_tracking_tolerance_; ///< The number of digits required for tracking to given tolerance, condition number notwithstanding.
			mpfr_float tracking_tolerance_; ///< The tracking tolerance.
			mpfr_float path_truncation_threshold_; ///< The threshold for path truncation.

			// permanent temporaries
			mpfr_float next_stepsize_; /// The next stepsize
			unsigned next_precision_; ///< The next precision
			SuccessCode step_success_code_; ///< The code for step success.

			mpfr current_time_; ///< The current time.
			mpfr delta_t_; ///< The current delta_t.
			mpfr_float current_stepsize_; ///< The current stepsize.

			std::tuple< Vec<dbl>, Vec<mpfr> > current_space_; ///< The current space value. 
			std::tuple< Vec<dbl>, Vec<mpfr> > tentative_space_; ///< After correction, the tentative next space value
			std::tuple< Vec<dbl>, Vec<mpfr> > temporary_space_; ///< After prediction, the tentative next space value.


			

			std::tuple< double, mpfr_float > condition_number_estimate_; ///< An estimate on the condition number of the Jacobian		
			std::tuple< double, mpfr_float > error_estimate_; ///< An estimate on the error of a step.
			std::tuple< double, mpfr_float > norm_J_; ///< An estimate on the norm of the Jacobian
			std::tuple< double, mpfr_float > norm_J_inverse_;///< An estimate on the norm of the inverse of the Jacobian
			std::tuple< double, mpfr_float > norm_delta_z_; ///< The norm of the change in space resulting from a step.
			std::tuple< double, mpfr_float > size_proportion_; ///< The proportion of the space step size, taking into account the order of the predictor.
 
			config::Stepping stepping_config_; ///< The stepping configuration.
			config::Newton newton_config_; ///< The newton configuration.
			config::Security security_config_; ///< The security settings.
			config::AdaptiveMultiplePrecisionConfig AMP_config_; ///< The Adaptive Multiple Precision settings.
		}; // re: class Tracker


	} // namespace tracking
} // namespace bertini





#endif

