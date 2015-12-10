//This file is part of Bertini 2.0.
//
//base_tracker.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//base_tracker.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with base_tracker.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  base_tracker.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015

/**
\file base_tracker.hpp

\brief 

\brief Contains the abstract base Tracker type, from which all other Trackers inherit.
*/

#ifndef BERTINI_BASE_TRACKER_HPP
#define BERTINI_BASE_TRACKER_HPP

#include <algorithm>
#include "tracking/step.hpp"
#include "limbo.hpp"
#include "logging.hpp"


namespace bertini{

	namespace tracking{

		/**
		\class Tracker

		\brief Base tracker class for trackers offered in Bertini2.
	
		\see AMPTracker
		
		## Using a tracker

		Trackers in Bertini2 are the engine for tracking a path from one space-time pair to another.  The path is implicitly described by the system being tracked.

		1. Create a system.
		2. Create a tracker, associating it to the system.
		3. Set the config for the track.
		4. Track from your start point and time, to the target time.
		5. Profit.

		Specific examples are given with the implemented tracker types.  So far, these types are avaiable:

		* AMPTracker

		Please see their detailed documentation for description of how to use them correctly.


		## Purpose 

		Since the Bertini trackers have common functionality, and we want to be able to call arbitrary algorithms using and tracker type, we use inheritance.  That is, there is common functionality in all trackers, such as

		* Setup
		* TrackPath
		* Refine
		
		which we want all Trackers to be able to do.  However, the internal behaviour of a particular tracker varies -- which is why it is a different type.  In particular, the fixed precision trackers produce and work in a fixed precision, whereas the AMPTracker varies precision to get around or near branch points while tracking. 

		Hence, the use of trackers in Bertini2 is through pointers or references to Trackers, enabling the use of any kind of tracking in any algorithm, and further allowing the development of new tracker types as the theory and practice advance.


		## Creating a new tracker type

		To create a new Tracker type, inherit from this, and override the following functions:

		\code
		
		public:
		SuccessCode Refine(Vec<mpfr> & new_space,
							Vec<mpfr> const& start_point, mpfr const& current_time) override
		{}
		
		private:

		void TrackerLoopInitialization(mpfr const& start_time, Vec<mpfr> const& start_point) override
		{}

		SuccessCode InitialRefinement() override
		{}

		SuccessCode PreIterationCheck() const override
		{}

		SuccessCode TrackerIteration() override
		{}

		void CopyFinalSolution(Vec<mpfr> & solution_at_endtime) const override
		{}
		
		\endcode
	
		and optionally override the following functions

		\code
		void ResetCounters() override
		{}

		void PostTrackCleanup() override
		{}
		\endcode
		where you probably want to call this base function, which is why it is protected, not private.



		*/
		class Tracker
		{

		public:

			Tracker(System const& sys) : tracked_system_(sys)
			{
				Predictor(predict::DefaultPredictor());
			}



			/**
			\brief Get the tracker set up for tracking.

			Pass the tracker the configuration for tracking, to get it set up.
			*/
			void Setup(config::Predictor new_predictor_choice,
			           mpfr_float const& tracking_tolerance,
						mpfr_float const& path_truncation_threshold,
						config::Stepping const& stepping,
						config::Newton const& newton
						)
			{
				Predictor(new_predictor_choice);
				
				tracking_tolerance_ = tracking_tolerance;
				digits_tracking_tolerance_ = unsigned(ceil(-log10(tracking_tolerance)));

				path_truncation_threshold_ = path_truncation_threshold;

				stepping_config_ = stepping;
				newton_config_ = newton;
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
									) const
			{	

				BOOST_LOG_TRIVIAL(severity_level::debug) << "tracking path from t=" << start_time << " to t=" << endtime << " from x = " << start_point;

				if (start_point.size()!=tracked_system_.NumVariables())
					throw std::runtime_error("start point size must match the number of variables in the system to be tracked");

				
				TrackerLoopInitialization(start_time, start_point);
				

				SuccessCode initial_refinement_code = InitialRefinement();
				if (initial_refinement_code!=SuccessCode::Success)
					return initial_refinement_code;


				BOOST_LOG_TRIVIAL(severity_level::trace) << "starting while loop";
				// as precondition to this while loop, the correct container, either dbl or mpfr, must have the correct data.
				while (current_time_!=endtime)
				{	
					SuccessCode pre_iteration_code = PreIterationCheck();
					if (pre_iteration_code!=SuccessCode::Success)
					{
						PostTrackCleanup();
						return pre_iteration_code;
					}


					// compute the next delta_t
					if (abs(endtime-current_time_) < abs(current_stepsize_))
						delta_t_ = endtime-current_time_;
					else
						delta_t_ = current_stepsize_ * (endtime - current_time_)/abs(endtime - current_time_);


					step_success_code_ = TrackerIteration();


					if (step_success_code_==SuccessCode::Success)
					{	
						BOOST_LOG_TRIVIAL(severity_level::trace) << "tracker iteration successful";

						num_successful_steps_taken_++; 
						num_consecutive_successful_steps_++;
						current_time_ += delta_t_;
						num_consecutive_failed_steps_ = 0;
					}
					else if (step_success_code_==SuccessCode::GoingToInfinity)
					{
						PostTrackCleanup();
						return step_success_code_;
					}
					else
					{
						BOOST_LOG_TRIVIAL(severity_level::trace) << "tracker iteration unsuccessful";

						num_consecutive_successful_steps_=0;
						num_failed_steps_taken_++;
						num_consecutive_failed_steps_++;
					}
				}// re: while


				CopyFinalSolution(solution_at_endtime);
				PostTrackCleanup();
				return SuccessCode::Success;
			}


			/**
			\brief Refine a point given in multiprecision.
			
			Runs Newton's method using the current settings for tracking, including the min and max number of iterations allowed, the tracking tolerance, precision, etc.  YOU must ensure that the input point has the correct precision.

			\return The SuccessCode indicating whether the refinement completed.  

			\param new_space The result of refinement.
			\param start_point The seed for Newton's method for refinement.
			\param current_time The current time value for refinement.
			*/
			virtual
			SuccessCode Refine(Vec<mpfr> & new_space,
								Vec<mpfr> const& start_point, mpfr const& current_time) const = 0;


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
			config::Predictor Predictor() const
			{
				return predictor_choice_;
			}


			/**
			\brief get a const reference to the system.
			*/
			const class System& System() const
			{
				return tracked_system_;
			}

			/**
			\brief See how many steps have been taken.

			\return The total number of steps taken, including successes and fails.
			*/
			unsigned NumTotalStepsTaken () const
			{
				return num_failed_steps_taken_ + num_successful_steps_taken_;
			}

			virtual ~Tracker() = default;

		private:

			

			/**
			\brief Set up initialization of the internals for tracking a path.

			\param start_time The time at which to start tracking.
			\param start_point The point from which to start tracking.
			*/
			virtual
			void TrackerLoopInitialization(mpfr const& start_time, Vec<mpfr> const& start_point) const = 0;

			/**
			\brief Ensure that any pre-checks on precision or accuracy of start point pass.

			\return Anything but Success will terminate tracking.
			*/
			virtual 
			SuccessCode InitialRefinement() const = 0;

			/**
			\brief Check internal state for whether tracking should continue.  

			\return Code for whether to go on.  Tracking will terminate if the returned value is not Success.
			*/
			virtual 
			SuccessCode PreIterationCheck() const = 0;

			/**
			\brief A single iteration of the tracker loop.

			\return Whether the tracker loop was successful or not.  Incrementing of counters for the base class happens automatically.
			*/
			virtual 
			SuccessCode TrackerIteration() const = 0;

			/**
			\brief Copy the solution from whatever internal variable it is stored in, into the output variable.

			\param solution_at_endtime The output variable into which to copy the final solution.
			*/
			virtual
			void CopyFinalSolution(Vec<mpfr> & solution_at_endtime) const = 0;


		protected:

			/**
			\brief Function to be called before exiting the tracker loop.
			*/
			virtual
			void PostTrackCleanup() const 
			{}

			/**
			\brief Reset counters used during tracking.

			Your custom tracker type should almost certainly call this function.
			*/
			virtual
			void ResetCounters() const
			{
				// reset a bunch of counters to 0.
				num_consecutive_successful_steps_ = 0;
				num_successful_steps_taken_ = 0;
				num_failed_steps_taken_ = 0;
				num_consecutive_failed_steps_ = 0;
				num_total_steps_taken_ = 0;
			}


			const class System& tracked_system_; ///< The system being tracked.

			// tracking the numbers of things
			mutable unsigned num_total_steps_taken_; ///< The number of steps taken, including failures and successes.
			mutable unsigned num_successful_steps_taken_;  ///< The number of successful steps taken so far.
			mutable unsigned num_consecutive_successful_steps_; ///< The number of CONSECUTIVE successful steps taken in a row.
			mutable unsigned num_consecutive_failed_steps_; ///< The number of CONSECUTIVE failed steps taken in a row. 
			mutable unsigned num_failed_steps_taken_; ///< The total number of failed steps taken.

			
			// configuration for tracking
			config::Predictor predictor_choice_; ///< The predictor to use while tracking.
			unsigned predictor_order_; ///< The order of the predictor -- one less than the error estimate order.

			config::Stepping stepping_config_; ///< The stepping configuration.
			config::Newton newton_config_; ///< The newton configuration.
			config::Security security_config_; ///< The security settings.



			unsigned digits_final_ = 0; ///< The number of digits to track to, due to being in endgame zone.
			unsigned digits_tracking_tolerance_; ///< The number of digits required for tracking to given tolerance, condition number notwithstanding.
			mpfr_float tracking_tolerance_; ///< The tracking tolerance.
			mpfr_float path_truncation_threshold_; ///< The threshold for path truncation.


			mutable mpfr current_time_; ///< The current time.
			mutable mpfr delta_t_; ///< The current delta_t.
			mutable mpfr_float current_stepsize_; ///< The current stepsize.


			// permanent temporaries
			mutable mpfr_float next_stepsize_; /// The next stepsize
			mutable SuccessCode step_success_code_; ///< The code for step success.
		};



	} // re: namespace tracking
} // re: namespace bertini


#endif


