//This file is part of Bertini 2.
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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file base_tracker.hpp

\brief 

\brief Contains the abstract base Tracker type, from which all other Trackers inherit.
*/

#ifndef BERTINI_BASE_TRACKER_HPP
#define BERTINI_BASE_TRACKER_HPP

#include <algorithm>
//#include "bertini2/tracking/step.hpp"
#include "bertini2/trackers/ode_predictors.hpp"
#include "bertini2/trackers/newton_corrector.hpp"
#include "bertini2/limbo.hpp"
#include "bertini2/logging.hpp"
#include "bertini2/detail/visitable.hpp"

// Must be at the end of the include list
#include "bertini2/trackers/events.hpp"

#include "bertini2/detail/is_template_parameter.hpp"
#include "bertini2/detail/configured.hpp"

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
		template<class D>
		class Tracker : 
			public Observable<>,
			public detail::Configured<
				typename TrackerTraits< D >::NeededConfigs
					>
		{
			using NeededTypes = typename TrackerTraits< D >::NeededTypes;
			using BaseComplexType = typename TrackerTraits<D>::BaseComplexType;
			using BaseRealType = typename TrackerTraits<D>::BaseRealType;

			using CT = BaseComplexType;
			using RT = BaseRealType;
			
			
		public:
			using Config = detail::Configured< typename TrackerTraits< D >::NeededConfigs >;
			using Stepping = SteppingConfig;
			using Newton = NewtonConfig;
			using PrecConf = typename TrackerTraits< D >::PrecisionConfig;

			Tracker(System const& sys) : tracked_system_(std::ref(sys))
			{
				predictor_ = std::make_shared< predict::ExplicitRKPredictor >(predict::DefaultPredictor(), tracked_system_);
				corrector_ = std::make_shared< correct::NewtonCorrector >(tracked_system_);
				SetPredictor(predict::DefaultPredictor());
			}




			/**
			\brief Get the tracker set up for tracking.

			Pass the tracker the configuration for tracking, to get it set up.
			*/
			void Setup(Predictor new_predictor_choice,
			           double const& tracking_tolerance,
						double const& path_truncation_threshold,
						SteppingConfig const& stepping,
						NewtonConfig const& newton)
			{
				SetPredictor(new_predictor_choice);
				corrector_->Settings(newton);
				
				SetTrackingTolerance(tracking_tolerance);

				path_truncation_threshold_ = path_truncation_threshold;

				this->template Set(stepping);
				this->template Set(newton);

				current_stepsize_ = BaseRealType(stepping.initial_step_size);
			}


			using Config::Get;



			/**
			\brief Set how tightly to track the path.

			Smaller values tracks the path more tightly, at the expense of higher computational time.

			\param tracking_tolerance The new value.  Newton iterations are performed until the step length is less than this number, or the max number of iterations has been reached, in which case the overall predict-correct step is viewed as a failure, and the step is undone.  This number must be positive.
			*/
			void SetTrackingTolerance(double const& tracking_tolerance)
			{
				if (tracking_tolerance <= 0)
					throw std::runtime_error("tracking tolerance must be strictly positive");

				tracking_tolerance_ = tracking_tolerance;
				digits_tracking_tolerance_ = NumTraits<double>::TolToDigits(tracking_tolerance);
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
			SuccessCode TrackPath(Vec<CT> & solution_at_endtime,
									CT const& start_time, CT const& endtime,
									Vec<CT> const& start_point
									) const
			{	
				if (start_point.size()!=GetSystem().NumVariables())
					throw std::runtime_error("start point size must match the number of variables in the system to be tracked");

				
				
				SuccessCode initialization_code = TrackerLoopInitialization(start_time, endtime, start_point);
				if (initialization_code!=SuccessCode::Success)
				{
					PostTrackCleanup();
					return initialization_code;
				}

				// as precondition to this while loop, the correct container, either dbl or mpfr, must have the correct data.
				while (!IsSymmRelDiffSmall(current_time_,endtime_, Eigen::NumTraits<CT>::epsilon()))
				{	
					SuccessCode pre_iteration_code = PreIterationCheck();
					if (pre_iteration_code!=SuccessCode::Success)
					{
						PostTrackCleanup();
						return pre_iteration_code;
					}

					using std::abs;
					// compute the next delta_t
					if (abs(endtime_-current_time_) < abs(current_stepsize_))
						delta_t_ = endtime_-current_time_;
					else
						delta_t_ = current_stepsize_ * (endtime_ - current_time_)/abs(endtime_ - current_time_);


					step_success_code_ = TrackerIteration();

					if (infinite_path_truncation_ && (CheckGoingToInfinity()==SuccessCode::GoingToInfinity))
					{	
						OnInfiniteTruncation();
						PostTrackCleanup();
						return SuccessCode::GoingToInfinity;
					}
					else if (step_success_code_==SuccessCode::Success)
						OnStepSuccess();
					else
						OnStepFail();

				}// re: while


				CopyFinalSolution(solution_at_endtime);
				PostTrackCleanup();
				return SuccessCode::Success;
			}




			/**
			\brief Refine a point to tolerance implicitly internally set.
			
			Runs Newton's method using the current settings for tracking, including the min and max number of iterations allowed, the tracking tolerance, precision, etc.  YOU must ensure that the input point has the correct precision.

			\return The SuccessCode indicating whether the refinement completed.  

			\param[out] new_space The result of refinement.
			\param start_point The seed for Newton's method for refinement.
			\param current_time The current time value for refinement.
			*/
			template<typename C>
			SuccessCode Refine(Vec<C> & new_space,
								Vec<C> const& start_point, C const& current_time) const
			{

				static_assert(detail::IsTemplateParameter<C,NeededTypes>::value,"complex type for refinement must be a used type for the tracker");
				return this->AsDerived().RefineImpl(new_space, start_point, current_time);
			}




			/**
			\brief Refine a point to a given tolerance.
			
			Runs Newton's method using the current settings for tracking, including the min and max number of iterations allowed, precision, etc, EXCEPT for the tracking tolerance and max number of iterations, which you feed in here.  YOU must ensure that the input point has the correct precision.

			\return The SuccessCode indicating whether the refinement completed.  

			\param[out] new_space The result of refinement.
			\param start_point The seed for Newton's method for refinement.
			\param current_time The current time value for refinement.
			\param tolerance The tolerance to which to refine.
			\param max_iterations The maximum number of iterations to use to refine.
			*/
			template<typename C>
			SuccessCode Refine(Vec<C> & new_space,
								Vec<C> const& start_point, C const& current_time, double const& tolerance, unsigned max_iterations) const
			{
				static_assert(detail::IsTemplateParameter<C,NeededTypes>::value,"complex type for refinement must be a used type for the tracker");

				return this->AsDerived().RefineImpl(new_space, start_point, current_time, tolerance, max_iterations);
			}


			/**
			\brief Change tracker to use a predictor

			\param new_predictor_choice The new predictor to be used.

			\see Predictor
			*/
			void SetPredictor(Predictor new_predictor_choice)
			{
				predictor_->PredictorMethod(new_predictor_choice);
				predictor_order_ = predictor_->Order();
			}


			/**
			\brief Query the currently set predictor
			*/
			Predictor GetPredictor() const
			{
				return predictor_->PredictorMethod();
			}


			/**
			\brief get a const reference to the system.
			*/
			void SetSystem(const System & new_sys)
			{
				tracked_system_ = std::ref(new_sys);
				predictor_->ChangeSystem(tracked_system_);
				corrector_->ChangeSystem(tracked_system_);
			}

			/**
			\brief get a const reference to the system.
			*/
			const System& GetSystem() const
			{
				return tracked_system_.get();
			}

			/**
			\brief See how many steps have been taken.

			\return The total number of steps taken, including successes and fails.
			*/
			unsigned NumTotalStepsTaken () const
			{
				return num_failed_steps_taken_ + num_successful_steps_taken_;
			}

			/**
			\brief Set how large the stepsize should be.

			\param new_stepsize The new value.
			*/ 
			void SetStepSize(RT const& new_stepsize) const
			{
				current_stepsize_ = new_stepsize;
			}


			/**
			\brief Switch resetting of initial step size to that of the stepping settings.

			By default, initial step size is retrieved from the stepping settings at the start of each path track.  To turn this off, and re-use the previous step size from the previously tracked path, turn off by calling this function with false.
			*/
			void ReinitializeInitialStepSize(bool should_reinitialize_stepsize)
			{
				reinitialize_stepsize_ = should_reinitialize_stepsize;
			}
			
			virtual ~Tracker() = default;

			auto TrackingTolerance() const
			{
				return tracking_tolerance_;
			}

		private:

			// convert the base tracker into the derived type.
			const D& AsDerived() const
			{
				return static_cast<const D&>(*this);
			}

			/**
			\brief Set up initialization of the internals for tracking a path.

			\param start_time The time at which to start tracking.
			\param end_time The time to which to track.
			\param start_point The point from which to start tracking.
			*/
			virtual
			SuccessCode TrackerLoopInitialization(CT const& start_time, CT const& end_time, Vec<CT> const& start_point) const = 0;
			

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
			void CopyFinalSolution(Vec<CT> & solution_at_endtime) const = 0;

			// virtual
			// void CopyFinalSolution(Vec<dbl> & solution_at_endtime) const = 0;


		protected:


			


			template <typename ComplexType>
			SuccessCode CheckGoingToInfinity() const
			{
				if (GetSystem().DehomogenizePoint(std::get<Vec<ComplexType> >(current_space_)).norm() > path_truncation_threshold_)
					return SuccessCode::GoingToInfinity;
				else
					return SuccessCode::Success;
			}



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
			void ResetCounters() const = 0;




			void ResetCountersBase() const
			{
				// reset a bunch of counters to 0.
				num_consecutive_successful_steps_ = 0;
				num_successful_steps_taken_ = 0;
				num_failed_steps_taken_ = 0;
				num_consecutive_failed_steps_ = 0;
				num_total_steps_taken_ = 0;
			}


			/**
			\brief Increment and reset counters after a successful TrackerIteration()

			Your custom override, if provided, should almost certainly call this function.
			*/
			void IncrementBaseCountersSuccess() const
			{
				num_successful_steps_taken_++; 
				num_consecutive_successful_steps_++;
				current_time_ += delta_t_;
				num_consecutive_failed_steps_ = 0;
			}

			virtual
			void OnStepSuccess() const = 0;


			/**
			\brief Increment and reset counters after a failed TrackerIteration()

			Your custom override, if provided, should almost certainly call this function.
			*/
			void IncrementBaseCountersFail() const
			{
				num_consecutive_successful_steps_=0;
				num_failed_steps_taken_++;
				num_consecutive_failed_steps_++;
			}


			virtual
			void OnStepFail() const = 0;

			/**
			\brief Check whether the path is going to infinity, as it tracks.  

			This check is necessary because a homotopy may be malformed, or may have encountered a probability-0 event.  That it is a 0 probability event is why this check is disable-able via a toggle.
			*/
			virtual 
			SuccessCode CheckGoingToInfinity() const = 0;

			virtual 
			void OnInfiniteTruncation() const = 0;


			std::reference_wrapper<const System> tracked_system_; ///< Reference to the system being tracked.

			bool infinite_path_truncation_ = true; /// Whether should check if the path is going to infinity while tracking.  On by default.
			bool reinitialize_stepsize_ = true; ///< Whether should re-initialize the stepsize with each call to Trackpath.  On by default.

			// tracking the numbers of things
			mutable unsigned num_total_steps_taken_; ///< The number of steps taken, including failures and successes.
			mutable unsigned num_successful_steps_taken_;  ///< The number of successful steps taken so far.
			mutable unsigned num_consecutive_successful_steps_; ///< The number of CONSECUTIVE successful steps taken in a row.
			mutable unsigned num_consecutive_failed_steps_; ///< The number of CONSECUTIVE failed steps taken in a row. 
			mutable unsigned num_failed_steps_taken_; ///< The total number of failed steps taken.

			
			// configuration for tracking
			std::shared_ptr<predict::ExplicitRKPredictor > predictor_; // The predictor to use while tracking
			unsigned predictor_order_; ///< The order of the predictor -- one less than the error estimate order.

			std::shared_ptr<correct::NewtonCorrector> corrector_;



			unsigned digits_final_ = 0; ///< The number of digits to track to, due to being in endgame zone.
			unsigned digits_tracking_tolerance_; ///< The number of digits required for tracking to given tolerance, condition number notwithstanding.
			NumErrorT tracking_tolerance_; ///< The tracking tolerance.
			NumErrorT path_truncation_threshold_; ///< The threshold for path truncation.

			mutable CT endtime_; ///< The time we are tracking to.
			mutable CT current_time_; ///< The current time.
			mutable CT delta_t_; ///< The current delta_t.
			mutable RT current_stepsize_; ///< The current stepsize.


			// permanent temporaries
			mutable RT next_stepsize_; /// The next stepsize
			mutable SuccessCode step_success_code_; ///< The code for step success.



			mutable unsigned num_steps_since_last_condition_number_computation_; ///< How many steps have passed since the most recent condition number estimate.
			mutable unsigned num_successful_steps_since_stepsize_increase_; ///< How many successful steps have been taken since increased stepsize.
			mutable unsigned num_successful_steps_since_precision_decrease_; ///< The number of successful steps since decreased precision.

			using TupOfVec = typename NeededTypes::ToTupleOfVec;
			using TupOfReal = typename NeededTypes::ToTupleOfReal;
			
			mutable TupOfVec current_space_; ///< The current space value. 
			mutable TupOfVec tentative_space_; ///< After correction, the tentative next space value
			mutable TupOfVec temporary_space_; ///< After prediction, the tentative next space value.


			mutable NumErrorT condition_number_estimate_; ///< An estimate on the condition number of the Jacobian		
			mutable NumErrorT error_estimate_; ///< An estimate on the error of a step.
			mutable NumErrorT norm_J_; ///< An estimate on the norm of the Jacobian
			mutable NumErrorT norm_J_inverse_;///< An estimate on the norm of the inverse of the Jacobian
			mutable NumErrorT norm_delta_z_; ///< The norm of the change in space resulting from a step.
			mutable NumErrorT size_proportion_; ///< The proportion of the space step size, taking into account the order of the predictor.




			public: 

			
			NumErrorT LatestConditionNumber() const
			{
				return this->condition_number_estimate_;
			}

			
			NumErrorT LatestErrorEstimate() const
			{
				return this->error_estimate_;
			}

			
			NumErrorT LatestNormOfStep() const
			{
				return this->norm_delta_z_;
			}


			unsigned NumVariables() const
			{
				return GetSystem().NumVariables();
			}

			auto CurrentTime() const
			{
				return current_time_;
			}

			auto DeltaT() const
			{
				return delta_t_;
			}

			auto CurrentStepsize() const
			{
				return current_stepsize_;
			}


			virtual Vec<CT> CurrentPoint() const = 0;


			virtual unsigned CurrentPrecision() const = 0;
		};



	} // re: namespace tracking
} // re: namespace bertini


#endif


