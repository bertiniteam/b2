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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file base_tracker.hpp

\brief Contains the abstract base Tracker type, from which all other Trackers inherit.
*/

#ifndef BERTINI_BASE_TRACKER_HPP
#define BERTINI_BASE_TRACKER_HPP

#include <algorithm>
//#include "bertini2/tracking/step.hpp"
#include "bertini2/trackers/ode_predictors.hpp"
#include "bertini2/trackers/newton_corrector.hpp"
#include "bertini2/logging.hpp"

#include "bertini2/detail/observable.hpp"

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
			public Observable,
			public detail::Configured<
				typename TrackerTraits< D >::NeededConfigs
					>
		{
			using NeededTypes = typename TrackerTraits< D >::NeededTypes;
			using BaseComplexT = typename TrackerTraits<D>::BaseComplexT;
			using BaseRealT = typename TrackerTraits<D>::BaseRealT;

			using ComplexT = BaseComplexT;
			using RealT = BaseRealT;


		public:
			using Config = detail::Configured< typename TrackerTraits< D >::NeededConfigs >;  ///< The configuration-holding base, carrying the config types this tracker needs.
			using Stepping = SteppingConfig;  ///< The stepping (step-size adjustment) configuration type.
			using Newton = NewtonConfig;  ///< The Newton corrector configuration type.
			using PrecConf = typename TrackerTraits< D >::PrecisionConfig;  ///< The precision configuration type for this tracker.

			/// \brief Construct a tracker associated to a system to be tracked.
			Tracker(System const& sys) : tracked_system_(std::ref(sys)),
				predictor_(predict::DefaultPredictor(), sys),
				corrector_(sys)
			{
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
				corrector_.Settings(newton);

				SetTrackingTolerance(tracking_tolerance);

				path_truncation_threshold_ = path_truncation_threshold;

				this->template Set<SteppingConfig>(stepping);
				this->template Set<NewtonConfig>(newton);

				current_stepsize_ = BaseRealT(stepping.initial_step_size);
			}


			/// \brief Bring the configuration getter into scope, for querying config settings.
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


			/// \brief Set the threshold beyond which a path is judged to be going to infinity (and truncated).
			/// \param tol The new (strictly positive) truncation threshold.
			void SetInfiniteTruncationTolerance(double const& tol)
			{
				if (tol <= 0)
					throw std::runtime_error("truncation threshold must be strictly positive");

				path_truncation_threshold_ = tol;
			}


			/**
			\brief (Re)draw the random probe direction used for condition-number estimation.

			Call this once, from the algorithm, at the start of tracking a point -- after any per-path
			RNG reseed and before TrackPath -- so the direction is fixed for that whole point's track
			(pre-endgame tracking and the endgame both) and is reproducible across runs/ranks.  Do NOT
			call it per TrackPath: the endgame issues many TrackPath calls per point, and refreshing
			the direction mid-track perturbs condition estimates and churns the RNG.
			*/
			void RefreshConditionDirection()
			{
				corrector_.RefreshRandomDirection();
				// Share the SAME per-path probe with the predictor, so the predictor's and
				// corrector's ||J^{-1}|| estimates use one consistent direction (ADR-0024).
				predictor_.SetConditionProbe(corrector_.ConditionProbe());
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
			SuccessCode TrackPath(Vec<ComplexT> & solution_at_endtime,
									ComplexT const& start_time, ComplexT const& endtime,
									Vec<ComplexT> const& start_point
									) const
			{
				if (start_point.size()!=static_cast<Eigen::Index>(GetSystem().NumVariables()))
					throw std::runtime_error("start point size must match the number of variables in the system to be tracked");



				SuccessCode initialization_code = TrackerLoopInitialization(start_time, endtime, start_point);
				if (initialization_code!=SuccessCode::Success)
				{
					PostTrackCleanup();
					return initialization_code;
				}

				NotifyObservers(TrackingStarted<typename TrackerTraits<D>::EventEmitterType>(static_cast<typename TrackerTraits<D>::EventEmitterType const&>(*this)));

				// as precondition to this while loop, the correct container, either complex_dbl or mpfr, must have the correct data.
				while (!IsSymmRelDiffSmall(current_time_,endtime_, Eigen::NumTraits<ComplexT>::epsilon()))
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
				predictor_.PredictorMethod(new_predictor_choice);
				predictor_order_ = predictor_.Order();
			}


			/**
			\brief Query the currently set predictor
			*/
			Predictor GetPredictor() const
			{
				return predictor_.PredictorMethod();
			}


			/**
			\brief get a const reference to the system.
			*/
			void SetSystem(const System & new_sys)
			{
				tracked_system_ = std::ref(new_sys);
				predictor_.ChangeSystem(tracked_system_);
				corrector_.ChangeSystem(tracked_system_);
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
			void SetStepSize(RealT const& new_stepsize) const
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

			/// \brief Get the currently set tracking tolerance.
			auto TrackingTolerance() const
			{
				return tracking_tolerance_;
			}

			/// \brief Get the currently set infinite-truncation (path) threshold.
			auto InfiniteTruncationTolerance() const
			{
				return path_truncation_threshold_;
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
			SuccessCode TrackerLoopInitialization(ComplexT const& start_time, ComplexT const& end_time, Vec<ComplexT> const& start_point) const = 0;


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
			void CopyFinalSolution(Vec<ComplexT> & solution_at_endtime) const = 0;

			// virtual
			// void CopyFinalSolution(Vec<complex_dbl> & solution_at_endtime) const = 0;


		protected:





			/// \brief Check, at a given complex type, whether the current space value exceeds the truncation threshold.
			/// \tparam ComplexT The complex number type at which to perform the check.
			/// \return SuccessCode::GoingToInfinity if the dehomogenized norm exceeds the threshold, else SuccessCode::Success.
			template <typename ComplexT>
			SuccessCode CheckGoingToInfinity() const
			{
				if (GetSystem().DehomogenizePoint(std::get<Vec<ComplexT> >(current_space_)).norm() > path_truncation_threshold_)
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




			/// \brief Reset the base-class step counters (successes, failures, consecutive runs, total) to zero.
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

			/// \brief Hook invoked after a successful tracker iteration; derived types update their own state here.
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


			/// \brief Hook invoked after a failed tracker iteration; derived types update their own state here.
			virtual
			void OnStepFail() const = 0;

			/**
			\brief Check whether the path is going to infinity, as it tracks.

			This check is necessary because a homotopy may be malformed, or may have encountered a probability-0 event.  That it is a 0 probability event is why this check is disable-able via a toggle.
			*/
			virtual
			SuccessCode CheckGoingToInfinity() const = 0;

			/// \brief Hook invoked when a path is truncated for going to infinity; derived types react here.
			virtual
			void OnInfiniteTruncation() const = 0;


			std::reference_wrapper<const System> tracked_system_; ///< Reference to the system being tracked.

			bool infinite_path_truncation_ = true; ///< Whether should check if the path is going to infinity while tracking.  On by default.
			bool reinitialize_stepsize_ = true; ///< Whether should re-initialize the stepsize with each call to Trackpath.  On by default.

			// tracking the numbers of things
			mutable unsigned num_total_steps_taken_; ///< The number of steps taken, including failures and successes.
			mutable unsigned num_successful_steps_taken_;  ///< The number of successful steps taken so far.
			mutable unsigned num_consecutive_successful_steps_; ///< The number of CONSECUTIVE successful steps taken in a row.
			mutable unsigned num_consecutive_failed_steps_; ///< The number of CONSECUTIVE failed steps taken in a row.
			mutable unsigned num_failed_steps_taken_; ///< The total number of failed steps taken.


			// configuration for tracking
			//
			// predictor and corrector are held BY VALUE so that copying a tracker
			// deep-copies them.  (They were previously shared_ptr, which made every
			// tracker copy share one predictor/corrector with its source — unusable
			// from multiple threads.)  Both hold only work buffers and settings; they
			// store no reference to the System, so plain memberwise copy is correct.
			// They are mutable for the same reason as the state members above: they
			// hold scratch space mutated during the logically-const TrackPath.
			mutable predict::ExplicitRKPredictor predictor_; ///< The predictor to use while tracking.
			unsigned predictor_order_; ///< The order of the predictor -- one less than the error estimate order.

			mutable correct::NewtonCorrector corrector_;  ///< The Newton corrector used to refine predicted points.



			unsigned digits_final_ = 0; ///< The number of digits to track to, due to being in endgame zone.
			unsigned digits_tracking_tolerance_ = 5; ///< The number of digits required for tracking to given tolerance, condition number notwithstanding.
			NumErrorT tracking_tolerance_ = 1e-5; ///< The tracking tolerance.
			NumErrorT path_truncation_threshold_ = 1e5; ///< The threshold for path truncation.

			mutable ComplexT endtime_; ///< The time we are tracking to.
			mutable ComplexT current_time_; ///< The current time.
			mutable ComplexT delta_t_; ///< The current delta_t.
			mutable RealT current_stepsize_; ///< The current stepsize.


			// permanent temporaries
			mutable RealT next_stepsize_; ///< The next stepsize.
			mutable SuccessCode step_success_code_; ///< The code for step success.



			mutable unsigned num_steps_since_last_condition_number_computation_; ///< How many steps have passed since the most recent condition number estimate.
			mutable unsigned num_successful_steps_since_stepsize_increase_; ///< How many successful steps have been taken since increased stepsize.

			using TupOfVec = typename NeededTypes::ToTupleOfVec;  ///< A tuple of complex vectors, one per numeric type the tracker uses.
			using TupOfReal = typename NeededTypes::ToTupleOfReal;  ///< A tuple of real vectors, one per numeric type the tracker uses.

			mutable TupOfVec current_space_; ///< The current space value.
			mutable TupOfVec tentative_space_; ///< After correction, the tentative next space value
			mutable TupOfVec temporary_space_; ///< After prediction, the tentative next space value.


			/// Metadata from the most recent predict or correct step (norms, condition number,
			/// size proportion, error estimate, norm of the Newton step).  Filled by the
			/// predictor/corrector via their StepMetadata out-parameter.
			mutable StepMetadata last_step_;




			public:


			/// \brief Get the condition-number estimate from the most recent predict/correct step.
			NumErrorT LatestConditionNumber() const
			{
				return this->last_step_.condition_number_estimate;
			}


			/// \brief Get the error estimate from the most recent predict/correct step.
			NumErrorT LatestErrorEstimate() const
			{
				return this->last_step_.error_estimate;
			}


			/// \brief Get the norm of the Newton step from the most recent predict/correct step.
			NumErrorT LatestNormOfStep() const
			{
				return this->last_step_.norm_delta_z;
			}

			/// \brief Turn the going-to-infinity path truncation check on or off.
			void SetInfiniteTruncation(bool b)
			{
				infinite_path_truncation_ = b;
			}

			/// \brief Query whether the going-to-infinity path truncation check is enabled.
			auto InfiniteTruncation()
			{
				return infinite_path_truncation_;
			}

			/// \brief Get the number of variables in the system being tracked.
			unsigned NumVariables() const
			{
				return static_cast<unsigned>(GetSystem().NumVariables());
			}

			/// \brief Get the current time value of the track.
			auto CurrentTime() const
			{
				return current_time_;
			}

			/// \brief Get the current time-step increment (delta t).
			auto DeltaT() const
			{
				return delta_t_;
			}

			/// \brief Get the current step size.
			auto CurrentStepsize() const
			{
				return current_stepsize_;
			}


			/// \brief Get the current space point of the track (in the derived tracker's working type).
			virtual Vec<ComplexT> CurrentPoint() const = 0;


			/// \brief Get the current working precision of the tracker.
			virtual unsigned CurrentPrecision() const = 0;
		};



	} // re: namespace tracking
} // re: namespace bertini


#endif


