//This file is part of Bertini 2.
//
//fixed_precision_tracker.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fixed_precision_tracker.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fixed_precision_tracker.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire


/**
\file fixed_precision_tracker.hpp

*/

#ifndef BERTINI_FIXED_PRECISION_TRACKER_HPP
#define BERTINI_FIXED_PRECISION_TRACKER_HPP

#include "bertini2/trackers/base_tracker.hpp"


namespace bertini{

	namespace tracking{

		using std::max;
		using std::min;
		using std::pow;

		using bertini::max;


	




		/** 
		\class FixedPrecisionTracker<prec>

		\brief Functor-like class for tracking paths on a system
		*/
		template<class DerivedT>
		class FixedPrecisionTracker : public Tracker<FixedPrecisionTracker<DerivedT>>
		{	
		public:

			using BaseComplexType = typename TrackerTraits<DerivedT>::BaseComplexType;
			using BaseRealType = typename TrackerTraits<DerivedT>::BaseRealType;

			using CT = BaseComplexType;
			using RT = BaseRealType;

			virtual ~FixedPrecisionTracker() = default;

			using EmitterType = FixedPrecisionTracker<DerivedT>;
			using Base = Tracker<FixedPrecisionTracker<DerivedT>>;

			using Config =  typename Base::Config;
			FORWARD_GET_CONFIGURED
			using Stepping = typename Base::Stepping;
			using Newton = typename Base::Newton;

			FixedPrecisionTracker(System const& sys) : Base(sys){}

			/**
			\brief An additional no-op call, provided for conformity of interface with AMP tracker in generic code.
			*/
			void PrecisionSetup(FixedPrecisionConfig const&)
			{ }


			Vec<CT> CurrentPoint() const override
			{
				return std::get<Vec<CT>>(this->current_space_);
			}



			void ResetCounters() const override
			{
				Base::ResetCountersBase();

				this->num_successful_steps_since_stepsize_increase_ = 0;
				// initialize to the frequency so guaranteed to compute it the first try 	
				this->num_steps_since_last_condition_number_computation_ = this->Get<Stepping>().frequency_of_CN_estimation;
			}


			/**
			\brief Ensure that number of steps, stepsize, and precision still ok.

			\return Success if ok to keep going, and a different code otherwise. 
			*/
			SuccessCode PreIterationCheck() const override
			{
				if (this->num_successful_steps_taken_ >= Get<Stepping>().max_num_steps)
					return SuccessCode::MaxNumStepsTaken;
				if (this->current_stepsize_ < Get<Stepping>().min_step_size)
					return SuccessCode::MinStepSizeReached;

				return SuccessCode::Success;
			}






			
			


			void PostTrackCleanup() const override
			{
				this->NotifyObservers(TrackingEnded<EmitterType>(*this));
			}

			/**
			\brief Copy from the internally stored current solution into a final solution.
			
			If preservation of precision is on, this function first returns to the initial precision.

			\param[out] solution_at_endtime The solution at the end time
			*/
			void CopyFinalSolution(Vec<CT> & solution_at_endtime) const override
			{

				// the current precision is the precision of the output solution point.

				unsigned num_vars = this->GetSystem().NumVariables();
				solution_at_endtime.resize(num_vars);
				for (unsigned ii=0; ii<num_vars; ii++)
				{
					solution_at_endtime(ii) = std::get<Vec<CT> >(this->current_space_)(ii);
				}

			}




			/**
			\brief Run an iteration of the tracker loop.

			Predict and correct, adjusting precision and stepsize as necessary.

			\return Success if the step was successful, and a non-success code if something went wrong, such as a linear algebra failure or AMP Criterion violation.
			*/
			SuccessCode TrackerIteration() const override
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RT>::Real, 
			              				typename Eigen::NumTraits<CT>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				this->NotifyObservers(NewStep<EmitterType >(*this));

				Vec<CT>& predicted_space = std::get<Vec<CT> >(this->temporary_space_); // this will be populated in the Predict step
				Vec<CT>& current_space = std::get<Vec<CT> >(this->current_space_); // the thing we ultimately wish to update
				CT current_time = CT(this->current_time_);
				CT delta_t = CT(this->delta_t_);

				SuccessCode predictor_code = Predict(predicted_space, current_space, current_time, delta_t);

				if (predictor_code!=SuccessCode::Success)
				{
					this->NotifyObservers(FirstStepPredictorMatrixSolveFailure<EmitterType >(*this));

					this->next_stepsize_ = RT(Get<Stepping>().step_size_fail_factor)*this->current_stepsize_;

					UpdateStepsize();

					return predictor_code;
				}

				this->NotifyObservers(SuccessfulPredict<EmitterType , CT>(*this, predicted_space));

				Vec<CT>& tentative_next_space = std::get<Vec<CT> >(this->tentative_space_); // this will be populated in the Correct step

				CT tentative_next_time = current_time + delta_t;

				SuccessCode corrector_code = Correct(tentative_next_space,
													 predicted_space,
													 tentative_next_time);

				if (corrector_code == SuccessCode::GoingToInfinity)
				{
					// there is no corrective action possible...
					return corrector_code;
				}
				else if (corrector_code!=SuccessCode::Success)
				{
					this->NotifyObservers(CorrectorMatrixSolveFailure<EmitterType >(*this));

					this->next_stepsize_ = RT(Get<Stepping>().step_size_fail_factor)*this->current_stepsize_;
					UpdateStepsize();

					return corrector_code;
				}

				
				this->NotifyObservers(SuccessfulCorrect<EmitterType , CT>(*this, tentative_next_space));

				// copy the tentative vector into the current space vector;
				current_space = tentative_next_space;
				return SuccessCode::Success;
			}


			/**
			Check whether the path is going to infinity.
			*/
			SuccessCode CheckGoingToInfinity() const override
			{
				return Base::template CheckGoingToInfinity<CT>();
			}

			



			/**
			\brief Commit the next stepsize, and adjust internals.
			*/
			SuccessCode UpdateStepsize() const
			{
				this->SetStepSize(this->next_stepsize_);
				return SuccessCode::Success;
			}





			///////////
			//
			//  overrides for counter adjustment after a TrackerIteration()
			//
			////////////////

			/**
			\brief Increment and reset counters after a successful TrackerIteration()
			*/
			void OnStepSuccess() const override
			{
				Base::IncrementBaseCountersSuccess();
				this->NotifyObservers(SuccessfulStep<EmitterType >(*this));
			}

			/**
			\brief Increment and reset counters after a failed TrackerIteration()
			*/
			void OnStepFail() const override
			{
				Base::IncrementBaseCountersFail();
				this->num_successful_steps_since_stepsize_increase_ = 0;
				this->NotifyObservers(FailedStep<EmitterType >(*this));
			}



			void OnInfiniteTruncation() const override
			{
				this->NotifyObservers(InfinitePathTruncation<EmitterType>(*this));
			}

			//////////////
			//
			//


			/**
			\brief Wrapper function for calling the correct predictor.
			
			This function computes the next predicted space value, and sets some internals based on the prediction, such as the norm of the Jacobian.

			The real type and complex type must be commensurate.

			\param[out] predicted_space The result of the prediction
			\param current_space The current space point.
			\param current_time The current time value.
			\param delta_t The time differential for this step.  Allowed to be complex.
			*/
			SuccessCode Predict(Vec<CT> & predicted_space, 
								Vec<CT> const& current_space, 
								CT const& current_time, CT const& delta_t) const
			{

				return this->predictor_->Predict(
			                predicted_space,
							this->tracked_system_,
							current_space, current_time,
							delta_t,
							this->condition_number_estimate_,
							this->num_steps_since_last_condition_number_computation_,
							Get<Stepping>().frequency_of_CN_estimation,
							this->tracking_tolerance_);
			}



			/**
			\brief Run Newton's method.

			Wrapper function for calling Correct and getting the error estimates etc directly into the tracker object.

			\param corrected_space[out] The spatial result of the correction loop.
			\param current_space The start point in space for running the corrector loop.
			\param current_time The current time value.

			\return A SuccessCode indicating whether the loop was successful in converging in the max number of allowable newton steps, to the current path tolerance.
			*/
			SuccessCode Correct(Vec<CT> & corrected_space, 
								Vec<CT> const& current_space, 
								CT const& current_time) const
			{
				return this->corrector_->Correct(corrected_space,
												this->tracked_system_,
												current_space,
												current_time, 
												this->tracking_tolerance_,
												Get<Newton>().min_num_newton_iterations,
												Get<Newton>().max_num_newton_iterations);
			}



			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.  The tolerance is the tracking tolerance specified during Setup(...).


			\param[out] new_space The result of running the refinement.
			\param start_point The base point for running Newton's method.
			\param current_time The current time value.

			\return Code indicating whether was successful or not.  Regardless, the value of new_space is overwritten with the correction result.
			*/
			SuccessCode RefineImpl(Vec<CT> & new_space,
								Vec<CT> const& start_point, CT const& current_time) const
			{
				return this->corrector_->Correct(new_space,
							   this->tracked_system_,
							   start_point,
							   current_time, 
							   this->tracking_tolerance_,
							   Get<Newton>().min_num_newton_iterations,
							   Get<Newton>().max_num_newton_iterations);
			}






			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.


			\param[out] new_space The result of running the refinement.
			\param start_point The base point for running Newton's method.
			\param current_time The current time value.
			\param tolerance The tolerance for convergence.  This is a tolerance on \f$\Delta x\f$, not on function residuals.
			\param max_iterations The maximum number of permitted Newton iterations.  

			\return Code indicating whether was successful or not.  Regardless, the value of new_space is overwritten with the correction result.
			*/
			SuccessCode RefineImpl(Vec<CT> & new_space,
								Vec<CT> const& start_point, CT const& current_time,
								NumErrorT const& tolerance, unsigned max_iterations) const
			{
				return this->corrector_->Correct(new_space,
							   this->tracked_system_,
							   start_point,
							   current_time, 
							   tolerance,
							   1,
							   max_iterations);
			}


			/////////////////////////////////////////////
			//////////////////////////////////////
			/////////////////////////////
			////////////////////  data members stored in this class
			////////////
			//////
			//

			// no additional state variables needed for the FixedPrecision base tracker types

		};




		class DoublePrecisionTracker : public FixedPrecisionTracker<DoublePrecisionTracker>
		{
		public:
			using BaseComplexType = dbl;
			using BaseRealType = double;

			using EmitterType = typename TrackerTraits<DoublePrecisionTracker>::EventEmitterType;

			BERTINI_DEFAULT_VISITABLE()


			/**
			\brief Construct a tracker, associating to it a System.
			*/
			DoublePrecisionTracker(class System const& sys) : FixedPrecisionTracker<DoublePrecisionTracker>(sys)
			{	}


			DoublePrecisionTracker() = delete;

			virtual ~DoublePrecisionTracker() = default;


			unsigned CurrentPrecision() const override
			{
				return DoublePrecision();
			}


			/**
			\brief Set up the internals of the tracker for a fresh start.  

			Copies the start time, current stepsize, and start point.  Adjusts the current precision to match the precision of the start point.  Zeros counters.

			\param start_time The time at which to start tracking.
			\param end_time The time to which to track.
			\param start_point The space values from which to start tracking.
			*/
			SuccessCode TrackerLoopInitialization(BaseComplexType const& start_time,
			                               BaseComplexType const& end_time,
										   Vec<BaseComplexType> const& start_point) const override
			{
				this->NotifyObservers(Initializing<EmitterType,BaseComplexType>(*this,start_time, end_time, start_point));

				// set up the master current time and the current step size
				this->current_time_ = start_time;
				this->endtime_ = end_time;
				std::get<Vec<BaseComplexType> >(this->current_space_) = start_point;
				if (this->reinitialize_stepsize_)
					this->SetStepSize(min(BaseRealType(Get<Stepping>().initial_step_size),abs(start_time-end_time)/Get<Stepping>().min_num_steps));

				ResetCounters();

				return SuccessCode::Success;
			}


		private:

		}; // re: DoublePrecisionTracker


		class MultiplePrecisionTracker : public FixedPrecisionTracker<MultiplePrecisionTracker>
		{
		public:
			using BaseComplexType = mpfr;
			using BaseRealType = mpfr_float;

			using EmitterType = FixedPrecisionTracker<MultiplePrecisionTracker>;

			BERTINI_DEFAULT_VISITABLE()


			/**
			\brief Construct a tracker, associating to it a System.

			The precision of the tracker will be whatever the current default is.  The tracker cannot change its precision, and will require the default precision to be this precision whenever tracking is started.  That is, the precision is fixed.
			*/
			MultiplePrecisionTracker(class System const& sys) : FixedPrecisionTracker<MultiplePrecisionTracker>(sys), precision_(DefaultPrecision())
			{	}

			
			MultiplePrecisionTracker() = delete;

			virtual ~MultiplePrecisionTracker() = default;


			unsigned CurrentPrecision() const override
			{
				return precision_;
			}


			/**
			\brief Set up the internals of the tracker for a fresh start.  

			Copies the start time, current stepsize, and start point.  Adjusts the current precision to match the precision of the start point.  Zeros counters.

			\param start_time The time at which to start tracking.
			\param end_time The time to which to track.
			\param start_point The space values from which to start tracking.
			*/
			SuccessCode TrackerLoopInitialization(BaseComplexType const& start_time,
			                               BaseComplexType const& end_time,
										   Vec<BaseComplexType> const& start_point) const override
			{

				if (start_point(0).precision()!=DefaultPrecision())
				{
					std::stringstream err_msg;
					err_msg << "start point for fixed multiple precision tracker has differing precision from default (" << start_point(0).precision() << "!=" << DefaultPrecision() << "), tracking cannot start";
					throw std::runtime_error(err_msg.str());
				}

				if (start_point(0).precision()!=CurrentPrecision())
				{
					std::stringstream err_msg;
					err_msg << "start point for fixed multiple precision tracker has differing precision from tracker's precision (" << start_point(0).precision() << "!=" << CurrentPrecision() << "), tracking cannot start";
					throw std::runtime_error(err_msg.str());
				}

				if (DefaultPrecision()!=CurrentPrecision())
				{
					std::stringstream err_msg;
					err_msg << "current default precision differs from tracker's precision (" << DefaultPrecision() << "!=" << CurrentPrecision() << "), tracking cannot start";
					throw std::runtime_error(err_msg.str());
				}


				this->NotifyObservers(Initializing<EmitterType,BaseComplexType>(*this,start_time, end_time, start_point));

				// set up the master current time and the current step size
				this->current_time_ = start_time;
				this->endtime_ = end_time;
				std::get<Vec<BaseComplexType> >(this->current_space_) = start_point;
				if (this->reinitialize_stepsize_)
					this->SetStepSize(min(mpfr_float(Get<Stepping>().initial_step_size),mpfr_float(abs(start_time-end_time)/Get<Stepping>().min_num_steps)));

				ResetCounters();

				return SuccessCode::Success;
			}

			bool PrecisionSanityCheck() const
			{	
				return GetSystem().precision() == precision_ &&
						DefaultPrecision()==precision_ && 
						std::get<Vec<mpfr> >(current_space_)(0).precision() == precision_ &&
						std::get<Vec<mpfr> >(tentative_space_)(0).precision() == precision_ &&
						std::get<Vec<mpfr> >(temporary_space_)(0).precision() == precision_ &&
						Precision(this->endtime_)==precision_
						        ;				
			}
		private:

			unsigned precision_;
		}; // re: MultiplePrecisionTracker
	} // namespace tracking
} // namespace bertini


#endif



