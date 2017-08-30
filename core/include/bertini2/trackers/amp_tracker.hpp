//This file is part of Bertini 2.
//
//amp_tracker.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_tracker.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_tracker.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire


/**
\file amp_tracker.hpp

\brief 

\brief Contains the AMPTracker type, and necessary functions.
*/

#ifndef BERTINI_AMP_TRACKER_HPP
#define BERTINI_AMP_TRACKER_HPP

#pragma once 

#include "bertini2/trackers/base_tracker.hpp"


namespace bertini{

	namespace tracking{

		
		using std::max;
		using std::min;
		using std::pow;

		using bertini::max;

		/**
		The minimum time that can be represented effectively using a given precision
		*/
		inline
		mpfr_float MinTimeForCurrentPrecision(unsigned precision, mpfr_float const& time_to_go, int safety_digits = 3)
		{
			mpfr_float t = pow( mpfr_float(10), safety_digits-long(precision)) * time_to_go;
			if (precision==DoublePrecision() && t<1e-150)
				return mpfr_float("1e-150");
			else
				return t;
		}

		/**
		\brief Just another name for MinTimeForCurrentPrecision
		*/
		inline
		mpfr_float MinStepSizeForPrecision(unsigned precision, mpfr_float const& time_to_go, int safety_digits = 3)
		{
			return MinTimeForCurrentPrecision(precision, time_to_go, safety_digits);
		}


		/**
		 \brief Compute the cost function for arithmetic versus precision.

		 From \cite AMP2, \f$C(P)\f$.  As currently implemented, this is 
		 \f$ 10.35 + 0.13 P \f$, where P is the precision.
		
		 This function tells you the relative cost of arithmetic at a given precision. 

		 \todo Recompute this cost function for boost::multiprecision::mpfr_float
		*/
		inline
		double ArithmeticCost(unsigned precision)
		{
			if (precision==DoublePrecision())
				return 1;
			else
				return 10.35 + 0.13 * precision;
		}




		/**
		 \Compute a stepsize satisfying AMP Criterion B with a given precision
		*/
		inline
		mpfr_float StepsizeSatisfyingCriterionB(unsigned precision,
										unsigned digits_B,
										unsigned num_newton_iterations,
										unsigned predictor_order = 0)
		{
			return pow(mpfr_float(10), -( digits_B - precision )*num_newton_iterations/(predictor_order+1));
		}
		

		/**
		\brief The minimum number of digits based on a log10 of a stepsize.

		\param log_of_stepsize log10 of a stepsize
		\param time_to_go How much time you have left to track.
		*/
		inline
		unsigned MinDigitsForLogOfStepsize(mpfr_float const& log_of_stepsize, mpfr_float const& time_to_go, unsigned safety_digits = 3)
		{
			return (ceil(log_of_stepsize) + ceil(log10(abs(time_to_go))) + safety_digits).convert_to<unsigned>();
		}

		/**
		\todo this function assumes you are going to 0.
		*/
		inline
		unsigned MinDigitsForStepsizeInterval(mpfr_float const& min_stepsize, mpfr_float const& max_stepsize, mpfr_float const& time_to_go)
		{
			return max(MinDigitsForLogOfStepsize(-log10(min_stepsize),time_to_go),   
				       MinDigitsForLogOfStepsize(-log10(max_stepsize),time_to_go));
		}
		

		/**
		 \brief Compute precision and stepsize minimizing the ArithmeticCost() of tracking.
		 
		 For a given range of precisions, an old stepsize, and a maximum stepsize, the Cost of tracking is computed, and a minimizer found.  
		
		 This function is used in the AMPTracker tracking loop, both in case of successful steps and in Criterion errors.


		 \param[out] new_precision The minimizing precision.
		 \param[out] new_stepsize The minimizing stepsize.
		 \param[in] min_precision The minimum considered precision.
		 \param[in] min_stepsize The minimum permitted stepsize.
		 \param[in] max_precision The maximum considered precision.
		 \param[in] max_stepsize The maximum permitted stepsize.  
		 \param[in] criterion_B_rhs The right hand side of CriterionB from \cite AMP1, \cite AMP2
		 \param num_newton_iterations The number of allowed Newton corrector iterations.
		 \param AMP_config The configuration of AMP settings for tracking.
		 \param predictor_order The order of the predictor being used.  This is the order itself, not the order of the error estimate.

		 \see ArithmeticCost
		*/
		template<typename RealT>
		void MinimizeTrackingCost(unsigned & new_precision, RealT & new_stepsize, 
						  unsigned min_precision, RealT const& min_stepsize,
						  unsigned max_precision, RealT const& max_stepsize,
						  unsigned digits_B,
						  unsigned num_newton_iterations,
						  unsigned predictor_order = 0)
		{
			double min_cost = Eigen::NumTraits<double>::highest();
			new_precision = MaxPrecisionAllowed()+1; // initialize to an impossible value.
			new_stepsize = min_stepsize; // initialize to minimum permitted step size.

			auto minimizer_routine = 
				[&min_cost, &new_stepsize, &new_precision, &digits_B, num_newton_iterations, predictor_order, max_stepsize](unsigned p)
				{
					RealT candidate_stepsize = min(StepsizeSatisfyingCriterionB(p, digits_B, num_newton_iterations, predictor_order),
					                              max_stepsize);
					using std::abs;
					double current_cost = ArithmeticCost(p) / abs(double(candidate_stepsize));

					if (current_cost < min_cost)
					{
						min_cost = current_cost;
						new_stepsize = candidate_stepsize;
						new_precision = p;
					}
				};

			unsigned lowest_mp_precision_to_test = min_precision;

			if (min_precision<=DoublePrecision())
				minimizer_routine(DoublePrecision());			


			if (lowest_mp_precision_to_test < LowestMultiplePrecision())
				lowest_mp_precision_to_test = LowestMultiplePrecision();
			else
				lowest_mp_precision_to_test = (lowest_mp_precision_to_test/PrecisionIncrement()) * PrecisionIncrement(); // use integer arithmetic to round.


			if (max_precision < LowestMultiplePrecision())
				max_precision = LowestMultiplePrecision();
			else
				max_precision = (max_precision/PrecisionIncrement()) * PrecisionIncrement(); // use integer arithmetic to round.

			for (unsigned p = lowest_mp_precision_to_test; p <= max_precision; p+=PrecisionIncrement())
				minimizer_routine(p);
		}








		/** 
		\class AMPTracker

		\brief Functor-like class for tracking paths on a system
		
		
		## Explanation
		
		The bertini::AMPTracker class enables tracking using Adaptive Multiple Precision on an arbitrary square homotopy.  

		The intended usage is to:

		1. Create a system, and instantiate some settings.
		2. Create an AMPTracker, associating it to the system you are going to solve or track on.
		3. Run AMPTracker::Setup and AMPTracker::PrecisionSetup, getting the settings in line for tracking.
		4. Repeatedly, or as needed, call the AMPTracker::TrackPath function, feeding it a start point, and start and end times.  The initial precision is that of the start point.  

		Working precision and stepsize are adjusted automatically to get around nearby singularities to the path.  If the endpoint is singular, this may very well fail, as prediction and correction get more and more difficult with proximity to singularities.  

		The TrackPath method is intended to allow the user to track to nonsingular endpoints, or to an endgame boundary, from which an appropriate endgame will be called.
		
		## Some notes

		The AMPTracker has internal state.  That is, it stores the current state of the path being tracked as data members.  After tracking has concluded, these statistics may be extracted.  If you need additional accessors for these data, contact the software authors.

		The class additionally uses the Boost.Log library for logging.  At time of this writing, the trivial logger is being used.  As development continues we will move toward using a more sophisticated logging model.  Suggestions are welcome.
		
		This class, like the other Tracker classes, uses mutable members to store the current state of the tracker.  The tracking tolerances, settings, etc, remain constant throughout a track, but the internal state such as the current time or space values, will change.  You can also expect the precision of the tracker to differ after a certain calls, too.
		
		## Example Usage
		
		Below we demonstrate a basic usage of the AMPTracker class to track a single path.  

		The pattern is as described above: create an instance of the class, feeding it the system to be tracked, and some configuration.  Then, use the tracker to track paths of the system.

		\code{.cpp}
		DefaultPrecision(30); // set initial precision.  This is not strictly necessary.

		using namespace bertini::tracking;

		// 1. Create the system
		Var x = MakeVariable("x");
		Var y = MakeVariable("y");
		Var t = MakeVariable("t");

		System sys;

		VariableGroup v{x,y};

		sys.AddFunction(pow(x,2) + (1-t)*x - 1);
		sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
		sys.AddPathVariable(t);
		sys.AddVariableGroup(v);

		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		//  2. Create the Tracker object, associating the system to it.
		bertini::tracking::AMPTracker tracker(sys);

		SteppingConfig stepping_preferences;
		NewtonConfig newton_preferences;

		// 3. Get the settings into the tracker 
		tracker.Setup(Predictor::Euler,
		              	mpfr_float("1e-5"),
						mpfr_float("1e5"),
						stepping_preferences,
						newton_preferences);

		tracker.PrecisionSetup(AMP);
	
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
		
		## Testing

		* Test suite driving this class: AMP_tracker_basics.
		* File: test/tracking_basics/tracker_test.cpp
		* Functionality tested: Can use an AMPTracker to track in various situations, including tracking on nonsingular paths.  Also track to a singularity on the square root function.  Furthermore test that tracking fails to start from a singular start point, and tracking fails if a singularity is directly on the path being tracked.

		*/
		class AMPTracker : public Tracker<AMPTracker>
		{
			friend class Tracker<AMPTracker>;
		public:
			BERTINI_DEFAULT_VISITABLE()
			
			typedef Tracker<AMPTracker> Base;
			typedef typename TrackerTraits<AMPTracker>::EventEmitterType EmitterType;

			enum UpsampleRefinementOption
			{
			   upsample_refine_off  = 0,
			   upsample_refine_on   = 1
			};
		

			/**
			\brief Construct an Adaptive Precision tracker, associating to it a System.
			*/
			AMPTracker(class System const& sys) : Tracker(sys), current_precision_(DefaultPrecision())
			{	
				Set<PrecConf>(AMPConfigFrom(sys));
			}
			

			/**
			\brief Special additional setup call for the AMPTracker, selecting the config for adaptive precision.
			*/
			void PrecisionSetup(AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				Set<PrecConf>(AMP_config);
			}


			const unsigned GetCurrentPrecision() const
			{
				return current_precision_;
			}

			
			/**
			\brief Switch preservation of precision after tracking on / off

			By default, precision is preserved after tracking, so the precision of the ambient workspace is returned to its previous state once Adaptive Precision tracking is done.
			*/
			void PrecisionPreservation(bool should_preseve_precision)
			{
				preserve_precision_ = should_preseve_precision;
			}

			
			virtual ~AMPTracker() = default;


			Vec<mpfr> CurrentPoint() const override
			{
				if (this->CurrentPrecision()==DoublePrecision())
				{
					const auto& curr_vector = std::get<Vec<dbl>>(this->current_space_);
					Vec<mpfr> returnme(NumVariables());
					for (unsigned ii = 0; ii < NumVariables(); ++ii)
					{
						returnme(ii) = mpfr(curr_vector(ii));
					}
					return returnme;
				}
				else
					return std::get<Vec<mpfr>>(this->current_space_);
			}

			

		private:

			/**
			\brief Set up the internals of the tracker for a fresh start.  

			Copies the start time, current stepsize, and start point.  Adjusts the current precision to match the precision of the start point.  Zeros counters.

			\param start_time The time at which to start tracking.
			\param end_time The time to which to track.
			\param start_point The space values from which to start tracking.
			*/
			SuccessCode TrackerLoopInitialization(mpfr const& start_time,
			                               mpfr const& end_time,
										   Vec<mpfr> const& start_point) const override
			{
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(
				        (!preserve_precision_ 
				         || 
				         -log10(tracking_tolerance_) <= start_point(0).precision())
				         && "when tracking a path, either preservation of precision must be turned off (so precision can be higher at the end of tracking), or the initial precision must be high enough to support the resulting points to the desired tolerance"
				         );
				#endif

				NotifyObservers(Initializing<AMPTracker,mpfr>(*this,start_time, end_time, start_point));

				initial_precision_ = Precision(start_point(0));
				DefaultPrecision(initial_precision_);
				// set up the master current time and the current step size
				current_time_.precision(initial_precision_);
				current_time_ = start_time;

				endtime_highest_precision_.precision(initial_precision_);
				endtime_highest_precision_ = end_time;

				endtime_.precision(initial_precision_);
				endtime_ = end_time;

				current_stepsize_.precision(initial_precision_);
				if (reinitialize_stepsize_)
				{
					mpfr_float segment_length = abs(start_time-end_time)/Get<Stepping>().min_num_steps;
					SetStepSize(min(mpfr_float(Get<Stepping>().initial_step_size),segment_length));
				}

				// populate the current space value with the start point, in appropriate precision
				if (initial_precision_==DoublePrecision())
					MultipleToDouble( start_point);
				else
					MultipleToMultiple(initial_precision_, start_point);
				
				ChangePrecision<upsample_refine_off>(initial_precision_);
				
				ResetCounters();

				ChangePrecision<upsample_refine_off>(start_point(0).precision());

				return InitialRefinement();
			}





			void ResetCounters() const override
			{
				Tracker::ResetCountersBase();
				num_precision_decreases_ = 0;
				num_successful_steps_since_stepsize_increase_ = 0;
				num_successful_steps_since_precision_decrease_ = 0;
				// initialize to the frequency so guaranteed to compute it the first try 	
				num_steps_since_last_condition_number_computation_ = this->Get<Stepping>().frequency_of_CN_estimation;
			}

			/** 
			\brief Run an initial refinement of the start point, to ensure in high enough precision to start.

			\return Whether initial refinement was successful.  
			*/
			SuccessCode InitialRefinement() const
			{
				SuccessCode initial_refinement_code = RefineStoredPoint();
				if (initial_refinement_code!=SuccessCode::Success)
				{
					do {
						if (current_precision_ > Get<PrecConf>().maximum_precision)
						{
							NotifyObservers(SingularStartPoint<EmitterType>(*this));
							return SuccessCode::SingularStartPoint;
						}

						if (current_precision_==DoublePrecision())
							initial_refinement_code = ChangePrecision<upsample_refine_on>(LowestMultiplePrecision());
						else
							initial_refinement_code = ChangePrecision<upsample_refine_on>(current_precision_+PrecisionIncrement());
					}
					while (initial_refinement_code!=SuccessCode::Success);
				}
				return SuccessCode::Success;
			}		




			/**
			\brief Ensure that number of steps, stepsize, and precision still ok.

			\return Success if ok to keep going, and a different code otherwise. 
			*/
			SuccessCode PreIterationCheck() const override
			{
				if (num_successful_steps_taken_ >= Get<Stepping>().max_num_steps)
					return SuccessCode::MaxNumStepsTaken;
				if (current_stepsize_ < Get<Stepping>().min_step_size)
					return SuccessCode::MinStepSizeReached;
				if (current_precision_ > Get<PrecConf>().maximum_precision)
					return SuccessCode::MaxPrecisionReached;

				return SuccessCode::Success;
			}



			
			


			void PostTrackCleanup() const override
			{
				if (preserve_precision_)
					ChangePrecision(initial_precision_);
				NotifyObservers(TrackingEnded<EmitterType>(*this));
			}

			/**
			\brief Copy from the internally stored current solution into a final solution.
			
			If preservation of precision is on, this function first returns to the initial precision.

			\param[out] solution_at_endtime The solution at the end time
			*/
			void CopyFinalSolution(Vec<mpfr> & solution_at_endtime) const override
			{

				// the current precision is the precision of the output solution point.
				if (current_precision_==DoublePrecision())
				{
					unsigned num_vars = GetSystem().NumVariables();
					solution_at_endtime.resize(num_vars);
					for (unsigned ii=0; ii<num_vars; ii++)
						solution_at_endtime(ii) = mpfr(std::get<Vec<dbl> >(current_space_)(ii));
				}
				else
				{
					unsigned num_vars = GetSystem().NumVariables();
					solution_at_endtime.resize(num_vars);
					for (unsigned ii=0; ii<num_vars; ii++)
					{
						solution_at_endtime(ii).precision(current_precision_);
						solution_at_endtime(ii) = std::get<Vec<mpfr> >(current_space_)(ii);
					}
				}
			}




			/**
			\brief Run an iteration of the tracker loop.

			Predict and correct, adjusting precision and stepsize as necessary.

			\return Success if the step was successful, and a non-success code if something went wrong, such as a linear algebra failure or AMP Criterion violation.
			*/
			SuccessCode TrackerIteration() const override
			{
				if (current_precision_==DoublePrecision())
					return TrackerIteration<dbl>();
				else
					return TrackerIteration<mpfr>();
			}


			/**
			\brief Run an iteration of AMP tracking.

			\return SuccessCode indicating whether the iteration was successful.
			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template <typename ComplexType>
			SuccessCode TrackerIteration() const // not an override, because it is templated
			{	

				using RealType = typename Eigen::NumTraits<ComplexType>::Real;

				#ifndef BERTINI_DISABLE_ASSERTS
				assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
				#endif

				NotifyObservers(NewStep<EmitterType>(*this));

				Vec<ComplexType>& predicted_space = std::get<Vec<ComplexType> >(temporary_space_); // this will be populated in the Predict step
				Vec<ComplexType>& current_space = std::get<Vec<ComplexType> >(current_space_); // the thing we ultimately wish to update
				ComplexType current_time = ComplexType(current_time_);
				ComplexType delta_t = ComplexType(delta_t_);

				SuccessCode predictor_code = Predict<ComplexType, RealType>(predicted_space, current_space, current_time, delta_t);
				if (predictor_code==SuccessCode::MatrixSolveFailureFirstPartOfPrediction)
				{
					NotifyObservers(FirstStepPredictorMatrixSolveFailure<EmitterType>(*this));
					InitialMatrixSolveError();
					return predictor_code;
				}
				else if (predictor_code==SuccessCode::MatrixSolveFailure)
				{
					NotifyObservers(PredictorMatrixSolveFailure<EmitterType>(*this));
					NewtonConvergenceError();// decrease stepsize, and adjust precision as necessary
					return predictor_code;
				}	
				else if (predictor_code==SuccessCode::HigherPrecisionNecessary)
				{	
					NotifyObservers(PredictorHigherPrecisionNecessary<EmitterType>(*this));
					AMPCriterionError<ComplexType>();
					return predictor_code;
				}


				NotifyObservers(SuccessfulPredict<AMPTracker, ComplexType>(*this, predicted_space));

				Vec<ComplexType>& tentative_next_space = std::get<Vec<ComplexType> >(tentative_space_); // this will be populated in the Correct step

				ComplexType tentative_next_time = current_time + delta_t;

				SuccessCode corrector_code = Correct<ComplexType, RealType>(tentative_next_space,
													 predicted_space,
													 tentative_next_time);

				if (corrector_code==SuccessCode::MatrixSolveFailure || corrector_code==SuccessCode::FailedToConverge)
				{
					NotifyObservers(CorrectorMatrixSolveFailure<EmitterType>(*this));
					NewtonConvergenceError();
					return corrector_code;
				}
				else if (corrector_code == SuccessCode::HigherPrecisionNecessary)
				{
					NotifyObservers(CorrectorHigherPrecisionNecessary<EmitterType>(*this));
					AMPCriterionError<ComplexType>();
					return corrector_code;
				}
				else if (corrector_code == SuccessCode::GoingToInfinity)
				{
					// there is no corrective action possible...
					return corrector_code;
				}

				NotifyObservers(SuccessfulCorrect<AMPTracker, ComplexType>(*this, tentative_next_space));

				// copy the tentative vector into the current space vector;
				current_space = tentative_next_space;
				return AdjustAMPStepSuccess<ComplexType>();
			}


			/**
			Check whether the path is going to infinity.
			*/
			SuccessCode CheckGoingToInfinity() const override
			{
				if (current_precision_ == DoublePrecision())
					return Base::CheckGoingToInfinity<dbl>();
				else
					return Base::CheckGoingToInfinity<mpfr>();
			}

			/**
			\brief Commit the next precision and stepsize, and adjust internals.

			\tparam refine_if_necessary Flag indicating whether to refine if precision if increasing.

			\see UpsampleRefinementOption
			*/
			template <UpsampleRefinementOption refine_if_necessary = upsample_refine_off>
			SuccessCode UpdatePrecisionAndStepsize() const
			{
				SetStepSize(next_stepsize_);
				return ChangePrecision<refine_if_necessary>(next_precision_);
			}







			/**
			\brief Increase stepsize or decrease precision, because of consecutive successful steps.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			If the most recent step was successful, maybe adjust down precision and up stepsize.  

			The number of consecutive successful steps is recorded as state in this class, and if this number exceeds a user-determine threshold, the precision or stepsize are allowed to favorably change.  If not, then precision can only go up or remain the same.  Stepsize can only decrease.  These changes depend on the AMP criteria and current tracking tolerance.
			*/
			template <typename ComplexType>
			SuccessCode AdjustAMPStepSuccess() const
			{
				// TODO: think about why we consider reducing the stepsize?  this is despite documentation stating that it can only increase
				mpfr_float min_stepsize = current_stepsize_ * Get<Stepping>().step_size_fail_factor;
				mpfr_float max_stepsize = min( current_stepsize_ * Get<Stepping>().step_size_success_factor,  mpfr_float(Get<Stepping>().max_step_size));


				unsigned min_precision = MinRequiredPrecision_BCTol<ComplexType>();
				unsigned max_precision = max(min_precision,current_precision_);

				if (num_successful_steps_since_stepsize_increase_ < Get<Stepping>().consecutive_successful_steps_before_stepsize_increase)
					max_stepsize = current_stepsize_; // disallow stepsize changing 


				if ( (num_successful_steps_since_precision_decrease_ < Get<PrecConf>().consecutive_successful_steps_before_precision_decrease)
				    ||
				    (num_precision_decreases_ >= Get<PrecConf>().max_num_precision_decreases))
					min_precision = max(min_precision, current_precision_); // disallow precision changing 


				MinimizeTrackingCost(next_precision_, next_stepsize_, 
							min_precision, min_stepsize,
							max_precision, max_stepsize,
							DigitsB<ComplexType>(),
							Get<NewtonConfig>().max_num_newton_iterations,
							predictor_order_);


				if ( (next_stepsize_ > current_stepsize_) || (next_precision_ < current_precision_) )
					num_successful_steps_since_stepsize_increase_ = 0;
				else
					++num_successful_steps_since_stepsize_increase_;
				

				if (next_precision_ < current_precision_)
				{ 
					++num_precision_decreases_;
					num_successful_steps_since_precision_decrease_ = 0;
				}
				else
					++num_successful_steps_since_precision_decrease_;

				return UpdatePrecisionAndStepsize();
			}



			void InitialMatrixSolveError() const
			{
				next_stepsize_ = current_stepsize_;

				if (current_precision_==DoublePrecision())
					next_precision_ = LowestMultiplePrecision();
				else
					next_precision_ = current_precision_+(1+num_consecutive_failed_steps_) * PrecisionIncrement();

				UpdatePrecisionAndStepsize();
			}

			/**
			\brief The convergence_error function from \cite AMP2.  
	
			Adjust precision and stepsize due to Newton's method failure to converge in the maximum number of steps.

			Decrease stepsize, then increase precision until the new stepsize is greater than a precision-dependent minimum stepsize.  
			
			This function is used in the AMPTracker loop, when Newton's method fails to converge.

			\see MinStepSizeForPrecision
			*/
			void NewtonConvergenceError() const
			{
				next_precision_ = current_precision_;
				next_stepsize_ = Get<Stepping>().step_size_fail_factor*current_stepsize_;

				while (next_stepsize_ < MinStepSizeForPrecision(next_precision_, abs(current_time_ - endtime_)))
				{
					if (next_precision_==DoublePrecision())
						next_precision_=LowestMultiplePrecision();
					else
						next_precision_+=PrecisionIncrement();
				}

				UpdatePrecisionAndStepsize();
			}







			/**
			\brief Adjust step size and precision due to AMP Criterion violation
			
			This function adjusts internals to the tracker object.
			
			Precision is REQUIRED to increase at least one increment.  Stepsize is REQUIRED to decrease.

			\tparam ComplexType The complex number type.
			*/
			template<typename ComplexType>
			void AMPCriterionError() const
			{	
				using RealT = typename Eigen::NumTraits<ComplexType>::Real;

				unsigned min_next_precision; // sure, i could use a trigraph here, but it'd be terrible
				if (current_precision_==DoublePrecision())
					min_next_precision = LowestMultiplePrecision(); // precision increases
				else
					min_next_precision = current_precision_ + (1+num_consecutive_failed_steps_)*PrecisionIncrement(); // precision increases


				mpfr_float min_stepsize = MinStepSizeForPrecision(current_precision_, abs(current_time_ - endtime_));
				mpfr_float max_stepsize = current_stepsize_ * Get<Stepping>().step_size_fail_factor;  // Stepsize decreases.

				if (min_stepsize > max_stepsize)
				{
					// stepsizes are incompatible, must increase precision
					next_precision_ = min_next_precision;
					// decrease stepsize somewhat less than the fail factor
					next_stepsize_ = max(current_stepsize_ * (1+Get<Stepping>().step_size_fail_factor)/2, min_stepsize);
				}
				else
				{
					unsigned digits_B = DigitsB<ComplexType>();

					unsigned min_precision = max(min_next_precision,
					                             digits_B,
					                             DigitsC<ComplexType>(),
					                             MinDigitsForStepsizeInterval(min_stepsize, max_stepsize, abs(current_time_ - endtime_)),
					                             digits_final_
					                             );
					
						unsigned a = ceil(digits_B - (predictor_order_+1)* -log10(max_stepsize)/Get<NewtonConfig>().max_num_newton_iterations).convert_to<unsigned>();

					unsigned max_precision = max(min_precision, a);

					MinimizeTrackingCost(next_precision_, next_stepsize_, 
							min_precision, min_stepsize,
							Get<PrecConf>().maximum_precision, max_stepsize,
							digits_B,
							Get<NewtonConfig>().max_num_newton_iterations,
							predictor_order_);
				}

				UpdatePrecisionAndStepsize();
			}





			/**
			\brief Get the raw right-hand side of Criterion B based on current state.
			*/
			template<typename ComplexType>
			NumErrorT B_RHS() const
			{	
				return max(amp::CriterionBRHS(this->norm_J_, 
				           					  this->norm_J_inverse_, 
				           					  Get<NewtonConfig>().max_num_newton_iterations, 
				           					  tracking_tolerance_, 
				           					  this->size_proportion_, 
				           					  Get<PrecConf>()), NumErrorT(0));
			}




			/**
			\brief Get the right hand side of Criterion B based on current state.
			
			Returns the larger of the right hand side of B, or 1.

			Uses:

			* the most recent estimate on the norm of \f$J\f$,
			* the most recent estimate of the norm of \f$J^{-1}\f$,
			* the number of allowed newton iterations,
			* the tracking tolerance,
			* the size_proportion of the latest prediction,
			* the AMP configuration.

			*/
			template<typename ComplexType>
			unsigned DigitsB() const
			{	
				return unsigned(B_RHS<ComplexType>());
			}



			

			/**
			\brief Get the raw right-hand side of Criterion B based on current state.
			*/
			template<typename ComplexType>
			NumErrorT C_RHS() const
			{	
				return max(amp::CriterionCRHS(this->norm_J_inverse_, 
				                              NumErrorT(std::get<Vec<ComplexType> > (current_space_).norm()), 
				                              tracking_tolerance_, 
				                              Get<PrecConf>()), NumErrorT(0));
			}



			/**
			\brief Get the right hand side of Criterion C based on current state.
			
			Returns the larger of the right hand side of C, or 1.

			Uses:

			* the most recent estimate of the norm of \f$J^{-1}\f$,
			* the most recent current space point,
			* the current tracking tolerance,
			* the AMP configuration.

			*/
			template<typename ComplexType>
			unsigned DigitsC() const
			{	
				return unsigned(C_RHS<ComplexType>());
			}





			/**
			\brief Get the minimum required precision based on current state.

			The current state determining the minimum required precision is:

			* the norm of the Jacobian,
			* the norm of the inverse of the Jacobian,
			* the size_proportion, related to stepsize and predictor order,
			* the norm of the current solution. 
			* the tracking tolerance.

			The min digits needed is the maximum of 

			* Criterion B, 
			* Criterion C, and 
			* the digits required by the tracking tolerance.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template <typename ComplexType>
			unsigned MinRequiredPrecision_BCTol() const
			{
				return max(DigitsB<ComplexType>(), 
				           DigitsC<ComplexType>(), 
				           digits_tracking_tolerance_, 
				           DoublePrecision()); 
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
				Tracker::IncrementBaseCountersSuccess();
				NotifyObservers(SuccessfulStep<EmitterType>(*this));
			}

			/**
			\brief Increment and reset counters after a failed TrackerIteration()
			*/
			void OnStepFail() const override
			{
				Tracker::IncrementBaseCountersFail();
				num_successful_steps_since_precision_decrease_ = 0;
				num_successful_steps_since_stepsize_increase_ = 0;
				NotifyObservers(FailedStep<EmitterType>(*this));
			}



			void OnInfiniteTruncation() const override
			{
				NotifyObservers(InfinitePathTruncation<EmitterType>(*this));
			}


			////////////////
			//
			//       Predict and Correct functions
			//
			/////////////////


			


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
			template<typename ComplexType, typename RealType, typename Derived>
			SuccessCode Predict(Vec<ComplexType> & predicted_space, 
								const Eigen::MatrixBase<Derived>& current_space,
								ComplexType const& current_time, ComplexType const& delta_t) const
			{
				
								
				
				
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");
				static_assert(std::is_same<typename Derived::Scalar, ComplexType>::value, "scalar types must match");


				if (predictor_->HasErrorEstimate())
					return predictor_->Predict(predicted_space,
									this->error_estimate_,
									this->size_proportion_,
									this->norm_J_,
									this->norm_J_inverse_,
									tracked_system_,
									current_space, current_time, 
									delta_t,
									this->condition_number_estimate_,
									num_steps_since_last_condition_number_computation_, 
									Get<Stepping>().frequency_of_CN_estimation, 
									tracking_tolerance_,
									Get<PrecConf>());
				else
					return predictor_->Predict(predicted_space,
									this->size_proportion_,
									this->norm_J_,
									this->norm_J_inverse_,
									tracked_system_,
									current_space, current_time, 
									delta_t,
									this->condition_number_estimate_,
									num_steps_since_last_condition_number_computation_, 
									Get<Stepping>().frequency_of_CN_estimation, 
									tracking_tolerance_,
									Get<PrecConf>());
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
								ComplexType const& current_time) const
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");




				return corrector_->Correct(corrected_space,
									this->norm_delta_z_,
									this->norm_J_,
									this->norm_J_inverse_,
									this->condition_number_estimate_,
									tracked_system_,
									current_space,
									current_time,
									tracking_tolerance_,
									Get<NewtonConfig>().min_num_newton_iterations,
									Get<NewtonConfig>().max_num_newton_iterations,
									Get<PrecConf>());
			}



			/**
			\brief Run Newton's method from the currently stored time and space value, in current precision.

			This overwrites the value of the current space, if successful, with the refined value.

			\return Whether the refinement was successful.
			*/
			SuccessCode RefineStoredPoint() const
			{
				SuccessCode code;
				if (current_precision_==DoublePrecision())
				{
					code = RefineImpl<dbl>(std::get<Vec<dbl> >(temporary_space_),std::get<Vec<dbl> >(current_space_), dbl(current_time_));
					if (code == SuccessCode::Success)
						std::get<Vec<dbl> >(current_space_) = std::get<Vec<dbl> >(temporary_space_);
				}
				else
				{
					code = RefineImpl<mpfr>(std::get<Vec<mpfr> >(temporary_space_),std::get<Vec<mpfr> >(current_space_), current_time_);
					if (code == SuccessCode::Success)
						std::get<Vec<mpfr> >(current_space_) = std::get<Vec<mpfr> >(temporary_space_);
				}
				return code;
			}



			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.  The tolerance is the tracking tolerance specified during Setup(...).

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			\param[out] new_space The result of running the refinement.
			\param start_point The base point for running Newton's method.
			\param current_time The current time value.

			\return Code indicating whether was successful or not.  Regardless, the value of new_space is overwritten with the correction result.
			*/
			template <typename ComplexType>
			SuccessCode RefineImpl(Vec<ComplexType> & new_space,
								Vec<ComplexType> const& start_point, ComplexType const& current_time) const
			{
				using RealType = typename Eigen::NumTraits<ComplexType>::Real;
				
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				auto target_precision = Precision(current_time);
				assert(Precision(start_point)==target_precision);
				ChangePrecision(target_precision);
				Precision(new_space,target_precision);
				



				return corrector_->Correct(new_space,
										   this->norm_delta_z_,
										   this->norm_J_,
										   this->norm_J_inverse_,
										   this->condition_number_estimate_,
										   tracked_system_,
										   start_point,
										   current_time,
										   tracking_tolerance_,
										   Get<NewtonConfig>().min_num_newton_iterations,
										   Get<NewtonConfig>().max_num_newton_iterations,
										   Get<PrecConf>());
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
			\param max_iterations The maximum allowable number of iterations to perform.
			\return Code indicating whether was successful or not.  Regardless, the value of new_space is overwritten with the correction result.
			*/
			template <typename ComplexType>
			SuccessCode RefineImpl(Vec<ComplexType> & new_space,
								Vec<ComplexType> const& start_point, ComplexType const& current_time,
								NumErrorT const& tolerance, unsigned max_iterations) const
			{
				auto target_precision = Precision(current_time);
				assert(Precision(start_point)==target_precision);
				ChangePrecision(target_precision);
				Precision(new_space,target_precision);

				return corrector_->Correct(new_space,
										this->norm_delta_z_,
										this->norm_J_,
										this->norm_J_inverse_,
										this->condition_number_estimate_,
										tracked_system_,
										start_point,
										current_time, 
										tolerance,
										1,
										max_iterations,
										Get<PrecConf>());
			}














			/////////////////
			//
			//  Functions for converting between precision types
			//
			///////////////////////

		public:
			/**
			Change precision of tracker to next_precision.  Converts the internal temporaries, and adjusts precision of system. Then refines if necessary.

			If the new precision is higher than current precision, a refine step will be called, which runs Newton's method.  This may fail, leaving the tracker in a state with higher precision internals, but garbage digits after the previously known digits.

			\param new_precision The precision to change to.
			\return SuccessCode indicating whether the change was successful.  If the precision increases, and the refinement loop fails, this could be not Success.  Changing down is guaranteed to succeed.
			*/
			template <UpsampleRefinementOption refine_if_necessary = upsample_refine_off>
			SuccessCode ChangePrecision(unsigned new_precision) const
			{
				if (new_precision==current_precision_) // no op
					return SuccessCode::Success;

				NotifyObservers(PrecisionChanged<EmitterType>(*this,current_precision_,new_precision));
				

				bool upsampling_needed = new_precision > current_precision_;
				// reset the counter for estimating the condition number.  
				num_steps_since_last_condition_number_computation_ = this->Get<Stepping>().frequency_of_CN_estimation;

				if (new_precision==DoublePrecision() && current_precision_>DoublePrecision())
				{
					// convert from multiple precision to double precision
					MultipleToDouble();
				}
				else if(new_precision > DoublePrecision() && current_precision_ == DoublePrecision())
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


				if (refine_if_necessary && upsampling_needed)
					return RefineStoredPoint();
				else
					return SuccessCode::Success;
			}





			private:
			
			/**
			\brief Converts from double to double

			Copies a multiple-precision into the double storage vector, and changes precision of the time and delta_t.

			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void DoubleToDouble(Vec<dbl> const& source_point) const
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == GetSystem().NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				#endif

				current_precision_ = DoublePrecision();
				DefaultPrecision(DoublePrecision());

				GetSystem().precision(16);

				std::get<Vec<dbl> >(current_space_) = source_point;
			}

			/**
			\brief Converts from multiple to double

			Changes the precision of the internal temporaries to double precision
			*/
			void DoubleToDouble() const
			{
				DoubleToDouble(std::get<Vec<dbl> >(current_space_));
			}
			


			/**
			\brief Converts from multiple to double

			Copies a multiple-precision into the double storage vector, and changes precision of the time and delta_t.

			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void MultipleToDouble(Vec<mpfr> const& source_point) const
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == GetSystem().NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				#endif
				previous_precision_ = current_precision_;
				current_precision_ = DoublePrecision();
				DefaultPrecision(DoublePrecision());

				GetSystem().precision(DoublePrecision());

				if (std::get<Vec<dbl> >(current_space_).size()!=source_point.size())
					std::get<Vec<dbl> >(current_space_).resize(source_point.size());

				for (unsigned ii=0; ii<source_point.size(); ii++)
					std::get<Vec<dbl> >(current_space_)(ii) = dbl(source_point(ii));

				endtime_.precision(DoublePrecision());
			}

			/**
			\brief Converts from multiple to double

			Changes the precision of the internal temporaries to double precision
			*/
			void MultipleToDouble() const
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
			void DoubleToMultiple(unsigned new_precision, Vec<dbl> const& source_point) const
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == GetSystem().NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				assert(new_precision > DoublePrecision() && "must convert to precision higher than DoublePrecision when converting to multiple precision");
				#endif
				previous_precision_ = current_precision_;
				current_precision_ = new_precision;
				DefaultPrecision(new_precision);
				GetSystem().precision(new_precision);
				predictor_->ChangePrecision(new_precision);
				corrector_->ChangePrecision(new_precision);

				endtime_ = endtime_highest_precision_;
				endtime_.precision(new_precision);

				current_time_.precision(new_precision);

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
			void DoubleToMultiple(unsigned new_precision) const
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
			void MultipleToMultiple(unsigned new_precision, Vec<mpfr> const& source_point) const
			{	
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(source_point.size() == GetSystem().NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
				assert(new_precision > DoublePrecision() && "must convert to precision higher than DoublePrecision when converting to multiple precision");
				#endif
				previous_precision_ = current_precision_;
				current_precision_ = new_precision;
				DefaultPrecision(new_precision);
				GetSystem().precision(new_precision);
				predictor_->ChangePrecision(new_precision);
				corrector_->ChangePrecision(new_precision);

				endtime_ = endtime_highest_precision_;
				endtime_.precision(new_precision);

				current_time_.precision(new_precision);

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
			void MultipleToMultiple(unsigned new_precision) const
			{
				MultipleToMultiple( new_precision, std::get<Vec<mpfr> >(current_space_));
			}


			/**
			\brief Change precision of all temporary internal state variables.

			This excludes those which cannot be re-written without copying -- the current space point most notably.

			\brief new_precision The new precision to adjust to.
			*/
			void AdjustTemporariesPrecision(unsigned new_precision) const
			{
				unsigned num_vars = GetSystem().NumVariables();

				//  the current_space value is adjusted in the appropriate ChangePrecision function
				std::get<Vec<mpfr> >(tentative_space_).resize(num_vars);
				for (unsigned ii = 0; ii < num_vars; ++ii)
					std::get<Vec<mpfr> >(tentative_space_)(ii).precision(new_precision);

				std::get<Vec<mpfr> >(temporary_space_).resize(num_vars);
				for (unsigned ii = 0; ii < num_vars; ++ii)
					std::get<Vec<mpfr> >(temporary_space_)(ii).precision(new_precision);
			}



			
			/**
			\brief Ensure that all internal state is in uniform precision.

			\return True if all internal state variables are in the same precision as current_precision_, false otherwise.
			*/
			bool PrecisionSanityCheck() const
			{	
				if (current_precision_==DoublePrecision())
				{
					return true;
				}
				else
				{
					assert(DefaultPrecision()==current_precision_ && "current precision differs from the default precision");

					return GetSystem().precision() == current_precision_ &&
							predictor_->precision() == current_precision_ &&
							std::get<Vec<mpfr> >(current_space_)(0).precision() == current_precision_ &&
							std::get<Vec<mpfr> >(tentative_space_)(0).precision() == current_precision_ &&
							std::get<Vec<mpfr> >(temporary_space_)(0).precision() == current_precision_ &&
							Precision(endtime_) == current_precision_ && 
							Precision(current_time_) == current_precision_
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





			////////////
			// state variables
			/////////////
			bool preserve_precision_ = false; ///< Whether the tracker should change back to the initial precision after tracking paths.

			mutable unsigned previous_precision_; ///< The previous precision of the tracker.
			mutable unsigned current_precision_; ///< The current precision of the tracker, the system, and all temporaries.
			mutable unsigned next_precision_; ///< The next precision
			mutable unsigned num_precision_decreases_; ///< The number of times precision has decreased this track.
			mutable unsigned initial_precision_; ///< The precision at the start of tracking.
			mutable unsigned num_successful_steps_since_precision_decrease_; ///< The number of successful steps since decreased precision.

			mutable mpfr endtime_highest_precision_;

		public:

			unsigned CurrentPrecision() const override
			{
				return current_precision_;
			}
		}; // re: class Tracker

	} // namespace tracking
} // namespace bertini


#endif



