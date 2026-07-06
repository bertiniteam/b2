//This file is part of Bertini 2.
//
//trackers/observers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//trackers/observers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with trackers/observers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire


/**
\file include/bertini2/trackers/observers.hpp

\brief Contains the trackers/observers base types
*/

#pragma once

#include "bertini2/trackers/events.hpp"

#include "bertini2/detail/observer.hpp"

#include "bertini2/trackers/base_tracker.hpp"
#include "bertini2/logging.hpp"
#include <boost/type_index.hpp>

namespace bertini {

	namespace tracking{


		/// \brief Observer that records the starting precision and the first precision increase during a track.
		template<class TrackerT>
		class FirstPrecisionRecorder
			: public TypedObserver< FirstPrecisionRecorder<TrackerT>, TrackerT,
			        TrackingStarted<typename TrackerTraits<TrackerT>::EventEmitterType>,
			        PrecisionChanged<typename TrackerTraits<TrackerT>::EventEmitterType> >
		{ BOOST_TYPE_INDEX_REGISTER_CLASS

			using EmitterT = typename TrackerTraits<TrackerT>::EventEmitterType;

		public:

			/// \brief On tracking start, record the starting precision.
			ObserveResult OnEvent(TrackingStarted<EmitterT> const& e)
			{
				precision_increased_ = false;
				starting_precision_ = e.Get().CurrentPrecision();
				return ObserveResult::KeepObserving;
			}

			/// \brief On a precision increase, record it and unsubscribe (only the first increase is wanted).
			ObserveResult OnEvent(PrecisionChanged<EmitterT> const& e)
			{
				auto next = e.Next();
				if (next > e.Previous())
				{
					precision_increased_ = true;
					next_precision_ = next;
					time_of_first_increase_ = e.Get().CurrentTime();
					// done: ask to be dropped instead of mutating the list mid-dispatch.
					return ObserveResult::Unsubscribe;
				}
				return ObserveResult::KeepObserving;
			}

			/// \brief Get the precision tracking started at.
			unsigned  StartPrecision() const
			{
				return starting_precision_;
			}

			/// \brief Get the precision after the first increase.
			unsigned NextPrecision() const
			{
				return next_precision_;
			}

			/// \brief Query whether precision increased at all during the track.
			bool DidPrecisionIncrease() const
			{
				return precision_increased_;
			}

			/// \brief Get the time value at which precision first increased.
			typename TrackerTraits<TrackerT>::BaseComplexT TimeOfIncrease() const
			{
				return time_of_first_increase_;
			}

			virtual ~FirstPrecisionRecorder() = default;

		private:

			unsigned starting_precision_;
			unsigned next_precision_;
			bool precision_increased_;
			typename TrackerTraits<TrackerT>::BaseComplexT time_of_first_increase_;
		};


		/// \brief Observer that records the minimum and maximum precision reached during a track.
		template<class TrackerT>
		class MinMaxPrecisionRecorder
			: public TypedObserver< MinMaxPrecisionRecorder<TrackerT>, TrackerT,
			        TrackingStarted<typename TrackerTraits<TrackerT>::EventEmitterType>,
			        PrecisionChanged<typename TrackerTraits<TrackerT>::EventEmitterType> >
		{ BOOST_TYPE_INDEX_REGISTER_CLASS

			using EmitterT = typename TrackerTraits<TrackerT>::EventEmitterType;

		public:

			/// \brief On tracking start, seed the min and max with the starting precision.
			ObserveResult OnEvent(TrackingStarted<EmitterT> const& e)
			{
				min_precision_ = e.Get().CurrentPrecision();
				max_precision_ = e.Get().CurrentPrecision();
				return ObserveResult::KeepObserving;
			}

			/// \brief On a precision change, update the running min and max.
			ObserveResult OnEvent(PrecisionChanged<EmitterT> const& e)
			{
				auto next_precision = e.Next();
				if (next_precision < min_precision_)
					min_precision_ = next_precision;
				if (next_precision > max_precision_)
					max_precision_ = next_precision;
				return ObserveResult::KeepObserving;
			}

			/// \brief Get the minimum precision observed.
			unsigned MinPrecision() const
			{
				return min_precision_;
			}

			/// \brief Get the maximum precision observed.
			unsigned MaxPrecision() const
			{
				return max_precision_;
			}

			/// \brief Set the recorded minimum precision.
			void MinPrecision(unsigned m)
			{ min_precision_ = m;}

			/// \brief Set the recorded maximum precision.
			void MaxPrecision(unsigned m)
			{ max_precision_ = m;}

			virtual ~MinMaxPrecisionRecorder() = default;

		private:

			unsigned min_precision_ = std::numeric_limits<unsigned>::max();
			unsigned max_precision_ = 0;
		};


		/// \brief Observer that collects the precision at every tracking event into a vector.
		template<class TrackerT>
		class PrecisionAccumulator : public Observer<TrackerT>
		{ BOOST_TYPE_INDEX_REGISTER_CLASS

			using EmitterT = typename TrackerTraits<TrackerT>::EventEmitterType;

			virtual ObserveResult Observe(AnyEvent const& e) override
			{
				const TrackingEvent<EmitterT>* p = dynamic_cast<const TrackingEvent<EmitterT>*>(&e);
				if (p)
				{
					precisions_.push_back(p->Get().CurrentPrecision());
				}
				return ObserveResult::KeepObserving;
			}


		public:
			/// \brief Get the recorded sequence of precisions.
			const std::vector<unsigned>& Precisions() const
			{
				return precisions_;
			}

			virtual ~PrecisionAccumulator() = default;

		private:
			std::vector<unsigned> precisions_;
		};

		/**
		Example usage:
		PathAccumulator<AMPTracker> path_accumulator;
		*/
		/// \brief Observer that accumulates the current space point at each event of type EventT.
		template<class TrackerT, template<class> class EventT = SuccessfulStep>
		class AMPPathAccumulator
			: public TypedObserver< AMPPathAccumulator<TrackerT,EventT>, TrackerT,
			        EventT<typename TrackerTraits<TrackerT>::EventEmitterType> >
		{ BOOST_TYPE_INDEX_REGISTER_CLASS

			using EmitterT = typename TrackerTraits<TrackerT>::EventEmitterType;

		public:

			/// \brief On each observed event, append the tracker's current point to the path.
			ObserveResult OnEvent(EventT<EmitterT> const& e)
			{
				path_.push_back(e.Get().CurrentPoint());
				return ObserveResult::KeepObserving;
			}

			/// \brief Get the accumulated path of space points.
			const std::vector<Vec<complex_mp> >& Path() const
			{
				return path_;
			}

			virtual ~AMPPathAccumulator() = default;

		private:
			std::vector<Vec<complex_mp> > path_;
		};



		/// \brief Observer that logs verbose ("gory detail") diagnostics for every tracking event.
		template<class TrackerT>
		class GoryDetailLogger : public Observer<TrackerT>
		{ BOOST_TYPE_INDEX_REGISTER_CLASS
		public:

			using EmitterT = typename TrackerTraits<TrackerT>::EventEmitterType;  ///< The event-emitter type.

			virtual ~GoryDetailLogger() = default;

			/// \brief Log a human-readable description of the event, dispatching on its dynamic type.
			virtual ObserveResult Observe(AnyEvent const& e) override
			{


				if (auto p = dynamic_cast<const Initializing<EmitterT,complex_dbl>*>(&e))
				{
					BOOST_LOG_TRIVIAL(severity_level::debug) << std::setprecision(static_cast<int>(p->Get().GetSystem().precision()))
						<< "initializing in double, tracking path\nfrom\tt = "
						<< p->StartTime() << "\nto\tt = " << p->EndTime()
						<< "\n from\tx = \n" << p->StartPoint()
						<< "\n tracking system " << p->Get().GetSystem() << "\n\n";
				}
				else if (auto p = dynamic_cast<const Initializing<EmitterT,complex_mp>*>(&e))
				{
					BOOST_LOG_TRIVIAL(severity_level::debug) << std::setprecision(static_cast<int>(p->Get().GetSystem().precision()))
						 << "initializing in multiprecision, tracking path\nfrom\tt = " << p->StartTime() << "\nto\tt = " << p->EndTime() << "\n from\tx = \n" << p->StartPoint()
						<< "\n tracking system " << p->Get().GetSystem() << "\n\n";
				}

				else if(auto p = dynamic_cast<const TrackingEnded<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::trace) << "tracking ended";

				else if (auto p = dynamic_cast<const NewStep<EmitterT>*>(&e))
				{
					auto& t = p->Get();
					BOOST_LOG_TRIVIAL(severity_level::trace) << "Tracker iteration " << t.NumTotalStepsTaken() << "\ncurrent precision: " << t.CurrentPrecision();


					BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(static_cast<int>(t.CurrentPrecision()))
						<< "t = " << t.CurrentTime()
						<< "\ncurrent stepsize: " << t.CurrentStepsize()
						<< "\ndelta_t = " << t.DeltaT() 
						<< "\ncurrent x size = " << t.CurrentPoint().size()
						<< "\ncurrent x = " << t.CurrentPoint();
				}



				else if (auto p = dynamic_cast<const SingularStartPoint<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::trace) << "singular start point";
				else if (auto p = dynamic_cast<const InfinitePathTruncation<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::trace) << "tracker iteration indicated going to infinity, truncated path";




				else if (auto p = dynamic_cast<const SuccessfulStep<EmitterT>*>(&e))
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << "tracker iteration successful\n\n\n";
				}

				else if (auto p = dynamic_cast<const FailedStep<EmitterT>*>(&e))
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << "tracker iteration unsuccessful\n\n\n";
				}






				else if (auto p = dynamic_cast<const SuccessfulPredict<EmitterT,complex_mp>*>(&e))
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(static_cast<int>(Precision(p->ResultingPoint()))) << "prediction successful (complex_mp), result:\n" << p->ResultingPoint();
				}
				else if (auto p = dynamic_cast<const SuccessfulPredict<EmitterT,complex_dbl>*>(&e))
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(static_cast<int>(Precision(p->ResultingPoint()))) << "prediction successful (complex_dbl), result:\n" << p->ResultingPoint();
				}

				else if (auto p = dynamic_cast<const SuccessfulCorrect<EmitterT,complex_mp>*>(&e))
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(static_cast<int>(Precision(p->ResultingPoint()))) << "correction successful (complex_mp), result:\n" << p->ResultingPoint();
				}
				else if (auto p = dynamic_cast<const SuccessfulCorrect<EmitterT,complex_dbl>*>(&e))
				{
					BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(static_cast<int>(Precision(p->ResultingPoint()))) << "correction successful (complex_dbl), result:\n" << p->ResultingPoint();
				}


				else if (auto p = dynamic_cast<const PredictorHigherPrecisionNecessary<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::trace) << "Predictor, higher precision necessary";
				else if (auto p = dynamic_cast<const CorrectorHigherPrecisionNecessary<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::trace) << "corrector, higher precision necessary";



				else if (auto p = dynamic_cast<const CorrectorMatrixSolveFailure<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::trace) << "corrector, matrix solve failure or failure to converge";
				else if (auto p = dynamic_cast<const PredictorMatrixSolveFailure<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::trace) << "predictor, matrix solve failure or failure to converge";
				else if (auto p = dynamic_cast<const FirstStepPredictorMatrixSolveFailure<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::trace) << "Predictor, matrix solve failure in initial solve of prediction";


				else if (auto p = dynamic_cast<const PrecisionChanged<EmitterT>*>(&e))
					BOOST_LOG_TRIVIAL(severity_level::debug) << "changing precision from " << p->Previous() << " to " << p->Next();

				else
					BOOST_LOG_TRIVIAL(severity_level::debug) << "unlogged event, of type: " << boost::typeindex::type_id_runtime(e).pretty_name();

				return ObserveResult::KeepObserving;
			}

		};



		/// \brief Observer that prints a message to standard out whenever a step fails.
		template<class TrackerT>
		class StepFailScreenPrinter
			: public TypedObserver< StepFailScreenPrinter<TrackerT>, TrackerT,
			        FailedStep<typename TrackerTraits<TrackerT>::EventEmitterType> >
		{ BOOST_TYPE_INDEX_REGISTER_CLASS
		public:

			using EmitterT = typename TrackerTraits<TrackerT>::EventEmitterType;  ///< The event-emitter type.

			/// \brief On a failed step, print a notice to standard out.
			ObserveResult OnEvent(FailedStep<EmitterT> const& /*e*/)
			{
				std::cout << "observed step failure" << std::endl;
				return ObserveResult::KeepObserving;
			}

			virtual ~StepFailScreenPrinter() = default;
		};

	} //re: namespace tracking

}// re: namespace bertini
