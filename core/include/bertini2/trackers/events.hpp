//This file is part of Bertini 2.
//
//events.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//events.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with events.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire


/**
\file tracking/events.hpp

\brief Contains the tracking/events base types
*/

#pragma once
#include "bertini2/detail/events.hpp"
#include "bertini2/limbo.hpp"

namespace bertini {

	namespace tracking{


	
	/**
	\brief Generic event for Tracking
	*/
	ADD_BERTINI_EVENT_TYPE(TrackingEvent,ConstEvent);

	/**
	\brief A successful step occurred
	*/
	ADD_BERTINI_EVENT_TYPE(SuccessfulStep,TrackingEvent);

	/**
	\brief A failed step occurred
	*/
	ADD_BERTINI_EVENT_TYPE(FailedStep,TrackingEvent);

	/**
	\brief Taking a new step -- beginning of procedure for attempting to step forward
	*/
	ADD_BERTINI_EVENT_TYPE(NewStep,TrackingEvent);

	/**
	\brief The start point for the TrackPath call was singular

	This is evidenced by step size shrinking too far, or running out of precision.
	*/
	ADD_BERTINI_EVENT_TYPE(SingularStartPoint,TrackingEvent);

	/**
	\brief The predict part of tracking step was successful.
	*/
	template<class ObservedT, typename NumT>
	class SuccessfulPredict : public TrackingEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		/**
		\brief The constructor for a SuccessfulPredict Event.

		\param obs The observable emitting the event.
		\param resulting_point The space point which is the result of the prediction.
		*/
		SuccessfulPredict(const ObservedT & obs, 
		             Vec<NumT> const& resulting_point) : TrackingEvent<ObservedT>(obs),
													resulting_point_(resulting_point)
		{}


		virtual ~SuccessfulPredict() = default;
		SuccessfulPredict() = delete;

		/**
		\brief Get the resulting point of the prediction.
		*/
		const Vec<NumT>& ResultingPoint() const {return resulting_point_;}
	private:
		const Vec<NumT>& resulting_point_;
	};

	/**
	\brief The correct part of a time step was successful
	*/
	template<class ObservedT, typename NumT>
	class SuccessfulCorrect : public TrackingEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		/**
		\brief The constructor for a SuccessfulCorrect Event.

		\param obs The observable emitting the event.
		\param resulting_point The space point which is the result of the Newton correct solve.
		*/
		SuccessfulCorrect(const ObservedT & obs, 
		             Vec<NumT> const& resulting_point) : TrackingEvent<ObservedT>(obs),
													resulting_point_(resulting_point)
		{}


		virtual ~SuccessfulCorrect() = default;
		SuccessfulCorrect() = delete;
		
		/**
		\brief Get the resulting point of the correction.
		*/
		const Vec<NumT>& ResultingPoint() const {return resulting_point_;}
	private:
		const Vec<NumT>& resulting_point_;
	};

	////////////
	//
	//  Precision events

	/**
	\brief A generic event involving precision
	*/
	ADD_BERTINI_EVENT_TYPE(PrecisionEvent,TrackingEvent);

	template<class ObservedT>
	class PrecisionChanged : public PrecisionEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		/**
		\brief The constructor for a PrecisionChanged Event.

		\param obs The observable emitting the event.
		\param previous The precision before changing.
		\param next The precision after changing.
		*/
		PrecisionChanged(const ObservedT & obs, 
		             unsigned previous, unsigned next) : PrecisionEvent<ObservedT>(obs),
													prev_(previous),
													next_(next)
		{}


		virtual ~PrecisionChanged() = default;
		PrecisionChanged() = delete;
		
		/**
		\brief Get the previous precision.
		*/
		auto Previous() const {return prev_;}

		/**
		\brief Get the next precision, what it changed to.
		*/
		auto Next() const {return next_;}
	private:
		const unsigned prev_, next_;
	};

	/**
	\brief Precision increased during tracking
	*/
	template<class ObservedT>
	class PrecisionIncreased : public PrecisionChanged<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		/**
		\brief The constructor for a PrecisionIncreased Event.

		\param obs The observable emitting the event.
		\param previous The precision before changing.
		\param next The precision after changing.
		*/
		PrecisionIncreased(const ObservedT & obs, 
		             unsigned previous, unsigned next) : PrecisionChanged<ObservedT>(obs, previous, next)
		{}
		virtual ~PrecisionIncreased() = default;
		PrecisionIncreased() = delete;
	};

	/**
	\brief Precision decreased during tracking
	*/
	template<class ObservedT>
	class PrecisionDecreased : public PrecisionChanged<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		/**
		\brief The constructor for a PrecisionDecreased Event.

		\param obs The observable emitting the event.
		\param previous The precision before changing.
		\param next The precision after changing.
		*/
		PrecisionDecreased(const ObservedT & obs, 
		             unsigned previous, unsigned next) : PrecisionChanged<ObservedT>(obs, previous, next)
		{}
		virtual ~PrecisionDecreased() = default;
		PrecisionDecreased() = delete;
	};

	/**
	\brief A step failed, because precision needs to increase
	*/
	ADD_BERTINI_EVENT_TYPE(HigherPrecisionNecessary,PrecisionEvent);

	/**
	\brief Prediction failed, because precision needs to increase
	*/
	ADD_BERTINI_EVENT_TYPE(PredictorHigherPrecisionNecessary,HigherPrecisionNecessary);

	/**
	\brief Correction failed, because precision needs to increase
	*/
	ADD_BERTINI_EVENT_TYPE(CorrectorHigherPrecisionNecessary,HigherPrecisionNecessary);

	/**
	\brief Linear algebra declared a failure
	*/
	ADD_BERTINI_EVENT_TYPE(MatrixSolveFailure,PrecisionEvent);

	/**
	\brief Linear algebra declared a failure, during prediction
	*/
	ADD_BERTINI_EVENT_TYPE(PredictorMatrixSolveFailure,MatrixSolveFailure);

	/**
	\brief Linear algebra declared a failure, during the first step of a multi-point prediction.  
	*/
	ADD_BERTINI_EVENT_TYPE(FirstStepPredictorMatrixSolveFailure,MatrixSolveFailure);

	/**
	\brief Linear algebra declared a failure, during correction
	*/
	ADD_BERTINI_EVENT_TYPE(CorrectorMatrixSolveFailure,MatrixSolveFailure);
	////////////
	//
	//  Stepsize events

	/**
	\brief Stepsize is changing
	*/
	ADD_BERTINI_EVENT_TYPE(StepsizeEvent,TrackingEvent);

	/**
	\brief Stepsize decreased

	This means that the step failed, and decreasing step size was the best thing to do.
	*/
	ADD_BERTINI_EVENT_TYPE(StepsizeDecreased,StepsizeEvent);

	/**
	\brief Stepsize increased. 

	This means steps have been successful lately.
	*/
	ADD_BERTINI_EVENT_TYPE(StepsizeIncreased,StepsizeEvent);


	///////////
	//
	//  beginning and end events

	/**
	\brief Made a call to TrackPath
	*/
	ADD_BERTINI_EVENT_TYPE(TrackingStarted,TrackingEvent);

	/**
	\brief Tracking is stopping for whatever reason.
	*/
	ADD_BERTINI_EVENT_TYPE(TrackingEnded, TrackingEvent);

	/**
	\brief Tracking terminated because the path was going to infinity.
	*/
	ADD_BERTINI_EVENT_TYPE(InfinitePathTruncation, TrackingEvent);

	/**
	\brief TrackPath is initializing the tracker.
	*/
	template<class ObservedT, typename NumT>
	class Initializing : public TrackingEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:

		/**
		\brief Constructor for an Initializing Event

		\param obs The observed object, the tracker.  `*this` probably.
		\param start_time The time \f$t_0\f$ at which tracking is starting.
		\param end_time The target time for tracking.
		\param start_point The space point \f$x_0\f$ for starting tracking.
		*/
		Initializing(const ObservedT & obs, 
		             NumT const& start_time, 
		             NumT const& endtime, 
		             Vec<NumT> const& start_point) : TrackingEvent<ObservedT>(obs),
													start_time_(start_time),
													end_time_(endtime),
													start_point_(start_point)
		{}


		virtual ~Initializing() = default;
		Initializing() = delete;

		/**
		\brief Get the time (tracking wise, not clock wise) at which tracking is starting/
		*/
		const NumT& StartTime() const {return start_time_;}

		/**
		\brief Get the target time for tracking.
		*/
		const NumT& EndTime() const {return end_time_;}

		/**
		\brief Get the start point, the input point for TrackPath
		*/
		const Vec<NumT>& StartPoint() const {return start_point_;}
	private:
		const NumT& start_time_;
		const NumT& end_time_;
		const Vec<NumT>& start_point_;
	};
}// re: namespace tracking
}// re: namespace bertini

