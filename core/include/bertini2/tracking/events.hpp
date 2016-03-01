//This file is part of Bertini 2.0.
//
//tracking/events.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracking/events.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/events.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  tracking/events.hpp
//
//  copyright 2016
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

/**
\file tracking/events.hpp

\brief Contains the tracking/events base types
*/

#pragma once
#include "bertini2/detail/events.hpp"
#include "limbo.hpp"

namespace bertini {

	namespace tracking{


	

	ADD_BERTINI_EVENT_TYPE(TrackingEvent,Event)

	ADD_BERTINI_EVENT_TYPE(SuccessfulStep,TrackingEvent)

	ADD_BERTINI_EVENT_TYPE(FailedStep,TrackingEvent)

	ADD_BERTINI_EVENT_TYPE(NewStep,TrackingEvent)

	ADD_BERTINI_EVENT_TYPE(SingularStartPoint,TrackingEvent)

	template<class ObservedT, typename NumT>
	class SuccessfulPredict : public TrackingEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		SuccessfulPredict(const ObservedT & obs, 
		             Vec<NumT> const& resulting_point) : TrackingEvent<ObservedT>(obs),
													resulting_point_(resulting_point)
		{}


		virtual ~SuccessfulPredict() = default;
		SuccessfulPredict() = delete;

		const Vec<NumT>& ResultingPoint() const {return resulting_point_;}
	private:
		const Vec<NumT>& resulting_point_;
	};

	template<class ObservedT, typename NumT>
	class SuccessfulCorrect : public TrackingEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		SuccessfulCorrect(const ObservedT & obs, 
		             Vec<NumT> const& resulting_point) : TrackingEvent<ObservedT>(obs),
													resulting_point_(resulting_point)
		{}


		virtual ~SuccessfulCorrect() = default;
		SuccessfulCorrect() = delete;
		
		const Vec<NumT>& ResultingPoint() const {return resulting_point_;}
	private:
		const Vec<NumT>& resulting_point_;
	};

	////////////
	//
	//  Precision events

	ADD_BERTINI_EVENT_TYPE(PrecisionEvent,TrackingEvent)

	template<class ObservedT>
	class PrecisionChanged : public PrecisionEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		PrecisionChanged(const ObservedT & obs, 
		             unsigned previous, unsigned next) : PrecisionEvent<ObservedT>(obs),
													prev_(previous),
													next_(next)
		{}


		virtual ~PrecisionChanged() = default;
		PrecisionChanged() = delete;
		
		auto Previous() const {return prev_;}
		auto Next() const {return next_;}
	private:
		const unsigned prev_, next_;
	};

	template<class ObservedT>
	class PrecisionIncreased : public PrecisionEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		PrecisionIncreased(const ObservedT & obs, 
		             unsigned previous, unsigned next) : PrecisionChanged<ObservedT>(obs, previous, next)
		{}
		virtual ~PrecisionIncreased() = default;
		PrecisionIncreased() = delete;
	};

	template<class ObservedT>
	class PrecisionDecreased : public PrecisionEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		PrecisionDecreased(const ObservedT & obs, 
		             unsigned previous, unsigned next) : PrecisionChanged<ObservedT>(obs, previous, next)
		{}
		virtual ~PrecisionDecreased() = default;
		PrecisionDecreased() = delete;
	};

	ADD_BERTINI_EVENT_TYPE(HigherPrecisionNecessary,PrecisionEvent)
	ADD_BERTINI_EVENT_TYPE(PredictorHigherPrecisionNecessary,HigherPrecisionNecessary)
	ADD_BERTINI_EVENT_TYPE(CorrectorHigherPrecisionNecessary,HigherPrecisionNecessary)


	ADD_BERTINI_EVENT_TYPE(MatrixSolveFailure,PrecisionEvent)
	ADD_BERTINI_EVENT_TYPE(PredictorMatrixSolveFailure,MatrixSolveFailure)
	ADD_BERTINI_EVENT_TYPE(FirstStepPredictorMatrixSolveFailure,MatrixSolveFailure)
	ADD_BERTINI_EVENT_TYPE(CorrectorMatrixSolveFailure,MatrixSolveFailure)
	////////////
	//
	//  Stepsize events

	ADD_BERTINI_EVENT_TYPE(StepsizeEvent,TrackingEvent)

	ADD_BERTINI_EVENT_TYPE(StepsizeDecreased,StepsizeEvent)

	ADD_BERTINI_EVENT_TYPE(StepsizeIncreased,StepsizeEvent)


	///////////
	//
	//  beginning and end events

	ADD_BERTINI_EVENT_TYPE(TrackingStarted,TrackingEvent)

	ADD_BERTINI_EVENT_TYPE(TrackingEnded, TrackingEvent)

	ADD_BERTINI_EVENT_TYPE(InfinitePathTruncation, TrackingEvent)

	template<class ObservedT, typename NumT>
	class Initializing : public TrackingEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
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

		const NumT& StartTime() const {return start_time_;}
		const NumT& EndTime() const {return end_time_;}
		const Vec<NumT>& StartPoint() const {return start_point_;}
	private:
		const NumT& start_time_;
		const NumT& end_time_;
		const Vec<NumT>& start_point_;
	};
}// re: namespace tracking
}// re: namespace bertini

