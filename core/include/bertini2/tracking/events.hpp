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

namespace bertini {

	namespace tracking{


	template<class ObservedT>
	class TrackingEvent : public Event<ObservedT>
	{
	public:
		
		TrackingEvent(const ObservedT & obs) : Event<ObservedT>(obs)
		{}

		virtual ~TrackingEvent() = default;

		TrackingEvent() = delete;
	private:
		
	};

	#define BERTINI_EVENT_TYPE(event_name,event_subtype) \
	template<class ObservedT> \
	class event_name : public event_subtype<ObservedT> \
	{ \
	public: \
		event_name(const ObservedT & obs) : event_subtype<ObservedT>(obs){} \
		virtual ~event_name() = default; \
		event_name() = delete; \
	};

	BERTINI_EVENT_TYPE(HigherPrecisionNecessary,TrackingEvent)


	template<class ObservedT>
	class SuccessfulStep : public TrackingEvent<ObservedT>
	{
	public:
		SuccessfulStep(const ObservedT & obs) : TrackingEvent<ObservedT>(obs)
		{}

		virtual ~SuccessfulStep() = default;
		SuccessfulStep() = delete;
	};

	template<class ObservedT>
	class FailedStep : public TrackingEvent<ObservedT>
	{
	public:
		FailedStep(const ObservedT & obs) : TrackingEvent<ObservedT>(obs)
		{}

		virtual ~FailedStep() = default;
		FailedStep() = delete;
	};




	////////////
	//
	//  Precision events

	template<class ObservedT>
	class PrecisionChanged : public TrackingEvent<ObservedT>
	{
	public:
		PrecisionChanged(const ObservedT & obs) : TrackingEvent<ObservedT>(obs)
		{}

		virtual ~PrecisionChanged() = default;
		PrecisionChanged() = delete;
	};

	template<class ObservedT>
	class PrecisionIncreased : public PrecisionChanged<ObservedT>
	{
	public:
		PrecisionIncreased(const ObservedT & obs) : PrecisionChanged<ObservedT>(obs)
		{}

		virtual ~PrecisionIncreased() = default;
		PrecisionIncreased() = delete;
	};

	template<class ObservedT>
	class PrecisionDecreased : public PrecisionChanged<ObservedT>
	{
	public:
		PrecisionDecreased(const ObservedT & obs) : PrecisionChanged<ObservedT>(obs)
		{}

		virtual ~PrecisionDecreased() = default;
		PrecisionDecreased() = delete;
	};


	////////////
	//
	//  Stepsize events

	template<class ObservedT>
	class StepsizeChanged : public TrackingEvent<ObservedT>
	{
	public:
		StepsizeChanged(const ObservedT & obs) : TrackingEvent<ObservedT>(obs)
		{}

		virtual ~StepsizeChanged() = default;
		StepsizeChanged() = delete;
	};

	template<class ObservedT>
	class StepsizeIncreased : public StepsizeChanged<ObservedT>
	{
	public:
		StepsizeIncreased(const ObservedT & obs) : StepsizeChanged<ObservedT>(obs)
		{}

		virtual ~StepsizeIncreased() = default;
		StepsizeIncreased() = delete;
	};

	template<class ObservedT>
	class StepsizeDecreased : public StepsizeChanged<ObservedT>
	{
	public:
		StepsizeDecreased(const ObservedT & obs) : StepsizeChanged<ObservedT>(obs)
		{}

		virtual ~StepsizeDecreased() = default;
		StepsizeDecreased() = delete;
	};


	///////////
	//
	//  beginning and end events

	template<class ObservedT>
	class TrackingStarted : public TrackingEvent<ObservedT>
	{
	public:
		TrackingStarted(const ObservedT & obs) : TrackingEvent<ObservedT>(obs)
		{}

		virtual ~TrackingStarted() = default;
		TrackingStarted() = delete;
	};

	template<class ObservedT>
	class TrackingFinished : public TrackingEvent<ObservedT>
	{
	public:
		TrackingFinished(const ObservedT & obs) : TrackingEvent<ObservedT>(obs)
		{}

		virtual ~TrackingFinished() = default;
		TrackingFinished() = delete;
	};
}// re: namespace tracking
}// re: namespace bertini

