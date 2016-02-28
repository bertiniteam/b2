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


	

	ADD_BERTINI_EVENT_TYPE(TrackingEvent,Event)

	ADD_BERTINI_EVENT_TYPE(SuccessfulStep,TrackingEvent)

	ADD_BERTINI_EVENT_TYPE(FailedStep,TrackingEvent)




	////////////
	//
	//  Precision events

	ADD_BERTINI_EVENT_TYPE(PrecisionEvent,TrackingEvent)

	ADD_BERTINI_EVENT_TYPE(PrecisionIncreased,PrecisionEvent)

	ADD_BERTINI_EVENT_TYPE(PrecisionDecreased,PrecisionEvent)

	ADD_BERTINI_EVENT_TYPE(HigherPrecisionNecessary,PrecisionEvent)

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
	
}// re: namespace tracking
}// re: namespace bertini

