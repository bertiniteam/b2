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


	template<class T>
	class TrackingEvent : public Event<T>
	{
	public:
		
		TrackingEvent(const T & obs) : Event<T>(obs)
		{}

		virtual ~TrackingEvent() = default;

		TrackingEvent() = delete;
	private:
		
	};

	template<class T>
	class SuccessfulStep : public TrackingEvent<T>
	{
	public:
		SuccessfulStep(const T & obs) : TrackingEvent<T>(obs)
		{}

		virtual ~SuccessfulStep() = default;
		SuccessfulStep() = delete;
	};

}// re: namespace tracking
}// re: namespace bertini

