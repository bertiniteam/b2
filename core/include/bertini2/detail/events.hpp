//This file is part of Bertini 2.0.
//
//detail/events.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//detail/events.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with detail/events.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  detail/events.hpp
//
//  copyright 2016
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

/**
\file detail/events.hpp

\brief Contains the detail/events base types
*/


#pragma once

namespace bertini {


	class EventBase
	{
	public:
		virtual ~EventBase() = default;
	};

	template<class ObsT>
	class Event : public EventBase
	{
	public:
		Event(ObsT const& obs) : current_observable_(obs)
		{}

		virtual ~Event() = default;

		ObsT const& Get() const
		{return current_observable_;}

		Event() = delete;
	protected:
		const ObsT& current_observable_;

		
	};
} //re: namespace bertini


