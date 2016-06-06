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
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Daniel Brake
// University of Notre Dame
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


#ifndef BERTINI_DETAIL_EVENTS_HPP
#define BERTINI_DETAIL_EVENTS_HPP
#include <boost/type_index.hpp>
namespace bertini {


	class AnyEvent
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~AnyEvent() = default;
	};

	template<class ObsT>
	class Event : public AnyEvent
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
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


	#define ADD_BERTINI_EVENT_TYPE(event_name,event_parenttype) \
	template<class ObservedT> \
	class event_name : public event_parenttype<ObservedT> \
	{ BOOST_TYPE_INDEX_REGISTER_CLASS \
	public: \
		event_name(const ObservedT & obs) : event_parenttype<ObservedT>(obs){} \
		virtual ~event_name() = default; \
		event_name() = delete; \
	}
} //re: namespace bertini

#endif
