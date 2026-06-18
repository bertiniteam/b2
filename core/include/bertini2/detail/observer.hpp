//This file is part of Bertini 2.
//
//bertini2/detail/observer.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/detail/observer.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/detail/observer.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire
//

/**
\file bertini2/detail/observer.hpp

\brief Contains the observer base types.

\defgroup observer

*/

#ifndef BERTINI_DETAIL_OBSERVER_HPP
#define BERTINI_DETAIL_OBSERVER_HPP

#include <tuple>
#include <typeindex>
#include <utility>
#include <vector>

#include <boost/fusion/adapted/std_tuple.hpp>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>

#include "bertini2/detail/events.hpp"

namespace bertini{


	/**
	\brief The result an observer returns from Observe(), telling the observable
	whether to keep sending it events.

	Returning `Unsubscribe` is the clean, dispatch-safe way for an observer to
	stop receiving events: the observable removes it *after* the current
	notification loop finishes (no mid-iteration mutation of the observer list).
	This replaces the old pattern of calling `RemoveObserver(*this)` from inside
	`Observe()`.
	*/
	enum class ObserveResult { KeepObserving, Unsubscribe };


	/**
	\brief Strawman base class for Observer objects.

	\see Observer
	*/
	class AnyObserver
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~AnyObserver() = default;

		/**
		\brief Observe the observable object being observed.  This is probably in response to NotifyObservers.

		This virtual function must be overridden by actual observers, defining how they observe the observable they are observing, probably filtering events and doing something specific for different ones.

		\param e The event which was emitted by the observed object.
		\return `ObserveResult::KeepObserving` to keep receiving events, or
		        `ObserveResult::Unsubscribe` to ask the observable to drop this
		        observer once the current notification finishes.
		*/
		virtual ObserveResult Observe(AnyEvent const& e) = 0;

		/**
		\brief Declares which event types this observer wants to receive.

		Return a non-empty vector to opt into type-indexed dispatch: the observable will
		only call Observe() for events whose dynamic type exactly matches one of the
		returned type_index values.  Return an empty vector (the default) to receive
		every event — this is the correct choice for observers that handle many or all
		event types, and is the automatic behavior for Python-defined observers.
		*/
		virtual std::vector<std::type_index> SubscribedEventTypes() const { return {}; }
	};


	/**
	\brief Actual observer type, which you should derive from to extract custom information from observable types.

	\tparam ObservedT The type of object the observer observes.  
	\tparam RetT The type of object the observer returns when it visits.

	\see PrecisionAccumulator, GoryDetailLogger, MultiObserver
	*/
	template<class ObservedT>
	class Observer : public AnyObserver
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~Observer() = default;

		
	};



	
	/**
	\brief A class which can glob together observer types into a new, single observer type.

	If there are pre-existing observers for the object you wish to observe, rather than making one of each, and attaching each to the observable, you can make many things one.

	https://frinkiac.com/?q=many%20guns%20into%20five

	\tparam ObservedT The type of thing the observer types you are gluing together observe.  They must all observe the same type of object.
	\tparam ObserverTypes The already-existing observer types you are gluing together.  You can put as many of them together as you want!
	*/
	template<class ObservedT, template<class> class... ObserverTypes>
	class MultiObserver : public Observer<ObservedT>
	{	BOOST_TYPE_INDEX_REGISTER_CLASS
	public:

		/**
		\brief Observe override which calls the overrides for the types you glued together.

		The bundle keeps observing as long as at least one of its children still
		wants events; it asks to unsubscribe only once every glued-in observer has
		returned Unsubscribe.

		\param e The emitted event which caused observation.
		*/
		ObserveResult Observe(AnyEvent const& e) override
		{
		    using namespace boost::fusion;
		    bool keep = false;
		    auto f = [&e,&keep](auto &obs) {
		        if (obs.Observe(e) == ObserveResult::KeepObserving)
		            keep = true;
		    };
		    for_each(observers_, f);
		    return keep ? ObserveResult::KeepObserving : ObserveResult::Unsubscribe;
		}

		std::tuple<ObserverTypes<ObservedT>...> observers_;
		virtual ~MultiObserver() = default;
	};


}

#endif
