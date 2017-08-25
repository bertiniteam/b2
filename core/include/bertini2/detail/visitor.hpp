//This file is part of Bertini 2.
//
//visitor.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//visitor.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with visitor.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Dani Brake
// University of Notre Dame
//

/**
\file visitor.hpp

\brief Contains the visitor base types.
*/

#ifndef BERTINI_DETAIL_VISITOR_HPP
#define BERTINI_DETAIL_VISITOR_HPP

#include <tuple>
#include <utility>

#include <boost/fusion/adapted/std_tuple.hpp>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>

#include "bertini2/detail/events.hpp"

namespace bertini{

	/**
	\brief A strawman class for implementing the visitor pattern.  

	Deliberately empty, the non-base visitors must derive from it.  Don't directly derive from this class, use Visitor.

	\see Visitor, Observer, AnyObserver, MultiObserver
	*/
	class VisitorBase
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~VisitorBase() = default;
	};


	/**
	\brief The first non-trivial visitor class.  

	Derive from this when creating new visitor types.

	\tparam VisitedT The type of object the visitor visits.
	\tparam RetT The type of object to be returned when visiting.  Default is `void`.
	*/
	template<class VisitedT, typename RetT = void>
	class Visitor
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		typedef RetT ReturnType;
		virtual ReturnType Visit(VisitedT const &) = 0;
		virtual ~Visitor() = default;
	};

	
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
		*/
		virtual void Observe(AnyEvent const& e) = 0;
	};


	/**
	\brief Actual observer type, which you should derive from to extract custom information from observable types.

	\tparam ObservedT The type of object the observer observes.  
	\tparam RetT The type of object the observer returns when it visits.

	\see PrecisionAccumulator, GoryDetailLogger, MultiObserver
	*/
	template<class ObservedT, typename RetT = void>
	class Observer : public Visitor<ObservedT, RetT>, public AnyObserver
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

		\param e The emitted event which caused observation.
		*/
		void Observe(AnyEvent const& e) override
		{	
		    using namespace boost::fusion;
		    auto f = [&e](auto &obs) { obs.Observe(e); };
		    for_each(observers_, f);
		}

		/**
		Since Observers are also Visitors, this function calls Visit override for each of the observer types.
		*/
		void Visit(ObservedT const& t) override
		{
			using namespace boost::fusion;
		    auto f = [&t](auto &obs) { obs.Visit(t); };
		    for_each(observers_, f);
		}

		std::tuple<ObserverTypes<ObservedT>...> observers_;
		virtual ~MultiObserver() = default;
	};


}

#endif
