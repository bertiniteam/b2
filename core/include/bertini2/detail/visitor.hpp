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
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Daniel Brake
// University of Notre Dame
//

/**
\file visitor.hpp

\brief Contains the visitor base types
*/

#ifndef BERTINI_DETAIL_VISITOR_HPP
#define BERTINI_DETAIL_VISITOR_HPP

#include <tuple>
#include <utility>

#include <boost/fusion/include/std_tuple.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>

#include "bertini2/detail/events.hpp"

namespace bertini{

	/**
	A strawman class for implementing the visitor pattern.  

	Deliberately empty, the non-base visitors must derive from it.  Don't directly derive from this class, use Visitor.
	*/
	class VisitorBase
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~VisitorBase() = default;
	};


	/**
	The first non-trivial visitor class.  Derive from this.
	*/
	template<class VisitedT, typename RetT = void>
	class Visitor
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		typedef RetT ReturnType;
		virtual ReturnType Visit(VisitedT const &) = 0;
		virtual ~Visitor() = default;
	};

	

	class AnyObserver
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~AnyObserver() = default;
		virtual void Observe(AnyEvent const& e) = 0;
	};

	template<class ObservedT, typename RetT = void>
	class Observer : public Visitor<ObservedT, RetT>, public AnyObserver
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~Observer() = default;

		
	};



	

	template<class ObservedT, template<class> class... ObserverTypes>
	class MultiObserver : public Observer<ObservedT>
	{	BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		void Observe(AnyEvent const& e) override
		{	
		    using namespace boost::fusion;
		    auto f = [&e](auto &obs) { obs.Observe(e); };
		    for_each(observers_, f);
		}

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
