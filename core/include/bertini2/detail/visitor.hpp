//This file is part of Bertini 2.0.
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

//  visitor.hpp
//
//  copyright 2016
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

/**
\file visitor.hpp

\brief Contains the visitor base types
*/

#ifndef BERTINI_DETAIL_VISITOR_HPP
#define BERTINI_DETAIL_VISITOR_HPP

#include "bertini2/detail/events.hpp"

namespace bertini{

	/**
	A strawman class for implementing the visitor pattern.  

	Deliberately empty, the non-base visitors must derive from it.  Don't directly derive from this class, use Visitor.
	*/
	class VisitorBase
	{
	public:
		virtual ~VisitorBase() = default;
	};


	/**
	The first non-trivial visitor class.  Derive from this.
	*/
	template<class VisitedT, typename RetT = void>
	class Visitor
	{
	public:
		typedef RetT ReturnType;
		virtual ReturnType Visit(VisitedT const &) = 0;
		virtual ~Visitor() = default;
	};

	

	class AnyObserver
	{
	public:
		virtual ~AnyObserver() = default;
		virtual void Observe(AnyEvent const& e) = 0;
	};

	template<class ObservedT, typename RetT = void>
	class Observer : public Visitor<ObservedT, RetT>, public AnyObserver
	{
	public:
		virtual ~Observer() = default;

		
	};


	


}

#endif
