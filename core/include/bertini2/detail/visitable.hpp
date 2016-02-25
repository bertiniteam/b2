//This file is part of Bertini 2.0.
//
//visitable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//visitable.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with visitable.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  visitable.hpp
//
//  copyright 2016
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

/**
\file visitable.hpp

\brief Contains the visitable base types
*/

#pragma once

#include "detail/visitor.hpp"

namespace bertini{
	/**
	The base class for visitable types.  

	Implemented based on Alexandrescu, 2001, and Hythem Sidky's SAPHRON package, with his permission.
	*/
	template< typename RetT = void, template<class,class> class CatchAll = DefaultCatchAll>
	class VisitableBase
	{
	public:
		typedef RetT ReturnType;
		virtual ~VisitableBase() = default;
		virtual ReturnType Accept(VisitorBase&) = 0; // the implementation will either be provided by a macro, or by hand, for each class which is visitable.

	protected:
		template<typename T>
		static
		ReturnType AcceptBase(T& visited, VisitorBase& guest)
		{
			if (auto p = dynamic_cast<Visitor<T>*>(&guest))
				return p->Visit(visited);
			else
				return CatchAll<T,RetT>::OnUnknownVisitor(visited, guest);
		}
	};

	
	/** provide a macro for classes which want default Accept implementation, having nothing fancy to do when accepting.
	*/
	#define BERTINI_DEFAULT_VISITABLE() \
		virtual ReturnType Accept(VisitorBase& guest) override \
		{ return AcceptBase(*this, guest); }


	template< typename RetT = void, template<class,class> class CatchAll = DefaultCatchAll>
	class Observable : public VisitableBase<RetT, CatchAll>
	{	
	public:
		typedef RetT ReturnType;
		virtual ~Observable() = default;

		typedef Observer<Observable,RetT> ObsT;

		void AddObserver(ObsT& new_observer)
		{
			current_watchers_.push_back(new_observer);
		}

	private:
		template<class ObservedT>
		using ObserverContainer = std::vector<Observer<ObservedT, ReturnType>*>;

		ObserverContainer<Observable> current_watchers_;
	};

} // namespace bertini
