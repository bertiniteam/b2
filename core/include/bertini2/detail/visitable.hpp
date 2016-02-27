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

#ifndef BERTINI_DETAIL_VISITABLE_HPP
#define BERTINI_DETAIL_VISITABLE_HPP

#include "detail/visitor.hpp"
#include "detail/events.hpp"

namespace bertini{


	/**
	Defines the default behaviour when an unknown visitor is encountered during visitation.  

	Assumes the return type is default constructible.
	*/
	template<class VisitedT, typename RetT>
	struct DefaultConstruct
	{
		static RetT OnUnknownVisitor(VisitedT&, VisitorBase&)
		{
			return RetT();
		}
	};

	template<class VisitedT, typename RetT>
	using DefaultCatchAll = DefaultConstruct<VisitedT, RetT>;


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


	template<typename RetT = void, template<class,class> class CatchAll = DefaultCatchAll>
	class Observable : public VisitableBase<RetT, CatchAll>
	{	
	public:

		virtual ~Observable() = default;



		void AddObserver(AnyObserver* new_observer)
		{
			std::cout << "added observer\n";
			current_watchers_.push_back(new_observer);
		}

	protected:


		void NotifyObservers(AnyEvent const& e) const
		{

			for (auto& obs : current_watchers_)
				obs->Observe(e);

		}


	private:

		using ObserverContainer = std::vector<AnyObserver*>;

		ObserverContainer current_watchers_;
	};

} // namespace bertini


#endif

