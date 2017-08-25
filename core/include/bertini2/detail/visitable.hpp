//This file is part of Bertini 2.
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
\file visitable.hpp

\brief Contains the visitable base types
*/

#ifndef BERTINI_DETAIL_VISITABLE_HPP
#define BERTINI_DETAIL_VISITABLE_HPP

#include "bertini2/detail/visitor.hpp"
#include "bertini2/detail/events.hpp"

namespace bertini{

	namespace policy{

		/**
		\brief A policy for what to do when a visited type is encoutered by a visitor it doesn't know.

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

		/**
		The default policy for what to do when a visitable is visited by an unknown visitor.
		*/
		template<class VisitedT, typename RetT>
		using DefaultCatchAll = DefaultConstruct<VisitedT, RetT>;
	}

	


	/**
	\brief The base class for visitable types.  

	Implemented based on Alexandrescu, 2001, and Hythem Sidky's SAPHRON package, with his permission.

	\see BERTINI_DEFAULT_VISITABLE
	*/
	template< typename RetT = void, template<class,class> class CatchAll = policy::DefaultCatchAll>
	class VisitableBase
	{
	public:
		typedef RetT ReturnType;
		virtual ~VisitableBase() = default;
		virtual ReturnType Accept(VisitorBase&) = 0; // the implementation will either be provided by a macro, or by hand, for each class which is visitable.

	protected:

		/**
		\brief Abstract method for how to accept a visitor.  

		This function simply Forwards to the Visit method of the visitor, if the visitor's type is known.  If known, invokes the behaviour defined by the CatchAll template parameter.

		\tparam T The type of the visited object.  Should be inferred by the compiler.

		\param visited The visited visitable object.
		\param guest The visiting object.
		*/
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

	
	/** 
	\brief macro for classes which want default Accept implementation, having nothing fancy to do when accepting.
	*/
	#define BERTINI_DEFAULT_VISITABLE() \
		virtual ReturnType Accept(VisitorBase& guest) override \
		{ return AcceptBase(*this, guest); }


	/**
	\brief An abstract observable type, maintaining a list of observers, who can be notified in case of Events.
	
	Some known observable types are Tracker and Endgame.

	\tparam RetT The return type of the Visit method of the observer or visitor.  Default is `void`.
	\tparam CatchAll The policy to be invoked when the visited type doesn't know the visitor.  Default is the DefaultCatchAll.
	*/
	template<typename RetT = void, template<class,class> class CatchAll = policy::DefaultCatchAll>
	class Observable : public VisitableBase<RetT, CatchAll>
	{	
	public:

		virtual ~Observable() = default;


		/**
		\brief Add an observer, to observe this observable.
		*/
		void AddObserver(AnyObserver* new_observer) const
		{
			if (find(begin(current_watchers_), end(current_watchers_), new_observer)==end(current_watchers_))
				current_watchers_.push_back(new_observer);
		}

		/**
		\brief Remove an observer from this observable.
		*/
		void RemoveObserver(AnyObserver* observer) const
		{
			current_watchers_.erase(std::remove(current_watchers_.begin(), current_watchers_.end(), observer), current_watchers_.end());
		}

	protected:

		/**
		\brief Sends an Event (more particularly, AnyEvent) to all watching observers of this object.

		This function could potentially be improved by filtering on the observer's desired event types, if known at compile time.  This could potentially be a performance bottleneck (hopefully not!) since filtering can use `dynamic_cast`ing.  One hopes this cost is overwhelmed by things like linear algebra and system evaluation.

		\param e The event to emit.  Its type should be derived from AnyEvent.
		*/
		void NotifyObservers(AnyEvent const& e) const
		{

			for (auto& obs : current_watchers_)
				obs->Observe(e);

		}

		void NotifyObservers(AnyEvent & e) const
		{

			for (auto& obs : current_watchers_)
				obs->Observe(e);
		}


	private:

		using ObserverContainer = std::vector<AnyObserver*>;

		mutable ObserverContainer current_watchers_;
	};

} // namespace bertini


#endif

