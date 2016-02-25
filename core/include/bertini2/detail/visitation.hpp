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
\file visit.hpp

\brief Contains the Visit___ types
*/

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
		virtual ReturnType Visit(VisitedT&) = 0;
		virtual ~Visitor() = default;
	};



	/**
	Defines the default behaviour when an unknown visitor is encountered during visitation.  

	Assumes the return type is default constructible.
	*/
	template<class VisitedT, typename RetT>
	struct DefaultCatchAll
	{
		static RetT OnUnknownVisitor(VisitedT&, VisitorBase&)
		{
			return RetT();
		}
	};


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



	template<class ObservedT, typename RetT = void>
	class Observer : Visitor<ObservedT, RetT>
	{
	public:
		virtual ~Observer() = default;
	};



	namespace {
		template<class ObservedT, typename RetT>
		using ObserverContainer = std::vector<Observer<ObservedT, RetT>*>;
	}

	template< typename RetT = void, template<class,class> class CatchAll = DefaultCatchAll>
	class Observable : public VisitableBase<RetT, CatchAll>
	{
		ObserverContainer<Observable, RetT> watchers_;
	public:
		virtual ~Observable() = default;

		typedef Observer<Observable,RetT> ObsT;

		void AddObserver(ObsT& new_observer)
		{
			watchers_.push_back(new_observer);
		}
	};


}
