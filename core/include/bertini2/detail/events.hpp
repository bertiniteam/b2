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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire
//

//  detail/events.hpp


/**
\file detail/events.hpp

\brief Contains the detail/events base types
*/


#ifndef BERTINI_DETAIL_EVENTS_HPP
#define BERTINI_DETAIL_EVENTS_HPP
#include <boost/type_index.hpp>
namespace bertini {

	/**
	\brief Strawman Event type, enabling polymorphism.

	This class is abstract, and should never actually be created.

	\par Lifetime contract (important!)
	Events are short-lived stack temporaries.  An observable emits one with
	`NotifyObservers(SomeEvent<T>(*this))`, so the event lives only for the
	duration of that `NotifyObservers` call — i.e. only while `Observe()` is
	running.  Moreover, the concrete event subclasses hold their payload by
	`const&` into the emitter's internal buffers (e.g. `resulting_point_`,
	`start_point_`), which may themselves be reused on the next step.

	Therefore an observer **must not** store a reference or pointer to an event
	(or to anything it returns by reference) past the return of `Observe()`.
	If you need to keep data, copy it out by value inside `Observe()`.  This is
	exactly what the Python observer layer does: it reads values during the call
	and copies them into its own storage.
	*/
	class AnyEvent
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~AnyEvent() = default;
	};

	template<class ObsT, bool IsConst = true>
	class Event;

	/**
	\brief For emission of events from observables.
	
	An observable object probably wants to emit events to notify observers that things are happening.  The observers filter the events, at this time by dynamic casting (slow, I know, but it works for now.  Premature optimization, etc.)  A better method than dynamic casting would be to filter based on a union of bits or something.  If you, the reader, want to help make this better, please contact a developer!

	Say I am an observable object, and I want to emit an event.  Events attach the type of object emitting them, and in fact (a refence to) the emitter itself.  So if my type is `T`, I would do something like `NotifyObservers(Event<T>(*this))`.  Then an Observer can filter based on a heirarchy of event types, etc.  

	\tparam ObsT The Observed type.  When emitting an event, you pass in the type of object emitting the event, and the object itself.  Then the observer can `Get` the emitting object, and do (const) stuff to it.

	\see Observable, Observable::NotifyObservers, ADD_BERTINI_EVENT_TYPE, AMPPathAccumulator, AnyEvent
	*/
	template<class ObsT>
	class Event<ObsT, true> : public AnyEvent
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:

		/**
		\brief Constructor for an event.  

		\param obs The observable emitting the event to observers.  Passed by const reference, this permits arbitrary `const` function calls by the observer.
		*/
		Event(ObsT const& obs) : current_observable_(obs)
		{}

		virtual ~Event() = default;

		/**
		\brief Get the emitting object, by `const` reference.

		\return The observable who emitted the event.  This permits calls of arbitrary const functions, particularly getters.

		\see AMPPathAccumulator for simple example of filtering using `dynamic_cast`s
		*/
		ObsT const& Get() const
		{return current_observable_;}

		Event() = delete;

		using HeldT = const ObsT&;  ///< How the observable is held: by const reference.
	protected:
		const ObsT& current_observable_;  ///< The observable object this event carries.
	};


	/**
	\brief For emission of events from observables.
	
	An observable object probably wants to emit events to notify observers that things are happening.  The observers filter the events, at this time by dynamic casting (slow, I know, but it works for now.  Premature optimization, etc.)  A better method than dynamic casting would be to filter based on a union of bits or something.  If you, the reader, want to help make this better, please contact a developer!

	Say I am an observable object, and I want to emit an event.  Events attach the type of object emitting them, and in fact (a refence to) the emitter itself.  So if my type is `T`, I would do something like `NotifyObservers(Event<T>(*this))`.  Then an Observer can filter based on a heirarchy of event types, etc.  

	\tparam ObsT The Observed type.  When emitting an event, you pass in the type of object emitting the event, and the object itself.  Then the observer can `Get` the emitting object, and do (const) stuff to it.

	\see Observable, Observable::NotifyObservers, ADD_BERTINI_EVENT_TYPE, AMPPathAccumulator, AnyEvent
	*/
	template<class ObsT>
	class Event<ObsT,false> : public AnyEvent
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:

		/**
		\brief Constructor for an event.  

		\param obs The observable emitting the event to observers.  Passed by const reference, this permits arbitrary `const` function calls by the observer.
		*/
		Event(ObsT & obs) : current_observable_(obs)
		{}

		virtual ~Event() = default;

		/**
		\brief Get the emitting object, by `const` reference.

		\return The observable who emitted the event.  This permits calls of arbitrary const functions, particularly getters.

		\see AMPPathAccumulator for simple example of filtering using `dynamic_cast`s
		*/
		ObsT & Get()
		{return current_observable_;}

		Event() = delete;

		using HeldT = ObsT&;  ///< How the observable is held: by mutable reference.
	protected:
		ObsT& current_observable_;  ///< The observable object this event carries.
	};

	/// \brief An event holding a const reference to its observable.
	template<typename ObsT>
	using ConstEvent = Event<ObsT,true>;

	/// \brief An event holding a mutable reference to its observable.
	template<typename ObsT>
	using MutableEvent = Event<ObsT,false>;

	/**
	\brief Defines a new event type in a hierarchy

	\param event_name The name of the new Event type you are making.
	\param event_parenttype The name of the parent Event type in the heirarchy.
	*/
	#define ADD_BERTINI_EVENT_TYPE(event_name,event_parenttype) template<class ObservedT> \
	class event_name : public event_parenttype<ObservedT> \
	{ BOOST_TYPE_INDEX_REGISTER_CLASS \
	public: \
		using HeldT = typename event_parenttype<ObservedT>::HeldT; \
		event_name(HeldT obs) : event_parenttype<ObservedT>(obs){} \
		virtual ~event_name() = default; \
		event_name() = delete; }
	
} //re: namespace bertini

#endif
