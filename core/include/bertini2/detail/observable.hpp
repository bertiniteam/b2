//This file is part of Bertini 2.
//
//bertini2/detail/observable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/detail/observable.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/detail/observable.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
\file bertini2/detail/observable.hpp

\brief Contains the observable base types

\defgroup observable
*/

#ifndef BERTINI_DETAIL_OBSERVABLE_HPP
#define BERTINI_DETAIL_OBSERVABLE_HPP

#include "bertini2/detail/observer.hpp"
#include "bertini2/detail/events.hpp"

#include <typeindex>
#include <unordered_map>

namespace bertini{

	/**
	\brief An abstract observable type, maintaining a list of observers, who can be notified in case of Events.
	
	Some known observable types are Tracker and Endgame.
	*/
	class Observable
	{	
	public:

		virtual ~Observable() = default;


		/**
		\brief Add an observer, to observe this observable.

		Observers that override SubscribedEventTypes() with a non-empty list are
		registered in a type-indexed map so NotifyObservers only calls them for
		events they declared interest in.  Observers returning an empty list (the
		default) are placed in a catch-all list and receive every event.
		*/
		void AddObserver(AnyObserver& new_observer) const
		{
			auto types = new_observer.SubscribedEventTypes();
			if (types.empty())
			{
				if (find_if(begin(untyped_watchers_), end(untyped_watchers_), [&](const auto& held_obs)
				            { return &held_obs.get() == &new_observer; }) == end(untyped_watchers_))
					untyped_watchers_.push_back(std::ref(new_observer));
			}
			else
			{
				for (auto& ti : types)
				{
					auto& bucket = typed_watchers_[ti];
					if (find_if(begin(bucket), end(bucket), [&](const auto& held_obs)
					            { return &held_obs.get() == &new_observer; }) == end(bucket))
						bucket.push_back(std::ref(new_observer));
				}
			}
		}

		/**
		\brief Remove an observer from this observable.
		*/
		void RemoveObserver(AnyObserver& observer) const
		{
			auto erase_from = [&](auto& container) {
				auto new_end = std::remove_if(container.begin(), container.end(),
				                              [&](const auto& held_obs)
				                              { return &held_obs.get() == &observer; });
				container.erase(new_end, container.end());
			};

			erase_from(untyped_watchers_);
			for (auto& [ti, bucket] : typed_watchers_)
				erase_from(bucket);
		}

	protected:

		/**
		\brief Sends an event to observers that subscribed to its exact dynamic type,
		then to all catch-all (untyped) observers.

		Snapshots each observer list before iterating so that an observer calling
		RemoveObserver(*this) inside Observe() does not invalidate the loop iterator.
		*/
		void NotifyObservers(AnyEvent const& e) const
		{
			auto it = typed_watchers_.find(std::type_index(typeid(e)));
			if (it != typed_watchers_.end())
			{
				ObserverList snapshot = it->second;
				for (auto& obs : snapshot)
					obs.get().Observe(e);
			}
			ObserverList untyped_snapshot = untyped_watchers_;
			for (auto& obs : untyped_snapshot)
				obs.get().Observe(e);
		}

		void NotifyObservers(AnyEvent& e) const
		{
			auto it = typed_watchers_.find(std::type_index(typeid(e)));
			if (it != typed_watchers_.end())
			{
				ObserverList snapshot = it->second;
				for (auto& obs : snapshot)
					obs.get().Observe(e);
			}
			ObserverList untyped_snapshot = untyped_watchers_;
			for (auto& obs : untyped_snapshot)
				obs.get().Observe(e);
		}

	private:

		using ObserverList = std::vector<std::reference_wrapper<AnyObserver>>;

		mutable std::unordered_map<std::type_index, ObserverList> typed_watchers_;
		mutable ObserverList untyped_watchers_;
	};

} // namespace bertini


#endif

