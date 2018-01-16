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
\file bertini2/detail/observable.hpp

\brief Contains the observable base types
*/

#ifndef BERTINI_DETAIL_OBSERVABLE_HPP
#define BERTINI_DETAIL_OBSERVABLE_HPP

#include "bertini2/detail/observer.hpp"
#include "bertini2/detail/events.hpp"

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
		*/
		void AddObserver(AnyObserver& new_observer) const
		{
			if (find_if(begin(current_watchers_), end(current_watchers_), [&](const auto& held_obs)
			                              { return &held_obs.get() == &new_observer; })==end(current_watchers_))
				
				current_watchers_.push_back(std::ref(new_observer));
		}

		/**
		\brief Remove an observer from this observable.
		*/
		void RemoveObserver(AnyObserver& observer) const
		{

			auto new_end = std::remove_if(current_watchers_.begin(), current_watchers_.end(),
			                              [&](const auto& held_obs)
			                              { return &held_obs.get() == &observer; });

			current_watchers_.erase(new_end, current_watchers_.end());


			// current_watchers_.erase(std::remove(current_watchers_.begin(), current_watchers_.end(), std::ref(observer)), current_watchers_.end());
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
				obs.get().Observe(e);

		}

		void NotifyObservers(AnyEvent & e) const
		{

			for (auto& obs : current_watchers_)
				obs.get().Observe(e);
		}


	private:

		using ObserverContainer = std::vector<std::reference_wrapper<AnyObserver>>;

		mutable ObserverContainer current_watchers_;
	};

} // namespace bertini


#endif

