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

\defgroup observable Observables
*/

#ifndef BERTINI_DETAIL_OBSERVABLE_HPP
#define BERTINI_DETAIL_OBSERVABLE_HPP

#include "bertini2/detail/observer.hpp"
#include "bertini2/detail/events.hpp"

#include <algorithm>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/core/demangle.hpp>

namespace bertini{

	/**
	\brief Thrown when an observer is attached to an observable it cannot observe.

	The Python bindings translate this to a TypeError.
	*/
	struct IncompatibleObserver : public std::runtime_error
	{
		using std::runtime_error::runtime_error;
	};

	/**
	\brief An abstract observable type, maintaining a list of observers, who can be notified in case of Events.

	Some known observable types are Tracker and Endgame.
	*/
	class Observable
	{
	public:

		virtual ~Observable() = default;

		Observable() = default;

		/**
		Copies do NOT inherit the watcher list: observers subscribed to the source
		object, not to the copy.  This makes copies safe to use on other threads —
		a fresh copy notifies nobody until observers are explicitly added to it.
		Copy-assignment likewise leaves the target's own watchers untouched.
		*/
		Observable(Observable const&) {}
		/// \brief Copy-assignment: leaves the target's own watcher list untouched (see the copy ctor note).
		Observable& operator=(Observable const&) { return *this; }


		/**
		\brief Add an observer (non-owning) to observe this observable.

		Stores only a raw reference: the caller is responsible for keeping the
		observer alive while it is attached.  This is the right overload for
		stack-allocated C++ observers added and removed within a scope.  For an
		observer whose lifetime you would rather have tracked, see the
		shared_ptr overload.

		Observers that override SubscribedEventTypes() with a non-empty list are
		registered in a type-indexed map so NotifyObservers only calls them for
		events they declared interest in.  Observers returning an empty list (the
		default) are placed in a catch-all list and receive every event.

		Safe to call from inside an observer's Observe() (i.e. during a
		notification): the add is deferred and applied once the current
		notification loop finishes.  This is what lets a "meta-observer" attach
		new observers in response to events.

		\throws IncompatibleObserver if the observer's ObservedKind() is not one
		        this observable can be observed as.
		*/
		void AddObserver(AnyObserver& new_observer) const
		{
			std::lock_guard<std::recursive_mutex> lock(notify_mutex_);
			RejectIfIncompatible(new_observer);
			Watcher w;
			w.raw = &new_observer;        // non-owning: owned stays null
			EnqueueOrApplyAdd(std::move(w));
		}

		/**
		\brief Add an observer whose lifetime is co-owned by this observable.

		The observable keeps a shared_ptr, so the observer stays alive (and keeps
		receiving events) for as long as it is attached, even after the caller
		drops their own reference — "attach it and forget it", with no dangling.
		It is released when RemoveObserver is called or this observable is
		destroyed.  (Stack-allocated C++ observers can't be shared_ptr-owned; use
		the reference overload for those.)

		\throws IncompatibleObserver if incompatible (see the reference overload).
		*/
		void AddObserver(std::shared_ptr<AnyObserver> const& new_observer) const
		{
			std::lock_guard<std::recursive_mutex> lock(notify_mutex_);
			RejectIfIncompatible(*new_observer);
			Watcher w;
			w.owned = new_observer;          // owning: keeps the observer alive
			w.raw   = new_observer.get();
			EnqueueOrApplyAdd(std::move(w));
		}

		/**
		\brief Remove an observer from this observable.

		Like AddObserver, this is safe to call during a notification: the removal
		is deferred until the current notification loop finishes.
		*/
		void RemoveObserver(AnyObserver& observer) const
		{
			std::lock_guard<std::recursive_mutex> lock(notify_mutex_);
			if (dispatch_depth_ > 0)
			{
				PendingOp op;
				op.kind = PendingOp::Remove;
				op.watcher.raw = &observer;
				pending_ops_.push_back(std::move(op));
			}
			else
				RemoveWatcher(&observer);
		}

		/**
		\brief Whether this observable may be observed as the given type.

		Drives the compatibility check in AddObserver.  The default accepts the
		wildcard `typeid(void)` and this observable's own dynamic type.  An
		observable that emits its events templated on a base type (so a single
		observer type serves many concrete observables) should also accept that
		base type — e.g. ZeroDim accepts AnyZeroDim.
		*/
		virtual bool ObservableIsA(std::type_index t) const
		{
			return t == std::type_index(typeid(void))
			    || t == std::type_index(typeid(*this));
		}

	protected:

		/**
		\brief Sends an event to observers that subscribed to its exact dynamic type,
		then to all catch-all (untyped) observers.

		Mutations of the observer lists are deferred for the duration of a
		notification (see AddObserver/RemoveObserver), so iterating the live lists
		here is safe even when an observer subscribes/unsubscribes mid-dispatch.
		An observer that returns ObserveResult::Unsubscribe is dropped once the
		(possibly nested) notification completes.  The depth counter keeps this
		correct under re-entrant emission.
		*/
		void NotifyObservers(AnyEvent const& e) const
		{
			DispatchEvent(e);
		}

		/// \brief Notify all observers of a (mutable) event.
		void NotifyObservers(AnyEvent& e) const
		{
			DispatchEvent(e);
		}

	private:

		/**
		One entry in an observer list.  Either non-owning (`owned` is null; the
		caller guarantees the observer's lifetime) or owning (`owned` co-owns the
		observer, keeping it alive while attached).  `raw` is always the identity
		used for dedup and removal, and the pointer dispatched through.
		*/
		struct Watcher
		{
			std::shared_ptr<AnyObserver> owned;   // non-null iff owning
			AnyObserver* raw = nullptr;           // identity + the pointer events go to
		};

		using ObserverList = std::vector<Watcher>;

		struct PendingOp
		{
			enum Kind { Add, Remove } kind;
			Watcher watcher;   // Add: the watcher (carries `owned` for owning adds); Remove: identity in .raw
		};

		void RejectIfIncompatible(AnyObserver const& obs) const
		{
			if (!ObservableIsA(obs.ObservedKind()))
				throw IncompatibleObserver(
					"this observable (" + boost::core::demangle(typeid(*this).name())
					+ ") does not accept an observer for "
					+ boost::core::demangle(obs.ObservedKind().name()));
		}

		void EnqueueOrApplyAdd(Watcher w) const
		{
			if (dispatch_depth_ > 0)
			{
				PendingOp op;
				op.kind = PendingOp::Add;
				op.watcher = std::move(w);   // carries `owned`, so the observer stays alive until drained
				pending_ops_.push_back(std::move(op));
			}
			else
				AddWatcher(std::move(w));
		}

		void AddWatcher(Watcher const& w) const
		{
			auto present_in = [&](ObserverList const& c) {
				return std::find_if(c.begin(), c.end(),
				                    [&](Watcher const& x){ return x.raw == w.raw; }) != c.end();
			};

			auto types = w.raw->SubscribedEventTypes();
			if (types.empty())
			{
				if (!present_in(untyped_watchers_))
					untyped_watchers_.push_back(w);
			}
			else
			{
				for (auto& ti : types)
				{
					auto& bucket = typed_watchers_[ti];
					if (std::find_if(bucket.begin(), bucket.end(),
					                 [&](Watcher const& x){ return x.raw == w.raw; }) == bucket.end())
						bucket.push_back(w);
				}
			}
		}

		void RemoveWatcher(AnyObserver* who) const
		{
			auto erase_from = [&](ObserverList& c) {
				c.erase(std::remove_if(c.begin(), c.end(),
				                       [&](Watcher const& x){ return x.raw == who; }),
				        c.end());
			};
			erase_from(untyped_watchers_);
			for (auto& [ti, bucket] : typed_watchers_)
				erase_from(bucket);
		}

		void DrainPendingOps() const
		{
			// Apply in request order so that, within one notification, a remove
			// followed by a re-add (or vice versa) lands on the intended state.
			for (auto& op : pending_ops_)
			{
				if (op.kind == PendingOp::Add)
					AddWatcher(op.watcher);
				else
					RemoveWatcher(op.watcher.raw);
			}
			pending_ops_.clear();
		}

		void DispatchEvent(AnyEvent const& e) const
		{
			std::lock_guard<std::recursive_mutex> lock(notify_mutex_);
			++dispatch_depth_;

			auto run = [&](ObserverList& bucket) {
				for (auto& w : bucket)
				{
					// w.raw stays valid for the whole call: owning entries hold a
					// shared_ptr, non-owning ones are the caller's responsibility.
					if (w.raw->Observe(e) == ObserveResult::Unsubscribe)
					{
						PendingOp op;
						op.kind = PendingOp::Remove;
						op.watcher.raw = w.raw;
						pending_ops_.push_back(std::move(op));
					}
				}
			};

			auto it = typed_watchers_.find(std::type_index(typeid(e)));
			if (it != typed_watchers_.end())
				run(it->second);
			run(untyped_watchers_);

			--dispatch_depth_;
			if (dispatch_depth_ == 0)
				DrainPendingOps();
		}

		mutable std::unordered_map<std::type_index, ObserverList> typed_watchers_;
		mutable ObserverList untyped_watchers_;

		mutable unsigned dispatch_depth_ = 0;
		mutable std::vector<PendingOp> pending_ops_;

		// Serializes notification and watcher-list mutation so the same observable can be safely
		// notified from several worker threads at once (the MPI-less threaded solve fires per-path
		// events from each tracking thread onto one shared solver).  Recursive because an observer's
		// Observe() may itself Add/RemoveObserver on this same observable during a dispatch (the
		// "meta-observer" pattern) -- that re-entrant call re-locks on the same thread and is routed
		// through the deferred pending_ops_ path via dispatch_depth_.  Held across Observe() calls,
		// so observer callbacks are serialized; the heavy numerical work runs outside this lock.
		// Single-threaded callers pay one uncontended recursive lock per event -- negligible.
		mutable std::recursive_mutex notify_mutex_;
	};

} // namespace bertini


#endif

