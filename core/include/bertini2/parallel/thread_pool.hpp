//This file is part of Bertini 2.
//
//bertini2/parallel/thread_pool.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/thread_pool.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/thread_pool.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file bertini2/parallel/thread_pool.hpp

\brief Thread-safe queue and worker thread pool for MPI+thread hybrid parallelism.

ThreadSafeQueue<T>: blocking FIFO backed by std::deque + std::mutex + std::condition_variable.

WorkerThreadPool<TaskT, ResultT, StateFactory, TrackFn>: spawns N tracking threads. Each
thread calls StateFactory() once at startup to build its own self-contained tracking state
(typically a struct holding a System copy and a Tracker copy that references that System).
The MPI communication thread submits tasks via submit() and collects results via collect().
All MPI calls remain on the main thread (compatible with MPI_THREAD_FUNNELED).

The StateFactory-based design avoids the "tracker holds a reference to an external System"
problem: the factory returns a heap-owned state struct (unique_ptr), so the System copy and
the Tracker copy that references it live together — at a stable address — for the lifetime
of the thread.

Example:
  struct TrackState {
      System     sys;
      AMPTracker tracker;
  };
  auto factory = [&]() {
      auto s = std::unique_ptr<TrackState>(new TrackState{ base_system, base_tracker });
      s->tracker.SetSystem(s->sys);       // repoint the copy at its own System
      SetThreadPrecision(DefaultPrecision());
      return s;
  };
  auto track_fn = [&](std::unique_ptr<TrackState>& state, SolnIndT const& idx) -> BeforeResult {
      // ... track using state->tracker ...
  };
  WorkerThreadPool<SolnIndT, BeforeResult, decltype(factory), decltype(track_fn)>
      pool(n_threads, factory, track_fn);
*/

#pragma once

#ifdef BERTINI2_HAVE_MPI

#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>
#include <type_traits>
#include <variant>
#include <vector>

namespace bertini {
namespace parallel {


template<typename T>
class ThreadSafeQueue
{
	std::deque<T>           queue_;
	std::mutex              mutex_;
	std::condition_variable cv_;

public:
	void push(T item)
	{
		{
			std::lock_guard<std::mutex> lock(mutex_);
			queue_.push_back(std::move(item));
		}
		cv_.notify_one();
	}

	T pop()
	{
		std::unique_lock<std::mutex> lock(mutex_);
		cv_.wait(lock, [this]{ return !queue_.empty(); });
		T item = std::move(queue_.front());
		queue_.pop_front();
		return item;
	}

	std::optional<T> try_pop()
	{
		std::lock_guard<std::mutex> lock(mutex_);
		if (queue_.empty())
			return std::nullopt;
		T item = std::move(queue_.front());
		queue_.pop_front();
		return item;
	}
};


struct PoolShutdownSentinel {};


/**
\brief Thread pool for parallel path tracking within a single MPI worker rank.

\tparam TaskT         Task type received from the MPI manager.
\tparam ResultT       Result type sent back after tracking.
\tparam StateFactory  Callable () -> StateT: constructs a per-thread state struct
                      (e.g. containing a System copy + Tracker copy) once at thread startup.
\tparam TrackFn       Callable (StateT&, TaskT const&) -> ResultT: tracks one path using
                      the thread-local state.

The StateFactory approach sidesteps the "tracker holds a reference to an external System"
issue: each thread owns its state struct, so System and Tracker lifetimes are co-managed.
*/
template<typename TaskT, typename ResultT, typename StateFactory, typename TrackFn>
class WorkerThreadPool
{
	using StateT   = std::invoke_result_t<StateFactory>;
	using WorkItem = std::variant<TaskT, PoolShutdownSentinel>;

	ThreadSafeQueue<WorkItem>  task_queue_;
	ThreadSafeQueue<ResultT>   result_queue_;
	std::vector<std::thread>   threads_;
	int                        n_threads_;

public:
	WorkerThreadPool(int n_threads, StateFactory state_factory, TrackFn track_fn)
		: n_threads_(n_threads)
	{
		threads_.reserve(n_threads);
		for (int i = 0; i < n_threads; ++i)
		{
			threads_.emplace_back([this, state_factory, track_fn]() mutable
			{
				// Build thread-local state once (System copy + Tracker copy).
				StateT state = state_factory();

				while (true)
				{
					WorkItem item = task_queue_.pop();

					if (std::holds_alternative<PoolShutdownSentinel>(item))
						break;

					TaskT const& task = std::get<TaskT>(item);
					ResultT result = track_fn(state, task);
					result_queue_.push(std::move(result));
				}
			});
		}
	}

	void submit(TaskT task)
	{
		task_queue_.push(WorkItem{std::move(task)});
	}

	// Block until one result is available, then return it.
	ResultT collect()
	{
		return result_queue_.pop();
	}

	// Non-blocking: return a result if one is ready, otherwise std::nullopt.
	std::optional<ResultT> try_collect()
	{
		return result_queue_.try_pop();
	}

	// Signal all threads to exit and join. Call only after all submitted tasks
	// have been collected (result_queue_ empty).
	void shutdown()
	{
		for (int i = 0; i < n_threads_; ++i)
			task_queue_.push(WorkItem{PoolShutdownSentinel{}});
		for (auto& t : threads_)
			t.join();
		threads_.clear();
	}
};



} // namespace parallel
} // namespace bertini

#endif // BERTINI2_HAVE_MPI
