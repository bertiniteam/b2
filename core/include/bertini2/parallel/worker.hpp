//This file is part of Bertini 2.
//
//bertini2/parallel/worker.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/worker.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/worker.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file bertini2/parallel/worker.hpp

\brief Worker-side loops for dynamic-assignment manager-worker MPI parallelism.

RunWorkerLoop<TaskT, ResultT>: single-threaded worker. Receives tasks from rank 0,
calls track_fn, sends results back. Exits on sentinel.

RunWorkerLoopThreaded<TaskT, ResultT>: multi-threaded worker. Spawns a WorkerThreadPool
of n_threads tracking threads. The MPI communication thread (main) receives tasks one at
a time (manager protocol unchanged) and submits them to the pool. Results are collected
from the pool and sent back as they complete.

Thread safety: all MPI calls remain on the main thread (MPI_THREAD_FUNNELED).
Each tracking thread owns its state (System copy + Tracker copy) built by state_factory.
*/

#pragma once

#ifdef BERTINI2_HAVE_MPI

#include "bertini2/parallel/path_result.hpp"
#include "bertini2/parallel/mpi_utils.hpp"
#include "bertini2/parallel/thread_pool.hpp"

#include <mpi.h>

#include <chrono>
#include <cstdlib>   // std::getenv
#include <functional>
#include <limits>
#include <string>
#include <thread>

namespace bertini {
namespace parallel {


/**
\brief Run the worker side of a single tracking phase (single-threaded).

\param comm      MPI communicator.
\param track_fn  Called with each TaskT; performs tracking (stores results locally).
\param pack_fn   Called with each TaskT after track_fn; returns a ResultT to send.
*/
template<typename TaskT, typename ResultT>
void RunWorkerLoop(
	MPI_Comm comm,
	std::function<void(TaskT const&)> track_fn,
	std::function<ResultT(TaskT const&)> pack_fn)
{
	// Announce capacity: a serial worker can have exactly one task in flight.
	int capacity = 1;
	MPI_Send(&capacity, 1, MPI_INT, 0, TAG_CAPACITY, comm);

	while (true)
	{
		TaskT task;
		mpi_recv_serialized(comm, 0, TAG_WORK_ITEM, task);

		if (detail::is_sentinel(task))
			break;

		track_fn(task);

		ResultT result = pack_fn(task);
		mpi_send_serialized(comm, 0, TAG_RESULT, result);
	}
}


/**
\brief Run the worker side of a single tracking phase with a per-thread tracker pool.

Spawns n_threads tracking threads, announces a capacity of n_threads to the manager
(which keeps that many tasks outstanding on this rank — see RunManagerLoop), then
runs an event loop on the main (MPI) thread:

  - forward any completed results from the pool to the manager,
  - MPI_Iprobe for an incoming task; if present, receive and submit it to the pool,
  - sleep briefly when there was nothing to do.

The event loop never blocks in MPI_Recv while results are pending — a blocking
receive here would deadlock: the manager only sends the next task after receiving
a result, and the result would be stuck in the local queue.

The manager sends the sentinel only when this rank has zero tasks outstanding, so
after the sentinel arrives there is nothing left to drain (the loop condition
covers the defensive case regardless).

All MPI calls occur on this (main) thread — compatible with MPI_THREAD_FUNNELED.

\param comm           MPI communicator.
\param state_factory  Callable () -> StateT: called once per thread at startup to build
                      a self-contained tracking state (e.g. {System copy, Tracker copy}).
\param track_fn       Callable (StateT&, TaskT const&) -> ResultT: tracks one path.
\param n_threads      Number of tracking threads. Must be >= 1.
*/
template<typename TaskT, typename ResultT, typename StateFactory, typename TrackFn>
void RunWorkerLoopThreaded(
	MPI_Comm comm,
	StateFactory state_factory,
	TrackFn track_fn,
	int n_threads)
{
	int capacity = n_threads;
	MPI_Send(&capacity, 1, MPI_INT, 0, TAG_CAPACITY, comm);

	WorkerThreadPool<TaskT, ResultT, StateFactory, TrackFn> pool(n_threads, state_factory, track_fn);

	int  in_flight    = 0;      // tasks submitted to the pool but not yet collected
	bool got_sentinel = false;

	while (!got_sentinel || in_flight > 0)
	{
		bool did_work = false;

		// Forward all completed results to the manager.
		while (auto opt = pool.try_collect())
		{
			mpi_send_serialized(comm, 0, TAG_RESULT, *opt);
			--in_flight;
			did_work = true;
		}

		// Non-blocking check for an incoming task (or the sentinel).
		if (!got_sentinel)
		{
			int flag = 0;
			MPI_Status status;
			MPI_Iprobe(0, TAG_WORK_ITEM, comm, &flag, &status);
			if (flag)
			{
				TaskT task;
				mpi_recv_serialized(comm, 0, TAG_WORK_ITEM, task);

				if (detail::is_sentinel(task))
					got_sentinel = true;
				else
				{
					pool.submit(std::move(task));
					++in_flight;
				}
				did_work = true;
			}
		}

		// Idle: nothing arrived and nothing finished.  Tracking a path takes
		// milliseconds to seconds, so a short sleep costs nothing measurable.
		if (!did_work)
			std::this_thread::sleep_for(std::chrono::microseconds(200));
	}

	pool.shutdown();
}


/**
\brief Read the thread count for a worker rank from OMP_NUM_THREADS.

Returns the value of OMP_NUM_THREADS if set and >= 1, otherwise 1 (serial).
HPC schedulers (SLURM) set OMP_NUM_THREADS automatically from --cpus-per-task.
*/
inline int WorkerThreadCount()
{
	const char* env = std::getenv("OMP_NUM_THREADS");
	if (env)
	{
		int n = std::atoi(env);
		if (n >= 1)
			return n;
	}
	return 1;
}


} // namespace parallel
} // namespace bertini

#endif // BERTINI2_HAVE_MPI
