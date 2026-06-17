//This file is part of Bertini 2.
//
//bertini2/parallel/manager.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/manager.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/manager.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file bertini2/parallel/manager.hpp

\brief Generic manager (rank 0) loop for dynamic-assignment manager-worker MPI parallelism.

RunManagerLoop<TaskT, ResultT> uses credit-based dispatch: each worker announces
its capacity (number of tracking threads) at loop start via TAG_CAPACITY; the
manager keeps up to that many tasks outstanding per rank, dispatching a new task
each time a result arrives. Workers are told to stop by receiving a sentinel task,
sent only once that rank has zero tasks outstanding.

Requirements:
  - TaskT must be Boost.MPI-serializable
  - ResultT must be Boost.MPI-serializable
  - sentinel must be a TaskT value that workers can detect via is_sentinel()
    (or be a dedicated out-of-band value)

For Phase 1 (before-EG), TaskT = std::size_t (path index) with sentinel =
std::numeric_limits<std::size_t>::max().

For Phase 2 (during-EG), TaskT = Phase2Task<ComplexT> with sentinel constructed
via Phase2Task<ComplexT>::sentinel().
*/

#pragma once

#ifdef BERTINI2_HAVE_MPI

#include "bertini2/parallel/path_result.hpp"
#include "bertini2/parallel/mpi_utils.hpp"

#include "bertini2/parallel/mpi_include.hpp"

#include <functional>
#include <limits>
#include <map>
#include <queue>

namespace bertini {
namespace parallel {


/**
\brief Run the manager side of a single tracking phase.

\param comm     MPI communicator (WorldComm()).
\param tasks    Queue of TaskT items to distribute; consumed by this call.
\param store_fn Called on rank 0 with each received ResultT.

All worker ranks must simultaneously be running RunWorkerLoop<TaskT, ResultT>().
*/
template<typename TaskT, typename ResultT>
void RunManagerLoop(
	MPI_Comm comm,
	std::queue<TaskT>& tasks,
	std::function<void(ResultT const&)> store_fn)
{
	int world_size = 0;
	MPI_Comm_size(comm, &world_size);
	const int num_workers = world_size - 1;
	if (num_workers <= 0)
		return;

	TaskT sentinel_value = detail::make_sentinel(TaskT{});

	// Each worker announces how many tasks it can have in flight at once
	// (1 for a serial worker, n_threads for a threaded one).
	std::map<int, int> capacity;     // rank -> max tasks in flight
	for (int ii = 0; ii < num_workers; ++ii)
	{
		int cap = 0;
		MPI_Status status;
		MPI_Recv(&cap, 1, MPI_INT, MPI_ANY_SOURCE, TAG_CAPACITY, comm, &status);
		capacity[status.MPI_SOURCE] = (cap >= 1) ? cap : 1;
	}

	std::map<int, int> outstanding;  // rank -> tasks dispatched but not yet returned
	int active_workers = num_workers;

	// Seed each worker up to its capacity.  A rank seeded with zero tasks
	// (queue exhausted) is immediately sent the sentinel.
	for (auto const& [rank, cap] : capacity)
	{
		outstanding[rank] = 0;
		for (int ii = 0; ii < cap && !tasks.empty(); ++ii)
		{
			TaskT task = tasks.front(); tasks.pop();
			mpi_send_serialized(comm, rank, TAG_WORK_ITEM, task);
			++outstanding[rank];
		}
		if (outstanding[rank] == 0)
		{
			mpi_send_serialized(comm, rank, TAG_WORK_ITEM, sentinel_value);
			--active_workers;
		}
	}

	// Drain: refill a rank as each result arrives; retire it with a sentinel
	// once the queue is empty and its last outstanding task has come back.
	while (active_workers > 0)
	{
		ResultT result;
		int worker = mpi_recv_serialized_any(comm, TAG_RESULT, result);

		store_fn(result);
		--outstanding[worker];

		if (!tasks.empty())
		{
			TaskT task = tasks.front(); tasks.pop();
			mpi_send_serialized(comm, worker, TAG_WORK_ITEM, task);
			++outstanding[worker];
		}
		else if (outstanding[worker] == 0)
		{
			mpi_send_serialized(comm, worker, TAG_WORK_ITEM, sentinel_value);
			--active_workers;
		}
	}
}

} // namespace parallel
} // namespace bertini

#endif // BERTINI2_HAVE_MPI
