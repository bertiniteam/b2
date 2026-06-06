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

RunManagerLoop<TaskT, ResultT> seeds each worker with a task from the queue,
then collects results and dispatches new tasks until the queue is empty.
Workers are told to stop by receiving a sentinel task.

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

#include <mpi.h>

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

	std::map<int, bool> in_flight;  // rank -> currently tracking

	// Seed workers
	for (int rank = 1; rank <= num_workers; ++rank)
	{
		if (!tasks.empty())
		{
			TaskT task = tasks.front(); tasks.pop();
			mpi_send_serialized(comm, rank, TAG_WORK_ITEM, task);
			in_flight[rank] = true;
		}
		else
		{
			mpi_send_serialized(comm, rank, TAG_WORK_ITEM, sentinel_value);
		}
	}

	// Drain
	while (!in_flight.empty())
	{
		ResultT result;
		int worker = mpi_recv_serialized_any(comm, TAG_RESULT, result);

		store_fn(result);
		in_flight.erase(worker);

		if (!tasks.empty())
		{
			TaskT task = tasks.front(); tasks.pop();
			mpi_send_serialized(comm, worker, TAG_WORK_ITEM, task);
			in_flight[worker] = true;
		}
		else
		{
			mpi_send_serialized(comm, worker, TAG_WORK_ITEM, sentinel_value);
		}
	}
}

} // namespace parallel
} // namespace bertini

#endif // BERTINI2_HAVE_MPI
