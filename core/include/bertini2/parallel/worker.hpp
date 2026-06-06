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

\brief Generic worker loop for dynamic-assignment manager-worker MPI parallelism.

RunWorkerLoop<TaskT, ResultT> receives tasks from rank 0, calls track_fn with
the full task object, then calls pack_fn to assemble and send the result back.
Exits when a sentinel task is received.

Must be paired with RunManagerLoop<TaskT, ResultT>() running on rank 0.
*/

#pragma once

#ifdef BERTINI2_HAVE_MPI

#include "bertini2/parallel/path_result.hpp"
#include "bertini2/parallel/mpi_utils.hpp"

#include <mpi.h>

#include <functional>
#include <limits>

namespace bertini {
namespace parallel {


/**
\brief Run the worker side of a single tracking phase.

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

} // namespace parallel
} // namespace bertini

#endif // BERTINI2_HAVE_MPI
