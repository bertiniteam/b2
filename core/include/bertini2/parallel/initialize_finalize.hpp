//This file is part of Bertini 2.
//
//bertini2/parallel/initialize_finalize.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/initialize_finalize.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/initialize_finalize.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/parallel/initialize_finalize.hpp 

\brief Provides the free initialization and finalization routines for Bertini2.
*/


#pragma once

#include <iostream>

#include "bertini2/logging.hpp"
#include "bertini2/io/splash.hpp"

#ifdef BERTINI2_HAVE_MPI
#include "bertini2/parallel/mpi_include.hpp"
#endif

namespace bertini{

	namespace parallel{

	/// \brief Finalize the parallel (MPI) environment.
	void Finalize();
	/// \brief Initialize the parallel (MPI) environment.
	void Initialize();

#ifdef BERTINI2_HAVE_MPI
	MPI_Comm WorldComm();
#endif

	// These work in both serial and parallel builds.
	// In serial builds (no MPI) they return constants matching a single-rank world.
	/// \brief Get this process's rank (0 in serial builds).
	int  Rank();
	/// \brief Get the number of ranks (1 in serial builds).
	int  Size();
	/// \brief Query whether this rank is the manager (rank 0).
	bool IsManager();   // true iff Rank() == 0
	/// \brief Query whether this rank is a worker (rank != 0).
	bool IsWorker();    // true iff Rank() != 0

	}

	namespace serial{
		/// \brief Finalize the serial environment (no-op counterpart to the parallel version).
		void Finalize();
		/// \brief Initialize the serial environment (no-op counterpart to the parallel version).
		void Initialize();
	}
}
