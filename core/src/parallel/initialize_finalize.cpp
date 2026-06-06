//This file is part of Bertini 2.
//
//src/parallel/initialize_finalize.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/parallel/initialize_finalize.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/parallel/initialize_finalize.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file src/parallel/initialize_finalize.cpp 

\brief Provides the free initialization and finalization routines for Bertini2.
*/



#include "bertini2/parallel/initialize_finalize.hpp"


namespace bertini{

	namespace parallel{

#ifdef BERTINI2_HAVE_MPI
	MPI_Comm WorldComm()
	{
		return MPI_COMM_WORLD;
	}
#endif

		void Initialize()
		{
#ifdef BERTINI2_HAVE_MPI
			// Check first so mpi4py (or another library) that already called MPI_Init
			// is handled gracefully — we skip the init in that case.
			int already_initialized = 0;
			MPI_Initialized(&already_initialized);
			if (!already_initialized)
				MPI_Init(nullptr, nullptr);
#endif
		}

		void Finalize()
		{
#ifdef BERTINI2_HAVE_MPI
			int already_finalized = 0;
			MPI_Finalized(&already_finalized);
			if (!already_finalized)
				MPI_Finalize();
#endif
		}

		int Rank()
		{
#ifdef BERTINI2_HAVE_MPI
			int initialized = 0;
			MPI_Initialized(&initialized);
			if (!initialized) return 0;
			int r = 0;
			MPI_Comm_rank(MPI_COMM_WORLD, &r);
			return r;
#else
			return 0;
#endif
		}

		int Size()
		{
#ifdef BERTINI2_HAVE_MPI
			int initialized = 0;
			MPI_Initialized(&initialized);
			if (!initialized) return 1;
			int s = 1;
			MPI_Comm_size(MPI_COMM_WORLD, &s);
			return s;
#else
			return 1;
#endif
		}

		bool IsManager() { return Rank() == 0; }
		bool IsWorker()  { return Rank() != 0; }

	} // namespace parallel

	namespace serial{
		void Initialize()
		{
			// Splash screen only on rank 0 — guard handles both serial and parallel builds.
			if (parallel::IsManager())
			{
				std::cout << "\n\n\n" << SplashScreen() << "\n\n\n";
				std::cout << "\n\n" << DependencyVersions() << "\n\n";
			}

			logging::Logging::Init();
		}

		void Finalize()
		{

		}
	}

}






