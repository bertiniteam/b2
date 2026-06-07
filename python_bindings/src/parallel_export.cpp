//This file is part of Bertini 2.
//
//python_bindings/src/parallel_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python_bindings/src/parallel_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python_bindings/src/parallel_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

//  python_bindings/src/parallel_export.cpp:  source file for exposing MPI/parallel utilities to python.

#include "parallel_export.hpp"

namespace bertini{
	namespace python{

		void ExportParallel()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".parallel");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("parallel") = new_submodule;

			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "MPI parallelism support (requires MPI_Comm passed from mpi4py)";

			def("rank", &bertini::parallel::Rank,
				"Return the rank of this process. Returns 0 in serial (non-MPI) builds.");
			def("size", &bertini::parallel::Size,
				"Return the total number of processes. Returns 1 in serial (non-MPI) builds.");
			def("is_manager", &bertini::parallel::IsManager,
				"Return True if this process is the manager (rank 0). Always True in serial builds.");
			def("is_worker", &bertini::parallel::IsWorker,
				"Return True if this process is a worker (rank > 0). Always False in serial builds.");
			def("initialize", &bertini::parallel::Initialize,
				"Initialize MPI (or no-op if already initialized). Safe to call with mpi4py pre-init.");
			def("finalize", &bertini::parallel::Finalize,
				"Finalize MPI (or no-op if already finalized).");
		}
	}
}
