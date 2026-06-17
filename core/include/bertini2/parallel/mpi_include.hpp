//This file is part of Bertini 2.
//
//bertini2/parallel/mpi_include.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/mpi_include.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/mpi_include.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file bertini2/parallel/mpi_include.hpp

\brief Single include point for the MPI C API.

Bertini 2 uses only the MPI **C** API.  Including `<mpi.h>` directly from C++ otherwise drags in the
deprecated MPI **C++** bindings (the `MPI::` namespace, `mpicxx.h`).  MPICH compiles those in by
default, and they fail to build under modern C++ (notably clang) -- so the project built against
OpenMPI but not Homebrew MPICH on macOS.

Defining `MPICH_SKIP_MPICXX` / `OMPI_SKIP_MPICXX` before `<mpi.h>` pulls in only the C API.  These are
also defined on the compile command line by CMake (`ENABLE_MPI` block); having them here too is belt
and suspenders, so the guard holds even if the build-system definition fails to reach a translation
unit.  Always include this header instead of `<mpi.h>` directly.
*/

#pragma once

#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

#include <mpi.h>
