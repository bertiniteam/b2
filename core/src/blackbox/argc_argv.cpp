//This file is part of Bertini 2.
//
//src/blackbox/argc_argv.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/blackbox/argc_argv.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/blackbox/argc_argv.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file src/blackbox/argc_argv.cpp 

\brief Provides the methods for parsing the command-line arguments.
*/

#include "bertini2/blackbox/argc_argv.hpp"
#include "bertini2/io/splash.hpp"
#include <iostream>
#include <string>

namespace bertini{

ParsedArgs ParseArgcArgv(int argc, char** argv)
{
	ParsedArgs result;

	for (int i = 1; i < argc; ++i)
	{
		std::string arg{argv[i]};

		if (arg == "--help" || arg == "-h")
		{
			std::cout <<
				"Usage: bertini2 [options] [input_file]\n"
				"\n"
				"  input_file          path to Bertini classic input file (default: \"input\")\n"
				"  -f <file>           specify input file explicitly\n"
				"  --help, -h          print this message and exit\n"
				"  --version           print version and exit\n"
				"\n"
				"Parallelism:\n"
#ifdef BERTINI2_HAVE_MPI
				"  MPI ranks:    mpirun --bind-to none -n N bertini2 [input_file]\n"
				"  Threads/rank: OMP_NUM_THREADS=T mpirun --bind-to none -n N bertini2 [input_file]\n"
				"    Total capacity = N ranks x T threads. OMP_NUM_THREADS defaults to 1.\n"
				"    Note: threading uses std::thread, not OpenMP. OMP_NUM_THREADS is reused\n"
				"    as the thread-count variable because HPC schedulers (e.g. SLURM) set it\n"
				"    automatically from --cpus-per-task, so no extra configuration is needed.\n"
#else
				"  This build was compiled without MPI. Rebuild with MPI present for\n"
				"  multi-rank parallelism (it is auto-detected at configure time).\n"
				"  Thread count: set OMP_NUM_THREADS (default: 1). Uses std::thread, not OpenMP.\n"
#endif
				;
			std::exit(0);
		}
		else if (arg == "--version")
		{
			std::cout << "Bertini2 " << bertini::Version() << "\n";
			std::cout << bertini::DependencyVersions();
			std::exit(0);
		}
		else if (arg == "-f" && i + 1 < argc)
		{
			result.input_file = argv[++i];
		}
		else if (arg[0] != '-')
		{
			result.input_file = arg;
		}
		else
		{
			std::cerr << "warning: unrecognized argument '" << arg << "'\n";
		}
	}

	return result;
}

}

