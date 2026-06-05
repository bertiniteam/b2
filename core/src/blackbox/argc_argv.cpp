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
				"Usage: bertini [options] [input_file]\n"
				"\n"
				"  input_file          path to Bertini classic input file (default: \"input\")\n"
				"  -f <file>           specify input file explicitly\n"
				"  --help, -h          print this message and exit\n"
				"  --version           print version and exit\n";
			std::exit(0);
		}
		else if (arg == "--version")
		{
			std::cout << "Bertini2\n";
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

