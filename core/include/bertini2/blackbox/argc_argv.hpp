//This file is part of Bertini 2.
//
//bertini2/blackbox/argc_argv.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/blackbox/argc_argv.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/blackbox/argc_argv.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/blackbox/argc_argv.hpp 

\brief Provides the methods for parsing the command-line arguments.
*/

#pragma once

#include <boost/filesystem.hpp>

namespace bertini{

/**
\brief Parsed results of the command-line arguments.
*/
struct ParsedArgs {
	boost::filesystem::path input_file{"input"}; // Bertini1 default filename
};

/**
Parse command-line arguments and return a ParsedArgs struct.

Recognized flags:
  -f `<file>`   use `<file>` as the input file (overrides default "input")
  `<file>`      bare positional argument treated as input file path
  --help      print usage and exit
  --version   print version and exit

\param argc The number of arguments.
\param argv The argument strings.
*/
ParsedArgs ParseArgcArgv(int argc, char** argv);

} //namespace bertini
