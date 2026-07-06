//This file is part of Bertini 2.
//
//endgames_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//nag_algorithms_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with nag_algorithms_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst
//

/**
\file test/nag_algorithms/nag_algorithms_test.cpp:  main source file for the algorithm testing executable for Bertini2
*/




#define BOOST_ALL_DYN_LINK 1

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 NAG Algorithm Testing"
#include <boost/test/unit_test.hpp>

#define BERTINI_TEST_MODULE "nag_algorithms"
#include "bertini2/mpfr_extensions.hpp"
#include "test/utility/enable_logging.hpp"

#include <cstdlib>

// Records are ON BY DEFAULT for every solver (ambient ./bertini_output when
// BERTINI_RECORDS_DIR is unset).  Tests must be hermetic: without this, a bare
// solve would litter the test runner's cwd AND recall paths recorded by a
// PREVIOUS ctest run, silently changing what a rerun actually exercises.  The
// value `none` is the portable off switch (an empty value is POSIX-only: Windows
// deletes a variable assigned an empty string, which would flip "off" back to
// the default); records tests attach their own directories via RecordTo, which
// always wins over the ambient resolution.
struct RecordsOffByDefaultInTests
{
	/// Export the explicit records-off sentinel for the whole test module.
	RecordsOffByDefaultInTests()
	{
#ifdef _WIN32
		_putenv_s("BERTINI_RECORDS_DIR", "none");
#else
		setenv("BERTINI_RECORDS_DIR", "none", 1);
#endif
	}
};

BOOST_GLOBAL_FIXTURE( RecordsOffByDefaultInTests );

// deliberately left blank.  link other files with this one.

