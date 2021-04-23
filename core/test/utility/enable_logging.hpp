//This file is part of Bertini 2.
//
//test/utility/enable_logging.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/utility/enable_logging.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/utility/enable_logging.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire
// spring 2018

/**
\file test/utility/enable_logging.hpp 

\brief Logging in Bertini using Boost.Log
*/


#pragma once

#include "logging.hpp"

struct LogInitter
{
	LogInitter()
	{
		bertini::logging::Logging::Init("bertini2_tests_" + std::string(BERTINI_TEST_MODULE) + "_%N.log","%Message%",10*1024*1024, bertini::logging::severity_level::trace);
	}
};


BOOST_TEST_GLOBAL_FIXTURE( LogInitter );
