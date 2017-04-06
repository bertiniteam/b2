//This file is part of Bertini 2.
//
//test/blackbox/blackbox.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/blackbox/blackbox.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/blackbox/blackbox.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame


//test/blackbox/blackbox.cpp
//


 
#define BOOST_TEST_DYN_LINK 1

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Blackbox Testing"
#include <boost/test/unit_test.hpp>


#include "logging.hpp"


using sec_level = boost::log::trivial::severity_level;

using LoggingInit = bertini::LoggingInit;


#include "bertini2/blackbox/switches_zerodim.hpp"
#include "bertini2/system/precon.hpp"

BOOST_GLOBAL_FIXTURE( LoggingInit );

BOOST_AUTO_TEST_SUITE(blackbox_test)

BOOST_AUTO_TEST_SUITE(zero_dim)

using namespace bertini;

BOOST_AUTO_TEST_CASE(make_zero_dim)
{
	auto sys = system::Precon::GriewankOsborn();
	blackbox::ZeroDimRT my_runtime_type_options; // make defaults

	auto zd_ptr = blackbox::MakeZeroDim(my_runtime_type_options, sys);

	zd_ptr->Run();
}



BOOST_AUTO_TEST_SUITE_END() // end the zerodim sub-suite

BOOST_AUTO_TEST_SUITE_END() // end the blackbox suite




