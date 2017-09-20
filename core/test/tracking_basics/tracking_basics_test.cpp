//This file is part of Bertini 2.
//
//tracking_basics_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracking_basics_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking_basics_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire




 
#define BOOST_TEST_DYN_LINK
//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Tracking Basics Testing"
#include <boost/test/unit_test.hpp>


#include "logging.hpp"


using sec_level = boost::log::trivial::severity_level;

using LoggingInit = bertini::LoggingInit;


BOOST_GLOBAL_FIXTURE( LoggingInit );


#include "mpfr_extensions.hpp"

double threshold_clearance_d(1e-15);
bertini::mpfr_float threshold_clearance_mp("1e-28");
unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS(30);




// deliberately left blank.  link other files with this one.

