//This file is part of Bertini 2.
//
//b2/test/endgames/interpolation.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//b2/test/endgames/interpolation.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with b2/test/endgames/interpolation.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University


#include <iostream>
#include <boost/test/unit_test.hpp>

#include "bertini2/system/start_systems.hpp"
#include "bertini2/num_traits.hpp"

#include "bertini2/endgames/config.hpp"
#include "bertini2/endgames/interpolation.hpp"




BOOST_AUTO_TEST_SUITE(interpolation)



BOOST_AUTO_TEST_SUITE(generic_tests_double_16)
using BaseComplexType = bertini::dbl;
unsigned ambient_precision = bertini::DoublePrecision();
#include "test/endgames/generic_interpolation.hpp"
BOOST_AUTO_TEST_SUITE_END() // generic tests at some precision


BOOST_AUTO_TEST_SUITE(generic_tests_double_30)
using BaseComplexType = bertini::dbl;
unsigned ambient_precision = 30;
#include "test/endgames/generic_interpolation.hpp"
BOOST_AUTO_TEST_SUITE_END() // generic tests at some precision


BOOST_AUTO_TEST_SUITE(generic_tests_mpfr_16)
using BaseComplexType = bertini::mpfr;
unsigned ambient_precision = bertini::DoublePrecision();
#include "test/endgames/generic_interpolation.hpp"
BOOST_AUTO_TEST_SUITE_END() // generic tests at some precision


BOOST_AUTO_TEST_SUITE(generic_tests_mpfr_30)
using BaseComplexType = bertini::mpfr;
unsigned ambient_precision = 30;
#include "test/endgames/generic_interpolation.hpp"
BOOST_AUTO_TEST_SUITE_END() // generic tests at some precision

BOOST_AUTO_TEST_SUITE(generic_tests_mpfr_100)
using BaseComplexType = bertini::mpfr;
unsigned ambient_precision = 100;
#include "test/endgames/generic_interpolation.hpp"
BOOST_AUTO_TEST_SUITE_END() // generic tests at some precision


BOOST_AUTO_TEST_SUITE_END() // interpolation suite

