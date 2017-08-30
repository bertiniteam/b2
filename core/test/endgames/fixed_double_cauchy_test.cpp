//This file is part of Bertini 2.
//
//fixed_double_cauchy_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fixed_double_cauchy_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fixed_double_cauchy_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
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

#include "bertini2/endgames/fixed_prec_endgame.hpp"
#include "bertini2/endgames/cauchy.hpp"

//THIS NEEDS TO BE IMPLEMENTED!
// #include "bertini2/endgames/observers.hpp"
#include "bertini2/trackers/observers.hpp"


BOOST_AUTO_TEST_SUITE(fixed_double_cauchy_endgame)


// this first suite tests whether cauchy endgame works correctly when ambient precision of mpfr is essentially double precision
BOOST_AUTO_TEST_SUITE(generic_tests_precision_16)

using namespace bertini::tracking;
using namespace bertini::endgame;

using TrackerType = DoublePrecisionTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::Cauchy;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = bertini::DoublePrecision();

#include "test/endgames/generic_cauchy_test.hpp"

BOOST_AUTO_TEST_SUITE_END() // generic tests at some precision






// this first suite tests whether cauchy endgame works correctly when ambient precision of mpfr is higher than double precision
BOOST_AUTO_TEST_SUITE(generic_tests_precision_30)

using namespace bertini::tracking;
using namespace bertini::endgame;

using TrackerType = DoublePrecisionTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::Cauchy;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = 30;

#include "test/endgames/generic_cauchy_test.hpp"

BOOST_AUTO_TEST_SUITE_END() // generic tests at some precision


BOOST_AUTO_TEST_SUITE_END() // the overall test suite for fixed multiple
