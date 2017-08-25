//This file is part of Bertini 2.
//
//fixed_double_powerseries_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fixed_double_powerseries_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fixed_double_powerseries_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
#include "bertini2/endgames/powerseries.hpp"


BOOST_AUTO_TEST_SUITE(fixed_double_powerseries_endgame)

BOOST_AUTO_TEST_SUITE(generic_tests_ambient_precision_16)

using namespace bertini::endgame;

using TrackerType = bertini::tracking::DoublePrecisionTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::PSEG;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = 16;
#include "test/endgames/generic_pseg_test.hpp"

BOOST_AUTO_TEST_SUITE_END()







BOOST_AUTO_TEST_SUITE(generic_tests_ambient_precision_30)

using namespace bertini::endgame;

using TrackerType = bertini::tracking::DoublePrecisionTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::PSEG;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = 30;
#include "test/endgames/generic_pseg_test.hpp"

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

