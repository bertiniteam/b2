//This file is part of Bertini 2.
//
//amp_cauchy_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_cauchy_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_cauchy_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  silviana amethyst, university of wisconsin eau claire
//
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Fall 2015, Spring 2016





#include <iostream>
#include <boost/test/unit_test.hpp>

#include "bertini2/system/start_systems.hpp"
#include "bertini2/num_traits.hpp"

#include "bertini2/endgames/amp_endgame.hpp"
#include "bertini2/endgames/cauchy.hpp"


#include "bertini2/endgames/observers.hpp"
#include "bertini2/trackers/observers.hpp"

// top level test suite for adaptive precision cauchy test
BOOST_AUTO_TEST_SUITE(adaptive_precision_cauchy_endgame)



BOOST_AUTO_TEST_SUITE(generic_tests_ambient_precision_16)


using namespace bertini::endgame;

using TrackerType = bertini::tracking::AMPTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::Cauchy;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = bertini::DoublePrecision();

#include "test/endgames/generic_cauchy_test.hpp"

BOOST_AUTO_TEST_SUITE_END()


// repeat the tests at precision 30, higher than double precision
BOOST_AUTO_TEST_SUITE(generic_tests_ambient_precision_30)

using namespace bertini::endgame;

using TrackerType = bertini::tracking::AMPTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::Cauchy;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = 30;

#include "test/endgames/generic_cauchy_test.hpp"

BOOST_AUTO_TEST_SUITE_END()




// maybe this is overkill to test yet again at ambient 50 digits?
BOOST_AUTO_TEST_SUITE(generic_tests_ambient_precision_50)

using namespace bertini::endgame;

using TrackerType = bertini::tracking::AMPTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::Cauchy;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = 50;

#include "test/endgames/generic_cauchy_test.hpp"

BOOST_AUTO_TEST_SUITE_END()



// Exercises the complex_dbl -> complex_mp container MIGRATION that the well-conditioned tests (and
// cyclic-5, which escalates 0/70) never hit.  A well-conditioned path starts in the hardware-double
// fast lane, but the final tolerance is set tighter than double can deliver, so the in-double refine
// fails, the escalation trigger fires, the endgame migrates its containers to mpfr mid-run, and
// finishes in mpfr -- to a known-correct root.
BOOST_AUTO_TEST_SUITE(adaptive_numeric_type_migration)

using namespace bertini;
using namespace bertini::endgame;

using TrackerType   = bertini::tracking::AMPTracker;
using TestedEGType  = EndgameSelector<TrackerType>::Cauchy;
using PrecisionConfig = bertini::tracking::TrackerTraits<TrackerType>::PrecisionConfig;

BOOST_AUTO_TEST_CASE(tight_tolerance_forces_double_to_mpfr_migration)
{
	DefaultPrecision(DoublePrecision());

	System sys;
	auto x = node::Variable::Make("x");
	auto t = node::Variable::Make("t");
	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t );
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	bertini::tracking::SteppingConfig stepping_preferences;
	bertini::tracking::NewtonConfig newton_preferences;
	newton_preferences.max_num_newton_iterations = 2;
	newton_preferences.min_num_newton_iterations = 1;
	tracker.Setup(bertini::tracking::Predictor::HeunEuler, 1e-5, 1e5, stepping_preferences, newton_preferences);
	tracker.PrecisionSetup(precision_config);
	tracker.ReinitializeInitialStepSize(false);

	complex_mp time(real_mp("0.1"), real_mp("0.0"));
	Vec<complex_mp> sample(1);
	sample << complex_mp(real_mp("0.9000000000000001"), real_mp("0.4358898943540673"));
	Vec<complex_mp> x_origin(1);
	x_origin << complex_mp(1,0);

	TestedEGType eg(tracker);
	eg.SetBoundaryTime(time);
	eg.SetFinalTolerance(1e-25);   // double (~16 digits) cannot reach this -> the endgame must migrate to mpfr

	auto code = eg.Run(sample);

	BOOST_CHECK(code == SuccessCode::Success);
	// Converged to a tolerance only mpfr can deliver, so the endgame must have crossed above double.
	BOOST_CHECK_GT(Precision(eg.FinalApproximation<complex_mp>()), DoublePrecision());
	// ...and to the correct root.
	BOOST_CHECK((eg.FinalApproximation<complex_mp>() - x_origin).template lpNorm<Eigen::Infinity>() < 1e-22);
}

BOOST_AUTO_TEST_SUITE_END() // re: adaptive_numeric_type_migration


BOOST_AUTO_TEST_SUITE_END() // re:

