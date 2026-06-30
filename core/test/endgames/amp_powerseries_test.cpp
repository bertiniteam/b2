//This file is part of Bertini 2.
//
//amp_powerseries_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_powerseries_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_powerseries_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
// Tim Hodges, Colorado State University


#include <iostream>
#include <boost/test/unit_test.hpp>

#include "bertini2/system/start_systems.hpp"
#include "bertini2/num_traits.hpp"

#include "bertini2/endgames/amp_endgame.hpp"
#include "bertini2/endgames/powerseries.hpp"

#include "bertini2/endgames/observers.hpp"
#include "bertini2/trackers/observers.hpp"
#define B2_OBSERVE_TRACKERS

BOOST_AUTO_TEST_SUITE(amp_powerseries_endgame_generic_tests_ambient_precision_16)

using namespace bertini::endgame;

using TrackerType = bertini::tracking::AMPTracker; // select a tracker type
using TestedEGType = typename EndgameSelector<TrackerType>::PSEG;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = 16;
#include "test/endgames/generic_pseg_test.hpp"

BOOST_AUTO_TEST_SUITE_END()




BOOST_AUTO_TEST_SUITE(amp_powerseries_endgame_generic_tests_ambient_precision_30)

using namespace bertini::endgame;

using TrackerType = bertini::tracking::AMPTracker; // select a tracker type
using TestedEGType = typename EndgameSelector<TrackerType>::PSEG;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = 30;
#include "test/endgames/generic_pseg_test.hpp"

BOOST_AUTO_TEST_SUITE_END()





BOOST_AUTO_TEST_SUITE(amp_powerseries_endgame_generic_tests_ambient_precision_50)

using namespace bertini::endgame;

using TrackerType = bertini::tracking::AMPTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::PSEG;
auto TestedPredictor = bertini::tracking::Predictor::HeunEuler;
unsigned ambient_precision = 50;
#include "test/endgames/generic_pseg_test.hpp"

BOOST_AUTO_TEST_SUITE_END()






BOOST_AUTO_TEST_SUITE(amp_powerseries_endgame_AMPspecific_tests)

using namespace bertini::endgame;

using TrackerType = bertini::tracking::AMPTracker; // select a tracker type
using TestedEGType = EndgameSelector<TrackerType>::PSEG;
using mpfr = bertini::complex_mp;

using namespace bertini;
BOOST_AUTO_TEST_CASE(ensure_uniform_precision_16_30_40)
{
	TimeCont<mpfr> times;
	SampCont<mpfr> samples;

	DefaultPrecision(16);
	times.emplace_back(RandomUnit<mpfr>());

	DefaultPrecision(30);
	times.emplace_back(RandomUnit<mpfr>());

	DefaultPrecision(40);
	times.emplace_back(RandomUnit<mpfr>());



	DefaultPrecision(16);
	samples.emplace_back(RandomOfUnits<mpfr>(4));

	DefaultPrecision(30);
	samples.emplace_back(RandomOfUnits<mpfr>(4));

	DefaultPrecision(40);
	samples.emplace_back(RandomOfUnits<mpfr>(4));


	auto max_precision = TestedEGType::EnsureAtUniformPrecision(times, samples);

	BOOST_CHECK_EQUAL(max_precision,40);

	for (const auto& t : times)
		BOOST_CHECK_EQUAL(Precision(t),40);

	for (const auto& s : samples)
		BOOST_CHECK_EQUAL(Precision(s(0)),40);
}



BOOST_AUTO_TEST_CASE(ensure_uniform_precision_all_uniform_to_start)
{
	TimeCont<mpfr> times;
	SampCont<mpfr> samples;

	DefaultPrecision(30);
	times.emplace_back(RandomUnit<mpfr>());

	DefaultPrecision(30);
	times.emplace_back(RandomUnit<mpfr>());

	DefaultPrecision(30);
	times.emplace_back(RandomUnit<mpfr>());



	DefaultPrecision(30);
	samples.emplace_back(RandomOfUnits<mpfr>(4));

	DefaultPrecision(30);
	samples.emplace_back(RandomOfUnits<mpfr>(4));

	DefaultPrecision(30);
	samples.emplace_back(RandomOfUnits<mpfr>(4));


	auto max_precision = TestedEGType::EnsureAtUniformPrecision(times, samples);

	BOOST_CHECK_EQUAL(max_precision,30);

	for (const auto& t : times)
		BOOST_CHECK_EQUAL(Precision(t),30);

	for (const auto& s : samples)
		BOOST_CHECK_EQUAL(Precision(s(0)),30);
}


// Exercises the complex_dbl -> complex_mp container migration for the PowerSeries endgame, the same way
// the Cauchy test does: a well-conditioned path starts in the hardware-double fast lane, but the final
// tolerance is set tighter than double can deliver, so the in-double refine fails, the escalation trigger
// fires, the endgame migrates its containers (times/samples/derivatives/rand_vector) to mpfr mid-run, and
// finishes in mpfr -- to a known-correct root.
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

	auto precision_config = bertini::tracking::TrackerTraits<TrackerType>::PrecisionConfig(sys);

	TrackerType tracker(sys);
	bertini::tracking::SteppingConfig stepping_preferences;
	bertini::tracking::NewtonConfig newton_preferences;
	newton_preferences.max_num_newton_iterations = 2;
	newton_preferences.min_num_newton_iterations = 1;
	tracker.Setup(bertini::tracking::Predictor::HeunEuler, 1e-5, 1e5, stepping_preferences, newton_preferences);
	tracker.PrecisionSetup(precision_config);
	tracker.ReinitializeInitialStepSize(false);

	mpfr time(real_mp("0.1"), real_mp("0.0"));
	Vec<mpfr> sample(1);
	sample << mpfr(real_mp("0.9000000000000001"), real_mp("0.4358898943540673"));
	Vec<mpfr> x_origin(1);
	x_origin << mpfr(1,0);

	TestedEGType eg(tracker);
	eg.SetBoundaryTime(time);
	eg.SetFinalTolerance(1e-25);   // double (~16 digits) cannot reach this -> the endgame must migrate to mpfr

	auto code = eg.Run(sample);

	BOOST_CHECK(code == SuccessCode::Success);
	// Converged to a tolerance only mpfr can deliver, so the endgame must have crossed above double.
	BOOST_CHECK_GT(Precision(eg.FinalApproximation<mpfr>()), DoublePrecision());
	// ...and to the correct root.
	BOOST_CHECK((eg.FinalApproximation<mpfr>() - x_origin).template lpNorm<Eigen::Infinity>() < 1e-22);
}
BOOST_AUTO_TEST_SUITE_END()
