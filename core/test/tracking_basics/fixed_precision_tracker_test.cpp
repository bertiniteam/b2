//This file is part of Bertini 2.
//
//fixed_precision_tracker_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fixed_precision_tracker_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fixed_precision_tracker_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire



#include <boost/test/unit_test.hpp>
#include "bertini2/system/start_systems.hpp"
#include "bertini2/trackers/fixed_precision_tracker.hpp"
#include "bertini2/trackers/observers.hpp"



extern double threshold_clearance_d;
extern bertini::real_mp threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;



BOOST_AUTO_TEST_SUITE(fixed_precision_tracker_basics)

using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using complex_dbl = std::complex<double>;
using mpfr = bertini::complex_mp;
using real_mp = bertini::real_mp;


template<typename NumT> using Vec = bertini::Vec<NumT>;
template<typename NumT> using Mat = bertini::Mat<NumT>;
using bertini::DefaultPrecision;
BOOST_AUTO_TEST_CASE(double_tracker_track_linear)
{
	using bertini::operator<<;
	DefaultPrecision(100);
	using namespace bertini::tracking;

	Var y = Variable::Make("y");
	Var t = Variable::Make("t");

	System sys;

	VariableGroup v{y};

	sys.AddFunction(y-t);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);


	bertini::tracking::DoublePrecisionTracker tracker(sys);


	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;


	tracker.Setup(Predictor::Euler,
	              double(1e-5),
					double(1e5),
					stepping_preferences,
					newton_preferences);

	
	complex_dbl t_start(1);
	complex_dbl t_end(0);
	
	Vec<complex_dbl> y_start(1);
	y_start << complex_dbl(1);

	Vec<complex_dbl> y_end;

	auto obs = GoryDetailLogger<DoublePrecisionTracker>();
		tracker.AddObserver(obs);

	auto code = tracker.TrackPath(y_end, t_start, t_end, y_start);
	BOOST_CHECK(code==bertini::SuccessCode::Success);

	BOOST_CHECK_EQUAL(y_end.size(),1);
	BOOST_CHECK(abs(y_end(0)-complex_dbl(0)) < 1e-5);

}
	
BOOST_AUTO_TEST_CASE(multiple_100_tracker_track_linear)
{
	DefaultPrecision(100);
	using namespace bertini::tracking;

	Var y = Variable::Make("y");
	Var t = Variable::Make("t");

	System sys;

	VariableGroup v{y};

	sys.AddFunction(y-t);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);


	bertini::tracking::MultiplePrecisionTracker tracker(sys);


	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;


	tracker.Setup(Predictor::Euler,
	              1e-5,
					1e5,
					stepping_preferences,
					newton_preferences);

	GoryDetailLogger<MultiplePrecisionTracker> tons_of_detail;
	tracker.AddObserver(tons_of_detail);

	
	mpfr t_start(1);
	mpfr t_end(0);
	
	Vec<mpfr> y_start(1);
	y_start << mpfr(1);

	Vec<mpfr> y_end;

	tracker.TrackPath(y_end,
	                  t_start, t_end, y_start);

	BOOST_CHECK_EQUAL(y_end.size(),1);
	BOOST_CHECK(abs(y_end(0)-mpfr(0)) < 1e-5);

}


// REGRESSION: the predictor's condition-number refresh counter was passed to SetNormsCond
// BY VALUE.  The tracker initializes the counter TO frequency_of_CN_estimation, so with the
// by-value bug the member never advanced past that value and `counter >= frequency` was true
// on EVERY step -- i.e. the estimate silently refreshed every step and the frequency knob was
// inert (harmless at the default of 1, wrong for anything larger).  With the counter passed by
// reference, the estimate must hold constant between refreshes: on a path whose conditioning
// varies continuously, the recorded estimate may change at most every frequency-th step.
namespace {
struct CondNumberRecorder : public bertini::Observer<bertini::tracking::DoublePrecisionTracker>
{
	using Emitter = bertini::tracking::TrackerTraits<bertini::tracking::DoublePrecisionTracker>::EventEmitterType;

	std::vector<double> conds; ///< LatestConditionNumber after each successful step.

	bertini::ObserveResult Observe(bertini::AnyEvent const& e) override
	{
		if (auto p = dynamic_cast<const bertini::tracking::SuccessfulStep<Emitter>*>(&e))
			conds.push_back(static_cast<double>(p->Get().LatestConditionNumber()));
		return bertini::ObserveResult::KeepObserving;
	}
};
}

BOOST_AUTO_TEST_CASE(condition_number_refresh_honors_frequency)
{
	DefaultPrecision(100);
	using namespace bertini::tracking;

	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var t = Variable::Make("t");

	// TWO variables, deliberately: for a univariate system the estimate is constant BY
	// ALGEBRA -- ||J||*||J^{-1} r|| = |J|*|r|/|J| = |r|, the fixed probe's norm -- so a
	// 1-var version of this test can only "pass" through floating-point jitter (and on
	// macOS it doesn't).  Here J = [[2x, 2y], [1, -1]] genuinely varies along the path.
	System sys;
	VariableGroup v{x, y};
	sys.AddFunction(x*x + y*y - t - 1);
	sys.AddFunction(x - y - t/2);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	DoublePrecisionTracker tracker(sys);

	SteppingConfig stepping_preferences;
	stepping_preferences.frequency_of_CN_estimation = 3;
	stepping_preferences.max_step_size = 0.01;    // plenty of steps to observe the cadence
	NewtonConfig newton_preferences;
	tracker.Setup(Predictor::Euler,
	              double(1e-5),
	              double(1e5),
	              stepping_preferences,
	              newton_preferences);

	CondNumberRecorder recorder;
	tracker.AddObserver(recorder);

	// a real solution at t = 1:  x = y + 1/2,  2y^2 + y + 1/4 = 2
	const double y0 = (-1.0 + sqrt(15.0)) / 4.0;
	Vec<complex_dbl> start(2);
	start << complex_dbl(y0 + 0.5), complex_dbl(y0);
	Vec<complex_dbl> end_point;
	auto code = tracker.TrackPath(end_point, complex_dbl(1), complex_dbl(0), start);
	tracker.RemoveObserver(recorder);
	BOOST_CHECK(code==bertini::SuccessCode::Success);

	BOOST_REQUIRE_GE(recorder.conds.size(), 9u);
	unsigned changes = 0;
	for (size_t ii = 1; ii < recorder.conds.size(); ++ii)
		if (recorder.conds[ii] != recorder.conds[ii-1])
			++changes;

	BOOST_CHECK_GE(changes, 1u);                         // it does refresh...
	BOOST_CHECK_LT(2*changes, recorder.conds.size());    // ...but at most every 3rd step,
	                                                     // not every step (the by-value bug)
}




// Regression test for the step budget (b2 #410, second half).
//
// max_num_steps is documented as "the maximum number of steps allowed during tracking",
// but it was compared against num_successful_steps_taken_ only.  A path whose steps FAIL
// therefore never approached its budget -- it could fail indefinitely while sitting at 0%
// of its allowance.  Observed in the wild: 3.9 million failed steps against 343 successful
// ones, thirty minutes at 100% CPU, the enclosing Solve() never returning.
//
// Here every step fails on purpose: a tracking tolerance of 1e-30 is unreachable in double
// precision, so the corrector never converges and every iteration is a failed step.  With
// max_num_steps = 10 the budget must stop it.
//
// The case DISCRIMINATES.  Counting only successes, this path never reaches the budget
// (successes stay at 0) and instead halves its stepsize ~330 times down to the 1e-100
// min_step_size floor, returning MinStepSizeReached.  Counting every step, it returns
// MaxNumStepsTaken after 10.  So the asserted code distinguishes the fix from its absence.
BOOST_AUTO_TEST_CASE(step_budget_counts_failed_steps_not_only_successes)
{
	DefaultPrecision(30);
	using namespace bertini::tracking;

	Var y = Variable::Make("y");
	Var t = Variable::Make("t");

	System sys;
	VariableGroup v{y};
	sys.AddFunction(y*y - t);          // sqrt path: a genuine tracking problem
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	bertini::tracking::DoublePrecisionTracker tracker(sys);

	SteppingConfig stepping_preferences;
	stepping_preferences.max_num_steps = 10;      // the budget under test
	NewtonConfig newton_preferences;

	tracker.Setup(Predictor::Euler,
	              double(1e-30),                  // UNREACHABLE at double precision
	              double(1e5),
	              stepping_preferences,
	              newton_preferences);

	complex_dbl t_start(1);
	complex_dbl t_end(0);
	Vec<complex_dbl> start_point(1);
	start_point << complex_dbl(1);

	Vec<complex_dbl> end_point;
	auto code = tracker.TrackPath(end_point, t_start, t_end, start_point);

	BOOST_CHECK(code != bertini::SuccessCode::Success);
	BOOST_CHECK(code == bertini::SuccessCode::MaxNumStepsTaken);
	// and the budget is what stopped it: total steps did not exceed the allowance
	BOOST_CHECK_LE(tracker.NumTotalStepsTaken(), 11u);
}


BOOST_AUTO_TEST_SUITE_END()




