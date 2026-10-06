//This file is part of Bertini 2.
//
//tracker_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracker_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracker_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

#include <chrono>
#include <functional>
#include <set>
#include <vector>

#include <boost/test/unit_test.hpp>
#include "bertini2/random.hpp"
#include "bertini2/system/start_systems.hpp"
#include "bertini2/trackers/tracker.hpp"
#include "bertini2/trackers/observers.hpp"
#include "bertini2/records/solver_recording.hpp"   // CanonicalName, for readable test messages



extern double threshold_clearance_d;
extern bertini::real_mp threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;



BOOST_AUTO_TEST_SUITE(AMP_tracker_basics)

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
BOOST_AUTO_TEST_CASE(minstepsize)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;
    real_mp remaining_time("1e-10");

    BOOST_CHECK_EQUAL(MinStepSizeForPrecision(16, remaining_time),real_mp("1e-23"));

    BOOST_CHECK_CLOSE(MinStepSizeForPrecision(40, remaining_time),real_mp("1e-47"), real_mp("1e-28"));
}


BOOST_AUTO_TEST_CASE(mindigits)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    real_mp remaining_time("1e-30");
    real_mp min_stepsize("1e-35");
    real_mp max_stepsize("1e-33");

    auto digits = MinDigitsForStepsizeInterval(min_stepsize, max_stepsize, remaining_time);

    BOOST_CHECK_EQUAL(digits, 8);
}

BOOST_AUTO_TEST_CASE(AMP_tracker_track_linear)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{y};

    sys.AddFunction(y-t);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                  1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

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







BOOST_AUTO_TEST_CASE(AMP_tracker_track_quadratic)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{y};

    sys.AddFunction(y-pow(t,2));
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                  1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(-1);

    Vec<mpfr> y_start(1);
    y_start << mpfr(1);

    Vec<mpfr> y_end;

    tracker.TrackPath(y_end,
                      t_start, t_end, y_start);

    BOOST_CHECK_EQUAL(y_end.size(),1);
    BOOST_CHECK(abs(y_end(0)-mpfr(1)) < 1e-5);
}



BOOST_AUTO_TEST_CASE(AMP_tracker_track_decic)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{y};

    sys.AddFunction(y-pow(t,10));
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                    1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(-2);

    Vec<mpfr> y_start(1);
    y_start << mpfr(1);

    Vec<mpfr> y_end;

    tracker.TrackPath(y_end,
                      t_start, t_end, y_start);

    BOOST_CHECK_EQUAL(y_end.size(),1);
    BOOST_CHECK(abs(y_end(0)-mpfr("1024.0")) < 1e-5);

}




BOOST_AUTO_TEST_CASE(AMP_tracker_track_square_root)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddFunction(x-t);
    sys.AddFunction(pow(y,2)-x);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                    1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    Vec<mpfr> end_point;

    bertini::SuccessCode tracking_success;


    start_point << mpfr(1), mpfr(1);
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(),2);
    BOOST_CHECK(abs(end_point(0)-mpfr(0)) < 1e-5);
    BOOST_CHECK(abs(end_point(1)-mpfr(0)) < 1e-5);


    start_point << mpfr(1), mpfr(-1);
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(),2);
    BOOST_CHECK(abs(end_point(0)-mpfr(0)) < 1e-5);
    BOOST_CHECK(abs(end_point(1)-mpfr(0)) < 1e-5);

    t_start = mpfr(-1);
    start_point << mpfr(-1), mpfr(0,-1);
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(),2);
    BOOST_CHECK(abs(end_point(0)-mpfr(0)) < 1e-5);
    BOOST_CHECK(abs(end_point(1)-mpfr(0)) < 1e-5);

    start_point << mpfr(-1), mpfr(0,1);
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(),2);
    BOOST_CHECK(abs(end_point(0)-mpfr(0)) < 1e-5);
    BOOST_CHECK(abs(end_point(1)-mpfr(0)) < 1e-5);
}




/*
1.  Goal:  Handle a singular start point?
Expected Behavior:  Doesn't start tracking.
System:
  f = x^2 + (1-t)*x;
  g = y^2 + (1-t)*y;
Start t:   1
Start point:  (0,0)
End t:   0
End point:  N/A
*/
BOOST_AUTO_TEST_CASE(AMP_tracker_doesnt_start_from_singular_start_point)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddFunction(pow(x,2) + (1-t)*x);
    sys.AddFunction(pow(y,2) + (1-t)*y);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                    1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    Vec<mpfr> end_point;

    bertini::SuccessCode tracking_success;


    start_point << mpfr(0), mpfr(0);
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK(tracking_success==bertini::SuccessCode::SingularStartPoint);
    BOOST_CHECK_EQUAL(end_point.size(),0);
}



BOOST_AUTO_TEST_CASE(AMP_tracker_tracking_DOES_SOMETHING_PREDICTABLE_from_near_to_singular_start_point)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddFunction(pow(x,2) + (1-t)*x);
    sys.AddFunction(pow(y,2) + (1-t)*y);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;

    stepping_preferences.max_num_steps = 1e2;

    tracker.Setup(Predictor::Euler,
                    1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    Vec<mpfr> end_point;

    bertini::SuccessCode tracking_success;

    start_point << mpfr("1e-28"), mpfr("1e-28");
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK(tracking_success!=bertini::SuccessCode::Success && tracking_success!=bertini::SuccessCode::NeverStarted);
    BOOST_CHECK_EQUAL(end_point.size(),0);
}



/*
2.  Goal:  Make sure we can handle a simple, mixed (x,y in both polynomials), nonhomogeneous system.
Expected Behavior:  Success.
System:
  f = x^2 + (1-t)*x - 1;
 g = y^2 + (1-t)*x*y - 2;
Start t:  1
Start point:  (1, 1.414)
End t:  0
End point:  (6.180339887498949e-01, 1.138564265110173e+00)
(Using Bertini default tracking tolerances.)
*/
BOOST_AUTO_TEST_CASE(AMP_simple_nonhomogeneous_system_trackable_initialprecision16)
{
    DefaultPrecision(16);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddFunction(pow(x,2) + (1-t)*x - 1);
    sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                    1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    Vec<mpfr> end_point;

    bertini::SuccessCode tracking_success;


    start_point << mpfr(1), mpfr("1.41421356237309504880168872421");
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(),2);
    BOOST_CHECK(abs(end_point(0)-mpfr("6.180339887498949e-01")) < 1e-5);
    BOOST_CHECK(abs(end_point(1)-mpfr("1.138564265110173e+00")) < 1e-5);
}



BOOST_AUTO_TEST_CASE(AMP_simple_nonhomogeneous_system_trackable_initialprecision30_tighter_track_tol)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddFunction(pow(x,2) + (1-t)*x - 1);
    sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;

    newton_preferences.max_num_newton_iterations = 6;

    tracker.Setup(Predictor::Euler,
                    1e-30,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    DefaultPrecision(30);
    Vec<mpfr> start_point(2);
    Vec<mpfr> end_point;

    bertini::SuccessCode tracking_success;


    start_point << mpfr(1), mpfr("1.41421356237309504880168872421");
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    DefaultPrecision(40);

    Vec<mpfr> true_solution(2);
    true_solution <<  mpfr("0.61803398874989484820458683436563811772030918"), mpfr("1.13856426511017256414753784441721594451116198");


    BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(),2);
    BOOST_CHECK(abs(end_point(0)-true_solution(0)) < 1e-30);
    BOOST_CHECK(abs(end_point(1)-true_solution(1)) < 1e-30);

    BOOST_CHECK( (end_point - true_solution).norm() < 1e-30);
}



BOOST_AUTO_TEST_CASE(AMP_simple_nonhomogeneous_system_trackable_initialprecision30)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddFunction(pow(x,2) + (1-t)*x - 1);
    sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                    1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);

    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    Vec<mpfr> end_point;

    bertini::SuccessCode tracking_success;
    tracker.PrecisionPreservation(true);

    start_point << mpfr(1), mpfr("1.414");
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK_EQUAL(DefaultPrecision(),30);
    BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(),2);
    BOOST_CHECK(abs(end_point(0)-mpfr("6.180339887498949e-01")) < 1e-5);
    BOOST_CHECK(abs(end_point(1)-mpfr("1.138564265110173e+00")) < 1e-5);
}



BOOST_AUTO_TEST_CASE(AMP_simple_nonhomogeneous_system_trackable_initialprecision100)
{
    DefaultPrecision(100);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddFunction(pow(x,2) + (1-t)*x - 1);
    sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);



    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                    1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);
    tracker.PrecisionPreservation(true);
    auto AMP = bertini::tracking::AMPConfigFrom(sys);
    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    Vec<mpfr> end_point;

    bertini::SuccessCode tracking_success;


    start_point << mpfr(1), mpfr("1.414");
    tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK_EQUAL(DefaultPrecision(),100);
    BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(),2);
    BOOST_CHECK(abs(end_point(0)-mpfr("6.180339887498949e-01")) < 1e-5);
    BOOST_CHECK(abs(end_point(1)-mpfr("1.138564265110173e+00")) < 1e-5);
}







/*
5.  Goal:  Fail when running into a singularity.
Expected Behavior:  Path failure near t=0.5.
Note:  There's a parameter in this system, which depends on the path variable.  This implicitly checks that functionality.
System:
  s= -1*(1-t) + 1*t;
 f = x^2-s;
 g = y^2-s;
Start t:  1
Start point:  (1,1)
End t:  0
End point:  N/A (Path should fail at t=0.5)
*/
BOOST_AUTO_TEST_CASE(AMP_tracker_fails_with_singularity_on_path)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    auto s = -1*(1-t) + 1*t;

    System sys;

    VariableGroup v{x,y};

    sys.AddFunction(pow(x,2) - s);
    sys.AddFunction(pow(y,2) - s);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);



    bertini::tracking::AMPTracker tracker(sys);


    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;


    tracker.Setup(Predictor::Euler,
                    1e-5,
                    1e5,
                    stepping_preferences,
                    newton_preferences);
    tracker.PrecisionPreservation(true);
    auto AMP = bertini::tracking::AMPConfigFrom(sys);
    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    Vec<mpfr> end_point;

    start_point << mpfr(1), mpfr(1);
    bertini::SuccessCode tracking_success = tracker.TrackPath(end_point,
                      t_start, t_end, start_point);

    BOOST_CHECK(tracking_success!=bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(DefaultPrecision(),30);
}


// has two solutions:
//
// 1.
// x = -0.61803398874989484820458683
// y = 1.6180339887498948482045868
//
// 2.
// x = 0.16180339887498948482045868
// y = -0.6180339887498948482045868



BOOST_AUTO_TEST_CASE(AMP_track_total_degree_start_system)
{
    using namespace bertini::tracking;
    DefaultPrecision(30);

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddVariableGroup(v);

    sys.AddFunction(x*y+1);
    sys.AddFunction(x+y-1);
    sys.Homogenize();
    sys.AutoPatch();

    BOOST_CHECK(sys.IsHomogeneous());
    BOOST_CHECK(sys.IsPatched());



    auto TD = bertini::start_system::TotalDegreeBinomial(sys);
    TD.Homogenize();
    BOOST_CHECK(TD.IsHomogeneous());
    BOOST_CHECK(TD.IsPatched());

    auto final_system = (1-t)*sys + t*TD;
    final_system.AddPathVariable(t);

    auto tracker = AMPTracker(final_system);
    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;
    tracker.Setup(Predictor::Euler,
                    1e-5, 1e5,
                    stepping_preferences, newton_preferences);

    tracker.PrecisionSetup(bertini::tracking::AMPConfigFrom(final_system));
    tracker.PrecisionPreservation(true);
    mpfr t_start(1), t_end(0);
    std::vector<Vec<mpfr> > solutions;
    for (unsigned ii = 0; ii < TD.NumStartPoints(); ++ii)
    {
        auto start_point = TD.StartPoint<mpfr>(ii);

        Vec<mpfr> result;
        bertini::SuccessCode tracking_success = tracker.TrackPath(result,t_start,t_end,start_point);
        BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
        BOOST_CHECK_EQUAL(DefaultPrecision(),30);
        solutions.push_back(final_system.DehomogenizePoint(result));
    }

    Vec<mpfr> solution_1(2);
    solution_1 << mpfr("-0.61803398874989484820458683","0"), mpfr("1.6180339887498948482045868","0");

    Vec<mpfr> solution_2(2);
    solution_2 << mpfr("1.6180339887498948482045868","0"), mpfr("-0.6180339887498948482045868","0");

    unsigned num_occurences(0);
    for (auto s : solutions)
    {
        if ( (s-solution_1).norm() < real_mp("1e-5"))
            num_occurences++;
    }
    BOOST_CHECK_EQUAL(num_occurences,1);

    num_occurences = 0;
    for (auto s : solutions)
    {
        if ( (s-solution_2).norm() < real_mp("1e-5"))
            num_occurences++;
    }
    BOOST_CHECK_EQUAL(num_occurences,1);
}



std::vector<Vec<mpfr> > track_total_degree(bertini::tracking::AMPTracker const& tracker, bertini::start_system::TotalDegreeBinomial const& TD)
{
    [[maybe_unused]] auto initial_precision = DefaultPrecision();
    using namespace bertini::tracking;
    mpfr t_start(1), t_end(0);
    std::vector<Vec<mpfr> > solutions;
    for (unsigned ii = 0; ii < TD.NumStartPoints(); ++ii)
    {
        auto start_point = TD.StartPoint<mpfr>(ii);

        Vec<mpfr> result;
        bertini::SuccessCode tracking_success;

        tracking_success = tracker.TrackPath(result,t_start,t_end,start_point);
        BOOST_CHECK(tracking_success==bertini::SuccessCode::Success);
        solutions.push_back(tracker.GetSystem().DehomogenizePoint(result));
    }

    return solutions;
}

BOOST_AUTO_TEST_CASE(AMP_track_TD_functionalized)
{
    using namespace bertini::tracking;
    DefaultPrecision(30);

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;

    VariableGroup v{x,y};

    sys.AddVariableGroup(v);

    sys.AddFunction(x*y+1);
    sys.AddFunction(x+y-1);
    sys.Homogenize();
    sys.AutoPatch();

    BOOST_CHECK(sys.IsHomogeneous());
    BOOST_CHECK(sys.IsPatched());



    auto TD = bertini::start_system::TotalDegreeBinomial(sys);
    TD.Homogenize();
    BOOST_CHECK(TD.IsHomogeneous());
    BOOST_CHECK(TD.IsPatched());

    auto final_system = (1-t)*sys + t*TD;
    final_system.AddPathVariable(t);

    auto tracker = AMPTracker(final_system);
    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;
    tracker.Setup(Predictor::Euler,
                    1e-5, 1e5,
                    stepping_preferences, newton_preferences);

    auto AMP = bertini::tracking::AMPConfigFrom(final_system);
    tracker.PrecisionSetup(AMP);
    tracker.PrecisionPreservation(true);
    auto solutions = track_total_degree(tracker, TD);

    BOOST_CHECK_EQUAL(DefaultPrecision(),30);
    Vec<mpfr> solution_1(2);
    solution_1 << mpfr("-0.61803398874989484820458683","0"), mpfr("1.6180339887498948482045868","0");

    Vec<mpfr> solution_2(2);
    solution_2 << mpfr("1.6180339887498948482045868","0"), mpfr("-0.6180339887498948482045868","0");

    unsigned num_occurences(0);
    for (auto s : solutions)
    {
        if ( (s-solution_1).norm() < real_mp("1e-5"))
            num_occurences++;
    }
    BOOST_CHECK_EQUAL(num_occurences,1);

    num_occurences = 0;
    for (auto s : solutions)
    {
        if ( (s-solution_2).norm() < real_mp("1e-5"))
            num_occurences++;
    }
    BOOST_CHECK_EQUAL(num_occurences,1);
}


BOOST_AUTO_TEST_CASE(AMP_tracker_track_circle_line_RKF45)
{
    // Regression test: RKF45 + AMPTracker used to throw MinimizeTrackingCost
    // because AdjustAMPStepSuccess capped max_precision at current_precision_,
    // collapsing the scan window to a single value that failed criterion B.
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;
    VariableGroup v{x, y};
    sys.AddVariableGroup(v);
    sys.AddPathVariable(t);
    sys.AddFunction(t * (pow(x, 2) - 1) + (1 - t) * (pow(x, 2) + pow(y, 2) - 4));
    sys.AddFunction(t * (y - 1) + (1 - t) * (2 * x + 5 * y));

    auto AMP = bertini::tracking::AMPConfigFrom(sys);

    bertini::tracking::AMPTracker tracker(sys);
    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;

    tracker.Setup(Predictor::RKF45,
                  1e-5,
                  1e5,
                  stepping_preferences,
                  newton_preferences);
    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    start_point << mpfr(1), mpfr(1);

    Vec<mpfr> end_point;
    bertini::SuccessCode tracking_success;
    tracking_success = tracker.TrackPath(end_point, t_start, t_end, start_point);

    BOOST_CHECK(tracking_success == bertini::SuccessCode::Success);
    BOOST_CHECK_EQUAL(end_point.size(), 2);
}


BOOST_AUTO_TEST_CASE(arithmetic_cost_double_precision_is_one)
{
    using namespace bertini::tracking;
    BOOST_CHECK_EQUAL(ArithmeticCost(bertini::DoublePrecision()), 1.0);
}

BOOST_AUTO_TEST_CASE(arithmetic_cost_increases_with_precision)
{
    using namespace bertini::tracking;
    BOOST_CHECK_GT(ArithmeticCost(30), 1.0);
    BOOST_CHECK_GT(ArithmeticCost(100), ArithmeticCost(30));
}

// StepsizeSatisfyingCriterionB(p, digits_B=50, newton=2, predictor=0) = 10^(-(50-p)*2).
// At p<=40 this is <= 1e-20, well below min_stepsize=1e-5, so those precisions are
// skipped.  Only p=50 survives (stepsize=1), so that must be the result.
BOOST_AUTO_TEST_CASE(minimize_tracking_cost_skips_precision_below_min_stepsize)
{
    DefaultPrecision(50);
    using namespace bertini::tracking;

    unsigned new_precision = 0;
    real_mp new_stepsize("0");
    real_mp min_stepsize("1e-5");
    real_mp max_stepsize("1");

    MinimizeTrackingCost(new_precision, new_stepsize,
                         bertini::DoublePrecision(), min_stepsize,
                         50u, max_stepsize,
                         50u, 2u, 0u);

    BOOST_CHECK_EQUAL(new_precision, 50u);
    BOOST_CHECK(new_stepsize >= min_stepsize);
}

// With min_stepsize=2, max_precision=30, digits_B=30, newton=2:
// criterion B gives at most stepsize=1 (at p=30), which is < 2, so nothing passes.
BOOST_AUTO_TEST_CASE(minimize_tracking_cost_throws_when_no_precision_satisfies_min_stepsize)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    unsigned new_precision = 0;
    real_mp new_stepsize("0");
    real_mp min_stepsize("2");
    real_mp max_stepsize("10");

    BOOST_CHECK_THROW(
        MinimizeTrackingCost(new_precision, new_stepsize,
                             bertini::DoublePrecision(), min_stepsize,
                             30u, max_stepsize,
                             30u, 2u, 0u),
        std::runtime_error
    );
}

// Regression test for the AMP stepping deadlock (b2 #410).
//
// MinimizeTrackingCost throws when no (precision, stepsize) pair anywhere in the
// allowed window satisfies criterion B.  Both AMP adjusters catch that throw and
// report FailedToSelectPrecisionAndStepsize -- but the throw happens BEFORE either
// value is assigned, so current_precision_ and current_stepsize_ survive the failed
// step untouched.  The stepping loop used to treat that code like any other failed
// step and simply try again; with nothing changed, the retry was bit-identical, and
// TrackPath never returned.  Observed in the wild: 750k failed steps at a frozen
// precision of 20 and a frozen stepsize, one core pegged, no progress, forever.
//
// Here the window is emptied on purpose -- a tracking tolerance demanding ~40 digits
// against maximum_precision = 20 -- so no valid pair exists.  TrackPath must RETURN.
// The timeout is the regression guard: without the fix this case hangs rather than
// fails, which would wedge CI instead of reporting.
BOOST_AUTO_TEST_CASE(amp_tracker_terminates_when_no_precision_stepsize_pair_exists,
                     * boost::unit_test::timeout(60))
{
    DefaultPrecision(16);
    using namespace bertini::tracking;

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;
    VariableGroup v{x, y};
    sys.AddVariableGroup(v);
    sys.AddPathVariable(t);
    sys.AddFunction(t * (pow(x, 2) - 1) + (1 - t) * (pow(x, 2) + pow(y, 2) - 4));
    sys.AddFunction(t * (y - 1) + (1 - t) * (2 * x + 5 * y));

    auto AMP = bertini::tracking::AMPConfigFrom(sys);
    AMP.maximum_precision = 20;   // ceiling BELOW the digits the tolerance demands

    bertini::tracking::AMPTracker tracker(sys);
    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;

    tracker.Setup(Predictor::Euler,
                  1e-40,          // ~40 digits required; the ceiling is 20
                  1e5,
                  stepping_preferences,
                  newton_preferences);
    tracker.PrecisionSetup(AMP);

    mpfr t_start(1);
    mpfr t_end(0);

    Vec<mpfr> start_point(2);
    start_point << mpfr(1), mpfr(1);

    Vec<mpfr> end_point;
    auto tracking_success = tracker.TrackPath(end_point, t_start, t_end, start_point);

    // The point of the test: it comes back at all, and it comes back a failure.
    BOOST_CHECK(tracking_success != bertini::SuccessCode::Success);
    BOOST_CHECK(tracking_success == bertini::SuccessCode::FailedToSelectPrecisionAndStepsize);
}

// The predicate the stepping loop consults.  A terminal code is one for which the
// tracker took no corrective action, so retrying repeats the identical step.
BOOST_AUTO_TEST_CASE(terminal_step_codes_are_exactly_the_ones_that_adjust_nothing)
{
    using namespace bertini::tracking;

    BOOST_CHECK(StepFailureIsTerminal(bertini::SuccessCode::FailedToSelectPrecisionAndStepsize));

    // every one of these leaves the tracker adjusted, so the loop may retry
    BOOST_CHECK(!StepFailureIsTerminal(bertini::SuccessCode::Success));
    BOOST_CHECK(!StepFailureIsTerminal(bertini::SuccessCode::HigherPrecisionNecessary));
    BOOST_CHECK(!StepFailureIsTerminal(bertini::SuccessCode::FailedToConverge));
    BOOST_CHECK(!StepFailureIsTerminal(bertini::SuccessCode::MatrixSolveFailure));
    BOOST_CHECK(!StepFailureIsTerminal(bertini::SuccessCode::MatrixSolveFailureFirstPartOfPrediction));
}

BOOST_AUTO_TEST_CASE(set_start_precision_overrides_to_higher_precision)
{
    DefaultPrecision(16);
    using namespace bertini::tracking;

    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;
    VariableGroup v{y};
    sys.AddFunction(y - t);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);
    bertini::tracking::AMPTracker tracker(sys);

    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;
    tracker.Setup(Predictor::Euler, 1e-5, 1e5, stepping_preferences, newton_preferences);
    tracker.PrecisionSetup(AMP);
    tracker.SetStartPrecision(30);

    Vec<mpfr> y_start(1);
    y_start << mpfr(1);
    mpfr t_start(1), t_end(0);
    Vec<mpfr> y_end;

    auto code = tracker.TrackPath(y_end, t_start, t_end, y_start);
    BOOST_CHECK(code == bertini::SuccessCode::Success);
    BOOST_CHECK(abs(y_end(0) - mpfr(0)) < 1e-5);
}

namespace {

// Tracked from all ones at t = 1.  Cube root: y in (1-t)(y^3 - 2) + t(y^3 - 1), to the cube
// root of 2 at t = 0, a smooth path that drops to double on its own.  Square root: (x, y) in
// x^2 - t, y - x, toward the double root at t = 0; the Jacobian's condition grows like 1/x, so
// a track to t = 1e-20 at tracking tolerance 1e-12 ends in multiple precision (30 digits),
// where the condition estimate uses the multiprecision probe direction.  (One unknown would
// not do: a 1x1 Jacobian has condition number 1.  And at a loose tolerance the corrector
// accepts points well off the path near the root, so the track never needs the digits.)
struct TestHomotopy
{
    enum class Kind { CubeRoot, SquareRoot };

    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");
    System sys;
    double tracking_tolerance;
    mpfr t_end;

    explicit TestHomotopy(Kind kind = Kind::CubeRoot)
        : tracking_tolerance(kind == Kind::CubeRoot ? 1e-5 : 1e-12),
          t_end(kind == Kind::CubeRoot ? mpfr(0) : mpfr("1e-20"))
    {
        if (kind == Kind::CubeRoot)
        {
            sys.AddFunction((1 - t) * (pow(y, 3) - 2) + t * (pow(y, 3) - 1));
            sys.AddVariableGroup(VariableGroup{y});
        }
        else
        {
            sys.AddFunction(pow(x, 2) - t);
            sys.AddFunction(y - x);
            sys.AddVariableGroup(VariableGroup{x, y});
        }
        sys.AddPathVariable(t);
    }

    // Every tracker made here draws the same condition-number probe direction: the same seed,
    // at the same precision.
    bertini::tracking::AMPTracker MakeTracker() const
    {
        using namespace bertini::tracking;
        DefaultPrecision(30);
        bertini::ReseedThisThread(20261005);
        AMPTracker tracker(sys);
        tracker.Setup(Predictor::RKF45, tracking_tolerance, 1e5, SteppingConfig(), NewtonConfig());
        tracker.PrecisionSetup(AMPConfigFrom(sys));
        return tracker;
    }
};

// Everything a track leaves behind that a caller can read: the code, the step count, the
// precision, the step size, the latest condition-number estimate (which is computed from the
// probe direction) and the endpoint, compared bit for bit and precision for precision.
struct TrackOutcome
{
    bertini::SuccessCode code;
    unsigned steps;
    unsigned precision;
    real_mp stepsize;
    double condition_number;
    Vec<mpfr> endpoint;
};

TrackOutcome Track(bertini::tracking::AMPTracker const& tracker, unsigned start_digits,
                   mpfr const& t_end)
{
    DefaultPrecision(start_digits);
    Vec<mpfr> start(tracker.NumVariables());
    for (Eigen::Index ii = 0; ii < start.size(); ++ii)
        start(ii) = mpfr(1);
    Vec<mpfr> endpoint;
    auto code = tracker.TrackPath(endpoint, mpfr(1), t_end, start);
    return { code, tracker.NumTotalStepsTaken(), tracker.CurrentPrecision(),
             tracker.CurrentStepsize(), static_cast<double>(tracker.LatestConditionNumber()),
             endpoint };
}

void CheckIdentical(TrackOutcome const& a, TrackOutcome const& b)
{
    BOOST_CHECK(a.code == b.code);
    BOOST_CHECK_EQUAL(a.steps, b.steps);
    BOOST_CHECK_EQUAL(a.precision, b.precision);
    BOOST_CHECK_EQUAL(a.condition_number, b.condition_number);
    BOOST_CHECK_EQUAL(a.stepsize.precision(), b.stepsize.precision());
    BOOST_CHECK(a.stepsize == b.stepsize);
    BOOST_REQUIRE_EQUAL(a.endpoint.size(), b.endpoint.size());
    for (Eigen::Index ii = 0; ii < a.endpoint.size(); ++ii)
    {
        BOOST_CHECK_EQUAL(a.endpoint(ii).precision(), b.endpoint(ii).precision());
        BOOST_CHECK(a.endpoint(ii) == b.endpoint(ii));
    }
}

} // namespace


/**
A track's result does not depend on what the tracker tracked before it (#378).

A tracker is reused path after path -- by the zero-dim solver on each thread, by every
endgame, and by anyone driving one by hand.  Until 4.0 the first step of a track was built at
the precision the PREVIOUS track ended in, and kept that precision, so the same path tracked
by the same tracker gave different bits depending on its predecessor; a threaded solve, where
the scheduler picks each thread's predecessors, was irreproducible.  Every pairing of where the
previous track ended and where this one starts: above, at and below it, double and not.
*/
BOOST_AUTO_TEST_CASE(a_track_does_not_depend_on_the_track_before_it)
{
    using Kind = TestHomotopy::Kind;
    // (digits the previous track starts at, digits this track starts at)
    std::vector<std::pair<unsigned, unsigned>> const pairings{
        {60, 16}, {60, 30}, {16, 30}, {30, 16}, {30, 60}, {40, 30},
        // below the 30 digits the probe direction was drawn at: a drop that rounds it
        {20, 30}, {20, 60}, {20, 16}};

    // the cube root ends in double; the square root, stopped just short of its double root,
    // ends in multiple precision
    for (auto const kind : {Kind::CubeRoot, Kind::SquareRoot})
    {
        TestHomotopy h(kind);
        mpfr const& t_end = h.t_end;

        for (auto const& [before_digits, digits] : pairings)
        {
            BOOST_TEST_CONTEXT((kind == Kind::CubeRoot ? "cube root" : "square root")
                               << ", previous track starts at " << before_digits
                               << " digits, this one at " << digits)
            {
                auto used = h.MakeTracker();
                // preserved, the previous track ends where it started -- elsewhere than this one
                used.PrecisionPreservation(true);
                auto const before = Track(used, before_digits, mpfr(0.5));
                used.PrecisionPreservation(false);
                BOOST_REQUIRE(before.code == bertini::SuccessCode::Success);
                BOOST_REQUIRE_NE(before.precision, digits);

                auto const after_another = Track(used, digits, t_end);

                auto fresh = h.MakeTracker();
                auto const first_ever = Track(fresh, digits, t_end);

                BOOST_CHECK(first_ever.code == bertini::SuccessCode::Success);
                if (kind == Kind::SquareRoot)   // the premise of this half
                    BOOST_CHECK_GT(first_ever.precision, 16u);
                CheckIdentical(after_another, first_ever);
            }
        }
    }
}


/**
The step size is held at the tracker's working precision.

An assignment takes its source's precision (preserve_related_precision), so a step size
assigned from a number of another precision carries that precision into every product it
enters.  Checked after tracks that end in double and in multiple precision, each after a
predecessor at another precision.
*/
BOOST_AUTO_TEST_CASE(the_step_size_is_at_the_working_precision_after_a_track)
{
    using Kind = TestHomotopy::Kind;
    for (auto const kind : {Kind::CubeRoot, Kind::SquareRoot})
    {
        TestHomotopy h(kind);
        mpfr const& t_end = h.t_end;
        auto tracker = h.MakeTracker();
        for (unsigned digits : {60u, 30u, 16u, 40u, 16u, 20u})
        {
            BOOST_TEST_CONTEXT((kind == Kind::CubeRoot ? "cube root" : "square root")
                               << ", starting at " << digits << " digits")
            {
                auto const outcome = Track(tracker, digits, t_end);
                BOOST_REQUIRE(outcome.code == bertini::SuccessCode::Success);
                BOOST_CHECK_EQUAL(outcome.stepsize.precision(), outcome.precision);
            }
        }
    }
}


namespace {

// x^2 + (1-t)x, y^2 + (1-t)y.  At t = 1 the start point (0, 0) is a double root: its refinement
// raises the precision until it gives up with SingularStartPoint, so the track fails during its
// initialization.  At t = 1/2 the point (-1/2, -1/2) is a regular root, tracked to t = 1/4.
struct SingularAndRegularStarts
{
    Var x = Variable::Make("x");
    Var y = Variable::Make("y");
    Var t = Variable::Make("t");
    System sys;

    SingularAndRegularStarts()
    {
        sys.AddFunction(pow(x, 2) + (1 - t) * x);
        sys.AddFunction(pow(y, 2) + (1 - t) * y);
        sys.AddPathVariable(t);
        sys.AddVariableGroup(VariableGroup{x, y});
    }

    bertini::tracking::AMPTracker MakeTracker() const
    {
        using namespace bertini::tracking;
        AMPTracker tracker(sys);
        tracker.Setup(Predictor::RKF45, 1e-5, 1e5, SteppingConfig(), NewtonConfig());
        tracker.PrecisionSetup(AMPConfigFrom(sys));
        return tracker;
    }

    static bertini::SuccessCode TrackSingular(bertini::tracking::AMPTracker const& tracker)
    {
        DefaultPrecision(30);
        Vec<mpfr> start(2);
        start << mpfr(0), mpfr(0);
        Vec<mpfr> end;
        return tracker.TrackPath(end, mpfr(1), mpfr(0), start);
    }

    static bertini::SuccessCode TrackRegular(bertini::tracking::AMPTracker const& tracker)
    {
        DefaultPrecision(16);
        Vec<mpfr> start(2);
        start << mpfr("-0.5"), mpfr("-0.5");
        Vec<mpfr> end;
        return tracker.TrackPath(end, mpfr("0.5"), mpfr("0.25"), start);
    }
};

} // namespace


/**
A precision recorder attached once reports every track.

Each recorder starts afresh when a track begins and stays attached, so after a track whose
start-point refinement raised the precision, the next track, which raises nothing, says so.
*/
BOOST_AUTO_TEST_CASE(a_precision_recorder_attached_once_reports_every_track)
{
    using namespace bertini::tracking;
    SingularAndRegularStarts h;
    auto tracker = h.MakeTracker();
    FirstPrecisionRecorder<AMPTracker> first;
    MinMaxPrecisionRecorder<AMPTracker> min_max;
    tracker.AddObserver(first);
    tracker.AddObserver(min_max);

    BOOST_REQUIRE(SingularAndRegularStarts::TrackSingular(tracker) == bertini::SuccessCode::SingularStartPoint);
    BOOST_CHECK(first.DidPrecisionIncrease());
    BOOST_CHECK_EQUAL(first.StartPrecision(), 30u);
    BOOST_CHECK_EQUAL(min_max.MinPrecision(), 30u);
    BOOST_CHECK_GT(min_max.MaxPrecision(), 30u);

    BOOST_REQUIRE(SingularAndRegularStarts::TrackRegular(tracker) == bertini::SuccessCode::Success);
    BOOST_CHECK(!first.DidPrecisionIncrease());
    BOOST_CHECK_EQUAL(first.StartPrecision(), 16u);
    BOOST_CHECK_EQUAL(min_max.MinPrecision(), 16u);
    BOOST_CHECK_EQUAL(min_max.MaxPrecision(), 16u);
}


/**
A track whose initialization fails is recorded as itself, not as the track before it.

The start-point refinement comes before TrackingStarted, and a failed refinement ends the track
before TrackingStarted is ever sent.  The record starts at the track's Initializing event, so a
recorder that saw an earlier track reports the failed one exactly as a fresh recorder does.
*/
BOOST_AUTO_TEST_CASE(a_track_that_fails_to_initialize_is_recorded_as_itself)
{
    using namespace bertini::tracking;
    SingularAndRegularStarts h;

    auto used = h.MakeTracker();
    FirstPrecisionRecorder<AMPTracker> used_first;
    MinMaxPrecisionRecorder<AMPTracker> used_min_max;
    used.AddObserver(used_first);
    used.AddObserver(used_min_max);
    BOOST_REQUIRE(SingularAndRegularStarts::TrackRegular(used) == bertini::SuccessCode::Success);
    BOOST_REQUIRE(SingularAndRegularStarts::TrackSingular(used) == bertini::SuccessCode::SingularStartPoint);

    auto fresh = h.MakeTracker();
    FirstPrecisionRecorder<AMPTracker> fresh_first;
    MinMaxPrecisionRecorder<AMPTracker> fresh_min_max;
    fresh.AddObserver(fresh_first);
    fresh.AddObserver(fresh_min_max);
    BOOST_REQUIRE(SingularAndRegularStarts::TrackSingular(fresh) == bertini::SuccessCode::SingularStartPoint);

    BOOST_CHECK_EQUAL(used_first.DidPrecisionIncrease(), fresh_first.DidPrecisionIncrease());
    BOOST_CHECK_EQUAL(used_first.StartPrecision(), fresh_first.StartPrecision());
    BOOST_CHECK_EQUAL(used_first.NextPrecision(), fresh_first.NextPrecision());
    BOOST_CHECK(used_first.TimeOfIncrease() == fresh_first.TimeOfIncrease());
    BOOST_CHECK_EQUAL(used_min_max.MinPrecision(), fresh_min_max.MinPrecision());
    BOOST_CHECK_EQUAL(used_min_max.MaxPrecision(), fresh_min_max.MaxPrecision());
}


BOOST_AUTO_TEST_CASE(set_start_precision_clear_restores_default_behavior)
{
    DefaultPrecision(30);
    using namespace bertini::tracking;

    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;
    VariableGroup v{y};
    sys.AddFunction(y - t);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(v);

    auto AMP = bertini::tracking::AMPConfigFrom(sys);
    bertini::tracking::AMPTracker tracker(sys);

    SteppingConfig stepping_preferences;
    NewtonConfig newton_preferences;
    tracker.Setup(Predictor::Euler, 1e-5, 1e5, stepping_preferences, newton_preferences);
    tracker.PrecisionSetup(AMP);
    tracker.SetStartPrecision(50);
    tracker.SetStartPrecision(); // clear override

    Vec<mpfr> y_start(1);
    y_start << mpfr(1);
    mpfr t_start(1), t_end(0);
    Vec<mpfr> y_end;

    auto code = tracker.TrackPath(y_end, t_start, t_end, y_start);
    BOOST_CHECK(code == bertini::SuccessCode::Success);
    BOOST_CHECK(abs(y_end(0) - mpfr(0)) < 1e-5);
}


/**
A bare tracker honours a stop request, with no solver anywhere in sight.

The request is a fact about the process rather than about any one object, and tracking is
where it gets noticed -- between steps, which is the only place it can be noticed without
preempting mid-step.  A user driving a tracker directly is assumed to know what they are
doing; they get the same behaviour a solver gets, because it is the same check.
*/
BOOST_AUTO_TEST_CASE(a_stop_request_stops_a_bare_tracker)
{
    using namespace bertini::tracking;
    using Var = std::shared_ptr<bertini::node::Variable>;
    using Variable = bertini::node::Variable;
    using bertini::System;
    using bertini::VariableGroup;

    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;
    sys.AddFunction(y - t);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(VariableGroup{y});

    AMPTracker tracker(sys);
    tracker.Setup(Predictor::Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig());
    tracker.PrecisionSetup(bertini::tracking::AMPConfigFrom(sys));

    bertini::Vec<bertini::complex_mp> start(1), result;
    start << bertini::complex_mp(1);

    bertini::ScopedStopRequest guard;
    bertini::RequestStop();

    auto const code = tracker.TrackPath(result, bertini::complex_mp(1),
                                        bertini::complex_mp("0.1"), start);

    BOOST_CHECK(code == bertini::SuccessCode::ExternallyTerminated);

    // and with the request withdrawn, the identical track succeeds
    bertini::ClearStopRequest();
    auto const after = tracker.TrackPath(result, bertini::complex_mp(1),
                                         bertini::complex_mp("0.1"), start);
    BOOST_CHECK(after == bertini::SuccessCode::Success);
}

/**
A bare tracker honours a wall-clock deadline with no solver anywhere.  A deadline already in
the past stops the track at its first step boundary with WallClockLimitReached and the
tracker still says where it was -- at the start, having taken no steps.  Clearing the
deadline makes the identical track succeed, and a generous deadline does not bite.
*/
BOOST_AUTO_TEST_CASE(a_wall_clock_deadline_stops_a_bare_tracker)
{
    using namespace bertini::tracking;
    using Var = std::shared_ptr<bertini::node::Variable>;
    using Variable = bertini::node::Variable;
    using bertini::System;
    using bertini::VariableGroup;

    Var y = Variable::Make("y");
    Var t = Variable::Make("t");

    System sys;
    sys.AddFunction(y - t);
    sys.AddPathVariable(t);
    sys.AddVariableGroup(VariableGroup{y});

    AMPTracker tracker(sys);
    tracker.Setup(Predictor::Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig());
    tracker.PrecisionSetup(bertini::tracking::AMPConfigFrom(sys));

    bertini::Vec<bertini::complex_mp> start(1), result;
    start << bertini::complex_mp(1);

    BOOST_CHECK(!tracker.MaxWallClockTime().has_value());
    tracker.SetMaxWallClockTime(std::chrono::steady_clock::now() - std::chrono::seconds(1));
    BOOST_CHECK(tracker.MaxWallClockTime().has_value());

    auto const code = tracker.TrackPath(result, bertini::complex_mp(1),
                                        bertini::complex_mp("0.1"), start);
    BOOST_CHECK(code == bertini::SuccessCode::WallClockLimitReached);
    BOOST_CHECK_EQUAL(tracker.NumTotalStepsTaken(), 0u);
    BOOST_CHECK(abs(tracker.CurrentTime() - bertini::complex_mp(1)) < 1e-30);
    BOOST_CHECK((tracker.CurrentPoint() - start).norm() < 1e-30);

    tracker.ClearMaxWallClockTime();
    BOOST_CHECK(!tracker.MaxWallClockTime().has_value());
    auto const after = tracker.TrackPath(result, bertini::complex_mp(1),
                                         bertini::complex_mp("0.1"), start);
    BOOST_CHECK(after == bertini::SuccessCode::Success);

    tracker.SetMaxWallClockDuration(std::chrono::hours(1));
    auto const generous = tracker.TrackPath(result, bertini::complex_mp(1),
                                            bertini::complex_mp("0.1"), start);
    BOOST_CHECK(generous == bertini::SuccessCode::Success);
    tracker.ClearMaxWallClockTime();
}


/**
SuccessCode::NeverStarted is a DEFAULT, never a verdict, and callers depend on that: a solve
records nothing for a path whose code reads NeverStarted, counts it as unreached, and reports
itself cut short on account of it.  A tracker that returned it would make a path that WAS
attempted look like one nobody got to -- silently, since the value is a legal SuccessCode.

Nothing returns it today.  This pins that by ending a track every way a track can end, and
checking that no route produces it.  The codes themselves are not asserted one by one, which
would pin the tracker's choice of diagnosis rather than the invariant; what is asserted is the
invariant, plus that the scenarios really did end in several different ways, so the test cannot
pass by doing nothing.
*/
BOOST_AUTO_TEST_CASE(a_tracker_never_returns_never_started)
{
    using namespace bertini::tracking;
    using Variable = bertini::node::Variable;
    using bertini::System;
    using bertini::VariableGroup;
    using bertini::SuccessCode;
    using Cmp = bertini::complex_mp;

    // One homotopy in y and t, built fresh per scenario (a tracker holds a reference to its
    // system, and these systems differ).
    auto one_variable_system = [](int which)
    {
        auto y = Variable::Make("y");
        auto t = Variable::Make("t");
        System sys;
        if (which == 0)
            sys.AddFunction(y - t);              // the well-behaved path, y(t) = t
        else if (which == 1)
            sys.AddFunction(y*y - t + 1);        // y=0 at t=1 is a SINGULAR start (dy = 2y = 0)
        else
            sys.AddFunction(y*t - 1);            // y(t) = 1/t, off to infinity as t -> 0
        sys.AddPathVariable(t);
        sys.AddVariableGroup(VariableGroup{y});
        return sys;
    };

    auto track = [&one_variable_system](int which, Cmp const& start_value, Cmp const& endtime,
                                        std::function<void(AMPTracker&)> const& configure)
    {
        auto sys = one_variable_system(which);
        AMPTracker tracker(sys);
        tracker.Setup(Predictor::RKF45, 1e-5, 1e5, SteppingConfig(), NewtonConfig());
        tracker.PrecisionSetup(bertini::tracking::AMPConfigFrom(sys));
        if (configure)
            configure(tracker);

        bertini::Vec<Cmp> start(1), result;
        start << start_value;
        return tracker.TrackPath(result, Cmp(1), endtime, start);
    };

    std::vector<SuccessCode> observed;

    // a track that simply works
    observed.push_back(track(0, Cmp(1), Cmp("0.1"), nullptr));

    // out of steps
    observed.push_back(track(0, Cmp(1), Cmp("0.1"), [](AMPTracker& tr){
        auto stepping = tr.Get<SteppingConfig>();
        stepping.max_num_steps = 1;
        tr.Set<SteppingConfig>(stepping);
    }));

    // the step size floor is already above the step it wants to take
    observed.push_back(track(0, Cmp(1), Cmp("0.1"), [](AMPTracker& tr){
        auto stepping = tr.Get<SteppingConfig>();
        stepping.min_step_size = 1.0;
        tr.Set<SteppingConfig>(stepping);
    }));

    // asked for more digits than the precision ceiling allows
    observed.push_back(track(0, Cmp(1), Cmp("0.1"), [](AMPTracker& tr){
        tr.SetTrackingTolerance(1e-60);
        auto amp = tr.Get<AdaptiveMultiplePrecisionConfig>();
        amp.maximum_precision = 20;
        tr.Set<AdaptiveMultiplePrecisionConfig>(amp);
    }));

    // a singular start point, refused by the initial refinement
    observed.push_back(track(1, Cmp(0), Cmp("0.1"), nullptr));

    // a path running off to infinity, under ADAPTIVE precision: it does not get to the
    // truncation threshold, because chasing y = 1/t toward the pole escalates precision until
    // the ceiling stops it first.  Truncation itself is covered below, with a fixed-precision
    // tracker, which cannot escalate.
    observed.push_back(track(2, Cmp(1), Cmp(0), [](AMPTracker& tr){
        tr.SetInfiniteTruncationTolerance(10.0);
    }));

    // somebody asked us to stop
    {
        bertini::ScopedStopRequest guard;
        bertini::RequestStop();
        observed.push_back(track(0, Cmp(1), Cmp("0.1"), nullptr));
    }

    // the wall-clock deadline is already behind us
    observed.push_back(track(0, Cmp(1), Cmp("0.1"), [](AMPTracker& tr){
        tr.SetMaxWallClockTime(std::chrono::steady_clock::now() - std::chrono::seconds(1));
    }));

    // and the other initialization path: a fixed-precision tracker, working and out of steps
    for (unsigned max_steps : {100000u, 1u})
    {
        auto y = Variable::Make("y");
        auto t = Variable::Make("t");
        System sys;
        sys.AddFunction(y - t);
        sys.AddPathVariable(t);
        sys.AddVariableGroup(VariableGroup{y});

        DoublePrecisionTracker tracker(sys);
        tracker.Setup(Predictor::RKF45, double(1e-5), double(1e5), SteppingConfig(), NewtonConfig());
        auto stepping = tracker.Get<SteppingConfig>();
        stepping.max_num_steps = max_steps;
        tracker.Set<SteppingConfig>(stepping);

        bertini::Vec<bertini::complex_dbl> start(1), result;
        start << bertini::complex_dbl(1);
        observed.push_back(tracker.TrackPath(result, bertini::complex_dbl(1),
                                             bertini::complex_dbl(0.1), start));
    }

    // truncation: y = 1/t past the threshold, in fixed precision so nothing escalates
    {
        auto y = Variable::Make("y");
        auto t = Variable::Make("t");
        System sys;
        sys.AddFunction(y*t - 1);
        sys.AddPathVariable(t);
        sys.AddVariableGroup(VariableGroup{y});

        DoublePrecisionTracker tracker(sys);
        tracker.Setup(Predictor::RKF45, double(1e-5), double(1e5), SteppingConfig(), NewtonConfig());
        tracker.SetInfiniteTruncationTolerance(10.0);

        bertini::Vec<bertini::complex_dbl> start(1), result;
        start << bertini::complex_dbl(1);
        observed.push_back(tracker.TrackPath(result, bertini::complex_dbl(1),
                                             bertini::complex_dbl(0), start));
    }

    for (auto code : observed)
    {
        BOOST_TEST_MESSAGE("ended as " << bertini::records::CanonicalName(code));
        BOOST_CHECK(code != SuccessCode::NeverStarted);
    }

    // the scenarios really did end differently, so the checks above mean something
    std::set<SuccessCode> distinct(observed.begin(), observed.end());
    BOOST_CHECK_GE(distinct.size(), 4u);
    BOOST_CHECK(distinct.count(SuccessCode::Success) == 1);          // at least one worked
    BOOST_CHECK(distinct.size() > 1);                                // and at least one did not
}

BOOST_AUTO_TEST_SUITE_END()
