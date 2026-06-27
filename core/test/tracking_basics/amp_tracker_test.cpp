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

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire




#include <boost/test/unit_test.hpp>
#include "bertini2/system/start_systems.hpp"
#include "bertini2/trackers/tracker.hpp"



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

	

	auto TD = bertini::start_system::TotalDegree(sys);
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



std::vector<Vec<mpfr> > track_total_degree(bertini::tracking::AMPTracker const& tracker, bertini::start_system::TotalDegree const& TD)
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

	

	auto TD = bertini::start_system::TotalDegree(sys);
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


BOOST_AUTO_TEST_SUITE_END()




