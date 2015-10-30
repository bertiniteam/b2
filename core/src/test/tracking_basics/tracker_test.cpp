//This file is part of Bertini 2.0.
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

//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015

//tracker_test.cpp
//


#include <boost/test/unit_test.hpp>

#include "tracking/tracker.hpp"

using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;

extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;



BOOST_AUTO_TEST_SUITE(tracker_basics)



BOOST_AUTO_TEST_CASE(tracker_track_linear)
{
	mpfr_float::default_precision(30);
	using namespace bertini::tracking;

	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{y};

	sys.AddFunction(y-t);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("0.0");
	
	Vec<mpfr> y_start(1);
	y_start << mpfr("1.0");

	Vec<mpfr> y_end;

	tracker.TrackPath(y_end,
	                  t_start, t_end, y_start);

	BOOST_CHECK_EQUAL(y_end.size(),1);
	BOOST_CHECK(abs(y_end(0)-mpfr("0.0")) < 1e-5);

}







BOOST_AUTO_TEST_CASE(tracker_track_quadratic)
{
	mpfr_float::default_precision(30);
	using namespace bertini::tracking;

	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{y};

	sys.AddFunction(y-pow(t,2));
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("-1.0");
	
	Vec<mpfr> y_start(1);
	y_start << mpfr("1.0");

	Vec<mpfr> y_end;

	tracker.TrackPath(y_end,
	                  t_start, t_end, y_start);

	BOOST_CHECK_EQUAL(y_end.size(),1);
	BOOST_CHECK(abs(y_end(0)-mpfr("1.0")) < 1e-5);
}



BOOST_AUTO_TEST_CASE(tracker_track_decic)
{
	mpfr_float::default_precision(30);
	using namespace bertini::tracking;

	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{y};

	sys.AddFunction(y-pow(t,10));
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              	mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("-2.0");
	
	Vec<mpfr> y_start(1);
	y_start << mpfr("1.0");

	Vec<mpfr> y_end;

	tracker.TrackPath(y_end,
	                  t_start, t_end, y_start);

	BOOST_CHECK_EQUAL(y_end.size(),1);
	BOOST_CHECK(abs(y_end(0)-mpfr("1024.0")) < 1e-5);

}




BOOST_AUTO_TEST_CASE(tracker_track_square_root)
{
	mpfr_float::default_precision(30);
	using namespace bertini::tracking;

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(x-t);
	sys.AddFunction(pow(y,2)-x);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              	mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("0");
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	SuccessCode tracking_success;


	start_point << mpfr("1"), mpfr("1");
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success==SuccessCode::Success);
	BOOST_CHECK_EQUAL(end_point.size(),2);
	BOOST_CHECK(abs(end_point(0)-mpfr("0.0")) < 1e-5);
	BOOST_CHECK(abs(end_point(1)-mpfr("0.0")) < 1e-5);


	start_point << mpfr("1"), mpfr("-1");
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success==SuccessCode::Success);
	BOOST_CHECK_EQUAL(end_point.size(),2);
	BOOST_CHECK(abs(end_point(0)-mpfr("0.0")) < 1e-5);
	BOOST_CHECK(abs(end_point(1)-mpfr("0.0")) < 1e-5);

	t_start = mpfr("-1.0");
	start_point << mpfr("-1"), mpfr("0","-1");
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success==SuccessCode::Success);
	BOOST_CHECK_EQUAL(end_point.size(),2);
	BOOST_CHECK(abs(end_point(0)-mpfr("0.0")) < 1e-5);
	BOOST_CHECK(abs(end_point(1)-mpfr("0.0")) < 1e-5);

	start_point << mpfr("-1"), mpfr("0","1");
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success==SuccessCode::Success);
	BOOST_CHECK_EQUAL(end_point.size(),2);
	BOOST_CHECK(abs(end_point(0)-mpfr("0.0")) < 1e-5);
	BOOST_CHECK(abs(end_point(1)-mpfr("0.0")) < 1e-5);
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
BOOST_AUTO_TEST_CASE(tracker_doesnt_start_from_singular_start_point)
{
	mpfr_float::default_precision(30);
	using namespace bertini::tracking;

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(pow(x,2) + (1-t)*x);
	sys.AddFunction(pow(y,2) + (1-t)*y);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              	mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("0");
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	SuccessCode tracking_success;


	start_point << mpfr("0"), mpfr("0");
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success==SuccessCode::SingularStartPoint);
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
BOOST_AUTO_TEST_CASE(simple_nonhomogeneous_system_trackable_initialprecision16)
{
	mpfr_float::default_precision(16);
	using namespace bertini::tracking;

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(pow(x,2) + (1-t)*x - 1);
	sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              	mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("0");
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	SuccessCode tracking_success;


	start_point << mpfr("1"), mpfr("1.414");
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success==SuccessCode::Success);
	BOOST_CHECK_EQUAL(end_point.size(),2);
	BOOST_CHECK(abs(end_point(0)-mpfr("6.180339887498949e-01")) < 1e-5);
	BOOST_CHECK(abs(end_point(1)-mpfr("1.138564265110173e+00")) < 1e-5);
}




BOOST_AUTO_TEST_CASE(simple_nonhomogeneous_system_trackable_initialprecision30)
{
	mpfr_float::default_precision(30);
	using namespace bertini::tracking;

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(pow(x,2) + (1-t)*x - 1);
	sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              	mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("0");
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	SuccessCode tracking_success;


	start_point << mpfr("1"), mpfr("1.414");
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success==SuccessCode::Success);
	BOOST_CHECK_EQUAL(end_point.size(),2);
	BOOST_CHECK(abs(end_point(0)-mpfr("6.180339887498949e-01")) < 1e-5);
	BOOST_CHECK(abs(end_point(1)-mpfr("1.138564265110173e+00")) < 1e-5);
}



BOOST_AUTO_TEST_CASE(simple_nonhomogeneous_system_trackable_initialprecision100)
{
	mpfr_float::default_precision(100);
	using namespace bertini::tracking;

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(pow(x,2) + (1-t)*x - 1);
	sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              	mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1.0");
	mpfr t_end("0");
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	SuccessCode tracking_success;


	start_point << mpfr("1"), mpfr("1.414");
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success==SuccessCode::Success);
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
BOOST_AUTO_TEST_CASE(tracker_fails_with_singularity_on_path)
{
	mpfr_float::default_precision(30);
	using namespace bertini::tracking;

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	auto s = -1*(1-t) + 1*t;

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(pow(x,2) - s);
	sys.AddFunction(pow(y,2) - s);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	config::Stepping stepping_preferences;
	config::Newton newton_preferences;


	tracker.Setup(config::Predictor::Euler,
	              	mpfr_float("1e-5"),
					mpfr_float("1e5"),
					stepping_preferences,
					newton_preferences,
					AMP);

	mpfr t_start("1");
	mpfr t_end("0");
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	start_point << mpfr("1"), mpfr("1");
	SuccessCode tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	BOOST_CHECK(tracking_success!=SuccessCode::Success);
}



BOOST_AUTO_TEST_SUITE_END()




