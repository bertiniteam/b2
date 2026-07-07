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






BOOST_AUTO_TEST_SUITE_END()




