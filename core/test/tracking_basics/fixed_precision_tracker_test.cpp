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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire



#include <boost/test/unit_test.hpp>
#include "bertini2/system/start_systems.hpp"
#include "bertini2/trackers/fixed_precision_tracker.hpp"
#include "bertini2/trackers/observers.hpp"



extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;



BOOST_AUTO_TEST_SUITE(fixed_precision_tracker_basics)

using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;
using bertini::MakeVariable;

using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;
using bertini::DefaultPrecision;
BOOST_AUTO_TEST_CASE(double_tracker_track_linear)
{
	using bertini::operator<<;
	DefaultPrecision(100);
	using namespace bertini::tracking;

	Var y = MakeVariable("y");
	Var t = MakeVariable("t");

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

	
	dbl t_start(1);
	dbl t_end(0);
	
	Vec<dbl> y_start(1);
	y_start << dbl(1);

	Vec<dbl> y_end;

	auto obs = GoryDetailLogger<DoublePrecisionTracker>();
		tracker.AddObserver(&obs);

	auto code = tracker.TrackPath(y_end, t_start, t_end, y_start);
	BOOST_CHECK(code==bertini::SuccessCode::Success);

	BOOST_CHECK_EQUAL(y_end.size(),1);
	BOOST_CHECK(abs(y_end(0)-dbl(0)) < 1e-5);

}
	
BOOST_AUTO_TEST_CASE(multiple_100_tracker_track_linear)
{
	DefaultPrecision(100);
	using namespace bertini::tracking;

	Var y = MakeVariable("y");
	Var t = MakeVariable("t");

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
	tracker.AddObserver(&tons_of_detail);

	
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




