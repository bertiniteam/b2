//This file is part of Bertini 2.0.
//
//test/tracking_basics/path_observers.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/tracking_basics/path_observers.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/tracking_basics/path_observers.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

//test/tracking_basics/path_observers.cpp
//


#include <boost/test/unit_test.hpp>

#include "bertini2/tracking/amp_tracker.hpp"
#include "bertini2/tracking/observers.hpp"



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




BOOST_AUTO_TEST_SUITE(AMP_tracker_basics)




BOOST_AUTO_TEST_CASE(accumulate_single_path_square_root)
{
	mpfr_float::default_precision(16);
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
					newton_preferences);

	tracker.AMPSetup(AMP);

	mpfr t_start(1);
	mpfr t_end(0);
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	SuccessCode tracking_success;

	AMPPathAccumulator<AMPTracker> path_accumulator;
	PrecisionAccumulator<AMPTracker> precision_accumulator;

	tracker.AddObserver(&path_accumulator);
	tracker.AddObserver(&precision_accumulator);

	start_point << mpfr(1), mpfr(1);
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

	for (auto iter : path_accumulator.Path())
		std::cout << iter << "\n";

	for (auto iter : precision_accumulator.Precisions())
		std::cout << iter << " ";
	
}




BOOST_AUTO_TEST_CASE(some_other_thing_square_root)
{
	mpfr_float::default_precision(16);
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
					newton_preferences);

	tracker.AMPSetup(AMP);

	mpfr t_start(1);
	mpfr t_end(0);
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	SuccessCode tracking_success;

	GoryDetailLogger<AMPTracker> tons_of_detail;

	tracker.AddObserver(&tons_of_detail);

	start_point << mpfr(1), mpfr(1);
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

}


BOOST_AUTO_TEST_SUITE_END()



