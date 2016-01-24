//This file is part of Bertini 2.0.
//
//powerseries_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//powerseries_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with euler_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  powerseries_test.cpp
//
//  copyright 2015
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Fall 2015
#include <iostream>
#include <typeinfo>
#include <boost/test/unit_test.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include "bertini2/limbo.hpp"
#include "bertini2/mpfr_complex.hpp"
#include <deque>
#include "bertini2/tracking/cauchy_endgame.hpp"
#include "bertini2/start_system.hpp"



using Variable = bertini::node::Variable;


using Var = std::shared_ptr<Variable>;
 
using VariableGroup = bertini::VariableGroup;

using System = bertini::System;
// using dbl = std::complex<double>;
 using mpfr = bertini::complex;
// using mpfr_float = boost::multiprecision::mpfr_float;

using dbl = std::complex<double>;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;

extern double threshold_clearance_d;
extern boost::multiprecision::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;


BOOST_AUTO_TEST_SUITE(cauchy_endgame_class_basics)

BOOST_AUTO_TEST_CASE(circle_track_d_for_cauchy_class_test)
{
	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); //f(x) = (x-1)^3 + t*0, need t*0 for derivative calculation. 

	sys.AddFunction((x - mpfr(1))*(1-t) + (x + mpfr(1))*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);
	
	bertini::tracking::config::Stepping stepping_preferences;
	bertini::tracking::config::Newton newton_preferences;

	tracker.Setup(bertini::tracking::config::Predictor::Euler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.

	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr(".1");
	times.push_back(time);
	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
	samples.push_back(sample);

	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

	auto first_track_sample =  My_Endgame.CircleTrack(time,sample);

	// std::cout << "first track sample is " << first_track_sample << '\n';

	BOOST_CHECK((first_track_sample - sample).norm() < My_Endgame.GetTrackToleranceDuringEndgame());

} // end circle_track_d_for_cauchy_class_test

BOOST_AUTO_TEST_CASE(circle_track_mp_cycle_num_greater_than_1_for_cauchy_class_test)
{
	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x - mpfr("1"),2)*(1-t) + (pow(x,2) + mpfr("1"))*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);
	
	bertini::tracking::config::Stepping stepping_preferences;
	bertini::tracking::config::Newton newton_preferences;

	tracker.Setup(bertini::tracking::config::Predictor::Euler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.

	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr("0.1");
	times.push_back(time);
	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
	samples.push_back(sample);

	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);


	auto first_track_sample =  My_Endgame.CircleTrack(time,sample);

	// std::cout << "first track sample is " << first_track_sample << '\n';

	BOOST_CHECK((first_track_sample - sample).norm() > My_Endgame.GetTrackToleranceDuringEndgame());

	auto second_track_sample =  My_Endgame.CircleTrack(time,first_track_sample);

	// std::cout << "second track sample is " << second_track_sample << '\n';

	BOOST_CHECK((second_track_sample - sample).norm() < My_Endgame.GetTrackToleranceDuringEndgame());
	
} // end circle_track_d_for_cauchy_class_test



BOOST_AUTO_TEST_SUITE_END()
