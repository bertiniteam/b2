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
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Tim Hodges, Colorado State University


//  cauchy_class_test.cpp

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
#include "bertini2/tracking/amp_cauchy_endgame.hpp"
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


BOOST_AUTO_TEST_SUITE(amp_cauchy_endgame)

using namespace bertini::tracking;
using namespace bertini::tracking::endgame;

BOOST_AUTO_TEST_CASE(circle_track_d)
{
	/*
		In this test we take a univariate polynomial with one solution and ensure that our CircleTrack function returns to the same position. 
		This is tested by checking that our cauchy_samples has the front and end with the same value roughly.
	*/
	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
	SampCont<mpfr> cauchy_samples; //samples are space values that may be a vector of complex numbers.


	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr("0.1");
	cauchy_times.push_back(time);
	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
	cauchy_samples.push_back(sample);

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	auto first_track_success =  my_endgame.CircleTrack(time,sample);

	BOOST_CHECK((my_endgame.GetCauchySamples<mpfr>().back() - sample).norm() < my_endgame.Tolerances().newton_during_endgame);

} // end circle_track_d

BOOST_AUTO_TEST_CASE(circle_track_mp_cycle_num_greater_than_1)
{
	/*
		In this test we take a univariate polynomial with one solution of multiplicity 2. We need to ensure that our CircleTrack function 
		returns to the same position after two calls. The first call should land us on the second path, and the second call should land us back 
		where we started. This is tested by checking that cauchy_samples has the front and end with the same value roughly.
	*/
	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
	SampCont<mpfr> cauchy_samples; //samples are space values that may be a vector of complex numbers.



	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr("0.1");
	cauchy_times.push_back(time);
	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
	cauchy_samples.push_back(sample);

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);


	auto tracking_success =  my_endgame.CircleTrack(time,sample);

	auto first_track_sample = my_endgame.GetCauchySamples<mpfr>().back();
	// std::cout << "first track sample is " << first_track_sample << '\n';

	BOOST_CHECK((first_track_sample - sample).norm() > my_endgame.Tolerances().newton_during_endgame);

	tracking_success =  my_endgame.CircleTrack(time,first_track_sample);

	auto second_track_sample = my_endgame.GetCauchySamples<mpfr>().back();

	// std::cout << "second track sample is " << second_track_sample << '\n';

	BOOST_CHECK((second_track_sample - sample).norm() < my_endgame.Tolerances().newton_during_endgame);
	
} // end circle_track_mp_cycle_num_greater_than_1


BOOST_AUTO_TEST_CASE(compute_c_over_k_dbl_for_cauchy_class)
{
	/*
		In this test we take a univariate polynomial with one solution of multiplicity. For our first power series approximation we need to 
		make sure that we have a stablization of the cycle number. This is achieved by computing an approximation for C over K which is described 
		on page 53 of the Berini Book, Numerically Solving Polynomials Systems with Bertini. 
		This test case is ensuring no calculations are greatly altered, so our test condition is to check against a computed value that has been already known. 
	*/	

	mpfr_float::default_precision(16);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	sys.AddFunction(pow(x-1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> pseg_times;
	SampCont<mpfr> pseg_samples;
	

	mpfr time(1);
	Vec<mpfr> sample(1);

	time = mpfr(".1"); // x = .1
	pseg_times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	pseg_samples.push_back(sample);


	time = mpfr(".05"); // x = .1/2 = .05
	pseg_times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	pseg_samples.push_back(sample);


	time = mpfr(".025"); // x = .05/2 = .025
	pseg_times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	pseg_samples.push_back(sample);


	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> endgame_settings;

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker,security_settings,endgame_settings);
	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto first_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	// std::cout << "first c over k is " << first_c_over_k << '\n';
	// std::cout << "first diff " << abs(first_c_over_k - mpfr_float("1.12917")) << '\n';
	// std::cout << "Tol is " << my_endgame.GetTrackToleranceDuringEndgame() << '\n';

	 BOOST_CHECK( abs(first_c_over_k - mpfr_float("1.12917")) < 1e-5 ); 

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto second_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	// std::cout << "second_c_over_k is " << second_c_over_k << '\n';
	// std::cout << "second diff is " << abs(second_c_over_k - mpfr_float("1.05888")) << '\n';
	BOOST_CHECK( abs(second_c_over_k - mpfr_float("1.05888")) <  1e-5); 

} // end compute c over k dbl for cauchy class 

BOOST_AUTO_TEST_CASE(compute_c_over_k_mp_for_cauchy_class)
{
	/*
		In this test we take a univariate polynomial with one solution of multiplicity. For our first power series approximation we need to 
		make sure that we have a stablization of the cycle number. This is achieved by computing an approximation for C over K which is described 
		on page 53 of the Berini Book, Numerically Solving Polynomials Systems with Bertini. 
		This test case is ensuring no calculations are greatly altered, so our test condition is to check against a computed value that has been already known. 
	*/	
	
	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	sys.AddFunction(pow(x-1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> pseg_times;
	SampCont<mpfr> pseg_samples;
	

	mpfr time(1);
	Vec<mpfr> sample(1);

	time = mpfr(".1"); // x = .1
	pseg_times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	pseg_samples.push_back(sample);


	time = mpfr(".05"); // x = .1/2 = .05
	pseg_times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	pseg_samples.push_back(sample);


	time = mpfr(".025"); // x = .05/2 = .025
	pseg_times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	pseg_samples.push_back(sample);


	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> endgame_settings;

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker,security_settings,endgame_settings);
	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto first_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	BOOST_CHECK( abs(first_c_over_k - mpfr_float("1.12917")) < 1e-5 ); 

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto second_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	BOOST_CHECK( abs(second_c_over_k - mpfr_float("1.05888")) <  1e-5); 

} // end compute c over k mp for cauchy class 

BOOST_AUTO_TEST_CASE(checking_agreement_function_dbl_for_cauchy_class)
{

	/*
		After we have computed C over K approximations we have to make sure that we lie within some agreed terms. This means our agreements need
		to be with some ratio of eachother. (Minimum is usually .75), and we need to agree so many times (usually 3). This test case checks for this. 

	*/	

	mpfr_float::default_precision(16);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t");
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> pseg_times; 
	SampCont<mpfr> pseg_samples;
	TimeCont<mpfr> c_over_k_array;

	mpfr time(1);
	Vec<mpfr> sample(1);

	time = mpfr(".1"); // x = .1
	pseg_times.push_back(time);
	sample << mpfr("1.100000000000000e+00", "6.244997998398396e-01"); // f(.1) = 1.100000000000000e+00 6.244997998398396e-01
	pseg_samples.push_back(sample);


	time = mpfr(".05"); // x = .1/2 = .05
	pseg_times.push_back(time);
	sample << mpfr("1.125000000000000e+00", "4.841229182759271e-01"); //f(.05) = 1.125000000000000e+00 4.841229182759271e-01
	pseg_samples.push_back(sample);


	time = mpfr(".025"); // x = .05/2 = .025
	pseg_times.push_back(time);
	sample << mpfr("1.123854702920065e+00", "3.736283818773165e-01"); // f(.025) = 1.123854702920065e+00 3.736283818773165e-01
	pseg_samples.push_back(sample);


	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> endgame_settings;

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker,security_settings,endgame_settings);
	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto first_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	c_over_k_array.push_back(first_c_over_k);

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << mpfr("1.111721780135302e+00", "2.886596646224579e-01"); // f(.0125) = 1.111721780135302e+00 2.886596646224579e-01
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto second_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	c_over_k_array.push_back(second_c_over_k);



	//Setting up a new sample for approximation.
	time = mpfr(".00625"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << mpfr("1.096071609421043e+00", "2.237005761359081e-01"); // f(.00625) = 1.096071609421043e+00 2.237005761359081e-01
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto third_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	c_over_k_array.push_back(third_c_over_k);

	auto stabilized = my_endgame.CheckForCOverKStabilization(c_over_k_array);

	BOOST_CHECK(stabilized == true);

} // end check agreement function dbl for cauchy 

BOOST_AUTO_TEST_CASE(checking_agreement_function_mp_for_cauchy_class)
{
	/*
		After we have computed C over K approximations we have to make sure that we lie within some agreed terms. This means our agreements need
		to be with some ratio of eachother. (Minimum is usually .75), and we need to agree so many times (usually 3). This test case checks for this. 
		
	*/	
	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t");
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> pseg_times;
	SampCont<mpfr> pseg_samples;
	TimeCont<mpfr> c_over_k_array;

	mpfr time;
	Vec<mpfr> sample(1);

	time = mpfr(".1"); // x = .1
	pseg_times.push_back(time);
	sample << mpfr("1.100000000000000e+00", "6.244997998398396e-01"); // f(.1) = 1.100000000000000e+00 6.244997998398396e-01
	pseg_samples.push_back(sample);


	time = mpfr(".05"); // x = .1/2 = .05
	pseg_times.push_back(time);
	sample << mpfr("1.125000000000000e+00", "4.841229182759271e-01"); //f(.05) = 1.125000000000000e+00 4.841229182759271e-01
	pseg_samples.push_back(sample);


	time = mpfr(".025"); // x = .05/2 = .025
	pseg_times.push_back(time);
	sample << mpfr("1.123854702920065e+00", "3.736283818773165e-01"); // f(.025) = 1.123854702920065e+00 3.736283818773165e-01
	pseg_samples.push_back(sample);


	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> endgame_settings;

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker,security_settings,endgame_settings);
	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto first_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	c_over_k_array.push_back(first_c_over_k);

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << mpfr("1.111721780135302e+00", "2.886596646224579e-01"); // f(.0125) = 1.111721780135302e+00 2.886596646224579e-01
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto second_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	c_over_k_array.push_back(second_c_over_k);



	//Setting up a new sample for approximation.
	time = mpfr(".00625"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << mpfr("1.096071609421043e+00", "2.237005761359081e-01"); // f(.00625) = 1.096071609421043e+00 2.237005761359081e-01
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto third_c_over_k = my_endgame.ComputeCOverK<mpfr>();

	c_over_k_array.push_back(third_c_over_k);

	BOOST_CHECK(my_endgame.CheckForCOverKStabilization(c_over_k_array) == true);

} // end check agreement function mp for cauchy 



BOOST_AUTO_TEST_CASE(check_closed_loop_for_cycle_num_1)
{
	/*
		In this test case we have a solution we are trying to find with multiplicity 1. If we call CircleTrack we should have a closed loop. 
		This test case is checking to make sure our CheckClosedLoop function returns that we have actually closed. 
	*/
	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
	SampCont<mpfr> cauchy_samples; //samples are space values that may be a vector of complex numbers.


	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr(".1");
	cauchy_times.push_back(time);
	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
	cauchy_samples.push_back(sample);


	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	my_endgame.SetCauchySamples(cauchy_samples);
	my_endgame.SetCauchyTimes(cauchy_times);


	auto tracking_success =  my_endgame.CircleTrack(time,sample);
	BOOST_CHECK(my_endgame.CheckClosedLoop<mpfr>() == true);

} // end check closed loop if cycle num is 1 for cauchy class test

BOOST_AUTO_TEST_CASE(check_closed_loop_for_cycle_num_greater_than_1)
{
	/*
		In this test case we have a solution we are trying to find with multiplicity 2. If we call CircleTrack we should not have a closed loop. 
		
		This test case is checking to make sure our CheckClosedLoop function returns that we have not closed.
	*/

	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);


	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
	SampCont<mpfr> cauchy_samples; //samples are space values that may be a vector of complex numbers.



	mpfr time("0.1");
	cauchy_times.push_back(time);

	Vec<mpfr> sample(1);
	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
	cauchy_samples.push_back(sample);


	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	my_endgame.SetCauchySamples(cauchy_samples);
	my_endgame.SetCauchyTimes(cauchy_times);

	auto tracking_success =  my_endgame.CircleTrack(time,sample);
	BOOST_CHECK(my_endgame.CheckClosedLoop<mpfr>() == false);

	tracking_success =  my_endgame.CircleTrack(my_endgame.GetCauchyTimes<mpfr>().back(),my_endgame.GetCauchySamples<mpfr>().back());
	BOOST_CHECK(my_endgame.CheckClosedLoop<mpfr>() == true);
	
} // end check closed loop if cycle num is greater than 1 for cauchy class test


BOOST_AUTO_TEST_CASE(compare_cauchy_ratios)
{
	
	//  After we have stabilized our estimation of c over k, we can start to perform cauchy loops to see if the minimum and maximum norms
	// of sample points are within a tolerance. If so, we can start actually using the Cauchy Integral Formula to extrapolate approximations 
	// at the origin. 

	// In this example, since the cycle number is 1, we should see that our RatioEGOperatingZoneTest function says that we are good to proceed 
	// in using the Cauchy Integral Formula. 
	
	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);


	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
	SampCont<mpfr> cauchy_samples; //samples are space values that may be a vector of complex numbers.


	mpfr time(".1");
	Vec<mpfr> sample(1);

	cauchy_times.push_back(time);
	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
	cauchy_samples.push_back(sample);

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	auto tracking_success =  my_endgame.CircleTrack(time,sample);
	BOOST_CHECK(my_endgame.RatioEGOperatingZoneTest<mpfr>() == true);

} // end compare cauchy ratios for cauchy class test

BOOST_AUTO_TEST_CASE(compare_cauchy_ratios_cycle_num_greater_than_1)
{
	/*
		After we have stabilized our estimation of c over k, we can start to perform cauchy loops to see if the minimum and maximum norms
		of sample points are within a tolerance. If so, we can start actually using the Cauchy Integral Formula to extrapolate approximations 
		at the origin. 

	*/
	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
	SampCont<mpfr> cauchy_samples; //samples are space values that may be a vector of complex numbers.



	mpfr time("0.1");
	Vec<mpfr> sample(1);

	cauchy_times.push_back(time);
	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
	cauchy_samples.push_back(sample);


	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);
	auto tracking_success =  my_endgame.CircleTrack(time,sample);
	BOOST_CHECK(my_endgame.RatioEGOperatingZoneTest<mpfr>() == true);

} // end compare cauchy ratios for cycle num greater than 1 cauchy class test


BOOST_AUTO_TEST_CASE(pre_cauchy_loops)
{
	/*
		In the cauchy endgame we want to compute cauchy loops and make sure that we have that the max and min norms of the samples around the 
		origin are within some heuristic value. The function InitialCauchyLoops does this while holding onto the cycle number. 
		This test case is making sure we succeed and get the correct cycle number. 
	*/
	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> pseg_times;
	SampCont<mpfr> pseg_samples;

	mpfr time("0.1");
	Vec<mpfr> sample(1);

	pseg_times.push_back(time);
	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
	pseg_samples.push_back(sample);

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	my_endgame.SetPSEGSamples(pseg_samples);
	my_endgame.SetPSEGTimes(pseg_times);

	auto success_of_pre_cauchy_loops =  my_endgame.InitialCauchyLoops<mpfr>();

	BOOST_CHECK(success_of_pre_cauchy_loops == SuccessCode::Success);
	BOOST_CHECK(my_endgame.CycleNumber() == 1);

}//end pre_cauchy_loops

BOOST_AUTO_TEST_CASE(pre_cauchy_loops_cycle_num_greater_than_1)
{
	/*
		In the cauchy endgame we want to compute cauchy loops and make sure that we have that the max and min norms of the samples around the 
		origin are within some heuristic value. The function InitialCauchyLoops does this while holding onto the cycle number. 
		This test case is making sure we succeed and get the correct cycle number. 
	*/
	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> pseg_times;
	SampCont<mpfr> pseg_samples;



	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr("0.1");
	pseg_times.push_back(time);
	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
	pseg_samples.push_back(sample);


	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	my_endgame.SetPSEGSamples(pseg_samples);
	my_endgame.SetPSEGTimes(pseg_times);

	BOOST_CHECK(my_endgame.InitialCauchyLoops<mpfr>() == SuccessCode::Success);
	BOOST_CHECK(my_endgame.CycleNumber() == 2);
}// end pre_cauchy_loops_cycle_num_greater_than_1

BOOST_AUTO_TEST_CASE(first_approximation_using_pseg)
{
	/*
		The function tested is ComputeFirstApproximation. This function is the culmination of tracking samples until our estimation of 
		c over k is stabilized and we have our Cauchy loop ratios within a tolerance. When this is acheived we can use the power series approximation
		to compute our first extrapolant at the origin. After this we will use the cauchy integral formula to compute all further extrapolants.
	*/
	mpfr_float::default_precision(16);


	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3) + 1)*t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	mpfr origin = mpfr(0,0);
	Vec<mpfr> x_origin(1);
	x_origin << mpfr(1,0);

	Vec<mpfr> first_approx(1);
	mpfr time("0.1");
	Vec<mpfr> sample(1);
	sample << mpfr("0.5","0"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic


	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);



	auto first_approx_success = my_endgame.InitialPowerSeriesApproximation(time,sample,origin,first_approx);

	BOOST_CHECK((first_approx - x_origin).norm() < 1e-2);
	BOOST_CHECK(my_endgame.CycleNumber() == 3);

}// end first_approximation_using_pseg

BOOST_AUTO_TEST_CASE(compute_cauchy_approximation_cycle_num_1)
{
	/*
		This test case uses all the sample points collected by CircleTrack around a non-singular point to compute an extrapolant using the 
		trapezoidal rule for integration. 
	*/

	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
	SampCont<mpfr> cauchy_samples; //samples are space values that may be a vector of complex numbers.


	mpfr time("0.1");
	Vec<mpfr> sample(1);
	Vec<mpfr> x_origin(1);

	cauchy_times.push_back(time);
	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
	cauchy_samples.push_back(sample);
	x_origin << mpfr(1,0);

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	my_endgame.SetCauchyTimes(cauchy_times);
	my_endgame.SetCauchySamples(cauchy_samples);

	auto first_track_success =  my_endgame.CircleTrack(time,sample);

	my_endgame.CycleNumber(1); // manually set cycle number to 1 for this test

	Vec<mpfr> first_cauchy_approx;
	auto code = my_endgame.ComputeCauchyApproximationOfXAtT0<mpfr>(first_cauchy_approx);

	BOOST_CHECK((first_cauchy_approx - x_origin).norm() < my_endgame.Tolerances().newton_during_endgame);

}// end compute_cauchy_approximation_cycle_num_1

BOOST_AUTO_TEST_CASE(compute_cauchy_approximation_cycle_num_greater_than_1)
{

	/*
		This test case uses all the sample points collected by CircleTrack around a non-singular point to compute an extrapolant using the 
		trapezoidal rule for integration. CircleTrack is called twice to close up the loop. 
	*/

	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	TimeCont<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
	SampCont<mpfr> cauchy_samples; //samples are space values that may be a vector of complex numbers.



	mpfr time(1);
	Vec<mpfr> sample(1);
	Vec<mpfr> x_origin(1);

	time = mpfr("0.1");
	cauchy_times.push_back(time);
	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
	cauchy_samples.push_back(sample);
	x_origin << mpfr(1,0);


	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	my_endgame.SetCauchyTimes(cauchy_times);
	my_endgame.SetCauchySamples(cauchy_samples);

	auto first_track_success =  my_endgame.CircleTrack(time,sample);
	auto second_track_success = my_endgame.CircleTrack(my_endgame.GetCauchyTimes<mpfr>().back(),my_endgame.GetCauchySamples<mpfr>().back());

	my_endgame.CycleNumber(2);

	Vec<mpfr> first_cauchy_approx;
	auto code = my_endgame.ComputeCauchyApproximationOfXAtT0<mpfr>(first_cauchy_approx);

	BOOST_CHECK((first_cauchy_approx - x_origin).norm() < my_endgame.Tolerances().newton_during_endgame);
}// end compute_cauchy_approximation_cycle_num_greater_than_1

BOOST_AUTO_TEST_CASE(find_cauchy_samples_cycle_num_1)
{
	/*
		To actually use the Cauchy integral formula we must track around the origin and collect samples. This is done by ComputeCauchySamples. 
	*/

	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);



	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr(".1");

	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 


	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	auto finding_cauchy_samples_success = my_endgame.ComputeCauchySamples(time,sample);

	BOOST_CHECK((my_endgame.GetCauchySamples<mpfr>().back() - my_endgame.GetCauchySamples<mpfr>().front()).norm() < my_endgame.Tolerances().newton_during_endgame);
	BOOST_CHECK(my_endgame.GetCauchySamples<mpfr>().size() == 4);
	BOOST_CHECK(my_endgame.CycleNumber() == 1);

}// end find_cauchy_samples_cycle_num_1

BOOST_AUTO_TEST_CASE(find_cauchy_samples_cycle_num_greater_than_1)
{
	/*
		To actually use the Cauchy integral formula we must track around the origin and collect samples. This is done by ComputeCauchySamples. 
	*/
	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);


	auto time = mpfr("0.1");

	Vec<mpfr> sample(1);
	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 


	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	auto finding_cauchy_samples_success = my_endgame.ComputeCauchySamples(time,sample);

	BOOST_CHECK((my_endgame.GetCauchySamples<mpfr>().back() - my_endgame.GetCauchySamples<mpfr>().front()).norm() < my_endgame.Tolerances().newton_during_endgame);
	BOOST_CHECK_EQUAL(my_endgame.GetCauchySamples<mpfr>().size(), 7);
	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(), 2);

}// end find_cauchy_samples_cycle_num_greater_than_1

BOOST_AUTO_TEST_CASE(full_test_cycle_num_1)
{
	/*
		Full blown test to see if we can actually find the non singular point at the origin. 
	*/

	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);


	
	


	auto time = mpfr(".1");

	Vec<mpfr> sample(1);
	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 

	Vec<mpfr> solution(1);
	solution << mpfr(1,0);

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	auto cauchy_endgame_success = my_endgame.CauchyEG(time,sample);

	BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - solution).norm() < my_endgame.Tolerances().newton_during_endgame);
	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(), 1);

}// end full_test_cycle_num_1

BOOST_AUTO_TEST_CASE(full_test_cycle_num_greater_than_1)
{
	/*
		Full blown test to see if we can actually find the singular point at the origin. 
	*/
	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;
	newton_preferences.max_num_newton_iterations = 2;
	newton_preferences.min_num_newton_iterations = 1;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	auto time = mpfr("0.1");

	Vec<mpfr> sample(1);
	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01");

	Vec<mpfr> x_origin(1); 
	x_origin << mpfr(1,0);

	
	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	auto cauchy_endgame_success = my_endgame.CauchyEG(time,sample);


	BOOST_CHECK(cauchy_endgame_success==SuccessCode::Success);
	BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - x_origin).norm() < my_endgame.Tolerances().newton_during_endgame);
	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(), 2);
	BOOST_CHECK_EQUAL(mpfr_float::default_precision(),30);
}// end full_test_cycle_num_greater_than_1


BOOST_AUTO_TEST_CASE(cauchy_endgame_test_cycle_num_greater_than_1_base_precision_16)
{
	/*
		Full blown test to see if we can actually find the singular point at the origin. 
	*/
	mpfr_float::default_precision(16);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;
	newton_preferences.max_num_newton_iterations = 2;
	newton_preferences.min_num_newton_iterations = 1;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	tracker.ReinitializeInitialStepSize(false);


	mpfr time(1);
	Vec<mpfr> sample(1);
	Vec<mpfr> x_origin(1);


	time = mpfr("0.1");

	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
	x_origin << mpfr(1,0);

	
	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);
	std::cout.precision(8);
	auto cauchy_endgame_success = my_endgame.CauchyEG(time,sample);

	// std::cout << "first cauchy approx is " << my_endgame.FinalApproximation<mpfr>() << '\n';
	BOOST_CHECK(cauchy_endgame_success==SuccessCode::Success);
	BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - x_origin).norm() < my_endgame.Tolerances().newton_during_endgame);
	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(), 2);
	BOOST_CHECK_EQUAL(mpfr_float::default_precision(),bertini::DoublePrecision());
}// end cauchy_endgame_test_cycle_num_greater_than_1



BOOST_AUTO_TEST_CASE(cauchy_mp_for_cauchy_class_multiple_variables)
{
	/*
		Full blown test to see if we can actually find the non singular point at the origin. This example has multiple variables. 
	*/

	mpfr_float::default_precision(30);



	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t"), y = std::make_shared<Variable>("y");
	VariableGroup vars{x,y};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	sys.AddFunction((pow(x-1,3))*(1-t) + (pow(x,3) + 1)*t);
	sys.AddFunction((pow(y-1,2))*(1-t) + (pow(y,2) + 1)*t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);


	mpfr current_time(1);
	Vec<mpfr> current_space(2);
	current_time = mpfr(".1");
	current_space <<  mpfr("5.000000000000001e-01", "9.084258952712920e-17") ,mpfr("9.000000000000001e-01","4.358898943540673e-01");

	Vec<mpfr> correct(2);
	correct << mpfr(1,0),mpfr(1,0);

	config::Endgame<mpfr_float> endgame_settings;
	config::Cauchy<mpfr_float> cauchy_settings;
	config::Security<mpfr_float> security_settings;

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker,cauchy_settings,endgame_settings,security_settings);

	my_endgame.CauchyEG(current_time,current_space);

	BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - correct).norm() < my_endgame.Tolerances().newton_during_endgame);

}// end cauchy_mp_for_cauchy_class_multiple_variables


BOOST_AUTO_TEST_CASE(griewank_osborne)
{
	// Griewank Osborne is a very classic example. Here we allow x and y to mix. There are six paths to be tracked and we know there values 
	// at t = 0.1. 

	// Three of these paths will converge to origin and three will diverge to infinity triggering a SecurityMaxNorm issue. 

	// This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
	// the endgame. 
	

	mpfr_float::default_precision(30);

	bertini::System griewank_osborn_sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t"), y = std::make_shared<Variable>("y");
	VariableGroup vars{x,y};
	griewank_osborn_sys.AddVariableGroup(vars); 

	griewank_osborn_sys.AddFunction((mpfr("29")/mpfr("16"))*pow(x,3)-2*x*y);
	griewank_osborn_sys.AddFunction((y - pow(x,2)));
	griewank_osborn_sys.Homogenize();
	griewank_osborn_sys.AutoPatch();

	BOOST_CHECK(griewank_osborn_sys.IsHomogeneous());
	BOOST_CHECK(griewank_osborn_sys.IsPatched());	

	auto griewank_TD = bertini::start_system::TotalDegree(griewank_osborn_sys);
	griewank_TD.Homogenize();
	BOOST_CHECK(griewank_TD.IsHomogeneous());
	BOOST_CHECK(griewank_TD.IsPatched());

	auto final_griewank_osborn_system = (1-t)*griewank_osborn_sys + t*griewank_TD;
	final_griewank_osborn_system.AddPathVariable(t);


	auto griewank_AMP = config::AMPConfigFrom(final_griewank_osborn_system);

	AMPTracker tracker(final_griewank_osborn_system);
	
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(griewank_AMP);


	mpfr t_start(1), t_endgame_boundary(0.1);
	std::vector<Vec<mpfr> > griewank_solutions;
	std::vector<Vec<mpfr> > griewank_homogenized_solutions;
	for (unsigned ii = 0; ii < griewank_TD.NumStartPoints(); ++ii)
	{
		mpfr_float::default_precision(30);
		final_griewank_osborn_system.precision(30);
		auto start_point = griewank_TD.StartPoint<mpfr>(ii);

		Vec<mpfr> result;
		SuccessCode tracking_success;

		tracking_success = tracker.TrackPath(result,t_start,t_endgame_boundary,start_point);
		BOOST_CHECK(tracking_success==SuccessCode::Success);

		griewank_homogenized_solutions.push_back(result);
		griewank_solutions.push_back(final_griewank_osborn_system.DehomogenizePoint(result));
	}

	Vec<mpfr> correct(2);
	correct << mpfr(0,0),mpfr(0,0);

	config::Endgame<mpfr_float> endgame_settings;
	config::Cauchy<mpfr_float> cauchy_settings;
	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker,cauchy_settings,endgame_settings,security_settings,tolerances);


	unsigned num_paths_diverging = 0;
	unsigned num_paths_converging = 0;
	for(auto s : griewank_homogenized_solutions) //current_space_values)
	{
		std::cout << ".1 solution is " << final_griewank_osborn_system.DehomogenizePoint(s) << '\n';
		mpfr_float::default_precision(30);
		final_griewank_osborn_system.precision(Precision(s(0)));
		SuccessCode endgame_success = my_endgame.CauchyEG(t_endgame_boundary,s);

		std::cout << "endgame solution is " << final_griewank_osborn_system.DehomogenizePoint(my_endgame.FinalApproximation<mpfr>()) << '\n';

		if(endgame_success == SuccessCode::Success)
		{
			// std::cout << "CONVERGED" <<'\n';
			num_paths_converging++;
		}
		if(endgame_success == SuccessCode::SecurityMaxNormReached || endgame_success == SuccessCode::GoingToInfinity)
		{
			// std::cout << "DIVERGED" << '\n';
			num_paths_diverging++;
		}
	}
	BOOST_CHECK_EQUAL(num_paths_converging,3);
	BOOST_CHECK_EQUAL(num_paths_diverging,3);

}//end compute griewank osborne

// has six solutions at t = .1:
// 0
// 0.687592791426802019127961784761 0.0567041721787413764699348206477
// 1.11076734170908975052327605226  0.207914822575706710605647487

// 1
// 1.02264960658155701356264444257  0.520917033127216794197167359926
// 1.11076734170909913190783413484  0.207914822575691493611316218448

// 2
// 0.989757601991555768794484038153 -0.57762120530600610801563732366
// 1.11076734170909918741898536609  0.20791482257569138952790765984

// 3
// 0.687592791426887395278555459299 0.0567041721787893780032385748768 
// 0.689232658290901023523389312686 -0.207914822575691576878043065335

// 4
// 1.02264960658081959931948637747  0.520917033127118696605766956908
// 0.689232658292013194430512839668 -0.207914822576518617066069857076

// 5
// 0.989757601991599268196007422024 -0.577621205306094375358982164819
// 0.689232658290901310861758919599 -0.207914822575714712814276165999
BOOST_AUTO_TEST_CASE(total_degree_start_system_cauchy_class_used_with_AMP)
{
	/*
	In this example we take a decoupled system, homogenize and patch it. Track to endgame boundary and then run our endgame on the space
	values we have. 

	*/
	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);
	using namespace bertini::tracking;
	mpfr_float::default_precision(30);

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddVariableGroup(v);

	sys.AddFunction(pow(x-1,3));
	sys.AddFunction(pow(y-1,2));
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

	auto AMP = config::AMPConfigFrom(final_system);

	auto tracker = AMPTracker(final_system);
	config::Stepping<mpfr_float> stepping_preferences;
	config::Newton newton_preferences;
	tracker.Setup(config::Predictor::HeunEuler,
	              	mpfr_float("1e-5"), mpfr_float("1e5"),
					stepping_preferences, newton_preferences);

	tracker.AMPSetup(AMP);
	

	mpfr t_start(1), t_endgame_boundary(0.1);
	std::vector<Vec<mpfr> > solutions;
	std::vector<Vec<mpfr> > homogenized_solutions;
	for (unsigned ii = 0; ii < 1; ++ii)
	{
		mpfr_float::default_precision(30);
		final_system.precision(30);
		auto start_point = TD.StartPoint<mpfr>(ii);

		Vec<mpfr> result;
		SuccessCode tracking_success;

		tracking_success = tracker.TrackPath(result,t_start,t_endgame_boundary,start_point);
		BOOST_CHECK(tracking_success==SuccessCode::Success);

		homogenized_solutions.push_back(result);
		solutions.push_back(final_system.DehomogenizePoint(result));
	}

	Vec<mpfr> correct(2);
	correct << mpfr(1,0),mpfr(1,0);

	endgame::CauchyEndgame<AMPTracker> my_endgame(tracker);

	tracker.Setup(config::Predictor::HeunEuler,
	              	mpfr_float("1e-6"), mpfr_float("1e5"),
					stepping_preferences, newton_preferences);


	std::vector<Vec<mpfr> > endgame_solutions;

	unsigned num_successful_occurences = 0;
	unsigned num_min_track_time_reached = 0;
	for (auto s : homogenized_solutions)
	{
		mpfr_float::default_precision(30);
		final_system.precision(Precision(s(0)));
		SuccessCode endgame_success = my_endgame.CauchyEG(t_endgame_boundary,s);
		if(endgame_success == SuccessCode::Success)
		{
			if((tracker.GetSystem().DehomogenizePoint(my_endgame.FinalApproximation<mpfr>())-correct).norm() < my_endgame.Tolerances().newton_during_endgame)
			{
				num_successful_occurences++;
			}
		}
	}
 	BOOST_CHECK_EQUAL(num_successful_occurences,6);

}


BOOST_AUTO_TEST_SUITE_END()
