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

// BOOST_AUTO_TEST_CASE(circle_track_d_for_cauchy_class_test)
// {
// 	mpfr_float::default_precision(16);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); //f(x) = (x-1)^3 + t*0, need t*0 for derivative calculation. 

// 	sys.AddFunction((x - mpfr(1))*(1-t) + (x + mpfr(1))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > cauchy_samples; //samples are space values that may be a vector of complex numbers.


// 	mpfr time(1);
// 	Vec<mpfr> sample(1);


// 	time = mpfr(".1");
// 	cauchy_times.push_back(time);
// 	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
// 	cauchy_samples.push_back(sample);

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 	auto first_track_success =  My_Endgame.CircleTrack(time,sample);

// 	BOOST_CHECK((My_Endgame.cauchy_samples_.back() - sample).norm() < My_Endgame.GetTrackToleranceDuringEndgame());

// } // end circle_track_d_for_cauchy_class_test

// BOOST_AUTO_TEST_CASE(circle_track_mp_cycle_num_greater_than_1_for_cauchy_class_test)
// {
// 	mpfr_float::default_precision(30);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); 

// 	sys.AddFunction( pow(x - mpfr("1"),2)*(1-t) + (pow(x,2) + mpfr("1"))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > cauchy_samples; //samples are space values that may be a vector of complex numbers.



// 	mpfr time(1);
// 	Vec<mpfr> sample(1);


// 	time = mpfr("0.1");
// 	cauchy_times.push_back(time);
// 	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
// 	cauchy_samples.push_back(sample);

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);


// 	auto tracking_success =  My_Endgame.CircleTrack(time,sample);

// 	auto first_track_sample = My_Endgame.cauchy_samples_.back();
// 	// std::cout << "first track sample is " << first_track_sample << '\n';

// 	BOOST_CHECK((first_track_sample - sample).norm() > My_Endgame.GetTrackToleranceDuringEndgame());

// 	tracking_success =  My_Endgame.CircleTrack(time,first_track_sample);

// 	auto second_track_sample = My_Endgame.cauchy_samples_.back();

// 	// std::cout << "second track sample is " << second_track_sample << '\n';

// 	BOOST_CHECK((second_track_sample - sample).norm() < My_Endgame.GetTrackToleranceDuringEndgame());
	
// } // end circle_track_mp_cycle_num_greater_than_1_for_cauchy_class_test


// BOOST_AUTO_TEST_CASE(compute_c_over_k_dbl_for_cauchy_class)
// {
// 	/* In the power series endgame there is a bound on that is calculated that will be used as a higher
// 	bound for the exhaustive search of the best cycle number. 

// 	This is calcuated using a modified version of the cycle test in the bertini book, page 53. 

// 	The bound calculated is then compared against the user defined setting MaxCycleNum which is default at 6. 

// 	We take a function (x-1)^3 because we know it should have a cycle number of 3. (This is true but will not be true near 0 or near 0.1).

// 	This is because the path will not look cubic globally. 
// 	*/
	
// 	mpfr_float::default_precision(16);

// 	bertini::System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	sys.AddFunction(pow(x - 1,3));  //f(x) = (x-1)^3

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 

// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> pseg_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > pseg_samples; //samples are space values that may be a vector of complex numbers.
	

// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	// Vec<mpfr> derivative(1);

// 	time = mpfr(.1); // x = .1
// 	pseg_times.push_back(time);
// 	sample << mpfr("-0.729"); // f(.1) = -0.729
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".05"); // x = .1/2 = .05
// 	pseg_times.push_back(time);
// 	sample << mpfr("-0.857375"); //f(.05) = -0.857375
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".025"); // x = .05/2 = .025
// 	pseg_times.push_back(time);
// 	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
// 	pseg_samples.push_back(sample);


// 	bertini::tracking::config::Security endgame_security_struct;
// 	bertini::tracking::config::Tolerances endgame_tolerances_struct;

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_security_struct,endgame_tolerances_struct);
// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto first_c_over_k = My_Endgame.ComputeCOverK();

// 	// std::cout << "first c over k is " << first_c_over_k << '\n';
// 	// std::cout << "first diff " << abs(first_c_over_k - mpfr_float("1.12917")) << '\n';
// 	// std::cout << "Tol is " << My_Endgame.GetTrackToleranceDuringEndgame() << '\n';

// 	 BOOST_CHECK( abs(first_c_over_k - mpfr_float("1.12917")) < 1e-5 ); 

// 	//Setting up a new sample for approximation.
// 	time = mpfr(".0125"); //.025/2 = .0125
// 	pseg_times.push_back(time);
// 	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
// 	pseg_samples.push_back(sample);

// 	//Get rid of earliest sample. 
// 	pseg_times.pop_front();
// 	pseg_samples.pop_front();

// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto second_c_over_k = My_Endgame.ComputeCOverK();

// 	// std::cout << "second_c_over_k is " << second_c_over_k << '\n';
// 	// std::cout << "second diff is " << abs(second_c_over_k - mpfr_float("1.05888")) << '\n';
// 	BOOST_CHECK( abs(second_c_over_k - mpfr_float("1.05888")) <  1e-5); 

// } // end compute c over k dbl for cauchy class 

// BOOST_AUTO_TEST_CASE(compute_c_over_k_mp_for_cauchy_class)
// {
// 	/* In the power series endgame there is a bound on that is calculated that will be used as a higher
// 	bound for the exhaustive search of the best cycle number. 

// 	This is calcuated using a modified version of the cycle test in the bertini book, page 53. 

// 	The bound calculated is then compared against the user defined setting MaxCycleNum which is default at 6. 

// 	We take a function (x-1)^3 because we know it should have a cycle number of 3. (This is true but will not be true near 0 or near 0.1).

// 	This is because the path will not look cubic globally. 
// 	*/
	
// 	mpfr_float::default_precision(30);

// 	bertini::System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	sys.AddFunction(pow(x - 1,3));  //f(x) = (x-1)^3

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 

// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> pseg_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > pseg_samples; //samples are space values that may be a vector of complex numbers.
	

// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	// Vec<mpfr> derivative(1);

// 	time = mpfr(.1); // x = .1
// 	pseg_times.push_back(time);
// 	sample << mpfr("-0.729"); // f(.1) = -0.729
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".05"); // x = .1/2 = .05
// 	pseg_times.push_back(time);
// 	sample << mpfr("-0.857375"); //f(.05) = -0.857375
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".025"); // x = .05/2 = .025
// 	pseg_times.push_back(time);
// 	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
// 	pseg_samples.push_back(sample);


// 	bertini::tracking::config::Security endgame_security_struct;
// 	bertini::tracking::config::Tolerances endgame_tolerances_struct;

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_security_struct,endgame_tolerances_struct);
// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto first_c_over_k = My_Endgame.ComputeCOverK();

// 	// std::cout << "first c over k is " << first_c_over_k << '\n';
// 	// std::cout << "first diff " << abs(first_c_over_k - mpfr_float("1.12917")) << '\n';
// 	// std::cout << "Tol is " << My_Endgame.GetTrackToleranceDuringEndgame() << '\n';

// 	 BOOST_CHECK( abs(first_c_over_k - mpfr_float("1.12917")) < 1e-5 ); 

// 	//Setting up a new sample for approximation.
// 	time = mpfr(".0125"); //.025/2 = .0125
// 	pseg_times.push_back(time);
// 	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
// 	pseg_samples.push_back(sample);

// 	//Get rid of earliest sample. 
// 	pseg_times.pop_front();
// 	pseg_samples.pop_front();

// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto second_c_over_k = My_Endgame.ComputeCOverK();

// 	// std::cout << "second_c_over_k is " << second_c_over_k << '\n';
// 	// std::cout << "second diff is " << abs(second_c_over_k - mpfr_float("1.05888")) << '\n';
// 	BOOST_CHECK( abs(second_c_over_k - mpfr_float("1.05888")) <  1e-5); 

// } // end compute c over k mp for cauchy class 

// BOOST_AUTO_TEST_CASE(checking_agreement_function_dbl_for_cauchy_class)
// {
// 	// In the power series endgame there is a bound on that is calculated that will be used as a higher
// 	// bound for the exhaustive search of the best cycle number. 

// 	// This is calcuated using a modified version of the cycle test in the bertini book, page 53. 

// 	// The bound calculated is then compared against the user defined setting MaxCycleNum which is default at 6. 

// 	// We take a function (x-1)^3 because we know it should have a cycle number of 3. (This is true but will not be true near 0 or near 0.1).

// 	// This is because the path will not look cubic globally. 
	
	
// 	mpfr_float::default_precision(16);

// 	bertini::System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t");
// 	sys.AddFunction( pow(x - mpfr("1"),3)*(1-t) + (pow(x,3) + mpfr("1"))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 

// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> pseg_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > pseg_samples; //samples are space values that may be a vector of complex numbers.
// 	std::deque<mpfr_float> c_over_k_array;

// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	// Vec<mpfr> derivative(1);

// 	time = mpfr(".1"); // x = .1
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.100000000000000e+00", "6.244997998398396e-01"); // f(.1) = 1.100000000000000e+00 6.244997998398396e-01
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".05"); // x = .1/2 = .05
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.125000000000000e+00", "4.841229182759271e-01"); //f(.05) = 1.125000000000000e+00 4.841229182759271e-01
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".025"); // x = .05/2 = .025
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.123854702920065e+00", "3.736283818773165e-01"); // f(.025) = 1.123854702920065e+00 3.736283818773165e-01
// 	pseg_samples.push_back(sample);


// 	bertini::tracking::config::Security endgame_security_struct;
// 	bertini::tracking::config::Tolerances endgame_tolerances_struct;

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_security_struct,endgame_tolerances_struct);
// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto first_c_over_k = My_Endgame.ComputeCOverK();

// 	c_over_k_array.push_back(first_c_over_k);

// 	//Setting up a new sample for approximation.
// 	time = mpfr(".0125"); //.025/2 = .0125
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.111721780135302e+00", "2.886596646224579e-01"); // f(.0125) = 1.111721780135302e+00 2.886596646224579e-01
// 	pseg_samples.push_back(sample);

// 	//Get rid of earliest sample. 
// 	pseg_times.pop_front();
// 	pseg_samples.pop_front();

// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto second_c_over_k = My_Endgame.ComputeCOverK();

// 	c_over_k_array.push_back(second_c_over_k);



// 	//Setting up a new sample for approximation.
// 	time = mpfr(".00625"); //.025/2 = .0125
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.096071609421043e+00", "2.237005761359081e-01"); // f(.00625) = 1.096071609421043e+00 2.237005761359081e-01
// 	pseg_samples.push_back(sample);

// 	//Get rid of earliest sample. 
// 	pseg_times.pop_front();
// 	pseg_samples.pop_front();

// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto third_c_over_k = My_Endgame.ComputeCOverK();

// 	c_over_k_array.push_back(third_c_over_k);

// 	auto stabilized = My_Endgame.CheckForCOverKStabilization(c_over_k_array);

// 	BOOST_CHECK(stabilized == true);

// } // end check agreement function dbl for cauchy 

// BOOST_AUTO_TEST_CASE(checking_agreement_function_mp_for_cauchy_class)
// {
// 	/* In the power series endgame there is a bound on that is calculated that will be used as a higher
// 	bound for the exhaustive search of the best cycle number. 

// 	This is calcuated using a modified version of the cycle test in the bertini book, page 53. 

// 	The bound calculated is then compared against the user defined setting MaxCycleNum which is default at 6. 

// 	We take a function (x-1)^3 because we know it should have a cycle number of 3. (This is true but will not be true near 0 or near 0.1).

// 	This is because the path will not look cubic globally. 
// 	*/
	
// 	mpfr_float::default_precision(30);

// 	bertini::System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t");
// 	sys.AddFunction( pow(x - mpfr("1"),3)*(1-t) + (pow(x,3) + mpfr("1"))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 

// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> pseg_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > pseg_samples; //samples are space values that may be a vector of complex numbers.
// 	std::deque<mpfr_float> c_over_k_array;

// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	// Vec<mpfr> derivative(1);

// 	time = mpfr(".1"); // x = .1
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.100000000000000e+00", "6.244997998398396e-01"); // f(.1) = 1.100000000000000e+00 6.244997998398396e-01
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".05"); // x = .1/2 = .05
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.125000000000000e+00", "4.841229182759271e-01"); //f(.05) = 1.125000000000000e+00 4.841229182759271e-01
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".025"); // x = .05/2 = .025
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.123854702920065e+00", "3.736283818773165e-01"); // f(.025) = 1.123854702920065e+00 3.736283818773165e-01
// 	pseg_samples.push_back(sample);


// 	bertini::tracking::config::Security endgame_security_struct;
// 	bertini::tracking::config::Tolerances endgame_tolerances_struct;

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_security_struct,endgame_tolerances_struct);
// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto first_c_over_k = My_Endgame.ComputeCOverK();

// 	c_over_k_array.push_back(first_c_over_k);

// 	//Setting up a new sample for approximation.
// 	time = mpfr(".0125"); //.025/2 = .0125
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.111721780135302e+00", "2.886596646224579e-01"); // f(.0125) = 1.111721780135302e+00 2.886596646224579e-01
// 	pseg_samples.push_back(sample);

// 	//Get rid of earliest sample. 
// 	pseg_times.pop_front();
// 	pseg_samples.pop_front();

// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto second_c_over_k = My_Endgame.ComputeCOverK();

// 	c_over_k_array.push_back(second_c_over_k);



// 	//Setting up a new sample for approximation.
// 	time = mpfr(".00625"); //.025/2 = .0125
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.096071609421043e+00", "2.237005761359081e-01"); // f(.00625) = 1.096071609421043e+00 2.237005761359081e-01
// 	pseg_samples.push_back(sample);

// 	//Get rid of earliest sample. 
// 	pseg_times.pop_front();
// 	pseg_samples.pop_front();

// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto third_c_over_k = My_Endgame.ComputeCOverK();

// 	c_over_k_array.push_back(third_c_over_k);

// 	auto stabilized = My_Endgame.CheckForCOverKStabilization(c_over_k_array);

// 	BOOST_CHECK(stabilized == true);

// } // end check agreement function mp for cauchy 

// BOOST_AUTO_TEST_CASE(find_tolerance_for_closed_loop_dbl)
// {

// 	/* In the power series endgame there is a bound on that is calculated that will be used as a higher
// 	bound for the exhaustive search of the best cycle number. 

// 	This is calcuated using a modified version of the cycle test in the bertini book, page 53. 

// 	The bound calculated is then compared against the user defined setting MaxCycleNum which is default at 6. 

// 	We take a function (x-1)^3 because we know it should have a cycle number of 3. (This is true but will not be true near 0 or near 0.1).

// 	This is because the path will not look cubic globally. 
// 	*/
	
// 	mpfr_float::default_precision(30);

// 	bertini::System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t");
// 	sys.AddFunction( pow(x - mpfr("1"),3)*(1-t) + (pow(x,3) + mpfr("1"))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);

// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> pseg_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > pseg_samples; //samples are space values that may be a vector of complex numbers.

// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	// Vec<mpfr> derivative(1);

// 	time = mpfr(".1"); // x = .1
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.100000000000000e+00", "6.244997998398396e-01"); // f(.1) = 1.100000000000000e+00 6.244997998398396e-01
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".05"); // x = .1/2 = .05
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.125000000000000e+00", "4.841229182759271e-01"); //f(.05) = 1.125000000000000e+00 4.841229182759271e-01
// 	pseg_samples.push_back(sample);


// 	time = mpfr(".025"); // x = .05/2 = .025
// 	pseg_times.push_back(time);
// 	sample << mpfr("1.123854702920065e+00", "3.736283818773165e-01"); // f(.025) = 1.123854702920065e+00 3.736283818773165e-01
// 	pseg_samples.push_back(sample);


// 	bertini::tracking::config::Security endgame_security_struct;
// 	bertini::tracking::config::Tolerances endgame_tolerances_struct;

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_security_struct,endgame_tolerances_struct);
// 	My_Endgame.SetPSEGTimes(pseg_times);
// 	My_Endgame.SetPSEGSamples(pseg_samples);

// 	auto tol = My_Endgame.FindToleranceForClosedLoop(pseg_times.back(),pseg_samples.back());

// 	 //std::cout << "tol is " << tol << '\n';

// 	BOOST_CHECK(tol < 1e-5); // tol is 1e-06

// } // end find closed loop tol dbl

// BOOST_AUTO_TEST_CASE(check_closed_loop_for_cycle_num_1_cauchy_class_test)
// {
// 	mpfr_float::default_precision(16);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); //f(x) = (x-1)^3 + t*0, need t*0 for derivative calculation. 

// 	sys.AddFunction((x - mpfr(1))*(1-t) + (x + mpfr(1))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > cauchy_samples; //samples are space values that may be a vector of complex numbers.


// 	mpfr time(1);
// 	Vec<mpfr> sample(1);


// 	time = mpfr(".1");
// 	cauchy_times.push_back(time);
// 	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
// 	cauchy_samples.push_back(sample);


// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 	My_Endgame.SetCauchySamples(cauchy_samples);
// 	My_Endgame.SetCauchyTimes(cauchy_times);


// 	auto tracking_success =  My_Endgame.CircleTrack(time,sample);

// 	// std::cout << "first sample is " << cauchy_samples.front() << '\n';
// 	// std::cout << "last sample is " << cauchy_samples.back() << '\n';

// 	// for(int ii = 0; ii < My_Endgame.cauchy_samples_.size(); ii++)
// 	// {
// 	// 	std::cout << "samples at i is " << My_Endgame.cauchy_samples_[ii] << '\n';
// 	// }

// 	auto closed_loop_tolerance = My_Endgame.FindToleranceForClosedLoop(time,sample);

// 	// std::cout << "closed_loop_tolerance is " << closed_loop_tolerance << '\n';

// 	auto check_closed_loop = My_Endgame.CheckClosedLoop(closed_loop_tolerance);

// 	BOOST_CHECK(check_closed_loop == true);

// } // end check closed loop if cycle num is 1 for cauchy class test

// BOOST_AUTO_TEST_CASE(check_closed_loop_for_cycle_num_greater_than_1_cauchy_class_test)
// {

// 	mpfr_float::default_precision(30);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); 

// 	sys.AddFunction( pow(x - mpfr("1"),2)*(1-t) + (pow(x,2) + mpfr("1"))*t);


// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > cauchy_samples; //samples are space values that may be a vector of complex numbers.



// 	mpfr time(1);
// 	Vec<mpfr> sample(1);


// 	time = mpfr("0.1");
// 	cauchy_times.push_back(time);
// 	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
// 	cauchy_samples.push_back(sample);


// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 	My_Endgame.SetCauchySamples(cauchy_samples);
// 	My_Endgame.SetCauchyTimes(cauchy_times);

// 	auto tracking_success =  My_Endgame.CircleTrack(time,sample);

// 	// for(int ii = 0; ii < My_Endgame.cauchy_samples_.size(); ii++)
// 	// {
// 	// 	std::cout << "samples at i is " << My_Endgame.cauchy_samples_[ii] << '\n';
// 	// }


// 	auto closed_loop_tolerance = My_Endgame.FindToleranceForClosedLoop(time,sample);

// 	// std::cout << "closed_loop_tolerance is " << closed_loop_tolerance << '\n';

// 	auto check_closed_loop = My_Endgame.CheckClosedLoop(closed_loop_tolerance);

// 	BOOST_CHECK(check_closed_loop == false);
	
// } // end check closed loop if cycle num is 1 for cauchy class test


// BOOST_AUTO_TEST_CASE(compare_cauchy_ratios_for_cauchy_class_test)
// {
// 	mpfr_float::default_precision(16);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); //f(x) = (x-1)^3 + t*0, need t*0 for derivative calculation. 

// 	sys.AddFunction((x - mpfr(1))*(1-t) + (x + mpfr(1))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);


	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	// std::cout << "tracker made " << '\n';

// 	std::deque<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > cauchy_samples; //samples are space values that may be a vector of complex numbers.


// 	mpfr time(1);
// 	Vec<mpfr> sample(1);


// 	time = mpfr(".1");
// 	cauchy_times.push_back(time);
// 	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
// 	cauchy_samples.push_back(sample);

// 		// std::cout << "sample made " << '\n';

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 		// std::cout << "endgame made " << '\n';

// 	auto tracking_success =  My_Endgame.CircleTrack(time,sample);

// 		// std::cout << "circle track done " << tracking_success << '\n';

// 	auto comparison_of_cauchy_ratios = My_Endgame.CompareCauchyRatios();

// 		// std::cout << "comparison made  made " << comparison_of_cauchy_ratios << '\n';

// 	//std::cout << "comparison is " << comparison_of_cauchy_ratios << '\n';


// 	BOOST_CHECK(comparison_of_cauchy_ratios == true);

// } // end compare cauchy ratios for cauchy class test

// BOOST_AUTO_TEST_CASE(compare_cauchy_ratios_cycle_num_greater_than_1_for_cauchy_class_test)
// {
// 	mpfr_float::default_precision(30);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); 

// 	sys.AddFunction( pow(x - mpfr("1"),2)*(1-t) + (pow(x,2) + mpfr("1"))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > cauchy_samples; //samples are space values that may be a vector of complex numbers.



// 	mpfr time(1);
// 	Vec<mpfr> sample(1);


// 	time = mpfr("0.1");
// 	cauchy_times.push_back(time);
// 	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
// 	cauchy_samples.push_back(sample);


// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 	auto tracking_success =  My_Endgame.CircleTrack(time,sample);

// 	auto comparison_of_cauchy_ratios = My_Endgame.CompareCauchyRatios();

// 	//std::cout << "comparison is " << comparison_of_cauchy_ratios << '\n';


// 	BOOST_CHECK(comparison_of_cauchy_ratios == true);

// } // end compare cauchy ratios for cycle num greater than 1 cauchy class test


// BOOST_AUTO_TEST_CASE(pre_cauchy_loops_for_cauchy_class_test)
// {
// 	mpfr_float::default_precision(16);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); //f(x) = (x-1)^3 + t*0, need t*0 for derivative calculation. 

// 	sys.AddFunction((x - mpfr(1))*(1-t) + (x + mpfr(1))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> pseg_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > pseg_samples; //samples are space values that may be a vector of complex numbers.

// 	mpfr time(1);
// 	Vec<mpfr> sample(1);


// 	time = mpfr(".1");
// 	pseg_times.push_back(time);
// 	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
// 	pseg_samples.push_back(sample);

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 	My_Endgame.SetPSEGSamples(pseg_samples);
// 	My_Endgame.SetPSEGTimes(pseg_times);

// 	auto success_of_pre_cauchy_loops =  My_Endgame.PreCauchyLoops();

// 	// for(int ii = 0; ii < My_Endgame.cauchy_samples_.size(); ii++)
// 	// {
// 	// std::cout << "samples at i is " << My_Endgame.cauchy_samples_[ii] << '\n';
//  // 	}

// 	// std::cout << "success code is " << success_of_pre_cauchy_loops << '\n';

// 	BOOST_CHECK(success_of_pre_cauchy_loops == bertini::tracking::SuccessCode::Success);
// 	BOOST_CHECK(My_Endgame.endgame_settings_.cycle_number == 1);

// }

// BOOST_AUTO_TEST_CASE(pre_cauchy_loops_cycle_num_greater_than_1_for_cauchy_class_test)
// {
// 	mpfr_float::default_precision(30);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); 

// 	sys.AddFunction( pow(x - mpfr("1"),2)*(1-t) + (pow(x,2) + mpfr("1"))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> pseg_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > pseg_samples; //samples are space values that may be a vector of complex numbers.



// 	mpfr time(1);
// 	Vec<mpfr> sample(1);


// 	time = mpfr("0.1");
// 	pseg_times.push_back(time);
// 	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
// 	pseg_samples.push_back(sample);


// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 	My_Endgame.SetPSEGSamples(pseg_samples);
// 	My_Endgame.SetPSEGTimes(pseg_times);

// 	auto success_of_pre_cauchy_loops =  My_Endgame.PreCauchyLoops();

// 	// for(int ii = 0; ii < My_Endgame.cauchy_samples_.size(); ii++)
// 	// {
// 	// std::cout << "samples at i is " << My_Endgame.cauchy_samples_[ii] << '\n';
//  // 	}

// 	// std::cout << "success code is " << success_of_pre_cauchy_loops << '\n';

// 	BOOST_CHECK(success_of_pre_cauchy_loops == bertini::tracking::SuccessCode::Success);
// 	BOOST_CHECK(My_Endgame.endgame_settings_.cycle_number == 2);

// }

// BOOST_AUTO_TEST_CASE(first_approximation_using_pseg_for_cauchy_class_test)
// {

// 	/* Compute approximation at origin using three sample points. 

// 	// Then, find a new sample and remove earliest known sample. 

// 	//Compute a new approximation with the current three samples. 

// 	Do this till we have three approximations. Check to see if the third approximation is better than second. 
// 	*/
// 	mpfr_float::default_precision(16);


// 	bertini::System sys;
// 	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars);
// 	sys.AddPathVariable(t);
// 	// Define homotopy system
// 	sys.AddFunction( pow(x - mpfr(1),3)*(1-t) + (pow(x,3) + mpfr(1))*t);

// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	mpfr origin = mpfr("0","0");
// 	Vec<mpfr> x_origin(1);
// 	x_origin << mpfr("1","0");
// 	//std::deque<mpfr> times; //times are not vectors they are just complex numbers.
// 	//std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.


// 	Vec<mpfr> first_approx(1);
// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	// Vec<mpfr> derivative(1);

// 	time = mpfr(".1"); // x = .1
// 	// times.push_back(time);
// 	sample << mpfr(".500000000000000","0"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic
// 	// samples.push_back(sample);


// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);



// 	auto first_approx_success = My_Endgame.ComputeFirstApproximation(time,sample,origin,first_approx);

// 	BOOST_CHECK((first_approx - x_origin).norm() < 1e-2);
// 	BOOST_CHECK(My_Endgame.endgame_settings_.cycle_number == 3);

// }

// BOOST_AUTO_TEST_CASE(compute_cauchy_approximation_cycle_num_1_for_cauchy_class_test)
// {

// 	mpfr_float::default_precision(16);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); //f(x) = (x-1)^3 + t*0, need t*0 for derivative calculation. 

// 	sys.AddFunction((x - mpfr(1))*(1-t) + (x + mpfr(1))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > cauchy_samples; //samples are space values that may be a vector of complex numbers.


// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	Vec<mpfr> x_origin(1);

// 	time = mpfr(".1");
// 	cauchy_times.push_back(time);
// 	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 
// 	cauchy_samples.push_back(sample);
// 	x_origin << mpfr("1.0","0.0");

// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 	My_Endgame.SetCauchyTimes(cauchy_times);
// 	My_Endgame.SetCauchySamples(cauchy_samples);

// 	auto first_track_success =  My_Endgame.CircleTrack(time,sample);

// 	My_Endgame.endgame_settings_.cycle_number = 1;

// 	auto first_cauchy_approx = My_Endgame.ComputeCauchyApproximationOfXAtT0(My_Endgame.cauchy_samples_);

// 	// std::cout << "first cauchy approx is " << first_cauchy_approx << '\n';

// 	BOOST_CHECK((first_cauchy_approx - x_origin).norm() < My_Endgame.GetTrackToleranceDuringEndgame());

// }

// BOOST_AUTO_TEST_CASE(compute_cauchy_approximation_cycle_num_greater_than_1_for_cauchy_class_test)
// {
// 	mpfr_float::default_precision(30);

// 	System sys;
// 	Var x = std::make_shared<Variable>("x");
// 	Var t = std::make_shared<Variable>("t"); 

// 	sys.AddFunction( pow(x - mpfr("1"),2)*(1-t) + (pow(x,2) + mpfr("1"))*t);

// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars); 
// 	sys.AddPathVariable(t);


// 	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

// 	bertini::tracking::AMPTracker tracker(sys);
	
// 	bertini::tracking::config::Stepping stepping_preferences;
// 	bertini::tracking::config::Newton newton_preferences;

// 	tracker.Setup(bertini::tracking::config::Predictor::Euler,
//                 mpfr_float("1e-5"),
//                 mpfr_float("1e5"),
//                 stepping_preferences,
//                 newton_preferences);
	
// 	tracker.AMPSetup(AMP);

// 	std::deque<mpfr> cauchy_times; //times are not vectors they are just complex numbers.
// 	std::deque< Vec<mpfr> > cauchy_samples; //samples are space values that may be a vector of complex numbers.



// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	Vec<mpfr> x_origin(1);

// 	time = mpfr("0.1");
// 	cauchy_times.push_back(time);
// 	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 
// 	cauchy_samples.push_back(sample);
// 	x_origin << mpfr("1.0","0.0");


// 	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

// 	My_Endgame.SetCauchyTimes(cauchy_times);
// 	My_Endgame.SetCauchySamples(cauchy_samples);

// 	auto first_track_success =  My_Endgame.CircleTrack(time,sample);
// 	auto second_track_success = My_Endgame.CircleTrack(My_Endgame.cauchy_times_.back(),My_Endgame.cauchy_samples_.back());

// 	My_Endgame.endgame_settings_.cycle_number = 2;

// 	for(int i = 0; i < My_Endgame.cauchy_samples_.size(); ++i){
// 		std::cout << My_Endgame.cauchy_samples_[i] << '\n';
// 	}

// 	auto first_cauchy_approx = My_Endgame.ComputeCauchyApproximationOfXAtT0(My_Endgame.cauchy_samples_);

// 	// std::cout << "first cauchy approx is " << first_cauchy_approx << '\n';

// 	BOOST_CHECK((first_cauchy_approx - x_origin).norm() < My_Endgame.GetTrackToleranceDuringEndgame());


// }

BOOST_AUTO_TEST_CASE(find_cauchy_samples_cycle_num_1_for_cauchy_class_test)
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



	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr(".1");

	sample << mpfr("7.999999999999999e-01", "2.168404344971009e-19"); // 


	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

	My_Endgame.FindCauchySamples(time,sample);

	// std::cout << "first cauchy approx is " << first_cauchy_approx << '\n';

	BOOST_CHECK((My_Endgame.cauchy_samples_.back() - My_Endgame.cauchy_samples_.front()).norm() < My_Endgame.GetTrackToleranceDuringEndgame());
	BOOST_CHECK(My_Endgame.cauchy_samples_.size() == 4);
	BOOST_CHECK(My_Endgame.endgame_settings_.cycle_number == 1);
}

BOOST_AUTO_TEST_CASE(find_cauchy_samples_cycle_num_greater_than_1_for_cauchy_class_test)
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




	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr("0.1");

	sample << mpfr("9.000000000000001e-01", "4.358898943540673e-01"); // 


	bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);

	My_Endgame.FindCauchySamples(time,sample);

	BOOST_CHECK((My_Endgame.cauchy_samples_.back() - My_Endgame.cauchy_samples_.front()).norm() < My_Endgame.GetTrackToleranceDuringEndgame());
	BOOST_CHECK(My_Endgame.cauchy_samples_.size() == 7);
	BOOST_CHECK(My_Endgame.endgame_settings_.cycle_number == 2);
}

BOOST_AUTO_TEST_SUITE_END()
