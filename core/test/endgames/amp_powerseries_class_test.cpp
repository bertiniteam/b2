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
//  copyright 2015, 2015
//
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS


#include <iostream>
#include <boost/test/unit_test.hpp>

#include "bertini2/tracking/amp_powerseries_endgame.hpp"
#include "bertini2/start_system.hpp"



using System = bertini::System;
using Variable = bertini::node::Variable;


using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


// using dbl = std::complex<double>;
// using mpfr = bertini::complex;
// using mpfr_float = boost::multiprecision::mpfr_float;

using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;

extern double threshold_clearance_d;
extern boost::multiprecision::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;


BOOST_AUTO_TEST_SUITE(amp_powerseries_endgame)

using namespace bertini::tracking;

BOOST_AUTO_TEST_CASE( basic_hermite_test_case_against_matlab_mp )
{
	/*This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 
	The test case will construct 3 samples with derivative,time, and space values for the function x^8 + 1.
	After the three samples have been constructed there will be a hermite interpolation. 
	We check this against the tracking tolerance for the endgame. 
	*/
	mpfr_float::default_precision(30);



	mpfr target_time(0,0); //our target time is the origin.
	unsigned int num_samples = 3;

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time;
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("1.00000001"); // f(.1) = 1.00000001
	samples.push_back(sample);
	derivative << mpfr("8e-7"); //f'(.1) = 8e-7
	derivatives.push_back(derivative);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("1.0000000000390625"); //f(.05) = 1.0000000000390625
	samples.push_back(sample);
	derivative << mpfr("6.25e-9"); //f'(.05) = 6.25e-9
	derivatives.push_back(derivative);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("1.000000000000152587890625"); // f(.025) = 1.000000000000152587890625
	samples.push_back(sample);
	derivative << mpfr("4.8828125e-11"); //f'(.025) = 4.8828125e-11
	derivatives.push_back(derivative);

	Vec< mpfr > first_approx = endgame::HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);


	BOOST_CHECK( (first_approx(0) - mpfr("0.9999999767578209232082898114211261253459","0")).norm() < 1e-7); 
	// answer was found using matlab for a check. difference is diff is 2.32422e-08

}//end basic hermite test case mp against matlab


BOOST_AUTO_TEST_CASE(hermite_test_case_mp_for_powerseries_class)
{
	/*This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 

	The test case will construct 3 samples with derivative,time, and space values for the function x^8 + 1.

	After the three samples have been constructed there will be a hermite interpolation. 

	Next, a new sample is constructed and the earliest sample is discarded. Leaving us three samples that are "nearer"
	to the target at the origin. 

	A new approximation is made, with the a new sample and approximation done afterwards. 

	We then check to make sure our approximations are getting better by checking the distance from the correct answer and the 
	various approximations made. 

	*/
	mpfr_float::default_precision(30);

	mpfr target_time(0,0);
	unsigned int num_samples = 3;

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time;
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("1.00000001"); // f(.1) = 1.00000001
	samples.push_back(sample);
	derivative << mpfr("8e-7"); //f'(.1) = 8e-7
	derivatives.push_back(derivative);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("1.0000000000390625"); //f(.05) = 1.0000000000390625
	samples.push_back(sample);
	derivative << mpfr("6.25e-9"); //f'(.05) = 6.25e-9
	derivatives.push_back(derivative);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("1.000000000000152587890625"); // f(.025) = 1.000000000000152587890625
	samples.push_back(sample);
	derivative << mpfr("4.8828125e-11"); //f'(.025) = 4.8828125e-11
	derivatives.push_back(derivative);


	 Vec< mpfr > first_approx = endgame::HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);
	 Vec< mpfr > correct(1);
	 correct << mpfr("1");


	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("1.00000000000000059604644775390625"); // f(.0125) = 1.00000000000000059604644775390625
	samples.push_back(sample);
	derivative << mpfr("3.814697265625e-13"); //f'(.0125) = 3.814697265625e-13
	derivatives.push_back(derivative);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();
	derivatives.pop_front();

	//Compute the second approximation.
	Vec< mpfr > second_approx = endgame::HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);


	// //Check to make sure we are doing better. 
	BOOST_CHECK(abs(second_approx(0)-correct(0)) < abs(first_approx(0)-correct(0)));

	//Setting up new sample for use in approximation.
	time = mpfr("0.00625"); //.0125/2 = 0.00625
	times.push_back(time);
	sample << mpfr("1.0000000000000000023283064365386962890625"); // f(.00625) = 1.0000000000000000023283064365386962890625
	samples.push_back(sample);
	derivative << mpfr("2.98023223876953125000000000000000e-15"); //f'(.00625) = 2.98023223876953125000000000000000×e-15
	derivatives.push_back(derivative);

	times.pop_front();
	samples.pop_front();
	derivatives.pop_front();


	Vec< mpfr > third_approx = endgame::HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);


	BOOST_CHECK(abs(third_approx(0)-correct(0)) < abs(second_approx(0)-correct(0)));

}//end hermite test case mp

BOOST_AUTO_TEST_CASE(hermite_test_case_dbl_for_powerseries_class)
{
	/*This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 

	The test case will construct 3 samples with derivative,time, and space values for the function x^8 + 1.

	After the three samples have been constructed there will be a hermite interpolation. 

	Next, a new sample is constructed and the earliest sample is discarded. Leaving us three samples that are "nearer"
	to the target at the origin. 

	A new approximation is made, with the a new sample and approximation done afterwards. 

	We then check to make sure our approximations are getting better by checking the distance from the correct answer and the 
	various approximations made. 

	*/

	mpfr_float::default_precision(16);



	mpfr target_time(0,0);
	unsigned int num_samples = 3;

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time;
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("1.00000001"); // f(.1) = 1.00000001
	samples.push_back(sample);
	derivative << mpfr("8e-7"); //f'(.1) = 8e-7
	derivatives.push_back(derivative);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("1.0000000000390625"); //f(.05) = 1.0000000000390625
	samples.push_back(sample);
	derivative << mpfr("6.25e-9"); //f'(.05) = 6.25e-9
	derivatives.push_back(derivative);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("1.000000000000152587890625"); // f(.025) = 1.000000000000152587890625
	samples.push_back(sample);
	derivative << mpfr("4.8828125e-11"); //f'(.025) = 4.8828125e-11
	derivatives.push_back(derivative);


	 Vec< mpfr > first_approx = endgame::HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);
	 Vec< mpfr > correct(1);
	 correct << mpfr("1");


	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("1.00000000000000059604644775390625"); // f(.0125) = 1.00000000000000059604644775390625
	samples.push_back(sample);
	derivative << mpfr("3.814697265625e-13"); //f'(.0125) = 3.814697265625e-13
	derivatives.push_back(derivative);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();
	derivatives.pop_front();

	//Compute the second approximation.
	Vec< mpfr > second_approx = endgame::HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);


	// //Check to make sure we are doing better. 
	BOOST_CHECK(abs(second_approx(0)-correct(0)) < abs(first_approx(0)-correct(0)));

	//Setting up new sample for use in approximation.
	time = mpfr("0.00625"); //.0125/2 = 0.00625
	times.push_back(time);
	sample << mpfr("1.0000000000000000023283064365386962890625"); // f(.00625) = 1.0000000000000000023283064365386962890625
	samples.push_back(sample);
	derivative << mpfr("2.98023223876953125000000000000000e-15"); //f'(.00625) = 2.98023223876953125000000000000000×e-15
	derivatives.push_back(derivative);

	times.pop_front();
	samples.pop_front();
	derivatives.pop_front();

	Vec< mpfr > third_approx = endgame::HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);

	// //Make sure we are doing better at approximating. 
	BOOST_CHECK(abs(third_approx(0)-correct(0)) < abs(second_approx(0)-correct(0)));

}//end hermite test case dbl

BOOST_AUTO_TEST_CASE(compute_bound_on_cycle_num_mp_for_powerseries_class)
{
	/* In the power series endgame there is a bound on that is calculated that will be used as a higher
	bound for the exhaustive search of the best cycle number. 

	This is calcuated using a modified version of the cycle test in the bertini book, page 53. 

	The bound calculated is then compared against the user defined setting MaxCycleNum which is default at 6. 

	We take a function (x-1)^3 because we know it should have a cycle number of 3. (This is true but will not be true near 0 or near 0.1).

	This is because the path will not look cubic globally. 
	*/
	
	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	sys.AddFunction(pow(x - 1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time;
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);


	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);


	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);

	config::Endgame<mpfr_float> endgame_settings;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.ComputeBoundOnCycleNumber<mpfr>();


	BOOST_CHECK(my_endgame.UpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

	auto first_upper_bound = my_endgame.UpperBoundOnCycleNumber();

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);


	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.ComputeBoundOnCycleNumber<mpfr>();


	BOOST_CHECK(my_endgame.UpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

} // end compute bound on cycle number mp

BOOST_AUTO_TEST_CASE(compute_bound_on_cycle_num_dbl_for_powerseries_class)
{
	/* In the power series endgame there is a bound on that is calculated that will be used as a higher
	bound for the exhaustive search of the best cycle number. 

	This is calcuated using a modified version of the cycle test in the bertini book, page 53. 

	The bound calculated is then compared against the user defined setting MaxCycleNum which is default at 6. 

	We take a function (x-1)^3 because we know it should have a cycle number of 3. (This is true but will not be true near 0 or near 0.1).

	This is because the path will not look cubic globally. 
	*/
	
	mpfr_float::default_precision(16);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	sys.AddFunction(pow(x - 1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.


	mpfr time;
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);


	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);


	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);


	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,security_settings,tolerances);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.ComputeBoundOnCycleNumber<mpfr>();


	BOOST_CHECK(my_endgame.UpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

	auto first_upper_bound = my_endgame.UpperBoundOnCycleNumber();

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);


	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();


	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.ComputeBoundOnCycleNumber<mpfr>();


	BOOST_CHECK(my_endgame.UpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

} // end compute bound on cycle number dbl


BOOST_AUTO_TEST_CASE(compute_cycle_number_test_mp_for_powerseries_class)
{
	/*
	Once we have an upper bound for the exhaustive search on the correct cycle number, see above test cases. 

	The next step is to take each possible cycle number, take dx/dt -> dx/ds and t -> s, where t = s^c where c is the proposed 
	cycle number. 

	With each conversion we do a Hermite interpolation with the first n-1 samples we have to predict the nth sample.

	The best approximation determines which proposed cycle number we will use to run a Hermite interpolation with the n samples to 
	approximate at the origin.
	*/
	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t");
	sys.AddFunction( pow(x-1,3) );  //f(x) = (x-1)^3 

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.


	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);


	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);


	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);



	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);


	config::PowerSeries power_series_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,power_series_settings,tolerances);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.ComputeCycleNumber<mpfr>();

	BOOST_CHECK(my_endgame.CycleNumber() == 1);

} // end compute cycle number mp

BOOST_AUTO_TEST_CASE(compute_cycle_number_test_dbl_for_powerseries_class)
{	
	/*
	Once we have an upper bound for the exhaustive search on the correct cycle number, see above test cases. 

	The next step is to take each possible cycle number, take dx/dt -> dx/ds and t -> s, where t = s^c where c is the proposed 
	cycle number. 

	With each conversion we do a Hermite interpolation with the first n-1 samples we have to predict the nth sample.

	The best approximation determines which proposed cycle number we will use to run a Hermite interpolation with the n samples to 
	approximate at the origin.
	*/

	mpfr_float::default_precision(16);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t");
	sys.AddFunction( pow(x-1,3) );  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.


	mpfr time;
	Vec<mpfr> sample(1);


	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);


	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);


	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);


	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);


	config::PowerSeries power_series_settings;
	config::Security<mpfr_float> security_settings;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,power_series_settings,security_settings);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);


	auto derivatives = my_endgame.ComputeCycleNumber<mpfr>();

	BOOST_CHECK(my_endgame.CycleNumber() == 1);

} // end compute cycle number mp

BOOST_AUTO_TEST_CASE(compute_approximation_of_x_at_t0_mp_for_powerseries_class)
{
	
	/*
	Compute approximation at origin using three sample points. 

	Then, find a new sample and remove earliest known sample. 
	
	Compute a new approximation with the current three samples. 
	
	Do this till we have three approximations. Check to see if the third approximation is better than second. 
	*/

	mpfr_float::default_precision(30);


	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x - mpfr(1),3)*(1-t) + (pow(x,3) + mpfr(1))*t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	mpfr origin = mpfr(0,0);
	Vec<mpfr> x_origin(1);
	x_origin << mpfr(1,0);
	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.



	mpfr time;
	Vec<mpfr> sample(1);
	// Vec<mpfr> derivative(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("0.50000000000000007812610562824908678293817200689660e0", "0.90818521453245102009102405377246269374768541575725e-16"); // f(.1) = 0.50000000000000007812610562824908678293817200689660e0 0.90818521453245102009102405377246269374768541575725e-16 from bertini classic
	samples.push_back(sample);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("0.60000000000000000073301140774132693211475245208482e0", "0.85317043225116681251211164764202768939258706233715e-18"); //f(.05) = 0.60000000000000000073301140774132693211475245208482e0 0.85317043225116681251211164764202768939258706233715e-18 from bertini classic.
	samples.push_back(sample);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("0.67729059415987117534436422700325955211292174181447e0", "0.38712848412230230944976856052769427012481164605204e-16"); // f(.025) = 0.67729059415987117534436422700325955211292174181447e0 0.38712848412230230944976856052769427012481164605204e-16 from bertini classic
	samples.push_back(sample);

	config::Endgame<mpfr_float> endgame_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings,tolerances);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	Vec<mpfr> first_approx;
	auto code = my_endgame.ComputeApproximationOfXAtT0(first_approx, origin);
	BOOST_CHECK(code==SuccessCode::Success);
	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("0.73905643972939615063207159129047946977060596291527e0", "0.44983246338361191567539019211879692583079653921201e-18"); // f(.0125) = 0.73905643972939615063207159129047946977060596291527e0 0.44983246338361191567539019211879692583079653921201e-18
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	Vec<mpfr> second_approx;
	code = my_endgame.ComputeApproximationOfXAtT0(second_approx,origin);
	BOOST_CHECK(code==SuccessCode::Success);
	//Setting up a new sample for approximation.
	time = mpfr(".00625"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("0.78910678115791459147153183840413839746566925123387e0", "0.22640341504967423865456128414532605156908222616485e-16"); // f(.00625) = 0.78910678115791459147153183840413839746566925123387e0 0.22640341504967423865456128414532605156908222616485e-16
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	Vec<mpfr> third_approx;
	code = my_endgame.ComputeApproximationOfXAtT0(third_approx, origin);
	BOOST_CHECK(code==SuccessCode::Success);
	BOOST_CHECK((third_approx - x_origin).norm() < (second_approx - x_origin).norm());

} // end compute approximation of x at t0 mp

BOOST_AUTO_TEST_CASE(compute_approximation_of_x_at_t0_dbl_for_powerseries_class)
{
	/* 
	Compute approximation at origin using three sample points. 

	Then, find a new sample and remove earliest known sample. 

	Compute a new approximation with the current three samples. 

	Do this till we have three approximations. Check to see if the third approximation is better than second. 
	*/
	mpfr_float::default_precision(16);


	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x - mpfr(1),3)*(1-t) + (pow(x,3) + mpfr(1))*t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	mpfr origin = mpfr(0,0);
	Vec<mpfr> x_origin(1);
	x_origin << mpfr(1,0);
	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.


	mpfr time(1);
	Vec<mpfr> sample(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr(".500000000000000","0"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic
	samples.push_back(sample);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr(".6000000000000000","0"); //f(.05) = 6.000000000000000e-01 8.165397611531455e-19 from bertini classic.
	samples.push_back(sample);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr(".6772905941598711", "0"); // f(.025) = 6.772905941598711e-01 3.869924129415447e-17 from bertini classic
	samples.push_back(sample);

	config::Endgame<mpfr_float> endgame_settings;
	config::Security<mpfr_float> security_settings;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings,security_settings);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);


	Vec<mpfr> first_approx;
	auto code = my_endgame.ComputeApproximationOfXAtT0(first_approx, origin);
	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr(".7390564397293962", "0"); // f(.0125) = 7.390564397293962e-01 4.641740550953566e-19
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	Vec<mpfr> second_approx;
	code = my_endgame.ComputeApproximationOfXAtT0(second_approx, origin);

	time = mpfr(".00625"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("0.78910678115791459147153183840413839746566925123387e0", "0.22640341504967423865456128414532605156908222616485e-16"); // f(.00625) = 0.78910678115791459147153183840413839746566925123387e0 0.22640341504967423865456128414532605156908222616485e-16
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	Vec<mpfr> third_approx;
	code = my_endgame.ComputeApproximationOfXAtT0(third_approx, origin);


	BOOST_CHECK((third_approx - x_origin).norm() < (second_approx - x_origin).norm());

} // end compute approximation of x at t0 mp

BOOST_AUTO_TEST_CASE(compute_initial_samples_mp_for_powerseries_class)
{
	/*
	Given a time value (usually 0.1) and a vector representing values for all other variables the first step in the power series 
	endgame is to get some initial samples to make the first hermite interpolation. 

	For this reason there exists a ComputeInitialSamples function. 

	This test case checks the computed samples vs. samples computed by hand with Matlab to see that they are within the 
	track tolerance during the endgame. 

	*/
	mpfr_float::default_precision(30);



	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x - 1,3)*(1-t) + (pow(x,3) + 1)*t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	mpfr sample_factor = mpfr(".5");
	mpfr origin = mpfr("0");
	Vec<mpfr> x_origin(1);
	x_origin << mpfr(1,0);
	std::deque<mpfr> correct_times; //times are not vectors they are just complex numbers, these are the correct times from bertini
	std::deque< Vec<mpfr> > correct_samples; //samples are space values that may be a vector of complex numbers, correct space values 
											// from bertini.

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque<Vec<mpfr> > samples;
	mpfr time(1);
	Vec<mpfr> sample(1);


	time = mpfr(".1"); // x = .1
	correct_times.push_back(time);
	sample << mpfr("5.000000000000001e-01", "9.084258952712920e-17"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic
	correct_samples.push_back(sample);

	time = mpfr(".05"); // x = .1/2 = .05
	correct_times.push_back(time);
	sample << mpfr("6.000000000000000e-01", "8.165397611531455e-19"); //f(.05) = 6.000000000000000e-01 8.165397611531455e-19 from bertini classic.
	correct_samples.push_back(sample);

	time = mpfr(".025"); // x = .05/2 = .025
	correct_times.push_back(time);
	sample << mpfr("6.772905941598711e-01", "3.869924129415447e-17"); // f(.025) = 6.772905941598711e-01 3.869924129415447e-17 from bertini classic
	correct_samples.push_back(sample);

	unsigned int num_samples = 3; 

	mpfr current_time(1);
	Vec<mpfr> current_space(1);
	current_time = mpfr(".1");
	current_space << mpfr("5.000000000000001e-01", "9.084258952712920e-17");

	config::Endgame<mpfr_float> endgame_settings;
	config::PowerSeries power_series_settings;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings,power_series_settings);

	my_endgame.ComputeInitialSamples(current_time, current_space, times, samples);

	

	for(unsigned ii = 0; ii < samples.size(); ++ii)
	{
		BOOST_CHECK((samples[ii] - correct_samples[ii]).norm() < my_endgame.Tolerances().newton_during_endgame);
	}

}//end compute initial samples mp

BOOST_AUTO_TEST_CASE(compute_initial_samples_dbl_for_powerseries_class)
{
	/*
	Given a time value (usually 0.1) and a vector representing values for all other variables the first step in the power series 
	endgame is to get some initial samples to make the first hermite interpolation. 

	For this reason there exists a ComputeInitialSamples function. 

	This test case checks the computed samples vs. samples computed by hand with Matlab to see that they are within the 
	track tolerance during the endgame. 

	*/
	mpfr_float::default_precision(16);



	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3)+1)*t);


	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	mpfr sample_factor = mpfr(".5");
	mpfr origin = mpfr("0");
	Vec<mpfr> x_origin(1);
	x_origin << mpfr(1,0);
	std::deque<mpfr> correct_times; //times are not vectors they are just complex numbers, these are the correct times from bertini
	std::deque< Vec<mpfr> > correct_samples; //samples are space values that may be a vector of complex numbers, correct space values 
											// from bertini.

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples;
	mpfr time(1);
	Vec<mpfr> sample(1);

	time = mpfr(".1"); // x = .1
	correct_times.push_back(time);
	sample << mpfr("5.000000000000001e-01", "9.084258952712920e-17"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic
	correct_samples.push_back(sample);

	time = mpfr(".05"); // x = .1/2 = .05
	correct_times.push_back(time);
	sample << mpfr("6.000000000000000e-01", "8.165397611531455e-19"); //f(.05) = 6.000000000000000e-01 8.165397611531455e-19 from bertini classic.
	correct_samples.push_back(sample);

	time = mpfr(".025"); // x = .05/2 = .025
	correct_times.push_back(time);
	sample << mpfr("6.772905941598711e-01", "3.869924129415447e-17"); // f(.025) = 6.772905941598711e-01 3.869924129415447e-17 from bertini classic
	correct_samples.push_back(sample);

	unsigned int num_samples = 3; 

	mpfr current_time(1);
	Vec<mpfr> current_space(1);
	current_time = mpfr(".1");
	current_space << mpfr("5.000000000000001e-01", "9.084258952712920e-17");

	config::PowerSeries power_series_settings;
	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,power_series_settings,security_settings,tolerances);


	my_endgame.ComputeInitialSamples(current_time, current_space, times, samples);



	for(unsigned ii = 0; ii < samples.size(); ++ii)
	{
		BOOST_CHECK((samples[ii] - correct_samples[ii]).norm() < my_endgame.Tolerances().newton_during_endgame);
	}

}//end compute initial samples dbl

BOOST_AUTO_TEST_CASE(pseg_mp_for_powerseries_class)
{
	/*
	The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
	one of the solutions at t = endgame_time. 

	This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
	the endgame. 
	*/

	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");

	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3)+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);


	mpfr current_time(1);
	Vec<mpfr> current_space(1);
	current_time = mpfr(".1");
	current_space << mpfr("5.000000000000001e-01", "9.084258952712920e-17");

	Vec<mpfr> correct(1);
	correct << mpfr(1,0);

	config::Endgame<mpfr_float> endgame_settings;
	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings,security_settings,tolerances);
	my_endgame.PSEG(current_time,current_space);


	BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - correct).norm() < 1e-11);//my_endgame.GetTrackToleranceDuringEndgame());

}//end pseg mp for power series class 

BOOST_AUTO_TEST_CASE(pseg_dbl_for_powerseries_class)
{
	/*
	The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
	one of the solutions at t = endgame_time. 

	This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
	the endgame. 
	*/

	mpfr_float::default_precision(16);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3)+1)*t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);

	mpfr current_time(1);
	Vec<mpfr> current_space(1);
	current_time = mpfr(".1");
	current_space << mpfr("5.000000000000001e-01", "9.084258952712920e-17");

	Vec<mpfr> correct(1);
	correct << mpfr(1,0);

	config::Endgame<mpfr_float> endgame_settings;
	config::PowerSeries power_series_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings,power_series_settings,tolerances);


	my_endgame.PSEG(current_time,current_space);

	BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - correct).norm() < 1e-11);//my_endgame.GetTrackToleranceDuringEndgame());
	
}//end pseg dbl for power series class






BOOST_AUTO_TEST_CASE(cycle_num_2_example)
{
	/*
	The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
	one of the solutions at t = endgame_time. 

	This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
	the endgame. 
	*/

	mpfr_float::default_precision(30);

	System sys;
	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2)-1)*t);


	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);


	mpfr start_time(1);
	Vec<mpfr> start_point(1); start_point << mpfr(1);
	mpfr t_endgame_boundary("0.1");
	
	Vec<mpfr> eg_boundary_point;
	auto init_success = tracker.TrackPath(eg_boundary_point, start_time, t_endgame_boundary, start_point);

	BOOST_CHECK(init_success==SuccessCode::Success);
	BOOST_CHECK(abs(eg_boundary_point(0) - 1)< 1e-5);

	Vec<mpfr> correct_root(1);
	correct_root << mpfr(1,0);

	config::Endgame<mpfr_float> endgame_settings;
	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings,security_settings,tolerances);
	my_endgame.PSEG(t_endgame_boundary,eg_boundary_point);

	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(),1);//my_endgame.
	BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - correct_root).norm() < 1e-11);//my_endgame.GetTrackToleranceDuringEndgame());

}//end pseg mp for power series class 








BOOST_AUTO_TEST_CASE(pseg_mp_for_powerseries_class_multiple_variables)
{
	/*
	The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
	one of the solutions at t = endgame_time. 

	In this test we do multiple variables decoupled, that has a high multiplicity (5) solution. 

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
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);


	mpfr current_time(1);
	Vec<mpfr> current_space(2);
	current_time = mpfr(".1");
	current_space <<  mpfr("5.000000000000001e-01", "9.084258952712920e-17") ,mpfr("9.000000000000001e-01","4.358898943540673e-01");

	Vec<mpfr> correct(2);
	correct << mpfr(1,0),mpfr(1,0);

	config::Endgame<mpfr_float> endgame_settings;
	config::PowerSeries power_series_settings;
	config::Security<mpfr_float> security_settings;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings,power_series_settings,security_settings);

	my_endgame.PSEG(current_time,current_space);

	BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - correct).norm() < 1e-10);//my_endgame.GetTrackToleranceDuringEndgame());

}//end pseg mp test case for power series class

BOOST_AUTO_TEST_CASE(griewank_osborne_for_powerseries_class_griewank_osborne)
{
	
	// The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
	// one of the solutions at t = endgame_time. 

	// Griewank Osborne is a very classic example. Here we allow x and y to mix. There are six paths to be tracked and we know there values 
	// at t = 0.1. 

	// Three of these paths will converge to origin and three will diverge to infinity triggering a SecurityMaxNorm issue. 

	// This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
	// the endgame. 
	

	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t"), y = std::make_shared<Variable>("y");
	VariableGroup vars{x,y};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	sys.AddFunction(((mpfr("29")/mpfr("16"))*pow(x,3)-2*x*y)*(1-t) + (pow(x,3) - 1)*t);
	sys.AddFunction((y - pow(x,2))*(1-t) + (pow(y,2) - 1)*t);

	auto AMP = config::AMPConfigFrom(sys);

	AMPTracker tracker(sys);
	
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;

	tracker.Setup(config::Predictor::HeunEuler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_settings,
                newton_settings);
	
	tracker.AMPSetup(AMP);


	mpfr endgame_time(1);
	Vec<mpfr> current_space_1(2);
	Vec<mpfr> current_space_2(2);
	Vec<mpfr> current_space_3(2);
	Vec<mpfr> current_space_4(2);
	Vec<mpfr> current_space_5(2);
	Vec<mpfr> current_space_6(2);

	endgame_time = mpfr(".1");

	current_space_1 <<  mpfr("1.028694756284462e-01", "-9.661822229074768e-01"), mpfr("-8.937287306314232e-01", "-2.480445481975048e-01");
	current_space_2 <<  mpfr("-5.219899550566304e-01", "-2.212788871407134e-17"), mpfr("3.684968541245056e-01", "2.471980953266950e-17");
	current_space_3 <<  mpfr("1.028694756284466e-01", "9.661822229074761e-01"), mpfr("-8.937287306314221e-01", "2.480445481975040e-01");
	current_space_4 <<  mpfr("6.098408897464429e-03", "1.058791184067875e-21"), mpfr("-9.109808533477256e+00", "-2.374402757743255e-17");
	current_space_5 <<  mpfr("1.220071827679809e+00", "7.657177843178875e-19"), mpfr("1.386185299689565e+00", "2.852806966352484e-18");
	current_space_6 <<  mpfr("-9.099192327775354e-01", "2.114194236346734e-17"), mpfr("8.573852849693505e-01", "-2.164338586824188e-17");

	std::vector<Vec<mpfr> > current_space_values;

	current_space_values.push_back(current_space_1);
	current_space_values.push_back(current_space_2);
	current_space_values.push_back(current_space_3);
	current_space_values.push_back(current_space_4);
	current_space_values.push_back(current_space_5);
	current_space_values.push_back(current_space_6);


	Vec<mpfr> correct(2);
	correct << mpfr(0,0),mpfr(0,0);

	config::Endgame<mpfr_float> endgame_settings;
	config::PowerSeries power_series_settings;
	config::Security<mpfr_float> security_settings;
	config::Tolerances<mpfr_float> tolerances;

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker,endgame_settings,power_series_settings,security_settings,tolerances);

	unsigned num_paths_diverging = 0;
	unsigned num_paths_converging = 0;
	for (const auto& s : current_space_values)
	{
		mpfr_float::default_precision(30);
		SuccessCode endgame_success = my_endgame.PSEG(endgame_time,s);
		if(endgame_success == SuccessCode::Success){
			BOOST_CHECK((my_endgame.FinalApproximation<mpfr>() - correct).norm() < 1e-11);// my_endgame.GetTrackToleranceDuringEndgame());
			num_paths_converging++;
		}
		if(endgame_success == SuccessCode::SecurityMaxNormReached){
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

BOOST_AUTO_TEST_CASE(total_degree_start_system_powerseries_class_used_with_AMP)
{
	/*
	In this example we take a decoupled system, homogenize and patch it. Track to endgame boundary and then run our endgame on the space
	values we have. 

	*/
	unsigned ambient_precision = 16;
	using namespace bertini::tracking;
	mpfr_float::default_precision(ambient_precision);

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
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;
	tracker.Setup(config::Predictor::HeunEuler,
	              	mpfr_float("1e-5"), mpfr_float("1e5"),
					stepping_settings, newton_settings);

	tracker.AMPSetup(AMP);
	

	mpfr t_start(1), t_endgame_boundary(0.1);
	std::vector<Vec<mpfr> > homogenized_solutions;
	for (unsigned ii = 0; ii < 1; ++ii)
	{
		mpfr_float::default_precision(ambient_precision);
		final_system.precision(ambient_precision);
		auto start_point = TD.StartPoint<mpfr>(ii);

		Vec<mpfr> result;
		SuccessCode tracking_success;

		tracking_success = tracker.TrackPath(result,t_start,t_endgame_boundary,start_point);
		BOOST_CHECK(tracking_success==SuccessCode::Success);

		homogenized_solutions.push_back(result);
	}

	Vec<mpfr> correct(2);
	correct << mpfr(1,0),mpfr(1,0);

	tracker.Setup(config::Predictor::HeunEuler,
	              	mpfr_float("1e-6"), mpfr_float("1e5"),
					stepping_settings, newton_settings);

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker);


	std::vector<Vec<mpfr> > endgame_solutions;

	unsigned num_successful_occurences = 0;
	unsigned num_min_track_time_reached = 0;
	for (auto const& s : homogenized_solutions)
	{
		SuccessCode endgame_success = my_endgame.PSEG(t_endgame_boundary,s);
		if(endgame_success == SuccessCode::Success)
		{
				num_successful_occurences++;

		}
		if(endgame_success == SuccessCode::MinStepSizeReached)
		{
			num_min_track_time_reached++;
		}
	}

 	BOOST_CHECK_EQUAL(num_successful_occurences,6);
 	BOOST_CHECK_EQUAL(num_min_track_time_reached,0);

}


BOOST_AUTO_TEST_CASE(parabola)
{
	/*
	In this example we take a decoupled system, homogenize and patch it. Track to endgame boundary and then run our endgame on the space
	values we have. 

	*/
	unsigned ambient_precision = 30;
	using namespace bertini::tracking;
	mpfr_float::default_precision(ambient_precision);

	Var x = std::make_shared<Variable>("x");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{x};

	sys.AddVariableGroup(v);

	sys.AddFunction(pow(x,2) - t);
	sys.AddPathVariable(t);
	Vec<mpfr> start_point(1);
	start_point << mpfr(1);

	auto AMP = config::AMPConfigFrom(sys);



	auto tracker = AMPTracker(sys);
	config::Stepping<mpfr_float> stepping_settings;
	config::Newton newton_settings;
	tracker.Setup(config::Predictor::HeunEuler,
	              	mpfr_float("1e-5"), mpfr_float("1e5"),
					stepping_settings, newton_settings);

	tracker.AMPSetup(AMP);
	

	mpfr t_start(1), t_endgame_boundary("0.1");

	Vec<mpfr> soln_at_EG_bdry;

	auto tracking_success = tracker.TrackPath(soln_at_EG_bdry,t_start,t_endgame_boundary,start_point);
	BOOST_CHECK(tracking_success==SuccessCode::Success);



	Vec<mpfr> correct_eg_soln(1);
	correct_eg_soln << mpfr(0,0);

	tracker.Setup(config::Predictor::HeunEuler,
	              	mpfr_float("1e-6"), mpfr_float("1e5"),
					stepping_settings, newton_settings);

	endgame::PowerSeriesEndgame<AMPTracker> my_endgame(tracker);


	auto endgame_success = my_endgame.PSEG(t_endgame_boundary,soln_at_EG_bdry);

	auto endgame_solution = my_endgame.FinalApproximation<mpfr>();
	BOOST_CHECK(abs(endgame_solution(0) - correct_eg_soln(0)) < 1e-11);
}

BOOST_AUTO_TEST_SUITE_END()
