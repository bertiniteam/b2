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
#include "bertini2/tracking/powerseries_endgame.hpp"
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


BOOST_AUTO_TEST_SUITE(powerseries_endgame_class_basics)

BOOST_AUTO_TEST_CASE( basic_hermite_test_case_against_matlab_mp )
{
	/*This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 
	The test case will construct 3 samples with derivative,time, and space values for the function x^8 + 1.
	After the three samples have been constructed there will be a hermite interpolation. 
	We check this against the tracking tolerance for the endgame. 
	*/
	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddFunction( pow(x,8) + 1 );  //f(x) = x^8 + 1


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

	mpfr target_time(1);
	target_time = mpfr("0","0"); //our target time is the origin.
	unsigned int num_samples = 3;

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time(1);
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

	bertini::tracking::config::Tolerances endgame_tolerances_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_tolerances_struct);
	My_Endgame.SetSTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetSDerivatives(derivatives);

	Vec< mpfr > first_approx = My_Endgame.HermiteInterpolateAndSolve(target_time,num_samples);
	Vec< mpfr > matlab_first_approx(1);
	matlab_first_approx << mpfr("1");



	BOOST_CHECK(abs(first_approx.norm() - mpfr("0.9999999767578209232082898114211261253459","0").norm()) < My_Endgame.GetTrackToleranceDuringEndgame()); // answer was found using matlab for a check. difference is diff is 2.32422e-08

}//end basic hermite test case mp


BOOST_AUTO_TEST_CASE(hermite_test_case_mp_for_powerseries_class)
{
	/*This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 

	The test case will construct 3 samples with derivative,time, and space values for the function x^8e-7 + 1.

	After the three samples have been constructed there will be a hermite interpolation. 

	Next, a new sample is constructed and the earliest sample is discarded. Leaving us three samples that are "nearer"
	to the target at the origin. 

	A new approximation is made, with the a new sample and approximation done afterwards. 

	We then check to make sure our approximations are getting better by checking the distance from the correct answer and the 
	various approximations made. 

	*/
	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	sys.AddFunction( pow(x,8) + 1 );  //f(x) = x^8 + 1

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

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

	mpfr target_time(1);
	target_time = mpfr("0","0"); //our target time is the origin.
	unsigned int num_samples = 3;

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time(1);
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

	bertini::tracking::config::Security endgame_security_struct;

	//First Approximation at the origin.
	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_security_struct);
	My_Endgame.SetSTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetSDerivatives(derivatives);

	 Vec< mpfr > first_approx = My_Endgame.HermiteInterpolateAndSolve(target_time,num_samples);
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

	My_Endgame.SetSTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetSDerivatives(derivatives);


	//Compute the second approximation.
	Vec< mpfr > second_approx = My_Endgame.HermiteInterpolateAndSolve(target_time,num_samples);


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

	My_Endgame.SetSTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetSDerivatives(derivatives);

	Vec< mpfr > third_approx = My_Endgame.HermiteInterpolateAndSolve(target_time,num_samples);


	BOOST_CHECK(abs(third_approx(0)-correct(0)) < abs(second_approx(0)-correct(0)));

}//end hermite test case mp

BOOST_AUTO_TEST_CASE(hermite_test_case_dbl_for_powerseries_class)
{
	// 	/*This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 

	// 	The test case will construct 3 samples with derivative,time, and space values for the function x^8e-7 + 1.

	// 	After the three samples have been constructed there will be a hermite interpolation. 

	// 	Next, a new sample is constructed and the earliest sample is discarded. Leaving us three samples that are "nearer"
	// 	to the target at the origin. 

	// 	A new approximation is made, with the a new sample and approximation done afterwards. 

	// 	We then check to make sure our approximations are getting better by checking the distance from the correct answer and the 
	// 	various approximations made. 

	// 	*/
	// 	//This is used to print out the vectors in maximum precision. 
 	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);
	mpfr_float::default_precision(16);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	sys.AddFunction( pow(x,8) + 1 );  //f(x) = x^8 + 1

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

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

	mpfr target_time(1);
	target_time = mpfr("0","0"); //our target time is the origin.
	unsigned int num_samples = 3;

	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time(1);
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

	bertini::tracking::config::PowerSeries power_series_struct;

	//First Approximation at the origin.
	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,power_series_struct);
	My_Endgame.SetSTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetSDerivatives(derivatives);

	 Vec< mpfr > first_approx = My_Endgame.HermiteInterpolateAndSolve(target_time,num_samples);
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

	My_Endgame.SetSTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetSDerivatives(derivatives);


	//Compute the second approximation.
	Vec< mpfr > second_approx = My_Endgame.HermiteInterpolateAndSolve(target_time,num_samples);


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

	My_Endgame.SetSTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetSDerivatives(derivatives);

	Vec< mpfr > third_approx = My_Endgame.HermiteInterpolateAndSolve(target_time,num_samples);

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
	
	//This is used to print out the vectors in maximum precision. 
	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);
	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x");
	sys.AddFunction(pow(x - 1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

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
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time(1);
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(.1); // x = .1
	times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);
	derivative << mpfr("2.43"); //f'(.1) = 2.43
	derivatives.push_back(derivative);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);
	derivative << mpfr("2.7075"); //f'(.05) = 2.7075
	derivatives.push_back(derivative);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);
	derivative << mpfr("2.851875"); //f'(.025) = 2.851875
	derivatives.push_back(derivative);

	bertini::tracking::config::EndGame endgame_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_struct);
	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetDerivatives(derivatives);

	My_Endgame.BoundOnCycleNumber();


	BOOST_CHECK(My_Endgame.GetUpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

	auto first_upper_bound = My_Endgame.GetUpperBoundOnCycleNumber();

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);
	derivative << mpfr("2.92546875"); //f'(.0125) = 2.92546875 
	derivatives.push_back(derivative);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();
	derivatives.pop_front();

	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetDerivatives(derivatives);

	My_Endgame.BoundOnCycleNumber();


	BOOST_CHECK(My_Endgame.GetUpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

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
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time(1);
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(.1); // x = .1
	times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);
	derivative << mpfr("2.43"); //f'(.1) = 2.43
	derivatives.push_back(derivative);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);
	derivative << mpfr("2.7075"); //f'(.05) = 2.7075
	derivatives.push_back(derivative);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);
	derivative << mpfr("2.851875"); //f'(.025) = 2.851875
	derivatives.push_back(derivative);

	bertini::tracking::config::Security endgame_security_struct;
	bertini::tracking::config::Tolerances endgame_tolerances_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_security_struct,endgame_tolerances_struct);
	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetDerivatives(derivatives);

	My_Endgame.BoundOnCycleNumber();


	BOOST_CHECK(My_Endgame.GetUpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

	auto first_upper_bound = My_Endgame.GetUpperBoundOnCycleNumber();

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);
	derivative << mpfr("2.92546875"); //f'(.0125) = 2.92546875 
	derivatives.push_back(derivative);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();
	derivatives.pop_front();

	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetDerivatives(derivatives);

	My_Endgame.BoundOnCycleNumber();

	BOOST_CHECK(My_Endgame.GetUpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

} // end compute bound on cycle number mp


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
	sys.AddFunction(pow(x - 1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

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
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time(1);
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);
	derivative << mpfr("2.43"); //f'(.1) = 2.43
	derivatives.push_back(derivative);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);
	derivative << mpfr("2.7075"); //f'(.05) = 2.7075
	derivatives.push_back(derivative);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);
	derivative << mpfr("2.851875"); //f'(.025) = 2.851875
	derivatives.push_back(derivative);


	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);
	derivative << mpfr("2.92546875"); //f'(.0125) = 2.92546875
	derivatives.push_back(derivative);

	bertini::tracking::config::PowerSeries power_series_struct;
	bertini::tracking::config::Tolerances endgame_tolerances_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,power_series_struct,endgame_tolerances_struct);
	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetDerivatives(derivatives);
	//My_Endgame.SetUpperBoundOnCycleNumber(6); // this is found from test cases above.

	 My_Endgame.ComputeCycleNumber(time,sample);

	BOOST_CHECK(My_Endgame.GetCycleNumber() == 1);

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
	sys.AddFunction(pow(x - 1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

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
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.

	mpfr time(1);
	Vec<mpfr> sample(1);
	Vec<mpfr> derivative(1);

	time = mpfr(".1"); // x = .1
	times.push_back(time);
	sample << mpfr("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);
	derivative << mpfr("2.43"); //f'(.1) = 2.43
	derivatives.push_back(derivative);

	time = mpfr(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << mpfr("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);
	derivative << mpfr("2.7075"); //f'(.05) = 2.7075
	derivatives.push_back(derivative);

	time = mpfr(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << mpfr("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);
	derivative << mpfr("2.851875"); //f'(.025) = 2.851875
	derivatives.push_back(derivative);

	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);
	derivative << mpfr("2.92546875"); //f'(.0125) = 2.92546875
	derivatives.push_back(derivative);

	bertini::tracking::config::PowerSeries power_series_struct;
	bertini::tracking::config::Security endgame_security_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,power_series_struct,endgame_security_struct);
	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);
	My_Endgame.SetDerivatives(derivatives);
	//My_Endgame.SetUpperBoundOnCycleNumber(6); // this is found from test cases above.

	 My_Endgame.ComputeCycleNumber(time,sample);

	BOOST_CHECK(My_Endgame.GetCycleNumber() == 1);

} // end compute cycle number mp

BOOST_AUTO_TEST_CASE(compute_approximation_of_x_at_t0_mp_for_powerseries_class)
{
	
	// Compute approximation at origin using three sample points. 

	// Then, find a new sample and remove earliest known sample. 

	// Compute a new approximation with the current three samples. 

	//  Do this till we have three approximations. Check to see if the third approximation is better than second. 
	

	mpfr_float::default_precision(30);


	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x - mpfr(1),3)*(1-t) + (pow(x,3) + mpfr(1))*t);

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

	mpfr origin = mpfr("0","0");
	Vec<mpfr> x_origin(1);
	x_origin << mpfr("1","0");
	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.


	mpfr time(1);
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

	bertini::tracking::config::EndGame endgame_struct;
	bertini::tracking::config::Tolerances endgame_tolerances_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_struct,endgame_tolerances_struct);
	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);

	auto first_approx = My_Endgame.ComputeApproximationOfXAtT0(origin);
	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("0.73905643972939615063207159129047946977060596291527e0", "0.44983246338361191567539019211879692583079653921201e-18"); // f(.0125) = 0.73905643972939615063207159129047946977060596291527e0 0.44983246338361191567539019211879692583079653921201e-18
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);

	auto second_approx = My_Endgame.ComputeApproximationOfXAtT0(origin);

	//Setting up a new sample for approximation.
	time = mpfr(".00625"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("0.78910678115791459147153183840413839746566925123387e0", "0.22640341504967423865456128414532605156908222616485e-16"); // f(.00625) = 0.78910678115791459147153183840413839746566925123387e0 0.22640341504967423865456128414532605156908222616485e-16
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);

	auto third_approx = My_Endgame.ComputeApproximationOfXAtT0(origin);

	BOOST_CHECK((third_approx - x_origin).norm() < (second_approx - x_origin).norm());

} // end compute approximation of x at t0 mp

BOOST_AUTO_TEST_CASE(compute_approximation_of_x_at_t0_dbl_for_powerseries_class)
{
	/* Compute approximation at origin using three sample points. 

	// Then, find a new sample and remove earliest known sample. 

	//Compute a new approximation with the current three samples. 

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

	mpfr origin = mpfr("0","0");
	Vec<mpfr> x_origin(1);
	x_origin << mpfr("1","0");
	std::deque<mpfr> times; //times are not vectors they are just complex numbers.
	std::deque< Vec<mpfr> > samples; //samples are space values that may be a vector of complex numbers.
	std::deque< Vec<mpfr> > derivatives; //derivatives are also a vector of complex numbers.


	mpfr time(1);
	Vec<mpfr> sample(1);
	// Vec<mpfr> derivative(1);

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

	bertini::tracking::config::EndGame endgame_struct;
	bertini::tracking::config::Security endgame_security_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_struct,endgame_security_struct);
	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);


	auto first_approx = My_Endgame.ComputeApproximationOfXAtT0(origin);
	//Setting up a new sample for approximation.
	time = mpfr(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr(".7390564397293962", "0"); // f(.0125) = 7.390564397293962e-01 4.641740550953566e-19
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);

	auto second_approx = My_Endgame.ComputeApproximationOfXAtT0(origin);

	time = mpfr(".00625"); //.025/2 = .0125
	times.push_back(time);
	sample << mpfr("0.78910678115791459147153183840413839746566925123387e0", "0.22640341504967423865456128414532605156908222616485e-16"); // f(.00625) = 0.78910678115791459147153183840413839746566925123387e0 0.22640341504967423865456128414532605156908222616485e-16
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	My_Endgame.SetTimes(times);
	My_Endgame.SetSamples(samples);

	auto third_approx = My_Endgame.ComputeApproximationOfXAtT0(origin);


	BOOST_CHECK((third_approx - x_origin).norm() < (second_approx - x_origin).norm());

} // end compute approximation of x at t0 mp

// BOOST_AUTO_TEST_CASE(compute_initial_samples_mp_for_powerseries_class)
// {
// 	/*
// 	Given a time value (usually 0.1) and a vector representing values for all other variables the first step in the power series 
// 	endgame is to get some initial samples to make the first hermite interpolation. 

// 	For this reason there exists a ComputeInitialSamples function. 

// 	This test case checks the computed samples vs. samples computed by hand with Matlab to see that they are within the 
// 	track tolerance during the endgame. 

// 	*/
// 	mpfr_float::default_precision(50);



// 	bertini::System sys;
// 	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars);
// 	sys.AddPathVariable(t);
// 	// Define homotopy system
// 	sys.AddFunction( pow(x - mpfr("1"),3)*(1-t) + (pow(x,3) + mpfr("1"))*t);


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

// 	mpfr sample_factor = mpfr(".5");
// 	mpfr origin = mpfr("0");
// 	Vec<mpfr> x_origin(1);
// 	x_origin << mpfr("1","0");
// 	std::deque<mpfr> correct_times; //times are not vectors they are just complex numbers, these are the correct times from bertini
// 	std::deque< Vec<mpfr> > correct_samples; //samples are space values that may be a vector of complex numbers, correct space values 
// 											// from bertini.

// 	std::deque<mpfr> times; //times are not vectors they are just complex numbers.

// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	// Vec<dbl> derivative(1);

// 	time = mpfr(".1"); // x = .1
// 	correct_times.push_back(time);
// 	sample << mpfr("5.000000000000001e-01", "9.084258952712920e-17"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic
// 	correct_samples.push_back(sample);

// 	time = mpfr(".05"); // x = .1/2 = .05
// 	correct_times.push_back(time);
// 	sample << mpfr("6.000000000000000e-01", "8.165397611531455e-19"); //f(.05) = 6.000000000000000e-01 8.165397611531455e-19 from bertini classic.
// 	correct_samples.push_back(sample);

// 	time = mpfr(".025"); // x = .05/2 = .025
// 	correct_times.push_back(time);
// 	sample << mpfr("6.772905941598711e-01", "3.869924129415447e-17"); // f(.025) = 6.772905941598711e-01 3.869924129415447e-17 from bertini classic
// 	correct_samples.push_back(sample);

// 	unsigned int num_samples = 3; 

// 	mpfr current_time(1);
// 	Vec<mpfr> current_space(1);
// 	current_time = mpfr(".1");
// 	current_space << mpfr("5.000000000000001e-01", "9.084258952712920e-17");

// 	bertini::tracking::config::EndGame endgame_struct;
// 	bertini::tracking::config::PowerSeries power_series_struct;

// 	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_struct,power_series_struct);

// 	My_Endgame.ComputeInitialSamples(current_time,current_space);

// 	auto samples = My_Endgame.GetSamples();

// 	for(unsigned ii = 0; ii < samples.size(); ++ii)
// 	{
// 		BOOST_CHECK((samples[ii] - correct_samples[ii]).norm() < My_Endgame.GetTrackToleranceDuringEndgame());
// 	}

// }//end compute initial samples mp

// BOOST_AUTO_TEST_CASE(compute_initial_samples_dbl_for_powerseries_class)
// {
// 	/*
// 	Given a time value (usually 0.1) and a vector representing values for all other variables the first step in the power series 
// 	endgame is to get some initial samples to make the first hermite interpolation. 

// 	For this reason there exists a ComputeInitialSamples function. 

// 	This test case checks the computed samples vs. samples computed by hand with Matlab to see that they are within the 
// 	track tolerance during the endgame. 

// 	*/
// 	mpfr_float::default_precision(16);



// 	bertini::System sys;
// 	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t");
// 	VariableGroup vars{x};
// 	sys.AddVariableGroup(vars);
// 	sys.AddPathVariable(t);
// 	// Define homotopy system
// 	sys.AddFunction( pow(x - mpfr("1"),3)*(1-t) + (pow(x,3) + mpfr("1"))*t);


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

// 	mpfr sample_factor = mpfr(".5");
// 	mpfr origin = mpfr("0");
// 	Vec<mpfr> x_origin(1);
// 	x_origin << mpfr("1","0");
// 	std::deque<mpfr> correct_times; //times are not vectors they are just complex numbers, these are the correct times from bertini
// 	std::deque< Vec<mpfr> > correct_samples; //samples are space values that may be a vector of complex numbers, correct space values 
// 											// from bertini.

// 	std::deque<mpfr> times; //times are not vectors they are just complex numbers.

// 	mpfr time(1);
// 	Vec<mpfr> sample(1);
// 	// Vec<dbl> derivative(1);

// 	time = mpfr(".1"); // x = .1
// 	correct_times.push_back(time);
// 	sample << mpfr("5.000000000000001e-01", "9.084258952712920e-17"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic
// 	correct_samples.push_back(sample);

// 	time = mpfr(".05"); // x = .1/2 = .05
// 	correct_times.push_back(time);
// 	sample << mpfr("6.000000000000000e-01", "8.165397611531455e-19"); //f(.05) = 6.000000000000000e-01 8.165397611531455e-19 from bertini classic.
// 	correct_samples.push_back(sample);

// 	time = mpfr(".025"); // x = .05/2 = .025
// 	correct_times.push_back(time);
// 	sample << mpfr("6.772905941598711e-01", "3.869924129415447e-17"); // f(.025) = 6.772905941598711e-01 3.869924129415447e-17 from bertini classic
// 	correct_samples.push_back(sample);


// 	mpfr current_time(1);
// 	Vec<mpfr> current_space(1);
// 	current_time = mpfr(".1");
// 	current_space << mpfr("5.000000000000001e-01", "9.084258952712920e-17");

// 	bertini::tracking::config::PowerSeries power_series_struct;
// 	bertini::tracking::config::Security endgame_security_struct;
// 	bertini::tracking::config::Tolerances endgame_tolernaces_struct;

// 	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,power_series_struct,endgame_security_struct,endgame_tolernaces_struct);


// 	My_Endgame.ComputeInitialSamples(current_time,current_space);

// 	auto samples = My_Endgame.GetSamples();

// 	for(unsigned ii = 0; ii < samples.size(); ++ii)
// 	{
// 		BOOST_CHECK((samples[ii] - correct_samples[ii]).norm() < My_Endgame.GetTrackToleranceDuringEndgame());
// 	}

// }//end compute initial samples dbl

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

	sys.AddFunction( pow(x - mpfr("1"),3)*(1-t) + (pow(x,3) + mpfr("1"))*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);
	
	bertini::tracking::config::Stepping stepping_preferences;
	bertini::tracking::config::Newton newton_preferences;

	tracker.Setup(bertini::tracking::config::Predictor::Euler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);


	mpfr current_time(1);
	Vec<mpfr> current_space(1);
	current_time = mpfr(".1");
	current_space << mpfr("5.000000000000001e-01", "9.084258952712920e-17");

	Vec<mpfr> correct(1);
	correct << mpfr("1","0");

	bertini::tracking::config::EndGame endgame_struct;
	bertini::tracking::config::Security endgame_security_struct;
	bertini::tracking::config::Tolerances endgame_tolernaces_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_struct,endgame_security_struct,endgame_tolernaces_struct);
	My_Endgame.PSEG(current_time,current_space);

	// std::cout << "norm is "<< (My_Endgame.GetFinalApproximation() - correct).norm() << '\n';
	// std::cout << "Tolerance is " << My_Endgame.GetTrackToleranceDuringEndgame() << '\n';
	BOOST_CHECK((My_Endgame.GetFinalApproximation() - correct).norm() < 1e-4);//My_Endgame.GetTrackToleranceDuringEndgame());

}//end compute initial samples mp

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

	sys.AddFunction( pow(x - mpfr("1"),3)*(1-t) + (pow(x,3) + mpfr("1"))*t);

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);
	
	bertini::tracking::config::Stepping stepping_preferences;
	bertini::tracking::config::Newton newton_preferences;

	tracker.Setup(bertini::tracking::config::Predictor::Euler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);

	mpfr current_time(1);
	Vec<mpfr> current_space(1);
	current_time = mpfr(".1");
	current_space << mpfr("5.000000000000001e-01", "9.084258952712920e-17");

	Vec<mpfr> correct(1);
	correct << mpfr("1","0");

	bertini::tracking::config::EndGame endgame_struct;
	bertini::tracking::config::PowerSeries power_series_struct;
	bertini::tracking::config::Tolerances endgame_tolernaces_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_struct,power_series_struct,endgame_tolernaces_struct);


	My_Endgame.PSEG(current_time,current_space);
	// std::cout << "norm is "<< (My_Endgame.GetFinalApproximation() - correct).norm() << '\n';
	// std::cout << "Tolerance is " << My_Endgame.GetTrackToleranceDuringEndgame() << '\n';

	BOOST_CHECK((My_Endgame.GetFinalApproximation() - correct).norm() < 1e-4);//My_Endgame.GetTrackToleranceDuringEndgame());
	
}//end pseg dbl for power series class


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

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);
	
	bertini::tracking::config::Stepping stepping_preferences;
	bertini::tracking::config::Newton newton_preferences;

	tracker.Setup(bertini::tracking::config::Predictor::Euler,
                mpfr_float("1e-6"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
	tracker.AMPSetup(AMP);


	mpfr current_time(1);
	Vec<mpfr> current_space(2);
	current_time = mpfr(".1");
	current_space <<  mpfr("5.000000000000001e-01", "9.084258952712920e-17") ,mpfr("9.000000000000001e-01","4.358898943540673e-01");

	Vec<mpfr> correct(2);
	correct << mpfr("1","0"),mpfr("1","0");

	bertini::tracking::config::EndGame endgame_struct;
	bertini::tracking::config::PowerSeries power_series_struct;
	bertini::tracking::config::Security endgame_security_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_struct,power_series_struct,endgame_security_struct);

	My_Endgame.PSEG(current_time,current_space);

	BOOST_CHECK((My_Endgame.GetFinalApproximation() - correct).norm() < 1e-4);//My_Endgame.GetTrackToleranceDuringEndgame());

}//end pseg mp test case for power series class

BOOST_AUTO_TEST_CASE(pseg_mp_for_powerseries_class_griewank_osborne)
{
	/*
	The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
	one of the solutions at t = endgame_time. 

	Griewank Osborne is a very classic example. Here we allow x and y to mix. There are six paths to be tracked and we know there values 
	at t = 0.1. 

	Three of these paths will converge to origin and three will diverge to infinity triggering a SecurityMaxNorm issue. 

	This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
	the endgame. 
	*/

	mpfr_float::default_precision(30);

	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t"), y = std::make_shared<Variable>("y");
	VariableGroup vars{x,y};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	sys.AddFunction(((mpfr("29")/mpfr("16"))*pow(x,3)-2*x*y)*(1-t) + (pow(x,3) - 1)*t);
	sys.AddFunction((y - pow(x,2))*(1-t) + (pow(y,2) - 1)*t);

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
	correct << mpfr("0","0"),mpfr("0","0");

	bertini::tracking::config::EndGame endgame_struct;
	bertini::tracking::config::PowerSeries power_series_struct;
	bertini::tracking::config::Security endgame_security_struct;
	bertini::tracking::config::Tolerances endgame_tolerances_struct;

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,endgame_struct,power_series_struct,endgame_security_struct,endgame_tolerances_struct);

	unsigned num_paths_diverging = 0;
	unsigned num_paths_converging = 0;
	for(auto s : current_space_values)
	{
		bertini::tracking::SuccessCode endgame_success = My_Endgame.PSEG(endgame_time,s);
		if(endgame_success == bertini::tracking::SuccessCode::Success){
			BOOST_CHECK((My_Endgame.GetFinalApproximation() - correct).norm() < 1e-4);// My_Endgame.GetTrackToleranceDuringEndgame());
			num_paths_converging++;
		}
		if(endgame_success == bertini::tracking::SuccessCode::SecurityMaxNormReached){
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

	Take note that in this example we have two successes while the other hit a min track time. This is accounted for. 
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

	auto AMP = bertini::tracking::config::AMPConfigFrom(final_system);

	// std::cout << "AMP is " << AMP << '\n';
	// std::cout << "final system is " << final_system << '\n';
	// std::cout << "coefficient bound is " << final_system.CoefficientBound() << '\n';



	auto tracker = AMPTracker(final_system);
	config::Stepping stepping_preferences;
	config::Newton newton_preferences;
	tracker.Setup(config::Predictor::Euler,
	              	mpfr_float("1e-7"), mpfr_float("1e5"),
					stepping_preferences, newton_preferences);

	tracker.AMPSetup(AMP);
	

	mpfr t_start(1), t_endgame_boundary(0.1);
	std::vector<Vec<mpfr> > solutions;
	std::vector<Vec<mpfr> > homogenized_solutions;
	for (unsigned ii = 0; ii < TD.NumStartPoints(); ++ii)
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

	Vec<mpfr> solution_1(2);
	solution_1 << mpfr("0.687592791426802019127961784761","0.0567041721787413764699348206477"), mpfr("1.11076734170908975052327605226","0.207914822575706710605647487");

	Vec<mpfr> solution_2(2);
	solution_2 << mpfr("1.02264960658155701356264444257","0.520917033127216794197167359926"), mpfr("1.11076734170909913190783413484","0.207914822575691493611316218448");

	Vec<mpfr> solution_3(2);
	solution_3 << mpfr("0.989757601991555768794484038153","-0.57762120530600610801563732366"), mpfr("1.11076734170909918741898536609","0.20791482257569138952790765984");

	Vec<mpfr> solution_4(2);
	solution_4 << mpfr("0.687592791426887395278555459299","0.0567041721787893780032385748768"), mpfr("0.689232658290901023523389312686", "-0.207914822575691576878043065335");

	Vec<mpfr> solution_5(2);
	solution_5 << mpfr("1.02264960658081959931948637747","0.520917033127118696605766956908"), mpfr("0.689232658292013194430512839668","-0.207914822576518617066069857076");

	Vec<mpfr> solution_6(2);
	solution_6 << mpfr("0.989757601991599268196007422024","-0.577621205306094375358982164819"), mpfr("0.689232658290901310861758919599","-0.207914822575714712814276165999");


	unsigned num_occurences(0);
	for (auto s : solutions)
	{
		if ( (s-solution_1).norm() < mpfr_float("1e-5"))
			num_occurences++;
	}
	BOOST_CHECK_EQUAL(num_occurences,1);

	num_occurences = 0;
	for (auto s : solutions)
	{
		if ( (s-solution_2).norm() < mpfr_float("1e-5"))
			num_occurences++;
	}
	BOOST_CHECK_EQUAL(num_occurences,1);

	num_occurences = 0;
	for (auto s : solutions)
	{
		if ( (s-solution_3).norm() < mpfr_float("1e-5"))
			num_occurences++;
	}
	BOOST_CHECK_EQUAL(num_occurences,1);

	num_occurences = 0;
	for (auto s : solutions)
	{
		if ( (s-solution_4).norm() < mpfr_float("1e-5"))
			num_occurences++;
	}
	BOOST_CHECK_EQUAL(num_occurences,1);

	num_occurences = 0;
	for (auto s : solutions)
	{
		if ( (s-solution_5).norm() < mpfr_float("1e-5"))
			num_occurences++;
	}
	BOOST_CHECK_EQUAL(num_occurences,1);

	num_occurences = 0;
	for (auto s : solutions)
	{
		if ( (s-solution_6).norm() < mpfr_float("1e-5"))
			num_occurences++;
	}
	BOOST_CHECK_EQUAL(num_occurences,1);

	Vec<mpfr> correct(2);
	correct << mpfr("1","0"),mpfr("1","0");

	bertini::tracking::endgame::PowerSeriesEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker);


	std::vector<Vec<mpfr> > endgame_solutions;

	unsigned num_successful_occurences = 0;
	unsigned num_min_track_time_reached = 0;
	for (auto s : homogenized_solutions)
	{
		bertini::tracking::SuccessCode endgame_success = My_Endgame.PSEG(t_endgame_boundary,s);
		if(endgame_success == bertini::tracking::SuccessCode::Success)
		{
			// std::cout << "difference is " << (tracker.GetSystem().DehomogenizePoint(My_Endgame.GetFinalApproximation())-correct).norm() << '\n';
			// if((tracker.GetSystem().DehomogenizePoint(My_Endgame.GetFinalApproximation())-correct).norm() < My_Endgame.GetTrackToleranceDuringEndgame())
			// {
				num_successful_occurences++;
			// }
		}
		if(endgame_success == bertini::tracking::SuccessCode::MinStepSizeReached)
		{
			num_min_track_time_reached++;
		}
		// std::cout << "solution is " << My_Endgame.GetFinalApproximation() << '\n';
	}

	// std::cout << " num_occurences is " << num_occurences << '\n';
 	BOOST_CHECK_EQUAL(num_successful_occurences,6);
 	BOOST_CHECK_EQUAL(num_min_track_time_reached,0);

}

BOOST_AUTO_TEST_SUITE_END()
