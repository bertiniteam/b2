//This file is part of Bertini 2.
//
//generic_pseg_test.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//generic_pseg_test.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with generic_pseg_test.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University

/**
\file generic_pseg_test.hpp Defines tests that all PSEG's, combined with all tracker types,  must pass. 

This file in intended for inclusion into another file, which declares TrackerType and TestedEGType
*/


// there is deliberately a missing #pragma once here, because want to allow inclusion in multiple suites in same file.


using System = bertini::System;
using Variable = bertini::node::Variable;


using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;
using bertini::MakeVariable;

using SuccessCode = bertini::SuccessCode;


using dbl = bertini::dbl;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;
using mpq_rational = bertini::mpq_rational;

using bertini::Precision;
using bertini::DefaultPrecision;

template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;


using PrecisionConfig = bertini::tracking::TrackerTraits<TrackerType>::PrecisionConfig;

using BRT = bertini::tracking::TrackerTraits<TrackerType>::BaseRealType;
using BCT = bertini::tracking::TrackerTraits<TrackerType>::BaseComplexType;

template<typename ...T>
BCT ComplexFromString(T... s)
{return bertini::NumTraits<BCT>::FromString(s...);}

template<typename ...T>
BRT RealFromString(T... s)
{return bertini::NumTraits<BRT>::FromString(s...);}

/**
This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 
The test case will construct 3 samples with derivative,time, and space values for the function x^8 + 1.
After the three samples have been constructed there will be a hermite interpolation. 
We check this against the tracking tolerance for the endgame. 
*/
BOOST_AUTO_TEST_CASE( basic_hermite_test_case_against_matlab )
{
	DefaultPrecision(ambient_precision);



	BCT target_time(0,0); //our target time is the origin.
	unsigned int num_samples = 3;

	bertini::TimeCont<BCT> times; 
	bertini::SampCont<BCT> samples, derivatives;

	BCT time;
	Vec<BCT> sample(1), derivative(1);

	time = ComplexFromString(".1"); // x = .1
	times.push_back(time);
	sample << ComplexFromString("1.00000001"); // f(.1) = 1.00000001
	samples.push_back(sample);
	derivative << ComplexFromString("8e-7"); //f'(.1) = 8e-7
	derivatives.push_back(derivative);

	time = ComplexFromString(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << ComplexFromString("1.0000000000390625"); //f(.05) = 1.0000000000390625
	samples.push_back(sample);
	derivative << ComplexFromString("6.25e-9"); //f'(.05) = 6.25e-9
	derivatives.push_back(derivative);

	time = ComplexFromString(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << ComplexFromString("1.000000000000152587890625"); // f(.025) = 1.000000000000152587890625
	samples.push_back(sample);
	derivative << ComplexFromString("4.8828125e-11"); //f'(.025) = 4.8828125e-11
	derivatives.push_back(derivative);

	Vec< BCT > first_approx = HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);


	BOOST_CHECK( norm(first_approx(0) - ComplexFromString("0.9999999767578209232082898114211261253459","0")) < 1e-7); 
	// answer was found using matlab for a check. difference is diff is 2.32422e-08

}//end basic hermite test case mp against matlab







/**

This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 

The test case will construct 3 samples with derivative,time, and space values for the function x^8 + 1.

After the three samples have been constructed there will be a hermite interpolation. 

Next, a new sample is constructed and the earliest sample is discarded. Leaving us three samples that are "nearer"
to the target at the origin. 

A new approximation is made, with the a new sample and approximation done afterwards. 

We then check to make sure our approximations are getting better by checking the distance from the correct answer and the 
various approximations made. 

*/
BOOST_AUTO_TEST_CASE(hermite_interpolation)
{
	DefaultPrecision(ambient_precision);

	BCT target_time(0,0);
	unsigned int num_samples = 3;

	bertini::TimeCont<BCT> times; 
	bertini::SampCont<BCT> samples, derivatives;

	BCT time;
	Vec<BCT> sample(1), derivative(1);

	time = ComplexFromString(".1"); // x = .1
	times.push_back(time);
	sample << ComplexFromString("1.00000001"); // f(.1) = 1.00000001
	samples.push_back(sample);
	derivative << ComplexFromString("8e-7"); //f'(.1) = 8e-7
	derivatives.push_back(derivative);

	time = ComplexFromString(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << ComplexFromString("1.0000000000390625"); //f(.05) = 1.0000000000390625
	samples.push_back(sample);
	derivative << ComplexFromString("6.25e-9"); //f'(.05) = 6.25e-9
	derivatives.push_back(derivative);

	time = ComplexFromString(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << ComplexFromString("1.000000000000152587890625"); // f(.025) = 1.000000000000152587890625
	samples.push_back(sample);
	derivative << ComplexFromString("4.8828125e-11"); //f'(.025) = 4.8828125e-11
	derivatives.push_back(derivative);


	 Vec< BCT > first_approx = HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);
	 Vec< BCT > correct(1);
	 correct << ComplexFromString("1");


	//Setting up a new sample for approximation.
	time = ComplexFromString(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << ComplexFromString("1.00000000000000059604644775390625"); // f(.0125) = 1.00000000000000059604644775390625
	samples.push_back(sample);
	derivative << ComplexFromString("3.814697265625e-13"); //f'(.0125) = 3.814697265625e-13
	derivatives.push_back(derivative);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();
	derivatives.pop_front();

	//Compute the second approximation.
	Vec< BCT > second_approx = HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);


	// //Check to make sure we are doing better. 
	BOOST_CHECK(abs(second_approx(0)-correct(0)) < abs(first_approx(0)-correct(0)));

	//Setting up new sample for use in approximation.
	time = ComplexFromString("0.00625"); //.0125/2 = 0.00625
	times.push_back(time);
	sample << ComplexFromString("1.0000000000000000023283064365386962890625"); // f(.00625) = 1.0000000000000000023283064365386962890625
	samples.push_back(sample);
	derivative << ComplexFromString("2.98023223876953125000000000000000e-15"); //f'(.00625) = 2.98023223876953125000000000000000×e-15
	derivatives.push_back(derivative);

	times.pop_front();
	samples.pop_front();
	derivatives.pop_front();


	Vec< BCT > third_approx = HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);

	BOOST_CHECK((first_approx - correct).norm() < 1e-10);
	BOOST_CHECK((second_approx - correct).norm() < 1e-10);	
	BOOST_CHECK((third_approx - correct).norm() < 1e-10);

}//end hermite test case






/** In the power series endgame there is a bound on that is calculated that will be used as a higher
bound for the exhaustive search of the best cycle number. 

This is calcuated using a modified version of the cycle test in the bertini book, page 53. 

The bound calculated is then compared against the user defined setting MaxCycleNum which is default at 6. 

We take a function (x-1)^3 because we know it should have a cycle number of 3. (This is true but will not be true near 0 or near 0.1).

This is because the path will not look cubic globally. 
*/
BOOST_AUTO_TEST_CASE(compute_bound_on_cycle_num)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x");
	sys.AddFunction(pow(x-1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> times; 
	bertini::SampCont<BCT> samples, derivatives;

	BCT time;
	Vec<BCT> sample(1), derivative(1);

	time = ComplexFromString(".1"); // x = .1
	times.push_back(time);
	sample << ComplexFromString("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);


	time = ComplexFromString(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << ComplexFromString("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);


	time = ComplexFromString(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << ComplexFromString("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);

	bertini::endgame::EndgameConfig endgame_settings;

	TestedEGType my_endgame(tracker,endgame_settings);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.SetRandVec<BCT>(samples.back().size());
	my_endgame.ComputeBoundOnCycleNumber<BCT>();


	BOOST_CHECK(my_endgame.UpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

	auto first_upper_bound = my_endgame.UpperBoundOnCycleNumber();

	//Setting up a new sample for approximation.
	time = ComplexFromString(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << ComplexFromString("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);


	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.ComputeBoundOnCycleNumber<BCT>();


	BOOST_CHECK(my_endgame.UpperBoundOnCycleNumber() == 6); // max_cycle_num implemented max(5,6) = 6

} // end compute bound on cycle number 








/**
Once we have an upper bound for the exhaustive search on the correct cycle number, see above test cases. 

The next step is to take each possible cycle number, take dx/dt -> dx/ds and t -> s, where t = s^c where c is the proposed 
cycle number. 

With each conversion we do a Hermite interpolation with the first n-1 samples we have to predict the nth sample.

The best approximation determines which proposed cycle number we will use to run a Hermite interpolation with the n samples to 
approximate at the origin.
*/	
BOOST_AUTO_TEST_CASE(compute_cycle_number)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t");
	sys.AddFunction( pow(x-1,3) );  //f(x) = (x-1)^3 

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> times; 
	bertini::SampCont<BCT> samples; 


	BCT time(1);
	Vec<BCT> sample(1);


	time = ComplexFromString(".1"); // x = .1
	times.push_back(time);
	sample << ComplexFromString("-0.729"); // f(.1) = -0.729
	samples.push_back(sample);


	time = ComplexFromString(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << ComplexFromString("-0.857375"); //f(.05) = -0.857375
	samples.push_back(sample);


	time = ComplexFromString(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << ComplexFromString("-0.926859375"); // f(.025) = -0.926859375
	samples.push_back(sample);



	//Setting up a new sample for approximation.
	time = ComplexFromString(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << ComplexFromString("-0.962966796875"); // f(.0125) = -0.962966796875
	samples.push_back(sample);


	bertini::endgame::PowerSeriesConfig power_series_settings;

	TestedEGType my_endgame(tracker,power_series_settings);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.ComputeAllDerivatives<BCT>();

	my_endgame.SetRandVec<BCT>(samples.back().size());
	my_endgame.ComputeCycleNumber<BCT>(BCT(0));

	BOOST_CHECK(my_endgame.CycleNumber() == 1);

} // end compute cycle number 




/**
Compute approximation at origin using three sample points. 

Then, find a new sample and remove earliest known sample. 

Compute a new approximation with the current three samples. 

Do this till we have three approximations. Check to see if the third approximation is better than second. 
*/
BOOST_AUTO_TEST_CASE(compute_approximation_of_x_at_t0)
{
	DefaultPrecision(ambient_precision);


	bertini::System sys;
	Var x = MakeVariable("x"), t = MakeVariable("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3)+1)*t);

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);

	auto origin = BCT(0,0);
	Vec<BCT> x_origin(1);
	x_origin << BCT(1);
	bertini::TimeCont<BCT> times; 
	bertini::SampCont<BCT> samples; 



	BCT time;
	Vec<BCT> sample(1);

	time = ComplexFromString(".1"); // x = .1
	times.push_back(time);
	sample << ComplexFromString("0.50000000000000007812610562824908678293817200689660e0", "0.90818521453245102009102405377246269374768541575725e-16"); // f(.1) = 0.50000000000000007812610562824908678293817200689660e0 0.90818521453245102009102405377246269374768541575725e-16 from bertini classic
	samples.push_back(sample);

	time = ComplexFromString(".05"); // x = .1/2 = .05
	times.push_back(time);
	sample << ComplexFromString("0.60000000000000000073301140774132693211475245208482e0", "0.85317043225116681251211164764202768939258706233715e-18"); //f(.05) = 0.60000000000000000073301140774132693211475245208482e0 0.85317043225116681251211164764202768939258706233715e-18 from bertini classic.
	samples.push_back(sample);

	time = ComplexFromString(".025"); // x = .05/2 = .025
	times.push_back(time);
	sample << ComplexFromString("0.67729059415987117534436422700325955211292174181447e0", "0.38712848412230230944976856052769427012481164605204e-16"); // f(.025) = 0.67729059415987117534436422700325955211292174181447e0 0.38712848412230230944976856052769427012481164605204e-16 from bertini classic
	samples.push_back(sample);

	bertini::endgame::EndgameConfig endgame_settings;

	TestedEGType my_endgame(tracker,endgame_settings);
	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	Vec<BCT> first_approx;
	Vec<BCT> approx_1(1);
	my_endgame.SetRandVec<BCT>(samples.back().size());

	my_endgame.ComputeAllDerivatives<BCT>();

	auto code = my_endgame.ComputeApproximationOfXAtT0(first_approx, origin);
	approx_1 << ComplexFromString("1.04025","6.86403e-14");
	BOOST_CHECK(code==SuccessCode::Success);
	//Setting up a new sample for approximation.
	time = ComplexFromString(".0125"); //.025/2 = .0125
	times.push_back(time);
	sample << ComplexFromString("0.73905643972939615063207159129047946977060596291527e0", "0.44983246338361191567539019211879692583079653921201e-18"); // f(.0125) = 0.73905643972939615063207159129047946977060596291527e0 0.44983246338361191567539019211879692583079653921201e-18
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	Vec<BCT> second_approx;
	Vec<BCT> approx_2(1);

	my_endgame.ComputeAllDerivatives<BCT>();

	code = my_endgame.ComputeApproximationOfXAtT0(second_approx,origin);
	approx_2 << ComplexFromString("0.69995","2.05044e-16");
	BOOST_CHECK(code==SuccessCode::Success);
	//Setting up a new sample for approximation.
	time = ComplexFromString(".00625"); //.0125/2 = .00625
	times.push_back(time);
	sample << ComplexFromString("0.78910678115791459147153183840413839746566925123387e0", "0.22640341504967423865456128414532605156908222616485e-16"); // f(.00625) = 0.78910678115791459147153183840413839746566925123387e0 0.22640341504967423865456128414532605156908222616485e-16
	samples.push_back(sample);

	//Get rid of earliest sample. 
	times.pop_front();
	samples.pop_front();

	my_endgame.SetTimes(times);
	my_endgame.SetSamples(samples);

	my_endgame.ComputeAllDerivatives<BCT>();
	
	Vec<BCT> third_approx;
	Vec<BCT> approx_3(1);
	code = my_endgame.ComputeApproximationOfXAtT0(third_approx, origin);
	approx_3 << ComplexFromString("0.683002","-4.87707e-17");
	BOOST_CHECK(code==SuccessCode::Success);
} // end compute approximation of x at t0







/**
Given a time value (usually 0.1) and a vector representing values for all other variables the first step in the power series 
endgame is to get some initial samples to make the first hermite interpolation. 

For this reason there exists a ComputeInitialSamples function. 

This test case checks the computed samples vs. samples computed by hand with Matlab to see that they are within the 
track tolerance during the endgame. 
*/
BOOST_AUTO_TEST_CASE(compute_initial_samples)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x"), t = MakeVariable("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars); sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3) + 1)*t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);

	BCT sample_factor = ComplexFromString(".5");
	BCT origin = BCT(0);
	Vec<BCT> x_origin(1);
	x_origin << BCT(1);
	bertini::TimeCont<BCT> correct_times; 
	bertini::SampCont<BCT> correct_samples;

	bertini::TimeCont<BCT> times; 
	std::deque<Vec<BCT> > samples;
	BCT time(1);
	Vec<BCT> sample(1);


	time = ComplexFromString(".1"); // x = .1
	correct_times.push_back(time);
	sample << ComplexFromString("5.000000000000001e-01", "9.084258952712920e-17"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic
	correct_samples.push_back(sample);

	time = ComplexFromString(".05"); // x = .1/2 = .05
	correct_times.push_back(time);
	sample << ComplexFromString("6.000000000000000e-01", "8.165397611531455e-19"); //f(.05) = 6.000000000000000e-01 8.165397611531455e-19 from bertini classic.
	correct_samples.push_back(sample);

	time = ComplexFromString(".025"); // x = .05/2 = .025
	correct_times.push_back(time);
	sample << ComplexFromString("6.772905941598711e-01", "3.869924129415447e-17"); // f(.025) = 6.772905941598711e-01 3.869924129415447e-17 from bertini classic
	correct_samples.push_back(sample);

	BCT current_time(1);
	Vec<BCT> current_space(1);
	current_time = ComplexFromString(".1");
	current_space << ComplexFromString("5.000000000000001e-01", "9.084258952712920e-17");

	bertini::endgame::EndgameConfig endgame_settings;
	bertini::endgame::PowerSeriesConfig power_series_settings;

	TestedEGType my_endgame(tracker,endgame_settings,power_series_settings);

	auto tracking_success = my_endgame.ComputeInitialSamples(current_time, origin, current_space, times, samples);

	BOOST_REQUIRE(tracking_success==SuccessCode::Success);
	

	for(unsigned ii = 0; ii < samples.size(); ++ii)
	{
		BOOST_CHECK_EQUAL(samples[ii].size(),1);
		BOOST_CHECK((samples[ii] - correct_samples[ii]).norm() < 1e-5);
	}

}//end compute initial samples


/**
This test will check to see if we can compute initial samples where the time we start and target_time are not 0.1 and 0 respectively.
target_time is .1 + .1I and we are going to start at the time value .2

*/
BOOST_AUTO_TEST_CASE(compute_initial_samples_non_zero_target_time)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x"), t = MakeVariable("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars); sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3) + 1)*t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);

	BCT target_time = ComplexFromString(".1",".1");
	Vec<BCT> x_origin(1);
	x_origin << BCT(1);
	bertini::TimeCont<BCT> correct_times; 
	bertini::SampCont<BCT> correct_samples;

	bertini::TimeCont<BCT> times; 
	std::deque<Vec<BCT> > samples;
	BCT time(1);
	Vec<BCT> sample(1);


	time = ComplexFromString(".2"); // x = .2
	correct_times.push_back(time);
	sample << ComplexFromString("3.603621541081173e-01", "2.859583229930518e-18"); // f(.2) = 3.603621541081173e-01 2.859583229930518e-18 from bertini classic
	correct_samples.push_back(sample);

	time = ComplexFromString(".15",".05"); // x = ((.2) + (.1 + .1I)) / 2 = .15 + .05I
	correct_times.push_back(time);
	sample << ComplexFromString("4.205197361710131e-01","-6.739509600163453e-02"); //f(.15 + .05I) = 4.205197361710131e-01 -6.739509600163453e-02  from bertini classic.
	correct_samples.push_back(sample);

	time = ComplexFromString(".125",".075"); // x = ((.15 + .05I) + (.1 + .1I))/2 = .125 + .075I
	correct_times.push_back(time);
	sample << ComplexFromString("4.469031847163714e-01", "-1.059855741137409e-01"); // f(.125 + .075I) =  4.469031847163714e-01 -1.059855741137409e-01 from bertini classic
	correct_samples.push_back(sample);

	BCT current_time(1);
	Vec<BCT> current_space(1);
	current_time = ComplexFromString(".2");
	current_space << ComplexFromString("3.603621541081173e-01", "2.859583229930518e-18");

	bertini::endgame::EndgameConfig endgame_settings;
	bertini::endgame::PowerSeriesConfig power_series_settings;

	TestedEGType my_endgame(tracker,endgame_settings,power_series_settings);

	auto tracking_success = my_endgame.ComputeInitialSamples(current_time, target_time, current_space, times, samples);

	BOOST_REQUIRE(tracking_success==SuccessCode::Success);
	

	for(unsigned ii = 0; ii < samples.size(); ++ii)
	{
		BOOST_CHECK_EQUAL(samples[ii].size(),1);
		BOOST_CHECK((samples[ii] - correct_samples[ii]).norm() < 1e-5);
	}

}//end compute initial samples nonzero target time







/**
The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
one of the solutions at t = endgame_time. 

This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
the endgame. 
*/
BOOST_AUTO_TEST_CASE(pseg_full_run)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x"), t = MakeVariable("t");

	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3)+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-6,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);


	BCT current_time(1);
	Vec<BCT> current_space(1);
	current_time = ComplexFromString(".1");
	current_space << ComplexFromString("5.000000000000001e-01", "9.084258952712920e-17");

	Vec<BCT> correct(1);
	correct << BCT(1);

	bertini::endgame::EndgameConfig endgame_settings;
	bertini::endgame::SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,endgame_settings,security_settings);
	my_endgame.Run(current_time,current_space);


	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - correct).norm() < 1e-11);

}//end pseg mp for power series class 


/**
The function that runs the power series endgame is called Run. Run takes an endgame_time value and and endgame_space value that is
one of the solutions at t = endgame_time, it can also set a target_time other than t = 0. 

This test will start at t = 0.2 and use the endgame to find the solution at t = .1 + .1*I. This is checking to make sure the powerseries
endgame is general enough to move from t = a to t = b for a,b being generic complex numbers. 
*/
BOOST_AUTO_TEST_CASE(pseg_full_run_non_zero_target_time)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x"), t = MakeVariable("t");

	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3)+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-6,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);


	BCT current_time(1);
	BCT target_time(1);
	Vec<BCT> current_space(1);
	current_time = ComplexFromString(".2");
	target_time = ComplexFromString(".1", ".1");

	current_space << ComplexFromString("3.603621541081173e-01", "2.859583229930518e-18");

	Vec<BCT> correct(1);
	correct << ComplexFromString("4.680740395503238e-01", "-1.470429372721208e-01");

	bertini::endgame::EndgameConfig endgame_settings;
	bertini::endgame::SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,endgame_settings,security_settings);
	my_endgame.Run(current_time,current_space,target_time);


	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - correct).norm() < 1e-11);

}//end pseg mp for power series class non zero target time







/**
The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
one of the solutions at t = endgame_time. 

This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
the endgame. 
*/
BOOST_AUTO_TEST_CASE(full_run_cycle_num_2)
{
	
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2)-1)*t);


	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-6,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);


	BCT start_time(1);
	Vec<BCT> start_point(1); start_point << BCT(1);
	auto t_endgame_boundary = ComplexFromString("0.1");
	
	Vec<BCT> eg_boundary_point;
	auto init_success = tracker.TrackPath(eg_boundary_point, start_time, t_endgame_boundary, start_point);

	BOOST_CHECK(init_success==SuccessCode::Success);
	BOOST_CHECK(abs(eg_boundary_point(0) - BRT(1))< 1e-5);

	Vec<BCT> correct_root(1);
	correct_root << BCT(1);

	bertini::endgame::EndgameConfig endgame_settings;
	bertini::endgame::SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,endgame_settings,security_settings);
	my_endgame.Run(t_endgame_boundary,eg_boundary_point);

	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(),1);
	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - correct_root).norm() < 1e-11);

}//end pseg for power series class 







/**
The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
one of the solutions at t = endgame_time. 

In this test we do multiple variables decoupled, that has a high multiplicity (5) solution. 
*/
BOOST_AUTO_TEST_CASE(full_run_multiple_variables)
{
	DefaultPrecision(ambient_precision);


	bertini::System sys;
	Var x = MakeVariable("x"), t = MakeVariable("t"), y = MakeVariable("y");
	VariableGroup vars{x,y};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	sys.AddFunction((pow(x-1,3))*(1-t) + (pow(x,3) + 1)*t);
	sys.AddFunction((pow(y-1,2))*(1-t) + (pow(y,2) + 1)*t);

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-6,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);


	BCT current_time(1);
	Vec<BCT> current_space(2);
	current_time = ComplexFromString(".1");
	current_space <<  ComplexFromString("5.000000000000001e-01", "9.084258952712920e-17") ,ComplexFromString("9.000000000000001e-01","4.358898943540673e-01");

	Vec<BCT> correct(2);
	correct << BCT(1),BCT(1);

	bertini::endgame::EndgameConfig endgame_settings;
	bertini::endgame::PowerSeriesConfig power_series_settings;
	bertini::endgame::SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,endgame_settings,power_series_settings,security_settings);

	my_endgame.Run(current_time,current_space);

	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - correct).norm() < 1e-10);//my_endgame.GetTrackToleranceDuringEndgame());

}//end pseg mp test case for power series class




/**
The function that runs the power series endgame is called PSEG. PSEG takes an endgame_time value and and endgame_space value that is
one of the solutions at t = endgame_time. 

Griewank Osborne is a very classic example. Here we allow x and y to mix. There are six paths to be tracked and we know there values 
at t = 0.1. 

Three of these paths will converge to origin and three will diverge to infinity triggering a SecurityMaxNorm issue. 

This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
the endgame. 


has six solutions at t = .1:
0
0.687592791426802019127961784761 0.0567041721787413764699348206477
1.11076734170908975052327605226  0.207914822575706710605647487

1
1.02264960658155701356264444257  0.520917033127216794197167359926
1.11076734170909913190783413484  0.207914822575691493611316218448

2
0.989757601991555768794484038153 -0.57762120530600610801563732366
1.11076734170909918741898536609  0.20791482257569138952790765984

3
0.687592791426887395278555459299 0.0567041721787893780032385748768 
0.689232658290901023523389312686 -0.207914822575691576878043065335

4
1.02264960658081959931948637747  0.520917033127118696605766956908
0.689232658292013194430512839668 -0.207914822576518617066069857076

5
0.989757601991599268196007422024 -0.577621205306094375358982164819
0.689232658290901310861758919599 -0.207914822575714712814276165999
*/
BOOST_AUTO_TEST_CASE(griewank_osborne)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x"), t = MakeVariable("t"), y = MakeVariable("y");
	VariableGroup vars{x,y};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);

	sys.AddFunction((mpq_rational(29,16)*pow(x,3)-2*x*y)*(1-t) + (pow(x,3) - 1)*t);
	sys.AddFunction((y - pow(x,2))*(1-t) + (pow(y,2) - 1)*t);

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;

	tracker.Setup(TestedPredictor,
                1e-6,
                1e5,
                stepping_settings,
                newton_settings);
	
	tracker.PrecisionSetup(precision_config);


	BCT endgame_time(1);
	Vec<BCT> current_space_1(2);
	Vec<BCT> current_space_2(2);
	Vec<BCT> current_space_3(2);
	Vec<BCT> current_space_4(2);
	Vec<BCT> current_space_5(2);
	Vec<BCT> current_space_6(2);

	endgame_time = ComplexFromString(".1");

	current_space_1 <<  ComplexFromString("1.028694756284462e-01", "-9.661822229074768e-01"), ComplexFromString("-8.937287306314232e-01", "-2.480445481975048e-01");
	current_space_2 <<  ComplexFromString("-5.219899550566304e-01", "-2.212788871407134e-17"), ComplexFromString("3.684968541245056e-01", "2.471980953266950e-17");
	current_space_3 <<  ComplexFromString("1.028694756284466e-01", "9.661822229074761e-01"), ComplexFromString("-8.937287306314221e-01", "2.480445481975040e-01");
	current_space_4 <<  ComplexFromString("6.098408897464429e-03", "1.058791184067875e-21"), ComplexFromString("-9.109808533477256e+00", "-2.374402757743255e-17");
	current_space_5 <<  ComplexFromString("1.220071827679809e+00", "7.657177843178875e-19"), ComplexFromString("1.386185299689565e+00", "2.852806966352484e-18");
	current_space_6 <<  ComplexFromString("-9.099192327775354e-01", "2.114194236346734e-17"), ComplexFromString("8.573852849693505e-01", "-2.164338586824188e-17");

	std::vector<Vec<BCT> > current_space_values;

	current_space_values.push_back(current_space_1);
	current_space_values.push_back(current_space_2);
	current_space_values.push_back(current_space_3);
	current_space_values.push_back(current_space_4);
	current_space_values.push_back(current_space_5);
	current_space_values.push_back(current_space_6);


	Vec<BCT> correct(2);
	correct << BCT(0), BCT(0);

	bertini::endgame::EndgameConfig endgame_settings;
	bertini::endgame::PowerSeriesConfig power_series_settings;
	bertini::endgame::SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,endgame_settings,power_series_settings,security_settings);

	unsigned num_paths_diverging = 0;
	unsigned num_paths_converging = 0;
	for (const auto& s : current_space_values)
	{
		DefaultPrecision(ambient_precision);
		SuccessCode endgame_success = my_endgame.Run(endgame_time,s);
		if(endgame_success == SuccessCode::Success){
			BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - correct).norm() < 1e-11);// my_endgame.GetTrackToleranceDuringEndgame());
			num_paths_converging++;
		}
		else if(endgame_success == SuccessCode::SecurityMaxNormReached || endgame_success == SuccessCode::GoingToInfinity){
			num_paths_diverging++;
		}
		else
		{
			num_paths_diverging++;
		}

	}
	BOOST_CHECK_EQUAL(num_paths_converging,3);
	BOOST_CHECK_EQUAL(num_paths_diverging,3);

}//end compute griewank osborne





/**
In this example we take a decoupled system, homogenize and patch it. Track to endgame boundary and then run our endgame on the space
values we have. 
*/
BOOST_AUTO_TEST_CASE(total_degree_start_system)
{
	using namespace bertini::tracking;
	DefaultPrecision(ambient_precision);

	Var x = MakeVariable("x");
	Var y = MakeVariable("y");
	Var t = MakeVariable("t");

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

	//auto gamma = bertini::MakeRational(bertini::node::Rational::Rand());
	//gamma*
	
	auto final_system = (1-t)*sys + t*TD;
	final_system.AddPathVariable(t);

	auto precision_config = PrecisionConfig(final_system);



	auto tracker = TrackerType(final_system);
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;
	tracker.Setup(TestedPredictor,
	              	1e-5, 1e5,
					stepping_settings, newton_settings);

	tracker.PrecisionSetup(precision_config);

#ifdef B2_OBSERVE_TRACKERS
			GoryDetailLogger<TrackerType> tons_of_detail;
			tracker.AddObserver(&tons_of_detail);
#endif

	unsigned num_paths_to_run = 1;
	BCT t_start(1), t_endgame_boundary(0.1);
	std::vector<Vec<BCT> > homogenized_solutions;
	for (unsigned ii = 0; ii < num_paths_to_run; ++ii)
	{
		DefaultPrecision(ambient_precision);
		final_system.precision(ambient_precision);
		auto start_point = TD.StartPoint<BCT>(ii);

		Vec<BCT> result;
		SuccessCode tracking_success;

		tracking_success = tracker.TrackPath(result,t_start,t_endgame_boundary,start_point);
		BOOST_CHECK(tracking_success==SuccessCode::Success);

		homogenized_solutions.push_back(result);
	}

	Vec<BCT> correct(2);
	correct << BCT(1),BCT(1);

	tracker.Setup(TestedPredictor,
	              	1e-6, 1e5,
					stepping_settings, newton_settings);

	TestedEGType my_endgame(tracker);



	std::vector<Vec<BCT> > endgame_solutions;

	unsigned num_successful_occurences = 0;
	for (auto const& s : homogenized_solutions)
	{
		SuccessCode endgame_success = my_endgame.Run(t_endgame_boundary,s);
		if(endgame_success == SuccessCode::Success)
		{
			BOOST_CHECK_EQUAL(Precision(my_endgame.FinalApproximation<BCT>()), tracker.CurrentPrecision());

				num_successful_occurences++;
		}
	}

 	BOOST_CHECK_EQUAL(num_successful_occurences,num_paths_to_run);
}




/**
In this example we take a decoupled system, homogenize and patch it. Track to endgame boundary and then run our endgame on the space
values we have. 

*/
BOOST_AUTO_TEST_CASE(parabola)
{
	using namespace bertini::tracking;
	DefaultPrecision(ambient_precision);

	Var x = MakeVariable("x");
	Var t = MakeVariable("t");

	System sys;

	VariableGroup v{x};

	sys.AddVariableGroup(v);

	sys.AddFunction(pow(x,2) - t);
	sys.AddPathVariable(t);
	Vec<BCT> start_point(1);
	start_point << BCT(1);

	auto precision_config = PrecisionConfig(sys);



	auto tracker = TrackerType(sys);
	bertini::tracking::SteppingConfig stepping_settings;
	bertini::tracking::NewtonConfig newton_settings;
	tracker.Setup(TestedPredictor,
	              	1e-5, 1e5,
					stepping_settings, newton_settings);

	tracker.PrecisionSetup(precision_config);
	

	BCT t_start(1);
	auto t_endgame_boundary = ComplexFromString("0.1");

	Vec<BCT> soln_at_EG_bdry;

	auto tracking_success = tracker.TrackPath(soln_at_EG_bdry,t_start,t_endgame_boundary,start_point);
	BOOST_CHECK(tracking_success==SuccessCode::Success);



	Vec<BCT> correct_eg_soln(1);
	correct_eg_soln << BCT(0);

	tracker.Setup(TestedPredictor,
	              	1e-6, 1e5,
					stepping_settings, newton_settings);

	TestedEGType my_endgame(tracker);


	auto endgame_success = my_endgame.Run(t_endgame_boundary,soln_at_EG_bdry);
	BOOST_CHECK(endgame_success==SuccessCode::Success);

	auto endgame_solution = my_endgame.FinalApproximation<BCT>();

	BOOST_CHECK_EQUAL(Precision(my_endgame.FinalApproximation<BCT>()), tracker.CurrentPrecision());


	BOOST_CHECK_SMALL( abs(endgame_solution(0)-correct_eg_soln(0)), BRT(1e-10) );
}

/**
	Full blown test to see if we can actually track using an endgame to a nonzero target time. 
*/
BOOST_AUTO_TEST_CASE(pseg_full_run_nonzero_target_time)
{

	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction(pow(x-1,3)*(1-t) + (pow(x,3)+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
		
	bertini::tracking::SteppingConfig stepping_preferences;
	bertini::tracking::NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
	    1e-5,
	    1e5,
	    stepping_preferences,
	    newton_preferences);
		
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> pseg_times;
		bertini::SampCont<BCT> pseg_samples;

	auto start_time = ComplexFromString("0.2");
	Vec<BCT> start_sample(1);
	auto target_time = ComplexFromString(".15","-.01");
	Vec<BCT> first_approx(1);
	Vec<BCT> x_to_check_against(1);

	start_sample << ComplexFromString("3.603621541081173e-01", "2.859583229930518e-18"); 
	x_to_check_against << ComplexFromString("4.248924277564006e-01", "1.369835558109531e-02");

	bertini::endgame::EndgameConfig endgame_settings;
	bertini::endgame::SecurityConfig security_settings;
	bertini::endgame::PowerSeriesConfig power_series_settings;
	TestedEGType my_endgame(tracker,endgame_settings,power_series_settings, security_settings);

	auto endgame_success = my_endgame.Run(start_time,start_sample,target_time);
	BOOST_CHECK(endgame_success == SuccessCode::Success);

	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - x_to_check_against).norm() < 1e-10);

}// end cauchy_full_run_nonzero_target_time


