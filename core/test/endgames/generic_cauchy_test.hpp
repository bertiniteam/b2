//This file is part of Bertini 2.
//
//generic_cauchy_test.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//generic_cauchy_test.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with generic_cauchy_test.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
\file generic_cauchy_test.hpp Defines tests that all cauchy endgames's, combined with all tracker types,  must pass. 

This file in intended for inclusion into another test file, which declares TrackerType and TestedEGType
*/



// there is deliberately a missing #pragma once here, because want to allow inclusion in multiple suites in same file.


using System = bertini::System;
using Variable = bertini::node::Variable;


using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;



using dbl = bertini::dbl;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;
using mpq_rational = bertini::mpq_rational;

using bertini::MakeVariable;
using bertini::Precision;

template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;


using PrecisionConfig = bertini::tracking::TrackerTraits<TrackerType>::PrecisionConfig;

using BRT = bertini::tracking::TrackerTraits<TrackerType>::BaseRealType;
using BCT = bertini::tracking::TrackerTraits<TrackerType>::BaseComplexType;

using SuccessCode = bertini::SuccessCode;

template<typename ...T>
BCT ComplexFromString(T... s)
{return bertini::NumTraits<BCT>::FromString(s...);}

template<typename ...T>
BRT RealFromString(T... s)
{return bertini::NumTraits<BRT>::FromString(s...);}


using bertini::DefaultPrecision;



using namespace bertini::tracking;
using namespace bertini::endgame;

/**
	In this test we take a univariate polynomial with one solution and ensure that our CircleTrack function returns to the same position. 
	This is tested by checking that our cauchy_samples has the front and end with the same value roughly.
*/
BOOST_AUTO_TEST_CASE(circle_track_cycle_num_1)
{
	
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	auto origin = BCT(0,0);
	

	TestedEGType my_endgame(tracker);

	bertini::SampCont<BCT>& cauchy_samples = my_endgame.GetCauchySamples<BCT>();
	bertini::TimeCont<BCT>& cauchy_times = my_endgame.GetCauchyTimes<BCT>();



	cauchy_times.push_back(ComplexFromString("0.1"));
	cauchy_samples.push_back(Vec<BCT>(1));
	cauchy_samples.back() << ComplexFromString("7.999999999999999e-01", "2.168404344971009e-19"); // 

	const auto& time = cauchy_times.back();
	const auto& sample = cauchy_samples.back();

	auto first_track_success =  my_endgame.CircleTrack(time,origin,sample);

	BOOST_CHECK((my_endgame.GetCauchySamples<BCT>().back() - sample).template lpNorm<Eigen::Infinity>() < 1e-5);

}





/**
	In this test we take a univariate polynomial with one solution of multiplicity 2. We need to ensure that our CircleTrack function 
	returns to the same position after two calls. The first call should land us on the second path, and the second call should land us back 
	where we started. This is tested by checking that cauchy_samples has the front and end with the same value roughly.
*/
BOOST_AUTO_TEST_CASE(circle_track_cycle_num_greater_than_1)
{
	
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> cauchy_times; 
	bertini::SampCont<BCT> cauchy_samples; 



	BCT time(1);
	Vec<BCT> sample(1);
	auto origin = BCT(0,0);


	time = ComplexFromString("0.1");
	cauchy_times.push_back(time);
	sample << ComplexFromString("9.000000000000001e-01", "4.358898943540673e-01"); // 
	cauchy_samples.push_back(sample);

	TestedEGType my_endgame(tracker);

	my_endgame.SetCauchySamples(cauchy_samples);
	my_endgame.SetCauchyTimes(cauchy_times);
	
	auto tracking_success =  my_endgame.CircleTrack(time,origin,sample);

	const auto& first_track_sample = my_endgame.GetCauchySamples<BCT>().back();

	BOOST_CHECK((first_track_sample - sample).template lpNorm<Eigen::Infinity>() > 1e-5);

	tracking_success =  my_endgame.CircleTrack(time,origin,first_track_sample);

	const auto& second_track_sample = my_endgame.GetCauchySamples<BCT>().back();

	BOOST_CHECK((second_track_sample - sample).template lpNorm<Eigen::Infinity>() < 1e-5);
	
} // end circle_track_mp_cycle_num_greater_than_1


/**
	We are going to track around a nonzero target time. For a random nonzero time value we expect to not encircle any branch points. 
	Making the cycle number 1.
*/
BOOST_AUTO_TEST_CASE(circle_track__nonzero_target_time)
{
	
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> cauchy_times; 
	bertini::SampCont<BCT> cauchy_samples; 



	BCT time(1);
	Vec<BCT> sample(1);
	auto center = ComplexFromString(".19","-.01");


	time = ComplexFromString("0.2");
	cauchy_times.push_back(time);
	sample << ComplexFromString("3.603621541081173e-01", "2.859583229930518e-18"); // 
	cauchy_samples.push_back(sample);

	TestedEGType my_endgame(tracker);

	my_endgame.SetCauchySamples(cauchy_samples);
	my_endgame.SetCauchyTimes(cauchy_times);
	
	auto tracking_success =  my_endgame.CircleTrack(time,center,sample);


	const auto& first_track_sample = my_endgame.GetCauchySamples<BCT>().back();

	BOOST_CHECK((first_track_sample - sample).template lpNorm<Eigen::Infinity>() < 1e-5);

} // end circle_track_nonzero_target_time



/**
	In this test we take a univariate polynomial with one solution of multiplicity. For our first power series approximation we need to 
	make sure that we have a stablization of the cycle number. This is achieved by computing an approximation for C over K which is described 
	on page 53 of the Berini Book, Numerically Solving Polynomials Systems with Bertini \cite bertinibook.

	This test case is ensuring no calculations are greatly altered, so our test condition is to check against a computed value that has been already known. 
*/
BOOST_AUTO_TEST_CASE(compute_c_over_k_for_cauchy_class)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x");
	sys.AddFunction(pow(x-1,3));  //f(x) = (x-1)^3

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> pseg_times;
	bertini::SampCont<BCT> pseg_samples;
	

	BCT time(1);
	Vec<BCT> sample(1);

	time = ComplexFromString(".1"); // x = .1
	pseg_times.push_back(time);
	sample << ComplexFromString("-0.729"); // f(.1) = -0.729
	pseg_samples.push_back(sample);


	time = ComplexFromString(".05"); // x = .1/2 = .05
	pseg_times.push_back(time);
	sample << ComplexFromString("-0.857375"); //f(.05) = -0.857375
	pseg_samples.push_back(sample);


	time = ComplexFromString(".025"); // x = .05/2 = .025
	pseg_times.push_back(time);
	sample << ComplexFromString("-0.926859375"); // f(.025) = -0.926859375
	pseg_samples.push_back(sample);


	SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,security_settings);
	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto first_c_over_k = my_endgame.ComputeCOverK<BCT>();

	using std::abs;
	BOOST_CHECK( abs(first_c_over_k - RealFromString("1.12917")) < 1e-5 ); 

	//Setting up a new sample for approximation.
	time = ComplexFromString(".0125"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << ComplexFromString("-0.962966796875"); // f(.0125) = -0.962966796875
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto second_c_over_k = my_endgame.ComputeCOverK<BCT>();

	BOOST_CHECK( abs(second_c_over_k - RealFromString("1.05888")) <  1e-5); 

} // end compute c over k for cauchy class 







/**
	After we have computed C over K approximations we have to make sure that we lie within some agreed terms. This means our agreements need
	to be with some ratio of each other. (Minimum is usually .75), and we need to agree so many times (usually 3). This test case checks for this. 
	
*/	
BOOST_AUTO_TEST_CASE(stabilization_of_C_over_K)
{
	DefaultPrecision(ambient_precision);

	bertini::System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t");
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> pseg_times;
	bertini::SampCont<BCT> pseg_samples;
	bertini::TimeCont<BCT> c_over_k_array;

	BCT time;
	Vec<BCT> sample(1);

	time = ComplexFromString(".1"); // x = .1
	pseg_times.push_back(time);
	sample << ComplexFromString("1.100000000000000e+00", "6.244997998398396e-01"); // f(.1) = 1.100000000000000e+00 6.244997998398396e-01
	pseg_samples.push_back(sample);


	time = ComplexFromString(".05"); // x = .1/2 = .05
	pseg_times.push_back(time);
	sample << ComplexFromString("1.125000000000000e+00", "4.841229182759271e-01"); //f(.05) = 1.125000000000000e+00 4.841229182759271e-01
	pseg_samples.push_back(sample);


	time = ComplexFromString(".025"); // x = .05/2 = .025
	pseg_times.push_back(time);
	sample << ComplexFromString("1.123854702920065e+00", "3.736283818773165e-01"); // f(.025) = 1.123854702920065e+00 3.736283818773165e-01
	pseg_samples.push_back(sample);


	SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,security_settings);
	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto first_c_over_k = my_endgame.ComputeCOverK<BCT>();

	c_over_k_array.push_back(first_c_over_k);

	//Setting up a new sample for approximation.
	time = ComplexFromString(".0125"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << ComplexFromString("1.111721780135302e+00", "2.886596646224579e-01"); // f(.0125) = 1.111721780135302e+00 2.886596646224579e-01
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto second_c_over_k = my_endgame.ComputeCOverK<BCT>();

	c_over_k_array.push_back(second_c_over_k);



	//Setting up a new sample for approximation.
	time = ComplexFromString(".00625"); //.025/2 = .0125
	pseg_times.push_back(time);
	sample << ComplexFromString("1.096071609421043e+00", "2.237005761359081e-01"); // f(.00625) = 1.096071609421043e+00 2.237005761359081e-01
	pseg_samples.push_back(sample);

	//Get rid of earliest sample. 
	pseg_times.pop_front();
	pseg_samples.pop_front();

	my_endgame.SetPSEGTimes(pseg_times);
	my_endgame.SetPSEGSamples(pseg_samples);

	auto third_c_over_k = my_endgame.ComputeCOverK<BCT>();

	c_over_k_array.push_back(third_c_over_k);

	BOOST_CHECK(my_endgame.CheckForCOverKStabilization(c_over_k_array) == true);

} // end check agreement function mp for cauchy 




/**
	In this test case we have a solution we are trying to find with multiplicity 1. If we call CircleTrack we should have a closed loop. 
	This test case is checking to make sure our CheckClosedLoop function returns that we have actually closed. 
*/
BOOST_AUTO_TEST_CASE(check_closed_loop_for_cycle_num_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> cauchy_times; 
	bertini::SampCont<BCT> cauchy_samples; 


	BCT time(1);
	Vec<BCT> sample(1);
	auto origin = BCT(0,0);


	time = ComplexFromString(".1");
	cauchy_times.push_back(time);
	sample << ComplexFromString("7.999999999999999e-01", "2.168404344971009e-19"); // 
	cauchy_samples.push_back(sample);


	TestedEGType my_endgame(tracker);

	my_endgame.SetCauchySamples(cauchy_samples);
	my_endgame.SetCauchyTimes(cauchy_times);


	auto tracking_success =  my_endgame.CircleTrack(time,origin,sample);
	BOOST_CHECK(my_endgame.CheckClosedLoop<BCT>() == true);

} // end check closed loop if cycle num is 1 for cauchy class test



/**
		In this test case we have a solution we are trying to find with multiplicity 2. If we call CircleTrack we should not have a closed loop. 
		
		This test case is checking to make sure our CheckClosedLoop function returns that we have not closed.
	*/
BOOST_AUTO_TEST_CASE(check_closed_loop_for_cycle_num_greater_than_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);


	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> cauchy_times; 
	bertini::SampCont<BCT> cauchy_samples; 
	auto origin = BCT(0,0);



	auto time = ComplexFromString("0.1");
	cauchy_times.push_back(time);

	Vec<BCT> sample(1);
	sample << ComplexFromString("9.000000000000001e-01", "4.358898943540673e-01"); // 
	cauchy_samples.push_back(sample);


	TestedEGType my_endgame(tracker);

	my_endgame.SetCauchySamples(cauchy_samples);
	my_endgame.SetCauchyTimes(cauchy_times);

	auto tracking_success =  my_endgame.CircleTrack(time,origin,sample);
	BOOST_CHECK(my_endgame.CheckClosedLoop<BCT>() == false);

	tracking_success =  my_endgame.CircleTrack(my_endgame.GetCauchyTimes<BCT>().back(),origin,my_endgame.GetCauchySamples<BCT>().back());
	BOOST_CHECK(my_endgame.CheckClosedLoop<BCT>() == true);
	
} // end check closed loop if cycle num is greater than 1 for cauchy class test





/**
 After we have stabilized our estimation of c over k, we can start to perform cauchy loops to see if the minimum and maximum norms
	of sample points are within a tolerance. If so, we can start actually using the Cauchy Integral Formula to extrapolate approximations 
	at the origin. 

	In this example, since the cycle number is 1, we should see that our RatioEGOperatingZoneTest function says that we are good to proceed 
	in using the Cauchy Integral Formula. 
*/
BOOST_AUTO_TEST_CASE(compare_cauchy_ratios)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);


	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> cauchy_times; 
	bertini::SampCont<BCT> cauchy_samples; 


	auto time = ComplexFromString(".1");
	Vec<BCT> sample(1);
	auto origin = BCT(0,0);

	cauchy_times.push_back(time);
	sample << ComplexFromString("7.999999999999999e-01", "2.168404344971009e-19"); // 
	cauchy_samples.push_back(sample);

	TestedEGType my_endgame(tracker);
	my_endgame.SetCauchyTimes(cauchy_times);
	my_endgame.SetCauchySamples(cauchy_samples);

	auto tracking_success =  my_endgame.CircleTrack(time,origin,sample);
	BOOST_CHECK(my_endgame.RatioEGOperatingZoneTest<BCT>(origin) == true);

} // end compare cauchy ratios for cauchy class test




/**
	After we have stabilized our estimation of c over k, we can start to perform cauchy loops to see if the minimum and maximum norms
	of sample points are within a tolerance. If so, we can start actually using the Cauchy Integral Formula to extrapolate approximations 
	at the origin. 
*/
BOOST_AUTO_TEST_CASE(compare_cauchy_ratios_cycle_num_greater_than_1)
{
	
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> cauchy_times; 
	bertini::SampCont<BCT> cauchy_samples; 
	auto origin = BCT(0,0);



	auto time = ComplexFromString("0.1");
	Vec<BCT> sample(1);

	cauchy_times.push_back(time);
	sample << ComplexFromString("9.000000000000001e-01", "4.358898943540673e-01"); // 
	cauchy_samples.push_back(sample);

	TestedEGType my_endgame(tracker);

	my_endgame.SetCauchyTimes(cauchy_times);
	my_endgame.SetCauchySamples(cauchy_samples);
	
	auto tracking_success =  my_endgame.CircleTrack(time,origin,sample);
	BOOST_CHECK(my_endgame.RatioEGOperatingZoneTest<BCT>(origin) == true);

} // end compare cauchy ratios for cycle num greater than 1 cauchy class test


/**
	In the cauchy endgame we want to compute cauchy loops and make sure that we have that the max and min norms of the samples around the 
	origin are within some heuristic value. The function InitialCauchyLoops does this while holding onto the cycle number. 
	This test case is making sure we succeed and get the correct cycle number. 
*/
BOOST_AUTO_TEST_CASE(initial_cauchy_loops)
{
	
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> pseg_times;
	bertini::SampCont<BCT> pseg_samples;

	auto time = ComplexFromString("0.1");
	Vec<BCT> sample(1);
	auto origin = BCT(0,0);

	pseg_times.push_back(time);
	sample << ComplexFromString("7.999999999999999e-01", "2.168404344971009e-19"); // 
	pseg_samples.push_back(sample);

	TestedEGType my_endgame(tracker);

	my_endgame.SetPSEGSamples(pseg_samples);
	my_endgame.SetPSEGTimes(pseg_times);

	auto success_of_initial_cauchy_loops =  my_endgame.InitialCauchyLoops<BCT>(origin);

	BOOST_CHECK(success_of_initial_cauchy_loops == SuccessCode::Success);
	BOOST_CHECK(my_endgame.CycleNumber() == 1);

}//end initial_cauchy_loops

/*
	This test case is to compute pre cauchy loops around a nonzero target time. Need to ensure we will succeed at closing the loops, 
	and ensure that our ratios around the target time are small enough. 
*/
BOOST_AUTO_TEST_CASE(initial_cauchy_loops_nonzero_target_time)
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
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

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
	auto center_for_loops = ComplexFromString(".15","-.01");

	pseg_times.push_back(start_time);
	start_sample << ComplexFromString("3.603621541081173e-01", "2.859583229930518e-18"); // 
	pseg_samples.push_back(start_sample);

	TestedEGType my_endgame(tracker);

	my_endgame.SetPSEGSamples(pseg_samples);
	my_endgame.SetPSEGTimes(pseg_times);

	auto success_of_initial_cauchy_loops =  my_endgame.InitialCauchyLoops<BCT>(center_for_loops);

	BOOST_CHECK(success_of_initial_cauchy_loops == SuccessCode::Success);


}//end inital_cauchy_loops nonzero target time.





/**
	In the cauchy endgame we want to compute cauchy loops and make sure that we have that the max and min norms of the samples around the 
	origin are within some heuristic value. The function InitialCauchyLoops does this while holding onto the cycle number. 
	This test case is making sure we succeed and get the correct cycle number. 
*/
BOOST_AUTO_TEST_CASE(initial_cauchy_loops_cycle_num_greater_than_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> pseg_times;
	bertini::SampCont<BCT> pseg_samples;



	BCT time(1);
	Vec<BCT> sample(1);
	auto origin = BCT(0,0);


	time = ComplexFromString("0.1");
	pseg_times.push_back(time);
	sample << ComplexFromString("9.000000000000001e-01", "4.358898943540673e-01"); // 
	pseg_samples.push_back(sample);


	TestedEGType my_endgame(tracker);

	my_endgame.SetPSEGSamples(pseg_samples);
	my_endgame.SetPSEGTimes(pseg_times);

	BOOST_CHECK(my_endgame.InitialCauchyLoops<BCT>(origin) == SuccessCode::Success);
	BOOST_CHECK(my_endgame.CycleNumber() == 2);
	
}// end initial_cauchy_loops_cycle_num_greater_than_1



/**
	The function tested is ComputeFirstApproximation. This function is the culmination of tracking samples until our estimation of 
	c over k is stabilized and we have our Cauchy loop ratios within a tolerance. When this is acheived we can use the power series approximation
	to compute our first extrapolant at the origin. After this we will use the cauchy integral formula to compute all further extrapolants.
*/
BOOST_AUTO_TEST_CASE(first_approximation_using_pseg)
{
	DefaultPrecision(ambient_precision);


	bertini::System sys;
	Var x = MakeVariable("x"), t = MakeVariable("t");
	VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	// Define homotopy system
	sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3) + 1)*t);

	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	auto origin = BCT(0,0);
	Vec<BCT> x_origin(1);
	x_origin << BCT(1,0);

	Vec<BCT> first_approx(1);
	auto time = ComplexFromString("0.1");
	Vec<BCT> sample(1);
	sample << ComplexFromString("0.5","0"); // f(.1) = 5.000000000000001e-01 9.084258952712920e-17 from bertini classic

#ifdef B2_OBSERVE_TRACKERS
	GoryDetailLogger<TrackerType> tons_of_detail;
	tracker.AddObserver(&tons_of_detail);
#endif

	TestedEGType my_endgame(tracker);



	auto first_approx_success = my_endgame.InitialPowerSeriesApproximation(time,sample,origin,first_approx);
	BOOST_CHECK(first_approx_success == SuccessCode::Success);
	BOOST_CHECK((first_approx - x_origin).template lpNorm<Eigen::Infinity>() < 1e-2);
	BOOST_CHECK(my_endgame.CycleNumber() == 3);

}// end first_approximation_using_pseg

/**
	The function tested is ComputeFirstApproximation. This function is the culmination of tracking samples until our estimation of 
	c over k is stabilized and we have our Cauchy loop ratios within a tolerance. When this is acheived we can use the power series approximation
	to compute our first extrapolant at the origin. After this we will use the cauchy integral formula to compute all further extrapolants.

	Specifically this test case will be attempting to do this with a target time that is not the origin. 
*/
BOOST_AUTO_TEST_CASE(first_approximation_using_pseg_nonzero_target_time)
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
		
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

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
	x_to_check_against << ComplexFromString("0.424892","0.0136983");

	pseg_times.push_back(start_time);
	start_sample << ComplexFromString("3.603621541081173e-01", "2.859583229930518e-18"); // 
	pseg_samples.push_back(start_sample);

	TestedEGType my_endgame(tracker);

	my_endgame.SetPSEGSamples(pseg_samples);
	my_endgame.SetPSEGTimes(pseg_times);

#ifdef B2_OBSERVE_TRACKERS
	GoryDetailLogger<TrackerType> tons_of_detail;
	tracker.AddObserver(&tons_of_detail);
#endif


	auto first_approx_success = my_endgame.InitialPowerSeriesApproximation(start_time,start_sample,target_time,first_approx);

	BOOST_CHECK((first_approx - x_to_check_against).template lpNorm<Eigen::Infinity>() < 1e-2);
	BOOST_CHECK(my_endgame.CycleNumber() == 1);

}// end first_approximation_using_pseg_nonzero_target_time


/**
	This test case uses all the sample points collected by CircleTrack around a non-singular point to compute an extrapolant using the 
	trapezoidal rule for integration. 
*/
BOOST_AUTO_TEST_CASE(compute_cauchy_approximation_cycle_num_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> cauchy_times; 
	bertini::SampCont<BCT> cauchy_samples; 


	auto time = ComplexFromString("0.1");
	Vec<BCT> sample(1);
	auto origin = BCT(0,0);
	Vec<BCT> x_origin(1);

	cauchy_times.push_back(time);
	sample << ComplexFromString("7.999999999999999e-01", "2.168404344971009e-19"); // 
	cauchy_samples.push_back(sample);
	x_origin << BCT(1,0);

	TestedEGType my_endgame(tracker);

	my_endgame.SetCauchyTimes(cauchy_times);
	my_endgame.SetCauchySamples(cauchy_samples);

	auto first_track_success =  my_endgame.CircleTrack(time,origin,sample);

	my_endgame.CycleNumber(1); // manually set cycle number to 1 for this test

	Vec<BCT> first_cauchy_approx;
	auto code = my_endgame.ComputeCauchyApproximationOfXAtT0<BCT>(first_cauchy_approx);

	BOOST_CHECK((first_cauchy_approx - x_origin).template lpNorm<Eigen::Infinity>() < 1e-5);

}// end compute_cauchy_approximation_cycle_num_1



/**
	This test case uses all the sample points collected by CircleTrack around a non-singular point to compute an extrapolant using the 
	trapezoidal rule for integration. CircleTrack is called twice to close up the loop. 
*/
BOOST_AUTO_TEST_CASE(compute_cauchy_approximation_cycle_num_greater_than_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	bertini::TimeCont<BCT> cauchy_times; 
	bertini::SampCont<BCT> cauchy_samples; 



	BCT time(1);
	Vec<BCT> sample(1);
	auto origin = BCT(0,0);
	Vec<BCT> x_origin(1);

	time = ComplexFromString("0.1");
	cauchy_times.push_back(time);
	sample << ComplexFromString("9.000000000000001e-01", "4.358898943540673e-01"); // 
	cauchy_samples.push_back(sample);
	x_origin << BCT(1,0);


	TestedEGType my_endgame(tracker);

	my_endgame.SetCauchyTimes(cauchy_times);
	my_endgame.SetCauchySamples(cauchy_samples);

	auto first_track_success =  my_endgame.CircleTrack(time,origin,sample);
	auto second_track_success = my_endgame.CircleTrack(my_endgame.GetCauchyTimes<BCT>().back(),origin,my_endgame.GetCauchySamples<BCT>().back());

	my_endgame.CycleNumber(2);

	Vec<BCT> first_cauchy_approx;
	auto code = my_endgame.ComputeCauchyApproximationOfXAtT0<BCT>(first_cauchy_approx);

	BOOST_CHECK((first_cauchy_approx - x_origin).template lpNorm<Eigen::Infinity>() < 1e-5);
}// end compute_cauchy_approximation_cycle_num_greater_than_1


/**
	To actually use the Cauchy integral formula we must track around the origin and collect samples. This is done by ComputeCauchySamples. 
*/
BOOST_AUTO_TEST_CASE(cauchy_samples_cycle_num_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);



	BCT time(1);
	Vec<BCT> sample(1);
	auto origin = BCT(0,0);


	time = ComplexFromString(".1");

	sample << ComplexFromString("7.999999999999999e-01", "2.168404344971009e-19"); // 


	TestedEGType my_endgame(tracker);

	auto finding_cauchy_samples_success = my_endgame.ComputeCauchySamples(time,origin,sample);

	BOOST_CHECK((my_endgame.GetCauchySamples<BCT>().back() - my_endgame.GetCauchySamples<BCT>().front()).template lpNorm<Eigen::Infinity>() < 1e-5);
	BOOST_CHECK(my_endgame.GetCauchySamples<BCT>().size() == 4);
	BOOST_CHECK(my_endgame.CycleNumber() == 1);

}// end find_cauchy_samples_cycle_num_1




/**
	To actually use the Cauchy integral formula we must track around the origin and collect samples. This is done by ComputeCauchySamples. 
*/
BOOST_AUTO_TEST_CASE(find_cauchy_samples_cycle_num_greater_than_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-6,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);


	auto time = ComplexFromString("0.1");

	Vec<BCT> sample(1);
	auto origin = BCT(0,0);
	sample << ComplexFromString("9.000000000000001e-01", "4.358898943540673e-01"); // 


	TestedEGType my_endgame(tracker);

	auto finding_cauchy_samples_success = my_endgame.ComputeCauchySamples(time,origin,sample);

	BOOST_CHECK((my_endgame.GetCauchySamples<BCT>().back() - my_endgame.GetCauchySamples<BCT>().front()).template lpNorm<Eigen::Infinity>() < 1e-6);
	BOOST_CHECK_EQUAL(my_endgame.GetCauchySamples<BCT>().size(), 7);
	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(), 2);

}// end find_cauchy_samples_cycle_num_greater_than_1





/**
	Full blown test to see if we can actually find the non singular point at the origin. 
*/
BOOST_AUTO_TEST_CASE(full_test_cycle_num_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction((x-1)*(1-t) + (x+1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);


	
	


	auto time = ComplexFromString(".1");

	Vec<BCT> sample(1);
	sample << ComplexFromString("7.999999999999999e-01", "2.168404344971009e-19"); // 

	Vec<BCT> solution(1);
	solution << BCT(1,0);

	TestedEGType my_endgame(tracker);

	auto cauchy_endgame_success = my_endgame.Run(time,sample);

	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - solution).template lpNorm<Eigen::Infinity>() < 1e-5);
	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(), 1);

}// end full_test_cycle_num_1




/**
	Full blown test to see if we can actually find the singular point at the origin. 
*/
BOOST_AUTO_TEST_CASE(full_test_cycle_num_greater_than_1)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;
	newton_preferences.max_num_newton_iterations = 2;
	newton_preferences.min_num_newton_iterations = 1;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	auto time = ComplexFromString("0.1");

	Vec<BCT> sample(1);
	sample << ComplexFromString("9.000000000000001e-01", "4.358898943540673e-01");

	Vec<BCT> x_origin(1); 
	x_origin << BCT(1,0);

	
	TestedEGType my_endgame(tracker);

	auto cauchy_endgame_success = my_endgame.Run(time,sample);


	BOOST_CHECK(cauchy_endgame_success==SuccessCode::Success);
	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - x_origin).template lpNorm<Eigen::Infinity>() < 1e-5);
	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(), 2);
}// end full_test_cycle_num_greater_than_1





/*
	Full blown test to see if we can actually find the singular point at the origin. 
*/
BOOST_AUTO_TEST_CASE(cauchy_endgame_test_cycle_num_greater_than_1_base)
{
	DefaultPrecision(ambient_precision);

	System sys;
	Var x = MakeVariable("x");
	Var t = MakeVariable("t"); 

	sys.AddFunction( pow(x-1,2)*(1-t) + (pow(x,2) + 1)*t);

	VariableGroup vars{x};
	sys.AddVariableGroup(vars); 
	sys.AddPathVariable(t);


	auto precision_config = PrecisionConfig(sys);

	TrackerType tracker(sys);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;
	newton_preferences.max_num_newton_iterations = 2;
	newton_preferences.min_num_newton_iterations = 1;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);

	tracker.ReinitializeInitialStepSize(false);


	BCT time(1);
	Vec<BCT> sample(1);
	Vec<BCT> x_origin(1);


	time = ComplexFromString("0.1");

	sample << ComplexFromString("9.000000000000001e-01", "4.358898943540673e-01"); // 
	x_origin << BCT(1,0);

	
	TestedEGType my_endgame(tracker);
	auto cauchy_endgame_success = my_endgame.Run(time,sample);

	BOOST_CHECK(cauchy_endgame_success==SuccessCode::Success);
	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - x_origin).template lpNorm<Eigen::Infinity>() < 1e-5);
	BOOST_CHECK_EQUAL(my_endgame.CycleNumber(), 2);
}// end cauchy_endgame_test_cycle_num_greater_than_1




/**
	Full blown test to see if we can actually find the non singular point at the origin. This example has multiple variables. 
*/
BOOST_AUTO_TEST_CASE(cauchy_multiple_variables)
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
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-5,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);


	BCT current_time(1);
	Vec<BCT> current_space(2);
	current_time = ComplexFromString(".1");
	current_space <<  ComplexFromString("5.000000000000001e-01", "9.084258952712920e-17") ,ComplexFromString("9.000000000000001e-01","4.358898943540673e-01");

	Vec<BCT> correct(2);
	correct << BCT(1,0),BCT(1,0);

	EndgameConfig endgame_settings;
	CauchyConfig cauchy_settings;
	SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,cauchy_settings,endgame_settings,security_settings);

	my_endgame.Run(current_time,current_space);

	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - correct).template lpNorm<Eigen::Infinity>() < 1e-11);

}// end cauchy_multiple_variables



/**
	Test to see if we can compute Cauchy samples around a nonzero target time.
*/
BOOST_AUTO_TEST_CASE(compute_cauchy_samples_nonzero_target_time)
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
		
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

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

	start_sample << ComplexFromString("3.603621541081173e-01", "2.859583229930518e-18"); // 



#ifdef B2_OBSERVE_TRACKERS
	GoryDetailLogger<TrackerType> tons_of_detail;
	tracker.AddObserver(&tons_of_detail);
#endif

	EndgameConfig endgame_settings;
	CauchyConfig cauchy_settings;
	SecurityConfig security_settings;
	TestedEGType my_endgame(tracker,cauchy_settings,endgame_settings,security_settings);

	auto cauchy_samples_success = my_endgame.ComputeCauchySamples(start_time,target_time,start_sample);


	BOOST_CHECK(cauchy_samples_success == SuccessCode::Success);
}// end compute_cauchy_samples_nonzero_target_time


/**
	Full blown test to see if we can actually track using an endgame to a nonzero target time. 
*/
BOOST_AUTO_TEST_CASE(cauchy_full_run_nonzero_target_time)
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
		
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

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





	GoryDetailLogger<TrackerType> tons_of_detail;
	tracker.AddObserver(&tons_of_detail);

	EndgameConfig endgame_settings;
	CauchyConfig cauchy_settings;
	SecurityConfig security_settings;
	TestedEGType my_endgame(tracker,cauchy_settings,endgame_settings,security_settings);

	auto endgame_success = my_endgame.Run(start_time,start_sample,target_time);
	BOOST_CHECK(endgame_success == SuccessCode::Success);

	BOOST_CHECK((my_endgame.FinalApproximation<BCT>() - x_to_check_against).template lpNorm<Eigen::Infinity>() < 1e-10);

}// end cauchy_full_run_nonzero_target_time





/**
Griewank Osborne is a classic example. Here we allow x and y to mix. There are six paths to be tracked and we know there values 
at t = 0.1. 

Three of these paths will converge to origin and three will diverge to infinity triggering a SecurityMaxNorm issue. 

This test will check to see if the answer that we converge on compared to the correct answer are withing the track tolerance during 
the endgame. 
*/
BOOST_AUTO_TEST_CASE(griewank_osborne)
{
	
	

	DefaultPrecision(ambient_precision);

	bertini::System griewank_osborn_sys;
	Var x = MakeVariable("x"), t = MakeVariable("t"), y = MakeVariable("y");
	VariableGroup vars{x,y};
	griewank_osborn_sys.AddVariableGroup(vars); 

	griewank_osborn_sys.AddFunction((mpq_rational(29,16))*pow(x,3)-2*x*y);
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


	auto precision_config = PrecisionConfig(final_griewank_osborn_system);

	TrackerType tracker(final_griewank_osborn_system);
	
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;

	tracker.Setup(TestedPredictor,
                1e-6,
                1e5,
                stepping_preferences,
                newton_preferences);
	
	tracker.PrecisionSetup(precision_config);


	BCT t_start(1), t_endgame_boundary(0.1);
	std::vector<Vec<BCT> > griewank_solutions;
	std::vector<Vec<BCT> > griewank_homogenized_solutions;
	for (unsigned ii = 0; ii < griewank_TD.NumStartPoints(); ++ii)
	{
		DefaultPrecision(ambient_precision);
		final_griewank_osborn_system.precision(ambient_precision);
		auto start_point = griewank_TD.StartPoint<BCT>(ii);

		Vec<BCT> result;
		SuccessCode tracking_success;

		tracking_success = tracker.TrackPath(result,t_start,t_endgame_boundary,start_point);
		BOOST_CHECK(tracking_success==SuccessCode::Success);

		griewank_homogenized_solutions.push_back(result);
		griewank_solutions.push_back(final_griewank_osborn_system.DehomogenizePoint(result));
	}

	Vec<BCT> correct(2);
	correct << BCT(0,0),BCT(0,0);

	EndgameConfig endgame_settings;
	CauchyConfig cauchy_settings;
	SecurityConfig security_settings;

	TestedEGType my_endgame(tracker,cauchy_settings,endgame_settings,security_settings);

#ifdef B2_OBSERVE_TRACKERS
	GoryDetailLogger<TrackerType> tons_of_detail;
	tracker.AddObserver(&tons_of_detail);
#endif

	unsigned num_paths_diverging = 0;
	unsigned num_paths_converging = 0;
	for (auto& s : griewank_homogenized_solutions) //current_space_values)
	{
		auto init_prec = Precision(s(0));
		DefaultPrecision(init_prec);
		final_griewank_osborn_system.precision(init_prec);

		SuccessCode endgame_success = my_endgame.Run(BCT(t_endgame_boundary),s);

		if(endgame_success == SuccessCode::Success)
		{
			BOOST_CHECK_EQUAL(Precision(my_endgame.FinalApproximation<BCT>()), tracker.CurrentPrecision());
			num_paths_converging++;
		}
		if(endgame_success == SuccessCode::SecurityMaxNormReached || endgame_success == SuccessCode::GoingToInfinity)
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
BOOST_AUTO_TEST_CASE(total_degree_start_system)
{
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


	auto final_system = (1-t)*sys + t*TD;
	final_system.AddPathVariable(t);

	auto precision_config = PrecisionConfig(final_system);

	auto tracker = TrackerType(final_system);
	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;
	tracker.Setup(TestedPredictor,
	              	1e-5, 1e5,
					stepping_preferences, newton_preferences);

	tracker.PrecisionSetup(precision_config);
	
	unsigned num_paths_to_track{1};
	BCT t_start(1);
	std::vector<Vec<BCT> > solutions;
	std::vector<Vec<BCT> > homogenized_solutions;
	for (unsigned ii = 0; ii < num_paths_to_track; ++ii)
	{
		DefaultPrecision(ambient_precision);
		BCT t_endgame_boundary{0.1};
		final_system.precision(ambient_precision);
		auto start_point = TD.StartPoint<BCT>(ii);

		Vec<BCT> result;
		SuccessCode tracking_success;

		tracking_success = tracker.TrackPath(result,t_start,t_endgame_boundary,start_point);
		BOOST_CHECK(tracking_success==SuccessCode::Success);

		homogenized_solutions.push_back(result);
		solutions.push_back(final_system.DehomogenizePoint(result));

		BOOST_CHECK_EQUAL(Precision(result), tracker.CurrentPrecision());
	}

	Vec<BCT> correct(2);
	correct << BCT(1,0),BCT(1,0);

	TestedEGType my_endgame(tracker);

	tracker.Setup(TestedPredictor,
	              	1e-6, 1e5,
					stepping_preferences, newton_preferences);

#ifdef B2_OBSERVE_TRACKERS
			GoryDetailLogger<TrackerType> tons_of_detail;
			tracker.AddObserver(&tons_of_detail);
#endif

	std::vector<Vec<BCT> > endgame_solutions;

	unsigned num_successful_occurences = 0;
	unsigned num_min_track_time_reached = 0;
	for (const auto& s : homogenized_solutions)
	{
		auto eg_prec = Precision(s);
		DefaultPrecision(eg_prec);
		BCT t_endgame_boundary{0.1};
		final_system.precision(eg_prec);
		SuccessCode endgame_success = my_endgame.Run(t_endgame_boundary,s);
		if(endgame_success == SuccessCode::Success)
		{
			BOOST_CHECK_EQUAL(Precision(my_endgame.FinalApproximation<BCT>()), tracker.CurrentPrecision());
			if((tracker.GetSystem().DehomogenizePoint(my_endgame.FinalApproximation<BCT>())-correct).template lpNorm<Eigen::Infinity>() < 1e-11)
			{
				num_successful_occurences++;
			}
		}
	}
 	BOOST_CHECK_EQUAL(num_successful_occurences,num_paths_to_track);

}






