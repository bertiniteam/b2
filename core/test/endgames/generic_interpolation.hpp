//This file is part of Bertini 2.
//
//b2/test/endgames/interpolation.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//b2/test/endgames/interpolation.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with b2/test/endgames/interpolation.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University

// this file in intended to be included into other test files, with blanks filled in above.


using namespace bertini;
using namespace bertini::endgame;

using bertini::DefaultPrecision;


template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;

using BCT = BaseComplexType;
using BRT = Eigen::NumTraits<BCT>::Real;

template<typename ...T>
BCT ComplexFromString(T... s)
{return bertini::NumTraits<BCT>::FromString(s...);}

template<typename ...T>
BRT RealFromString(T... s)
{return bertini::NumTraits<BRT>::FromString(s...);}


BOOST_AUTO_TEST_CASE( constant_four_variate )
{
	int num_samples = 3;
	DefaultPrecision(ambient_precision);
	TimeCont<BCT> times; 
	SampCont<BCT> samples, derivatives;

	times.emplace_back(ComplexFromString("1.890235648907826359017283461724"));
	times.emplace_back(ComplexFromString("2.6751239470926520645293652"));
	times.emplace_back(ComplexFromString("3.566754186725389412876498127534875"));

	Vec<BCT> sample(4);

	sample << ComplexFromString("6.4789162736409137056734534678679"), ComplexFromString("-1.5877816237549123614917624"), ComplexFromString("5.947461892534781890236417801"),ComplexFromString("-3.87746985236816238746178293");
	for (int ii=0; ii<3; ++ii)
		samples.push_back(sample);

	Vec<BCT> derivative(4);
	derivative << ComplexFromString("0"), ComplexFromString("0"), ComplexFromString("0"),ComplexFromString("0");
	for (int ii=0; ii<3; ++ii)
		derivatives.push_back(derivative);

	BCT target_time = ComplexFromString("0.9471925368945182312341234123");
	auto result = HermiteInterpolateAndSolve(target_time,num_samples,times,samples,derivatives);
	BOOST_CHECK( (result - sample).norm() < pow(RealFromString("10"),ambient_precision-2)); 
}


/**
This test case illustrates the convergent nature of the HemiteInterpolateAndSolve function. 
The test case will construct 3 samples with derivative,time, and space values for the function x^8 + 1.
After the three samples have been constructed there will be a hermite interpolation. 
We check this against the tracking tolerance for the endgame. 
*/
BOOST_AUTO_TEST_CASE( eight_degree_univariate )
{
	DefaultPrecision(ambient_precision);



	BCT target_time(0,0); //our target time is the origin.
	unsigned int num_samples = 3;

	TimeCont<BCT> times; 
	SampCont<BCT> samples, derivatives;

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
BOOST_AUTO_TEST_CASE(eight_degree_univariate_advanced_gets_better)
{
	DefaultPrecision(ambient_precision);

	BCT target_time(0,0);
	unsigned int num_samples = 3;

	TimeCont<BCT> times; 
	SampCont<BCT> samples, derivatives;

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