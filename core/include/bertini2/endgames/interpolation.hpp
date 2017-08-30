//This file is part of Bertini 2.
//
//interpolation.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//interpolation.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with interpolation.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
\file interpolation.hpp Contains interpolation and extrapolation functions used in the endgames for estimating singular (and nonsingular) roots.
*/

#pragma once

#include "bertini2/common/config.hpp"

namespace bertini{
	namespace endgame{

/**

\brief Estimates the root to interpolating polynomial.

Input: 
		target_time is the time value that we wish to interpolate at.
		samples are space values that correspond to the time values in times. 
	 	derivatives are the dx_dt or dx_ds values at the (time,sample) values.

Output: 
		Since we have a target_time the function returns the corrsponding space value at that time. 

Details: 
		We compare our approximations to the tracked value to come up with the cycle number. 
		Also, we use the Hermite interpolation to interpolate at the origin. Once two interpolants are withing FinalTol we 
		say we have converged. 

\param[out] endgame_tracker_ The tracker used to compute the samples we need to start an endgame. 
\param endgame_time The time value at which we start the endgame. 
\param x_endgame The current space point at endgame_time.
\param times A deque that will hold all the time values of the samples we are going to use to start the endgame. 
\param samples a deque that will hold all the samples corresponding to the time values in times. 

\tparam CT The complex number type.
*/			
template<typename CT>		
	Vec<CT> HermiteInterpolateAndSolve(CT const& target_time, const unsigned int num_sample_points, const TimeCont<CT> & times, const SampCont<CT> & samples, const SampCont<CT> & derivatives, ContStart shift_from = ContStart::Back)
{
	assert((times.size() >= num_sample_points) && "must have sufficient number of sample times");
	assert((samples.size() >= num_sample_points) && "must have sufficient number of sample points");
	assert((derivatives.size() >= num_sample_points) && "must have sufficient number of derivatives");

	
	unsigned num_t, num_s, num_d;
	if (shift_from == ContStart::Back)
	{
		num_t = times.size()-1;
		num_s = samples.size()-1;
		num_d = derivatives.size()-1;
	}
	else
	{
		num_t = num_s = num_d = num_sample_points-1;
	}

	Mat< Vec<CT> > space_differences(2*num_sample_points,2*num_sample_points);
	Vec<CT> time_differences(2*num_sample_points);


	for(unsigned int ii=0; ii<num_sample_points; ++ii)
	{ 
		space_differences(2*ii,0)   = samples[    num_s-ii];		/*  F[2*i][0]    = samples[i];    */
		space_differences(2*ii+1,0) = samples[    num_s-ii]; 	/*  F[2*i+1][0]  = samples[i];    */
		space_differences(2*ii+1,1) = derivatives[num_d-ii];		/*  F[2*i+1][1]  = derivatives[i]; */
		time_differences(2*ii)      = times[      num_t-ii];		/*  z[2*i]       = times[i];       */
		time_differences(2*ii+1)    = times[      num_t-ii];		/*  z[2*i+1]     = times[i];       */
	}

	//Add first round of finite differences to fill out rest of matrix. 
	for(unsigned int ii=1; ii< num_sample_points; ++ii)
	{
		space_differences(2*ii,1) = (space_differences(2*ii,0) - space_differences(2*ii-1,0)) / (time_differences(2*ii) - time_differences(2*ii-1));
	
	}

	//Filling out finite difference matrix to get the diagonal for hermite interpolation polyonomial.
	for(unsigned int ii=2; ii < 2*num_sample_points; ++ii)
	{
		for(unsigned int jj=2; jj <=ii; ++jj)
		{
			space_differences(ii,jj) = 
				(space_differences(ii,jj-1) - space_differences(ii-1,jj-1)) 
					/ 
				(time_differences(ii) - time_differences(ii-jj));						
		}
	}

	//Start of Result from Hermite polynomial, this is using the diagonal of the 
	//finite difference matrix.
	Vec<CT> Result = space_differences(2*num_sample_points - 1,2*num_sample_points - 1); 


	//This builds the hermite polynomial from the highest term down. 
	//As we multiply the previous result we will construct the highest term down to the last term.
	for (unsigned ii=num_sample_points-1; ii >= 1; --ii)
	{
		Result = ((Result*(target_time - time_differences(ii)) + space_differences(2*ii, 2*ii)) * (target_time - time_differences(ii-1)) + space_differences(2*ii-1, 2*ii-1)).eval();  
	}
	
	// Last term in hermite polynomial.
	return (Result * (target_time - time_differences(0)) + space_differences(0,0)).eval(); 
} //re: HermiteInterpolateAndSolve

}}  // re: namespaces
