//This file is part of Bertini 2.0.
//
//base_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//base_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with predict.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  base_endgame.hpp
//
//  copyright 2015
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Fall 2015



#ifndef BERTINI_TRACKING_BASE_ENDGAME_HPP
#define BERTINI_TRACKING_BASE_ENDGAME_HPP

/**
\file base_endgame.hpp

\brief 

\brief Contains parent class for all endgames.
*/
#include <typeinfo>
#include "tracking.hpp"
#include "system.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include "limbo.hpp"
#include <boost/multiprecision/mpfr.hpp>
#include "mpfr_complex.hpp"


namespace bertini{ 

	namespace tracking {
		/*!
 		*  \addtogroup endgame
 		*  @{
		 */
		namespace endgame {
		/**
		\class Endgame

		\brief Base endgame class for all endgames offered in Bertini2.
	
		\see powerseries
		
		## Using an endgame

		Endgames in Bertini2 are the engine for finishing homotopy continuation where we may encounter singular solutions.
		The path is implicitly described by the system being tracked.

		## Purpose 

		Since the Bertini Endgames have common functionality, and we want to be able to call arbitrary algorithms using and tracker type, we use inheritance.  That is, there is common functionality in all endgames, such as

		* ComputeApproximationOfXAtT0
	
		




		## Creating a new endgame type

		To create a new endgame type, inherit from this, and override the following function:

		\code
		
		virtual Vec<mpfr> ComputeApproximationOfXAtT0() {Vec<mpfr> base(1); base << mpfr("0","0"); return base;}

		*/
			class Endgame
			{
			public:
				/**
				\brief The structs used to hold all relevant data for any endgame. 
				*/
				config::EndGame endgame_settings_;
				config::Tolerances endgame_tolerances_;
				config::Security endgame_security_;

				void SetMinTrackTime(mpfr new_min_track_time){endgame_settings_.min_track_time = new_min_track_time;}
				mpfr GetMinTrackTime(){return endgame_settings_.min_track_time;}

				void SetTrackToleranceDuringEndgame(mpfr_float new_track_tolerance_during_endgame){endgame_tolerances_.track_tolerance_during_endgame = new_track_tolerance_during_endgame;}
				mpfr_float GetTrackToleranceDuringEndgame(){return endgame_tolerances_.track_tolerance_during_endgame;}

				void SetPathTruncationThreshold(mpfr_float new_path_truncation_threshold){endgame_tolerances_.path_truncation_threshold = new_path_truncation_threshold;}
				mpfr_float GetPathTruncationThreshold(){return endgame_tolerances_.path_truncation_threshold;}

			 	/**
				 \Every endgame must have a ComputeApproximationOfXAtT0 function, so we can approximate the value at the origin. 
				*/
				virtual Vec<mpfr> ComputeApproximationOfXAtT0() {Vec<mpfr> base(1); base << mpfr("0","0"); return base;}

				/*
				Input: 
					endgame_time is the time when we start the endgame process usually this is .1

			    	x_endgame is the space value at endgame_time

			    	upper_bound_on_cycle_number is the largest possible cycle number we will use to compute dx_ds and s = t^(1/c).

					times is a data struct_ure that holds the time values for the corresponding samples.
					num_sample_points is the size of times

				Output: 
					A data struct_ure holding the space values corresponding to the time values in times.


				Details:
					The first sample will be (x_endgame) and the first time is endgame_time.
					From there we do a geometric progression using the sample factor which by default is 1/2.
					The next_time = endgame_time * sample_factor.
					We can track then to the next_time and that construct_s the next_sample.
				*/
				template<typename ComplexType, typename TrackerType>
				void ComputeInitialSamples(const TrackerType & endgame_tracker_, const ComplexType endgame_time,const Vec<ComplexType> x_endgame, unsigned int num_samples , std::deque<ComplexType> & times, std::deque< Vec<ComplexType> > & samples) // passed by reference to allow times to be filled as well.
				{
					samples.push_back(x_endgame);
					times.push_back(endgame_time);
					Vec<ComplexType> next_sample;

					for(int ii=2; ii <= num_samples; ++ii)//start at 2 since first sample is at the endgame boundary.
					{ 
						ComplexType next_time = times.back() * endgame_settings_.sample_factor;	

						SuccessCode tracking_success = endgame_tracker_.TrackPath(next_sample,times.back(),next_time,samples.back());

						samples.push_back(next_sample);
						times.push_back(next_time);		
					}				
				}

			};


			/*
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
			*/
			template<typename ComplexType>		
				Vec<ComplexType> HermiteInterpolateAndSolve(ComplexType const& target_time, const unsigned int num_sample_points, const std::deque<ComplexType> & times, const std::deque< Vec<ComplexType> > & samples, const std::deque< Vec<ComplexType> > & derivatives)
			{
				Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);

				Mat< Vec<ComplexType> > finite_difference_matrix(2*num_sample_points,2*num_sample_points);
				Vec<ComplexType> array_of_times(2*num_sample_points);
				
				for(unsigned int ii=0; ii<num_sample_points; ++ii)
				{ 
					finite_difference_matrix(2*ii,0) = samples[ii];			/*  F[2*i][0]    = samples[i];    */
     				finite_difference_matrix(2*ii+1,0) = samples[ii]; 		/*  F[2*i+1][0]  = samples[i];    */
      				finite_difference_matrix(2*ii+1,1) = derivatives[ii];	/*  F[2*i+1][1]  = derivatives[i]; */
     				array_of_times(2*ii) = times[ii];						/*  z[2*i]       = times[i];       */
     				array_of_times(2*ii+1) =  times[ii];					/*  z[2*i+1]     = times[i];       */
				}

				//Add first round of finite differences to fill out rest of matrix. 
				for(unsigned int ii=1; ii< num_sample_points; ++ii)
				{
					finite_difference_matrix(2*ii,1) = (ComplexType(1)/(array_of_times(2*ii) - array_of_times(2*ii - 1))) * (finite_difference_matrix(2*ii,0) - finite_difference_matrix(2*ii - 1,0));
				
				}

				//Filliing out finite difference matrix to get the diagonal for hermite interpolation polyonomial.
				for(unsigned int ii=2; ii < 2*num_sample_points; ++ii)
				{
					for(unsigned int jj=2; jj <=ii; ++jj)
					{
						finite_difference_matrix(ii,jj) = (ComplexType(1)/(array_of_times(ii) - array_of_times(ii - jj)))*(finite_difference_matrix(ii,jj-1) - finite_difference_matrix(ii-1,jj-1));						
					}
				}

				 unsigned int ii = num_sample_points - 1 ;
				 auto Result = finite_difference_matrix(2*num_sample_points - 1,2*num_sample_points - 1); 
				 //Start of Result from Hermite polynomial, this is using the diagonal of the 
				 //finite difference matrix.

				 while(ii >= 1)
				{
					Result = (Result*(target_time - array_of_times(ii)) + finite_difference_matrix(2*ii,2*ii)) * (target_time - array_of_times(ii - 1)) + finite_difference_matrix(2*ii - 1,2*ii - 1);  
					ii--;
				}
				//This builds the hermite polynomial from the highest term down. 
				//As we multiply the previous result we will construct the highest term down to the last term.
				Result = Result * (target_time - array_of_times(0)) + finite_difference_matrix(0,0); // Last term in hermite polynomial.
				return Result;
			}

			
		}// end namespace endgame
	}// end namespace tracking 
}// end namespace bertini
				

#endif
