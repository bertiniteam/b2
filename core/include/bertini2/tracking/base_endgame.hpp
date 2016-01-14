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

		* ComputeInitialSamples
		
		which we want all Endgames to be able to do.  However, the internal behaviour of a particular endgame varies -- which is why it is a different type. In particular, the sample points needed for the power series endgame are a geometric progression time values with there space values whereas for Cauchy these samples create a polygon around the origin to use the Cauchy Integral Formula.
		
		Also, there are some setting that are needed for each endgame.

		* Number of sample points.

		* Final Tolerance
		




		## Creating a new tracker type

		To create a new Tracker type, inherit from this, and override the following functions:

		\code
		
		virtual Vec<mpfr> ComputeApproximationOfXAtT0() {Vec<mpfr> base(1); base << mpfr("0","0"); return base;}

		virtual void ComputeInitialSamples() {}





		*/
			class Endgame
			{
			public:
				/**
				\brief The structs used to hold all relevant data for any endgame. 
				*/
				config::EndGame endgame_struct_;
				config::Tolerances endgame_tolerances_;
				config::Security endgame_security_;

				void SetEndgameStruct(config::EndGame new_endgame_settings)
				{
					endgame_struct_ = new_endgame_settings;
				}

				config::EndGame GetEndgameStruct()
				{
					return endgame_struct_;
				}

				void SetTolerancesStruct(config::Tolerances new_tolerances_settings)
				{
					endgame_tolerances_ = new_tolerances_settings;
				}

				config::Tolerances GetTolerancesStruct()
				{
					return endgame_tolerances_;
				}

				void SetSecurityStruct(config::Security new_endgame_security_settings)
				{
					endgame_security_ = new_endgame_security_settings;
				}

				config::Security GetSecurityStruct()
				{
					return endgame_security_;
				}


				void SetNumSamplePoints(unsigned int new_sample_points)
				{
					endgame_struct_.num_sample_points = new_sample_points;
				}
				unsigned int GetNumSamplePoints()
				{
					return endgame_struct_.num_sample_points;
				}

				void SetFinalTolerance(mpfr_float new_final_tolerance)
				{
					endgame_tolerances_.final_tolerance = new_final_tolerance;
				}
				mpfr GetFinalTolerance()
				{
					return endgame_tolerances_.final_tolerance;
				}

				void SetMinTrackTime(mpfr new_min_track_time)
				{
					endgame_struct_.min_track_time = new_min_track_time;
				}

				mpfr GetMinTrackTime()
				{
					return endgame_struct_.min_track_time;
				}

				void SetTrackToleranceDuringEndgame(mpfr_float new_track_tolerance_during_endgame)
				{
					endgame_tolerances_.track_tolerance_during_endgame = new_track_tolerance_during_endgame;
				}

				mpfr_float GetTrackToleranceDuringEndgame()
				{
					return endgame_tolerances_.track_tolerance_during_endgame;
				}

				void SetPathTruncationThreshold(mpfr_float new_path_truncation_threshold)
				{
					endgame_tolerances_.path_truncation_threshold = new_path_truncation_threshold;
				}

				mpfr_float GetPathTruncationThreshold()
				{
					return endgame_tolerances_.path_truncation_threshold;
				}

				void SetSecurityLevel(unsigned int new_security_level)
				{
					endgame_security_.level = new_security_level;
				}

				unsigned int GetSecurityLevel()
				{
					return endgame_security_.level;
				}

				void SetSecurityMaxNorm(mpfr_float new_max_norm)
				{
					endgame_security_.max_norm = new_max_norm;
				}

				mpfr GetSecurityMaxNorm()
				{
					return endgame_security_.max_norm;
				}

				void SetFinalApproximationAtOrigin(Vec<mpfr> new_approximation){endgame_struct_.final_approximation_at_origin = new_approximation;}

				Vec<mpfr> GetFinalApproximationAtOrigin() {return endgame_struct_.final_approximation_at_origin;}

			 	/**
				 \Every endgame must have a ComputeApproximationOfXAtT0 function, so we can approximate the value at the origin. 
				*/
				virtual Vec<mpfr> ComputeApproximationOfXAtT0() {Vec<mpfr> base(1); base << mpfr("0","0"); return base;}
				 	/**
				 \Every endgame must have a ComputeInitialSamples function, that computes the samples for either a Hermite interpolation or Cauchy integral computation.
				*/
				virtual void ComputeInitialSamples() {}

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
					// std::cout << "ii is " << ii << '\n';
					for(unsigned int jj=2; jj <=ii; ++jj)
					{
						// std::cout << "jj is " << jj << '\n';
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
				//As we multiply the previous result we will construct_ the highest term down to the last term.

				Result = Result * (target_time - array_of_times(0)) + finite_difference_matrix(0,0); // Last term in hermite polynomial.


				return Result;
			}

			
		}// end namespace endgame
	}// end namespace tracking 
}// end namespace bertini
				

#endif
