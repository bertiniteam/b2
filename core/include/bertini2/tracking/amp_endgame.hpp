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

\brief Contains parent class, Endgame, the parent class for all endgames.
*/


#include <typeinfo>
#include "tracking.hpp"
#include "system.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include "limbo.hpp"
#include <boost/multiprecision/mpfr.hpp>
#include "mpfr_complex.hpp"

#include "bertini2/enable_permuted_arguments.hpp"

namespace bertini{ 

	namespace tracking {

		namespace endgame {
			template<typename T> using SampCont = std::deque<Vec<T> >;
			template<typename T> using TimeCont = std::deque<T>;


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

			/**
			\brief

			\param[out] endgame_tracker_ The tracker used to compute the samples we need to start an endgame. 
			\param endgame_time The time value at which we start the endgame. 
			\param x_endgame The current space point at endgame_time.
			\param times A deque that will hold all the time values of the samples we are going to use to start the endgame. 
			\param samples a deque that will hold all the samples corresponding to the time values in times. 

			\tparam CT The complex number type.
			*/			
			template<typename CT>		
				Vec<CT> HermiteInterpolateAndSolve(CT const& target_time, const unsigned int num_sample_points, const std::deque<CT> & times, const SampCont<CT> & samples, const SampCont<CT> & derivatives)
			{
				Mat< Vec<CT> > finite_difference_matrix(2*num_sample_points,2*num_sample_points);
				Vec<CT> array_of_times(2*num_sample_points);
				
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
					finite_difference_matrix(2*ii,1) = (CT(1)/(array_of_times(2*ii) - array_of_times(2*ii - 1))) * (finite_difference_matrix(2*ii,0) - finite_difference_matrix(2*ii - 1,0));
				
				}

				//Filliing out finite difference matrix to get the diagonal for hermite interpolation polyonomial.
				for(unsigned int ii=2; ii < 2*num_sample_points; ++ii)
				{
					for(unsigned int jj=2; jj <=ii; ++jj)
					{
						finite_difference_matrix(ii,jj) = (CT(1)/(array_of_times(ii) - array_of_times(ii - jj)))*(finite_difference_matrix(ii,jj-1) - finite_difference_matrix(ii-1,jj-1));						
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

		/**
		\class Endgame

		\brief Base endgame class for all endgames offered in Bertini2.
	
		\see PowerSeriesEndgame
		\see CauchyEndgame
		
		## Using an endgame

		Endgames in Bertini2 are the engine for finishing homotopy continuation where we may encounter singular solutions.
		The path is implicitly described by the system being tracked.

		## Purpose 

		Since the Bertini Endgames have common functionality, and we want to be able to call arbitrary algorithms using and tracker type, we use inheritance. That is, there is common functionality in all endgames, such as

		ComputeInitialSamples

		Also, there are settings that will be kept at this level to not duplicate code. 
			
		## Creating a new endgame type

		 To create a new endgame type, inherit from this class. 
		*/
			
			template<class TrackerType>
			class AMPEndgame
			{

			protected:
				// state variables
				mutable std::tuple<Vec<dbl>, Vec<mpfr> > final_approximation_at_origin_; 
				mutable unsigned int cycle_number_ = 0; 


				/**
				\brief Settings that will be used in every endgame. For example, the number of samples points, the cycle number, and the final approximation of that we compute 
				using the endgame. 
				*/
				config::Endgame endgame_settings_;
				/**
				\brief There are tolerances that are specific to the endgame. These settings are stored inside of this data member. 
				*/
				config::Tolerances tolerances_;
				/**
				During the endgame we may be checking that we are not computing when we have detected divergent paths or other undesirable behavior. The setttings for these checks are 
				in this data member. 
				*/
				config::Security security_;

				/**
				\brief A tracker tha must be passed into the endgame through a constructor. This tracker is what will be used to track to all time values in the Cauchy endgame. 
				*/
				const TrackerType& tracker_;

				unsigned precision_;
				bool preserve_precision_ = true; ///< Whether the endgame should change back to the initial precision after running 


				/**
				\brief Change precision of all temporary internal state variables.

				This excludes those which cannot be re-written without copying.

				\brief new_precision The new precision to adjust to.
				*/
				void ChangeTemporariesPrecision(unsigned new_precision) const
				{
					tracker_.precision(new_precision);

					TrackerType::ChangeTemporariesPrecisionImpl();
				}

				/**
				\brief Converts from multiple to different precision multiple precision

				Changes the precision of the internal temporaries to desired precision

				\param new_precision The new precision.
				*/
				void MultipleToMultiple(unsigned new_precision) const
				{
					mpfr_float::default_precision(new_precision);
					precision_ = new_precision;

					const auto& source_point = std::get<Vec<mpfr> >(final_approximation_at_origin_);
					auto& target_point = std::get<Vec<mpfr> >(final_approximation_at_origin_);
					target_point.resize(source_point.size());

					for (unsigned ii=0; ii<source_point.size(); ii++)
					{
						target_point(ii).precision(new_precision);
						target_point(ii) = mpfr(source_point(ii));
					}
					ChangeTemporariesPrecision(new_precision);

					TrackerType::MultipleToMultipleImpl();
				}

				/**
				\brief Converts from multiple to different precision multiple precision

				Changes the precision of the internal temporaries to desired precision

				\param new_precision The new precision.
				*/
				void DoubleToMultiple(unsigned new_precision) const
				{
					mpfr_float::default_precision(new_precision);
					precision_ = new_precision;
					const auto& source_point = std::get<Vec<dbl> >(final_approximation_at_origin_);
					auto& target_point = std::get<Vec<mpfr> >(final_approximation_at_origin_);

					target_point.resize(source_point.size());
					for (unsigned ii=0; ii<source_point.size(); ii++)
					{
						target_point(ii).precision(new_precision);
						target_point(ii) = mpfr(source_point(ii));
					}
					ChangeTemporariesPrecision(new_precision);

					TrackerType::DoubleToMultipleImpl();
				}

				/**
				\brief Converts from multiple to different precision multiple precision

				Changes the precision of the internal temporaries to desired precision

				\param new_precision The new precision.
				*/
				void MultipleToDouble() const
				{
					mpfr_float::default_precision(DoublePrecision());
					precision_ = DoublePrecision();

					const auto& source_point = std::get<Vec<mpfr> >(final_approximation_at_origin_);
					auto& target_point = std::get<Vec<dbl> >(final_approximation_at_origin_);
					target_point.resize(source_point.size());

					for (unsigned ii=0; ii<source_point.size(); ii++)
						target_point(ii) = dbl(source_point(ii));

					tracker_.precision(16);

					TrackerType::MultipleToDoubleImpl();
				}

				bool PrecisionSanityCheck()
				{
					return true;
				}
			public:	

				/**
				Change precision of endgame object to next_precision.  Converts the internal temporaries, and adjusts precision of system. 

				\param new_precision The precision to change to.
				\return SuccessCode indicating whether the change was successful.  If the precision increases, and the refinement loop fails, this could be not Success.  Changing down is guaranteed to succeed.
				*/
				SuccessCode ChangePrecision(unsigned new_precision) const
				{
					if (new_precision==precision_) // no op
						return SuccessCode::Success;

					if (new_precision==DoublePrecision() && precision_>DoublePrecision())
						MultipleToDouble();
					else if(new_precision > DoublePrecision() && precision_ == DoublePrecision())
						DoubleToMultiple(new_precision);
					else
						MultipleToMultiple(new_precision);


					#ifndef BERTINI_DISABLE_ASSERTS
					assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
					#endif

					return SuccessCode::Success;
				}


				explicit AMPEndgame(TrackerType const& tr, const std::tuple< const config::Endgame&, const config::Security&, const config::Tolerances& >& settings )
			      : tracker_(tr),
			      	endgame_settings_( std::get<0>(settings) ),
			        security_( std::get<1>(settings) ),
			        tolerances_( std::get<2>(settings) )
			   	{}

			    template< typename... Ts >
   				AMPEndgame(TrackerType const& tr, const Ts&... ts ) : AMPEndgame(tr, Unpermute< config::Endgame, config::Security, config::Tolerances >( ts... ) ) 
   				{}


				unsigned CycleNumber() const { return cycle_number_;}
				void CycleNumber(unsigned c) { cycle_number_ = c;}
				void IncrementCycleNumber(unsigned inc) { cycle_number_ += inc;}

				
				const config::Endgame& EndgameSettings()
				{
					return endgame_settings_;
				}

				const config::Tolerances& Tolerances()
				{
					return tolerances_;
				}

				const config::Security& SecuritySettings()
				{
					return security_;
				}



				/**
				\brief Setter for the general settings in tracking_conifg.hpp under Endgame.
				*/
				void SetEndgameSettings(config::Endgame new_endgame_settings){endgame_settings_ = new_endgame_settings;}


				

				/**
				\brief Setter for the security settings in tracking_conifg.hpp under Security.
				*/
				void SetSecuritySettings(config::Security new_endgame_security_settings){ security_ = new_endgame_security_settings;}


				/**
				\brief Setter for the tolerance settings in tracking_conifg.hpp under Tolerances.
				*/
				void SetToleranceSettings(config::Tolerances new_tolerances_settings){tolerances_ = new_tolerances_settings;}


				/**
				\brief Getter for the tracker used inside of an instance of the endgame. 
				*/
				const TrackerType & GetTracker(){return tracker_;}

				template<typename CT>
				const Vec<CT>& FinalApproximation() const {return std::get<Vec<CT> >(final_approximation_at_origin_);}

				const System& GetSystem() const { return tracker_.GetSystem();}
				/*
				Input:  endgame_tracker_: a template parameter for the endgame created. This tracker will be used to compute our samples. 
						endgame_time: is the time when we start the endgame process usually this is .1
						x_endgame: is the space value at endgame_time
						times: a deque of time values. These values will be templated to be CT 
						samples: a deque of sample values that are in correspondence with the values in times. These values will be vectors with entries of CT. 

				Output: Since times and samples are sent in by reference there is no output value.


				Details:
					The first sample will be (x_endgame) and the first time is start_time.
					From there we do a geometric progression using the sample factor which by default is 1/2.
					The next_time = start_time * sample_factor.
					We can track then to the next_time and that construct the next_sample. This is done for Power Series Endgame and the Cauchy endgame.
				*/
				/**
				\brief Whether we are using the PowerSeriesEndgame or the CauchyEndgame we need to have an initial set of samples to start the endgame. This function populates two deques 
				so that we are ready to start the endgame. 

				\param[out] endgame_tracker_ The tracker used to compute the samples we need to start an endgame. 
				\param start_time The time value at which we start the endgame. 
				\param x_endgame The current space point at start_time.
				\param times A deque that will hold all the time values of the samples we are going to use to start the endgame. 
				\param samples a deque that will hold all the samples corresponding to the time values in times. 

				\tparam CT The complex number type.
				\tparam TrackerType The tracker type. 
				*/	
				template<typename CT>
				SuccessCode ComputeInitialSamples(const CT & start_time,const Vec<CT> & x_endgame, TimeCont<CT> & times, SampCont<CT> & samples) // passed by reference to allow times to be filled as well.
				{
					samples.resize(endgame_settings_.num_sample_points);
					times.resize(endgame_settings_.num_sample_points);

					samples[0] = x_endgame;
					times[0] = start_time;

					for(int ii=1; ii < endgame_settings_.num_sample_points; ++ii)
					{ 
						times[ii] = times[ii-1] * endgame_settings_.sample_factor;	
						auto tracking_success = tracker_.TrackPath(samples[ii],times[ii-1],times[ii],samples[ii-1]);
						if (tracking_success!=SuccessCode::Success)
							return tracking_success;
					}				
					return SuccessCode::Success;
				}

			};

			/*
				Input: Two numbers d and n. These numbers will be used to compute d choose n. 

				Output: The value of d choose n. 

				Details: This is used to help find the closed loop tolerance for CauchyEndgame. After discussions it was deemed not necessary and that we should be using 
						the tracking tolerance during endgame as the tolerance for closing a loop. 
			
				mpfr_float ComputeCombination(mpfr_float d, mpfr_float n)
				{
					mpfr_float return_value;
					// error checking
 					if (d >= n && n >= 0)
 					{ // find the smallest of n & d - n since C(d,n) = C(d,d-n)
   						return_value = 1;
  						n = min(n, d - n);
  						for (unsigned int ii = 0; ii < n; ii++)
   						{ // update retVal
     						return_value = return_value / (ii + mpfr_float("1.0"));
     						return_value = return_value * (d - ii);
   						}
  					}
  					else
  					{ // no combinations
   						return_value = 0;
  					}
					return return_value;
				}
			*/

			

			
		}// end namespace endgame
	}// end namespace tracking 
}// end namespace bertini
				

#endif
