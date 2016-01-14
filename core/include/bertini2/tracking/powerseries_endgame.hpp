//This file is part of Bertini 2.0.
//
//powerseries_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//powerseries_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with predict.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  powerseries_endgame.hpp
//
//  copyright 2015
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Fall 2015



#ifndef BERTINI_TRACKING_POWERSERIES_HPP
#define BERTINI_TRACKING_POWERSERIES_HPP

/**
\file powerseries_endgame.hpp

\brief 

\brief Contains the power series endgame type, and necessary functions.
*/
#include <deque>
#include <typeinfo>
#include "tracking.hpp"
#include "system.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include "tracking/base_endgame.hpp"
#include <cstdio>
#include <new>

namespace bertini{ 

	namespace tracking {

		namespace endgame {

		// * 
		// \class PowerSeriesEndgame

		// \brief class used to finish tracking paths on a system
		
		
		// ## Explanation
		
		// The bertini::PowerSeriesEndgame class enables us to finish tracking on possibly singular paths on an arbitrary square homotopy.  

		// The intended usage is to:

		// 1. Create a system, and instantiate some settings.
		// 2. Track using a tracker of the users choice to the engame boundary. 
		// 3. Create a PowerSeriesEndgame, associating it to the system you are going to solve or track on.
		// 4. For each path being tracked send the PowerSeriesEndgame the time value and other variable values at that time. 
		// 5. The PowerSeriesEndgame, if successful, will store the target systems solution at $t = 0$.

		// ## Example Usage
		
		// Below we demonstrate a basic usage of the PowerSeriesEndgame class to find the singularity at $t = 0$. 

		// The pattern is as described above: create an instance of the class, feeding it the system to be used, and the endgame boundary time and other variable values at the endgame boundary. 

		// \code{.cpp}
		// 	mpfr_float::default_precision(30); // set initial precision.  This is not strictly necessary.
		// 	using namespace bertini::tracking;

		// 	// 1. Create the system
		// 	System sys;
		// 	Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t"), y = std::make_shared<Variable>("y");
		// 	VariableGroup vars{x,y};
		// 	sys.AddVariableGroup(vars); 
		// 	sys.AddPathVariable(t);
		// 	// Define homotopy system
		// 	sys.AddFunction((pow(x-1,3))*(1-t) + (pow(x,3) + 1)*t);
		// 	sys.AddFunction((pow(y-1,2))*(1-t) + (pow(y,2) + 1)*t);

		// 	// We make the assumption that we have tracked to t = 0.1 
		// 	mpfr current_time(1);
		// 	Vec<mpfr> current_space(2);
		// 	current_time = mpfr(".1");
		// 	current_space <<  mpfr("5.000000000000001e-01", "9.084258952712920e-17") ,mpfr("9.000000000000001e-01","4.358898943540673e-01");

		// 	//  2. Create the PowerSeriesEndgame object, associating the system to it.
		// 	endgame::PowerSeriesEndgame My_Endgame(sys);

		// 	//Calling the PSEG member function actually runs the endgame. 
		// 	My_Endgame.PSEG(current_time,current_space);

		// 	//Access solution at t = 0, using the .get_final_approximation_at_origin() member function. 
		// 	auto Answer  = My_Endgame.get_final_approximation_at_origin();
	
		
		// If this documentation is insufficient, please contact the authors with suggestions, or get involved!  Pull requests welcomed.
		
		// ## Testing

		// * Test suite driving this class: endgames_test.
		// * File: test/endgames/powerseries_class_test.cpp
		// * Functionality tested: All member functions of the PowerSeriesEndgame have been tested in Double and Multiple Precision. There is also, test running different systems to find singular solutions at t = 0.

		
			template<typename TrackerType> 
			class PowerSeriesEndgame : public Endgame{

				config::PowerSeries power_series_struct_; 

				std::deque<mpfr> times_;
				std::deque<mpfr> s_times_;

				std::deque< Vec<mpfr> > samples_;

				std::deque< Vec<mpfr> > derivatives_;
				std::deque< Vec<mpfr> > s_derivatives_;

				const TrackerType & endgame_tracker_;
			public:

				void ClearTimesDerivativesAndSamples(){times_.clear(); s_times_.clear(); samples_.clear(); derivatives_.clear(); s_derivatives_.clear();}

				void SetSampleFactor(mpfr_float new_sample_factor) {endgame_struct_.sample_factor = new_sample_factor;}
				mpfr_float GetSampleFactor(){return endgame_struct_.sample_factor;}

				void SetNumSamples(unsigned int new_num_samples) {Endgame::SetNumSamplePoints(new_num_samples);}
				unsigned int GetNumSamples(){return Endgame::GetNumSamplePoints();}

				void SetSecurityLevel(unsigned int new_security_level) {Endgame::SetSecurityLevel(new_security_level);}
				unsigned int GetSecurityLevel(){return Endgame::GetSecurityLevel();}

				void SetSecurityMaxNorm(mpfr new_max_norm){Endgame::SetSecurityMaxNorm(new_max_norm);}
				mpfr GetSecurityMaxNorm(){return Endgame::GetSecurityMaxNorm();}

				void SetCycleNumberAmplification(unsigned int new_cycle_number_amplification) { power_series_struct_.cycle_number_amplification = new_cycle_number_amplification;}
				unsigned int GetCycleNumberAmplification() { return power_series_struct_.cycle_number_amplification;}

				void SetMaxCycleNumber(unsigned int new_max_cycle_number) { power_series_struct_.max_cycle_number = new_max_cycle_number;}
				unsigned int GetMaxCycleNumber() { return power_series_struct_.max_cycle_number;}

				void SetCycleNumber(unsigned int new_cycle_number) {power_series_struct_.cycle_number = new_cycle_number;}
				unsigned int GetCycleNumber() { return power_series_struct_.cycle_number;}

				void SetUpperBoundOnCycleNumber(unsigned int new_upper_bound) { power_series_struct_.upper_bound_on_cycle_number = new_upper_bound;}
				unsigned int GetUpperBoundOnCycleNumber() { return power_series_struct_.upper_bound_on_cycle_number;}

				void SetTimes(std::deque<mpfr> times_to_set) { times_ = times_to_set;}
				std::deque< mpfr > GetTimes() {return times_;}

				void SetSTimes(std::deque<mpfr> s_times_to_set) { s_times_ = s_times_to_set;}
				std::deque< mpfr > GetSTimes() {return s_times_;}

				void SetSamples(std::deque< Vec<mpfr> > samples_to_set) { samples_ = samples_to_set;}
				std::deque< Vec<mpfr> > GetSamples() {return samples_;}

				void SetDerivatives(std::deque< Vec<mpfr> > derivatives_to_set) { derivatives_ = derivatives_to_set;}
				std::deque< Vec<mpfr> > GetDerivatives() {return derivatives_;}

				void SetSDerivatives(std::deque< Vec<mpfr> > s_derivatives_to_set) { s_derivatives_ = s_derivatives_to_set;}
				std::deque< Vec<mpfr> > GetSDerivatives() {return s_derivatives_;}

				void SetFinalTol(mpfr_float new_final_tolerance) {Endgame::SetFinalTolerance(new_final_tolerance);}
				mpfr_float GetFinalTol(){return Endgame::GetFinalTolerance();}

				// template<typename TrackingType>
	
				const TrackerType & GetEndgameTracker(){return endgame_tracker_;}

				void SetFinalApproximation(Vec<mpfr> new_final_approximation){Endgame::SetFinalApproximationAtOrigin(new_final_approximation);}

				Vec<mpfr> GetFinalApproximation() {return Endgame::GetFinalApproximationAtOrigin();}

				void SetEndgameSettings(config::EndGame new_endgame_settings){Endgame::SetEndgameStruct(new_endgame_settings);}

				config::EndGame GetEndgameSettings(){return Endgame::GetEndgameStruct();}

				void SetPowerSeriesSettings(config::PowerSeries new_power_series_settings){power_series_struct_ = new_power_series_settings;}

				config::PowerSeries GetPowerSeriesSettings(){return power_series_struct_;}

				void SetSecuritySettings(config::Security new_security_settings){Endgame::SetSecurityStruct(new_security_settings);}

				config::Security GetSecuritySettings(){return Endgame::GetSecurityStruct();}

				void SetToleranceSettings(config::Tolerances new_tolerances_settings){Endgame::SetTolerancesStruct(new_tolerances_settings);}

				config::Tolerances GetToleranceSettings(){return Endgame::GetTolerancesStruct();}

				PowerSeriesEndgame(TrackerType const& tracker) : endgame_tracker_(tracker){} //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				
				PowerSeriesEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::PowerSeries new_power_series_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances

					SetEndgameSettings(new_endgame_settings);
					SetPowerSeriesSettings(new_power_series_settings);
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::PowerSeries new_power_series_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetPowerSeriesSettings(new_power_series_settings);
					SetSecuritySettings(new_security_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::PowerSeries new_power_series_settings,config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetPowerSeriesSettings(new_power_series_settings);
					SetTolerancesStruct(new_tolerances_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);
					SetTolerancesStruct(new_tolerances_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker,config::PowerSeries new_power_series_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetPowerSeriesSettings(new_power_series_settings);
					SetSecuritySettings(new_security_settings);
					SetTolerancesStruct(new_tolerances_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::PowerSeries new_power_series_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetPowerSeriesSettings(new_power_series_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetTolerancesStruct(new_tolerances_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::PowerSeries new_power_series_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetPowerSeriesSettings(new_power_series_settings);
					SetSecuritySettings(new_security_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::PowerSeries new_power_series_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetPowerSeriesSettings(new_power_series_settings);
					SetTolerancesStruct(new_tolerances_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker,config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetSecuritySettings(new_security_settings);
					SetTolerancesStruct(new_tolerances_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetEndgameSettings(new_endgame_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::PowerSeries new_power_series_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetPowerSeriesSettings(new_power_series_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetSecuritySettings(new_security_settings);
				} 

				PowerSeriesEndgame(TrackerType const& tracker, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Endgame PowerSeries Security Tolerances
					SetTolerancesStruct(new_tolerances_settings);
				} 





				~PowerSeriesEndgame() {};

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
				Vec<ComplexType> HermiteInterpolateAndSolve(ComplexType const& target_time, const unsigned int num_sample_points)
			{
				Eigen::IOFormat HeavyFmt(Eigen::FullPrecision);

				Mat< Vec<ComplexType> > finite_difference_matrix(2*num_sample_points,2*num_sample_points);
				Vec<ComplexType> array_of_times(2*num_sample_points);
				
				for(unsigned int ii=0; ii<num_sample_points; ++ii)
				{ 
					finite_difference_matrix(2*ii,0) = samples_[ii];			/*  F[2*i][0]    = samples[i];    */
     				finite_difference_matrix(2*ii+1,0) = samples_[ii]; 		/*  F[2*i+1][0]  = samples[i];    */
      				finite_difference_matrix(2*ii+1,1) = s_derivatives_[ii];	/*  F[2*i+1][1]  = derivatives[i]; */
     				array_of_times(2*ii) = s_times_[ii];						/*  z[2*i]       = times[i];       */
     				array_of_times(2*ii+1) =  s_times_[ii];					/*  z[2*i+1]     = times[i];       */

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

			/*
			Input: 
					sample2, sample1, sample0 are space values. 
				 	derivatives are the dx_dt or dx_ds values at the (time,sample) values.

			Output: An upper bound on the cycle number. 

			Details:
					Using the formula for the cycle test outline in the Bertini Book pg. 53, we can compute an upper bound
					on the cycle number. This upper bound is used for an exhaustive search in ComputeCycleNumber for the actual cycle number. 
			*/

			void BoundOnCycleNumber()
			{ 
				const Vec<mpfr> & sample0 = samples_[0];
				const Vec<mpfr> & sample1 = samples_[1];
				const Vec<mpfr> & sample2 = samples_[2];

				Vec<mpfr> rand_vector = Vec<mpfr>::Random(sample0.size()); //should be a row vector for ease in multiplying.
				

				// //DO NOT USE Eigen .dot() it will do conjugate transpose which is not what we want.
				// //Also, the .transpose*rand_vector returns an expression template that we do .norm of since abs is not available for that expression type. 
				mpfr_float estimate = abs(log(abs((((sample2 - sample1).transpose()*rand_vector).norm())/(((sample1 - sample0).transpose()*rand_vector).norm()))));
				estimate = abs(log(endgame_struct_.sample_factor))/estimate;
				// std::cout << "estimate is " << estimate << '\n';
				if (estimate < 1)
				{
				  	power_series_struct_.upper_bound_on_cycle_number = 1;
				}
				else
				{
					//casting issues between auto and unsigned integer. TODO: Try to stream line this.
					unsigned int upper_bound;
					mpfr_float upper_bound_before_casting = round(floor(estimate + mpfr_float(.5))*power_series_struct_.cycle_number_amplification);
					upper_bound = unsigned (upper_bound_before_casting);
					power_series_struct_.upper_bound_on_cycle_number = max(upper_bound,power_series_struct_.max_cycle_number);
				}
				// Make sure to use Eigen and transpose to use Linear algebra. DO NOT USE Eigen .dot() it will do conjugate transpose which is not what we want. 
				// need to have cycle_number_amplification, this is the 5 used for 5 * estimate
			}


			/*
			Input:	 
			    	current_time is the time value that we wish to interpolate at and compare against the hermite value for each possible cycle number.
			    	x_current_time is the corresponding space value at current_time.
			    	upper_bound_on_cycle_number is the largest possible cycle number we will use to compute dx_ds and s = t^(1/c).

					samples are space values that correspond to the time values in times. 
					derivatives are the dx_dt or dx_ds values at the (time,sample) values.
					num_sample_points is the size of samples, times, derivatives, this is also the amount of information used in Hermite Interpolation.

			Output: 
					The actual cycle number that we find to be the best at approximating the tracked value (current_time,x_current_time)

			Details: 
					This is done by an exhaustive search from 1 to upper_bound_on_cycle_number. There is a conversion to the s-space from t-space in this function. 

			*/
			template<typename ComplexType>
			void ComputeCycleNumber(const ComplexType time, const Vec<ComplexType> x__at_time)
			{

				//Compute upper bound for cycle number.
				BoundOnCycleNumber();
				// std::cout << "upper_bound is " << upper_bound_on_cycle_number << '\n';

				Vec<ComplexType> x_current_time = samples_[samples_.size()-1];
				samples_.pop_back(); //use the last sample to compute the cycle number.
				ComplexType current_time = times_[times_.size()-1];
				times_.pop_back();

				 // std::cout << "x_current_time is " << x_current_time << '\n';
				 // std::cout << "current time is " << current_time << '\n';

				//Compute Cycle Number
				 //num_sample_points - 1 because we are using the most current sample to do an exhaustive search for the best cycle number. 
				SetNumSamples(GetNumSamples() - 1);

				mpfr_float min_found_difference = power_series_struct_.min_difference_in_approximations;
				// std::cout << "upper_bound_on_cycle_number is " << upper_bound_on_cycle_number << '\n';
				// std::cout << "num samples is " << num_sample_points << '\n';
				// std::cout << "min found difference " << min_found_difference << '\n';
				for(unsigned int cc = 1; cc <= power_series_struct_.upper_bound_on_cycle_number; ++cc)
				{
					 // std::cout << "cc is " << cc << '\n';				
					for(unsigned int ii=0; ii<GetNumSamples(); ++ii)// using the last sample to predict to. 
					{ 
						// std::cout << "ii is " << ii << '\n';
						// std::cout << "times[ii] is " << times[ii] << '\n';
						// std::cout << "samples[ii] is " << samples[ii] << '\n';
						s_times_.push_back(pow(times_[ii],ComplexType(1)/ComplexType(cc)));
						s_derivatives_.push_back(derivatives_[ii]*(ComplexType(cc)*pow(times_[ii],ComplexType((cc - 1))/ComplexType(cc))));
					}

					Vec<ComplexType> approx = HermiteInterpolateAndSolve(current_time,GetNumSamples());

					  // std::cout << "computed approx is " << approx << '\n';
					  // std::cout << "norm is " << (approx - x_current_time).norm() << '\n';

					if((approx - x_current_time).norm() < min_found_difference)
					{
						min_found_difference = (approx - x_current_time).norm();
						power_series_struct_.cycle_number = cc;
					}

					s_derivatives_.clear();
					s_times_.clear();

				}// end cc loop over cycle number possibilities

				SetNumSamples(GetNumSamples() + 1);
				// std::cout << "cycle number is " << cycle_number << '\n';
				samples_.push_back(x_current_time);
				times_.push_back(current_time); //push most recent sample back on. 
			}//end ComputeCycleNumber

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
			template<typename ComplexType>
			void ComputeInitialSamples(const ComplexType endgame_time,const Vec<ComplexType> x_endgame) // passed by reference to allow times to be filled as well.
			{

				samples_.push_back(x_endgame);
				times_.push_back(endgame_time);
				Vec<ComplexType> next_sample;

				// std::cout << "x_endgame is " << x_endgame << '\n';
				// std::cout << "endgame_time is " << endgame_time << '\n';

				for(int ii=2; ii <= GetNumSamples(); ++ii)//start at 2 since first sample is at the endgame boundary.
				{ 
					 // std::cout << "ii is " << ii <<  std::endl;
					ComplexType next_time = times_.back() * endgame_struct_.sample_factor;

					SuccessCode tracking_success = endgame_tracker_.TrackPath(next_sample,times_.back(),next_time,samples_.back());

					// std::cout << "next_sample is " << next_sample << '\n';
					// std::cout << "next_time is " << next_time << '\n';
					samples_.push_back(next_sample);
					times_.push_back(next_time);		
				}				
			}
			/*
			Input: 
					time_t0 is the time value that we wish to interpolate at.
					samples are space values that correspond to the time values in times. 
				 	derivatives are the dx_dt or dx_ds values at the (time,sample) values. EMPTY at this point.
				 	sys is the system that we are working with. 
				 	sample_factor is the ratio usually 1/2 that tells us what the next samples should be. 
				 	num_sample_points is the number of time values, sample values, and derivative values that we use for hermite interpolation.


			Output: A new sample point that was found by interpolation to time_t0.

			Details: 

			*/
			template<typename ComplexType>
			Vec<ComplexType> ComputeApproximationOfXAtT0(const ComplexType time_t0)
			{


				//Checking to make sure all samples are of the same precision.
				unsigned max_precision = 0; 

				for(unsigned ii = 0; ii < GetNumSamples();++ii)
				{
					if(Precision(samples_[ii](0)) > max_precision)
					{
						max_precision = Precision(samples_[ii](0));
					}
				}

				for(unsigned ii = 0; ii < GetNumSamples();++ii)
				{
					if(Precision(samples_[ii](0)) < max_precision)
					{
						for(unsigned jj = 0; jj < endgame_tracker_.GetSystem().NumVariables();++jj)
						{
							samples_[ii](jj).precision(max_precision);
						}
						//SuccessCode refine_success = endgame_tracker_.Refine(samples_[ii],samples_[ii],times_[ii]);	
					}
				}

				endgame_tracker_.GetSystem().precision(max_precision);
				//Compute dx_dt for each sample.
				for(unsigned ii = 0; ii < GetNumSamples(); ++ii)
				{	
				// 	std::cout << "ii is " << ii <<'\n';
				// 	std::cout << "Precision of sample is " << Precision(samples_[ii](0)) << '\n';
				// 	std::cout << "Precision of system is " << endgame_tracker_.System().precision() << '\n';
					// uses LU look at Eigen documentation on inverse in Eigen/LU.
				 	Vec<ComplexType> derivative = ComplexType(-1)*(endgame_tracker_.GetSystem().Jacobian(samples_[ii],times_[ii]).inverse())*(endgame_tracker_.GetSystem().TimeDerivative(samples_[ii],times_[ii]));
					derivatives_.push_back(derivative);
				}

				ComputeCycleNumber(times_[times_.size()-1],samples_[samples_.size()-1]); //sending in last element so that ComplexType can be known for templating.


				// //Conversion to S-plane.

				for(unsigned ii = 0; ii < samples_.size(); ++ii){
					s_derivatives_.push_back(derivatives_[ii]*(ComplexType(power_series_struct_.cycle_number)*pow(times_[ii],(ComplexType(power_series_struct_.cycle_number) - ComplexType(1))/ComplexType(power_series_struct_.cycle_number))));
					s_times_.push_back(pow(times_[ii],ComplexType(1)/ComplexType(power_series_struct_.cycle_number)));
				}
				Vec<ComplexType> Approx = HermiteInterpolateAndSolve(time_t0,GetNumSamples());
				s_times_.clear();
				s_derivatives_.clear();
				derivatives_.clear();
				 // std::cout << "approx is " << Approx <<  '\n';
				return Approx;


			}//end ComputeApproximationOfXAtT0


			/*
			Input: 
					endgame_time is the endgame boundary default set to .1
					x_endgame_time is the space value we are given at endgame_time
					sys is the system/homotopy that we currently working with. 


			Output: A space value at time = 0. 

			Details: 
					Using successive hermite interpolations with a geometric progression of time, space, and derivative values we attempt to find the value 
					of the homotopy at time t = 0.
			*/
			template<typename ComplexType>
			SuccessCode PSEG(const ComplexType endgame_time, const Vec<ComplexType> x_endgame_time)
			{
				//Set up for the endgame.
	 			ClearTimesDerivativesAndSamples();
			 	mpfr norm_of_latest_approximation(1);
			 	norm_of_latest_approximation = mpfr("0","0"); // settting up the norm for the latest approximation. 

			 	mpfr approx_error(1);  //setting up the error of successive approximations. 
			 	approx_error = mpfr("1","0");
			 	
			 	ComplexType origin(1);
			 	origin  = ComplexType("0","0");

				ComputeInitialSamples(endgame_time, x_endgame_time);
			 	
			 	Vec<ComplexType> prev_approx = ComputeApproximationOfXAtT0(origin);
			 	Vec<ComplexType> dehom_of_prev_approx = endgame_tracker_.GetSystem().DehomogenizePoint(prev_approx);

			 	ComplexType current_time = times_.back(); 

			  	Vec<ComplexType> next_sample;
			  	Vec<ComplexType> latest_approx;
			    Vec<ComplexType> dehom_of_latest_approx;

				while (approx_error.norm() > Endgame::GetFinalTolerance().abs())
				{
					SetFinalApproximation(prev_approx);
			  		 auto next_time = times_.back() * endgame_struct_.sample_factor; //setting up next time value.

			  		if (next_time.abs() < Endgame::GetMinTrackTime().abs())
			  		{
			  			std::cout << "Error current time norm is less than min track time." << '\n';
			  			SetFinalApproximation(prev_approx);
			  			return SuccessCode::MinTrackTimeReached;
			  		}

					SuccessCode tracking_success = endgame_tracker_.TrackPath(next_sample,current_time,next_time,samples_.back());
					if(tracking_success == SuccessCode::FailedToConverge)
					{
						std::cout << "Failed to converge. " << '\n';
						return SuccessCode::FailedToConverge;
					}	
					else if(tracking_success == SuccessCode::HigherPrecisionNecessary)
					{ 
						std::cout << "Higher precision necessary. " << '\n';
						return SuccessCode::FailedToConverge;
					}
					else if(tracking_success == SuccessCode::GoingToInfinity)
					{
						std::cout << "Going to infinity. " << '\n';
						return SuccessCode:: FailedToConverge;
					}
					else if(tracking_success == SuccessCode::MatrixSolveFailure)
					{
						std::cout << "Matrix solve failure. " << '\n';
						return SuccessCode::MatrixSolveFailure;
					}
					else if(tracking_success == SuccessCode::MaxNumStepsTaken)
					{
						std::cout << "Max num steps taken. " << '\n';
						return SuccessCode::MaxNumStepsTaken;
					}
					else if(tracking_success == SuccessCode::MaxPrecisionReached)
					{
						std::cout << "Max precision reached. " << '\n';
						return SuccessCode::MaxPrecisionReached;
					}
					else if(tracking_success == SuccessCode::MinStepSizeReached)
					{
						std::cout << "Min step size reached." << '\n';
						return SuccessCode::MinStepSizeReached;
					}
					else if(tracking_success == SuccessCode::Failure)
					{
						std::cout << "Failure. " << '\n';
						return SuccessCode::Failure;
					}
					else if(tracking_success == SuccessCode::SingularStartPoint)
					{
						std::cout << "Singular start point. " << '\n';
						// std::cout << "latest norm is " << prev_approx.norm() << '\n';
						return SuccessCode::SingularStartPoint;
					}
					else
					{

					}


					//push most current time into deque, and lose the least recent time.
			 		times_.push_back(next_time);
			 		current_time = times_.back();
			 		times_.pop_front(); 


			 		samples_.push_back(next_sample);
			 		samples_.pop_front();

			 		latest_approx = ComputeApproximationOfXAtT0(origin);
			 		dehom_of_latest_approx = endgame_tracker_.GetSystem().DehomogenizePoint(latest_approx);


			 		if(GetSecurityLevel() <= 0)
					{
				 		if(dehom_of_latest_approx.norm() > GetSecurityMaxNorm().abs() && dehom_of_prev_approx.norm() > GetSecurityMaxNorm().abs()){
			 				return SuccessCode::SecurityMaxNormReached;
				 		}
				 	}

			 		approx_error = (latest_approx - prev_approx).norm();

			 		if(approx_error.norm() < Endgame::GetFinalTolerance().abs())
			 		{
			 			//std::cout << "Power series endgame converged, check final approximation at origin." << '\n';
			 			SetFinalApproximation(latest_approx);
			 			// std::cout << "approx_error norm is " << approx_error.norm() << '\n';
			 			//std::cout << "success" << '\n';

			 			return SuccessCode::Success;
			 		}

			 		prev_approx = latest_approx;
				    dehom_of_prev_approx = dehom_of_latest_approx;
				} //end while	

			// in case if we get out of the for loop without setting. 
			SetFinalApproximation(latest_approx);
			// std::cout << "approx_error norm is " << approx_error.norm() << '\n';
			return SuccessCode::Success;

			} //end PSEG
		}; // end powerseries class


		}//namespace endgame

	}//namespace tracking

} // namespace bertini
#endif
