//This file is part of Bertini 2.0.
//
//cauchy_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//cauchy_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with predict.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  cauchy_endgame.hpp
//
//  copyright 2015
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Spring 2015


#ifndef BERTINI_TRACKING_CAUCHY_HPP
#define BERTINI_TRACKING_CAUCHY_HPP

/**
\file cauchy_endgame.hpp

\brief Contains the cauchy endgame type, and necessary functions.
*/


#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include "tracking/base_endgame.hpp"
#include <cstdio>
#include <typeinfo> 
#include <deque>
#include <algorithm> 



namespace bertini
{ 

	namespace tracking 
	{

		namespace endgame 
		{

		/** 
		\class CauchyEndgame

		\brief Class used to finish tracking paths during Homotopy Continuation.
		
		
		## Explanation
		
		The bertini::CauchyEngame class enables us to finish tracking on possibly singular paths on an arbitrary square homotopy.  

		The intended usage is to:

		1. Create a system, tracker, and instantiate some settings.
		2. Using the tracker created track to the engame boundary. 
		3. Create a CauchyEndgame, associating it to the system you are going to solve or track on.
		4. For each path being tracked send the CauchyEndgame the time value and other variable values at that time. 
		5. The CauchyEndgame, if successful, will store the target systems solution at $t = 0$.

		## Example Usage
		Below we demonstrate a basic usage of the CauchyEndgame class to find the singularity at $t = 0$. 

		The pattern is as described above: create an instance of the class, feeding it the system to be used, and the endgame boundary time and other variable values at the endgame boundary. 

		\code{.cpp}
		mpfr_float::default_precision(30); // set initial precision.  This is not strictly necessary.
		using namespace bertini::tracking;

		// 1. Create the system
		System sys;
		Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t"), y = std::make_shared<Variable>("y");
		VariableGroup vars{x,y};
		sys.AddVariableGroup(vars); 
		sys.AddPathVariable(t);
		// Define homotopy system
		sys.AddFunction((pow(x-1,3))*(1-t) + (pow(x,3) + 1)*t);
		sys.AddFunction((pow(y-1,2))*(1-t) + (pow(y,2) + 1)*t);

		//2. Setup a tracker. 
		bertini::tracking::AMPTracker tracker(sys);
	
		bertini::tracking::config::Stepping stepping_preferences;
		bertini::tracking::config::Newton newton_preferences;

		tracker.Setup(bertini::tracking::config::Predictor::Euler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
		tracker.AMPSetup(AMP);

		// We make the assumption that we have tracked to t = 0.1 (value at t = 0.1 known from Bertini 1.5)
		mpfr current_time(1);
		Vec<mpfr> current_space(2);
		current_time = mpfr(".1");
		current_space <<  mpfr("5.000000000000001e-01", "9.084258952712920e-17") ,mpfr("9.000000000000001e-01","4.358898943540673e-01");

		//  3. (Optional) configure settings for specific endgame. 
		bertini::tracking::config::Endgame endgame_settings;
		bertini::tracking::config::Cauchy cauchy_settings;
		bertini::tracking::config::Security endgame_security_settings;

		//4. Create the CauchyEngame object, by sending in the tracker and any other settings. Notice the endgame is templated by the tracker. 
		bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,cauchy_settings,endgame_settings,endgame_security_settings);

		//Calling the CauchyEG member function actually runs the endgame. 
		My_Endgame.CauchyEG(current_time,current_space);

		//Access solution at t = 0, using the .get_final_approximation_at_origin() member function. 
		auto Answer  = My_Endgame.get_final_approximation_at_origin();
		\endcode
		
		If this documentation is insufficient, please contact the authors with suggestions, or get involved!  Pull requests welcomed.
		
		## Testing
		Test suite driving this class: endgames_test.

		File: test/endgames/cauchy_class_test.cpp

		Functionality tested: All member functions of the CauchyEndgame have been tested in variable precision. There is also, test running different systems to find singular solutions at t = 0.


		*/



		
			template<typename TrackerType> 
			class CauchyEndgame : public BaseEndgame<TrackerType>
			{
			private:

				/**
				\brief Settings that are specific to the Cauchy endgame. 
				*/	
				config::Cauchy cauchy_settings_; 

				/**
				\brief A deque of times that are specifically used to compute the power series approximation for the Cauchy endgame. 
				*/
				std::deque<mpfr> pseg_times_; 
				/**
				\brief A deque of samples that are in correspondence with the pseg_times_. These samples are also used to compute the first power series approximation for the Cauchy endgame.
				*/
				std::deque< Vec<mpfr> > pseg_samples_; //samples used for the first approximation. 
				/**
				\brief A deque of times that are collected by CircleTrack. These samples are used to compute all approximations of the origin after the first. 
				*/
				std::deque<mpfr> cauchy_times_; 
				/**
				\brief A deque of samples collected by CircleTrack. Computed a mean of the values of this deque, after a loop has been closed, will give an approximation of the origin.
				*/
				std::deque< Vec<mpfr> > cauchy_samples_; 
				/**
				\brief A tolerance that is computed from settings passed into the Cauchy endgame class. 
				*/
				mpfr_float min_closed_loop_tolerance_;
				/**
				\brief A tolerance that is computed from settings passed into the Cauchy endgame class. 
				*/
				mpfr_float max_closed_loop_tolerance_;

				


			public:

				// some getters
				const std::deque< Vec<mpfr> >& CauchySamples() const {return cauchy_samples_;}
				const std::deque<mpfr>& CauchyTimes() const {return cauchy_times_;}

				/**
				\brief Function that clears all samples and times from data members for the Cauchy endgame
				*/
				void ClearTimesAndSamples(){pseg_times_.clear(); pseg_samples_.clear(); cauchy_times_.clear(); cauchy_samples_.clear();}
				/**
				\brief Setter for the deque holding time values for the power series approximation of the Cauchy endgame. 
				*/
				void SetPSEGTimes(std::deque<mpfr> pseg_times_to_set) { pseg_times_ = pseg_times_to_set;}
				/**
				\brief Getter for the deque holding time values for the power series approximation of the Cauchy endgame.
				*/
				std::deque< mpfr > GetPSEGTimes() {return pseg_times_;}

				/**
				\brief Setter for the deque holding space values for the power series approximation of the Cauchy endgame. 
				*/
				void SetPSEGSamples(std::deque< Vec<mpfr> > pseg_samples_to_set) { pseg_samples_ = pseg_samples_to_set;}
				/**
				\brief Getter for the deque holding space values for the power series approximation of the Cauchy endgame. 
				*/
				std::deque< Vec<mpfr> > GetPSEGSamples() {return pseg_samples_;}

				/**
				\brief Setter for the deque holding space values for the Cauchy endgame. 
				*/
				void SetCauchySamples(std::deque< Vec<mpfr> > cauchy_samples_to_set) { cauchy_samples_ = cauchy_samples_to_set;}

				/**
				\brief Getter for the deque holding space values for the Cauchy endgame. 
				*/
				std::deque< Vec<mpfr> > GetCauchySamples() {return cauchy_samples_;}

				/**
				\brief Setter for the deque holding time values for the Cauchy endgame. 
				*/
				void SetCauchyTimes(std::deque<mpfr> cauchy_times_to_set) { cauchy_times_ = cauchy_times_to_set;}

				/**
				\brief Getter for the deque holding time values for the Cauchy endgame. 
				*/
				std::deque< mpfr > GetCauchyTimes() {return cauchy_times_;}
	
				
				/**
				\brief Setter for the specific settings in tracking_conifg.hpp under Cauchy.
				*/
				void SetCauchySettings(config::Cauchy new_cauchy_settings){cauchy_settings_ = new_cauchy_settings;}

				/**
				\brief Getter for the specific settings in tracking_conifg.hpp under Cauchy.
				*/
				config::Cauchy GetCauchySettings(){return cauchy_settings_;}
				

				explicit CauchyEndgame(TrackerType const& tr, 
				                            const std::tuple< config::Cauchy const&, 
				                            					const config::Endgame&, 
				                            					const config::Security&, 
				                            					const config::Tolerances& >& settings )
			      : BaseEndgame<TrackerType>(tr, std::get<1>(settings), std::get<2>(settings), std::get<3>(settings) ), 
			          cauchy_settings_( std::get<0>(settings) )
			   	{
			   		min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),this->Tolerances().final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(this->Tolerances().newton_during_endgame,this->Tolerances().final_tolerance_times_final_tolerance_multiplier);
			   	}

			    template< typename... Ts >
   				CauchyEndgame(TrackerType const& tr, const Ts&... ts ) : CauchyEndgame(tr, Unpermute<config::Cauchy, config::Endgame, config::Security, config::Tolerances >( ts... ) ) 
   				{}


				~CauchyEndgame() {};


				/*
				Input: 
					A starting time and starting sample space to start tracking. 

				Output: A successcode to see if we have successfully tracked around the origin or if we have encountered an error. 

				Details:
					This function uses the number of sample points to track in a polygonal path around the origin. 
					This function should be called the number of times it takes to loop around the origin and get back the original 
					point we started with. In essence this function will help determine the cycle number. A deque of cauchy_samples and cauchy_times is populated in this function. 
				*/
				/**
				\brief A function that will track around the origin once. 

				\param[out] starting_time The time value at which we start trackin around the origin. 
				\param starting_sample The current variable values at staring_time.

				\tparam ComplexType The complex number type.
				*/

				template<typename ComplexType> 
				SuccessCode CircleTrack(ComplexType const& starting_time, Vec<ComplexType> const& starting_sample)
				{	std::cout << "CircleTrack()\n";

					if (this->EndgameSettings().num_sample_points < 3) // need to make sure we won't track right through the origin.
					{
						std::stringstream err_msg;
						err_msg << "ERROR: The number of sample points " << this->EndgameSettings().num_sample_points << " for circle tracking must be >= 3";
						throw std::runtime_error(err_msg.str());
					}	

					// the initial sample has already been added to the sample repo... so don't do that here, please
					
					mpfr_float radius = abs(starting_time), angle = arg(starting_time);
					mpfr_float angle_increment = 2*acos(mpfr_float(-1)) / (this->EndgameSettings().num_sample_points);

					auto number_of_steps = 0;
					ComplexType i(0,1); 

					Vec<ComplexType> current_sample = starting_sample;
					Vec<ComplexType> next_sample;
					ComplexType current_time = starting_time;
					ComplexType next_time;

					

					
 	

					for (unsigned ii = 0; ii < this->EndgameSettings().num_sample_points; ++ii)
					{
						//setting up the time value for the next sample. 
						if (ii!=this->EndgameSettings().num_sample_points-1)
						{
							angle += angle_increment;
							next_time = polar(radius, angle);
						}
						else
							next_time = starting_time;

						std::cout << "circle tracking from " << current_time << " to " << next_time << ", starting sample:\n" << current_sample << "\n\n\n";
						auto tracking_success = this->GetTracker().TrackPath(next_sample, current_time, next_time, current_sample);	
						std::cout << "tracked end point:\n" << next_sample << "\n\n";
						if (tracking_success != SuccessCode::Success)
						{
							std::cout << "tracker fail in circle track, type " << int(tracking_success) << std::endl;
							return tracking_success;
						}
						std::cout << "refining point to " << this->Tolerances().final_tolerance/100 << '\n';
						this->GetTracker().Refine(next_sample,next_sample,next_time,
						                          this->Tolerances().final_tolerance/100,
						                          this->EndgameSettings().max_num_newton_iterations);
						std::cout << "refined point:\n" << next_sample << "\n\n";
						

						current_sample = next_sample;
						current_time = next_time;
						cauchy_times_.push_back(current_time);
						cauchy_samples_.push_back(current_sample);
					}

					return SuccessCode::Success;

				}//end CircleTrack

				/*
				Input: 
					All input is available through class data members. 

				Output: An mpfr_float that is an estimate on the c/k ratio shown in the book under Cauchy Endgame.

				Details:
					Using the formula for the cycle test outline in the Bertini Book pg. 53, we can compute an estimate for c/k.
					This estimate is used to help find stabilization of the cycle number. 
				*/

				/**
				\brief A function that uses the assumption of being in the endgame operating zone to compute an approximation of the ratio c over k. 
				When the cycle number stabilizes we will see that the different approximations of c over k will stabilize. 
				Returns the computed value of c over k. 
				*/
				
				mpfr_float ComputeCOverK()
				{//Obtain samples for computing C over K.
					assert(pseg_samples_.size()>=3);
					const Vec<mpfr> & sample0 = pseg_samples_[0];
					const Vec<mpfr> & sample1 = pseg_samples_[1];
					const Vec<mpfr> & sample2 = pseg_samples_[2];

					Vec<mpfr> rand_vector = Vec<mpfr>::Random(sample0.size()); //should be a row vector for ease in multiplying.


					// //DO NOT USE Eigen .dot() it will do conjugate transpose which is not what we want.
					// //Also, the .transpose*rand_vector returns an expression template that we do .norm of since abs is not available for that expression type. 
					mpfr_float estimate = abs(log(abs((((sample2 - sample1).transpose()*rand_vector).norm())/(((sample1 - sample0).transpose()*rand_vector).norm()))));
					estimate = abs(log(this->EndgameSettings().sample_factor))/estimate;
					if (estimate < 1)
					  	return mpfr_float(1);
					else
						return estimate;

				}//end ComputeCOverK

				/*
				Input: An array of c/k estimates. If we have stabilization than these estimates will be withing some treshold. 

				Output: Boolean of true or false. 

				Details: We consider the ratio of consecutive c/k estimates. If we are less than the minimum needed for stabilization we return false.
				Otherwise, we return true. 
				*/

				/**
				\brief Function to determine if ratios of c/k estimates are withing a user defined threshold. 			

				\param c_over_k_array A deque of computed values for c/k.
				*/

				bool CheckForCOverKStabilization(std::deque<mpfr_float> const& c_over_k_array)
				{	std::cout << "CheckForCOverKStabilization()\n";
					assert(c_over_k_array.size()>=cauchy_settings_.num_needed_for_stabilization);
					for(unsigned ii = 1; ii < cauchy_settings_.num_needed_for_stabilization ; ++ii)
					{
						auto a = abs(c_over_k_array[ii-1]);
						auto b = abs(c_over_k_array[ii]);

						auto divide = a;

						if(a < b)
							divide = a/b;
						else
							divide = b/a;

						if(divide <  cauchy_settings_.minimum_for_c_over_k_stabilization)
							return false;
					}
					return true;

				}//end CheckForCOverKStabilization


				/*
				Input: A time value and the space value above that time.
					
				Output: An mpfr_float representing a tolerance threshold for declaring a loop to be closed. 


				Details: Used in Bertini 1 as a heuristic for computing separatedness of roots. Decided to not be used since assumptions for this tolerance are not usually met. 


				template<typename ComplexType>
				mpfr_float FindToleranceForClosedLoop(ComplexType x_time, Vec<ComplexType> x_sample)
				{
					auto degree_max = std::max(this->GetTracker().AMP_config_.degree_bound,mpfr_float("2.0")); 
					auto K = this->GetTracker().AMP_config_.coefficient_bound;
					mpfr_float N;
					mpfr_float M;
					mpfr_float L;

					if(max_closed_loop_tolerance_ < min_closed_loop_tolerance_)
					{
						max_closed_loop_tolerance_ = min_closed_loop_tolerance_;
					}

					auto error_tolerance = mpfr_float("1e-13");
					if(x_sample.size() <= 1)
					{
						N = degree_max;
					}
					else
					{
						N = ComputeCombination(degree_max + x_sample[0].precision() - 1, x_sample[0].precision() - 1);
					}
					M = degree_max * (degree_max - 1) * N;
					auto jacobian_at_current_time = this->GetSystem().Jacobian(x_sample,x_time);


					auto minimum_singular_value = Eigen::JacobiSVD< Mat<ComplexType> >(jacobian_at_current_time).singularValues()(this->GetSystem().NumVariables() - 1 );

					auto norm_of_sample = x_sample.norm();
					L = pow(norm_of_sample,degree_max - 2);
					auto tol = K * L * M;
					// std::cout << "TOL IS " << tol << '\n';
					if (tol == 0) // fail-safe
   							tol = minimum_singular_value;
 					else
   					{
   					tol = mpfr_float("2.0") / tol;
    				tol = tol * minimum_singular_value;
  					}
 					// make sure that tol is between min_tol & max_tol
					if (tol > max_closed_loop_tolerance_) // tol needs to be <= max_tol
 				    tol = max_closed_loop_tolerance_;
					if (tol < min_closed_loop_tolerance_) // tol needs to be >= min_tol
  					tol = min_closed_loop_tolerance_;
 					return tol;
				}// end FindToleranceForClosedLoop
				*/


				/*
				Input: No input, tolerance used is tracking tolerance during endgame. 
					
				Output: A boolean declaring if we have closed the loop or not. 

				Details: Take a point used to start CircleTrack along with last point tracked to in CircleTrack. We see if the difference between the two points
				is less than the tracking tolerance during endgame. If so, we declare these points to be the same. 	
				*/

				/**
				\brief Function that determines if we have closed a loop after calling CircleTrack().
				*/
				bool CheckClosedLoop()
				{	std::cout << "CheckClosedLoop()\n"; std::cout << "have " << cauchy_samples_.size() << " samples to work with\n";
					if((cauchy_samples_.front() - cauchy_samples_.back()).norm() < this->GetTracker().TrackingTolerance())
					{
						std::cout << "loop is closed, no refinement needed\n";
						return true;
					}
					std::cout << "closed loop check: pre-refinement " << (cauchy_samples_.front() - cauchy_samples_.back()).norm() << ", required tol " << this->GetTracker().TrackingTolerance() << std::endl;
					this->GetTracker().Refine(cauchy_samples_.front(),cauchy_samples_.front(),cauchy_times_.front(),this->Tolerances().final_tolerance,this->EndgameSettings().max_num_newton_iterations);
					this->GetTracker().Refine(cauchy_samples_.back(),cauchy_samples_.back(),cauchy_times_.back(),this->Tolerances().final_tolerance,this->EndgameSettings().max_num_newton_iterations);
					std::cout << "closed loop check: post-refinement " << (cauchy_samples_.front() - cauchy_samples_.back()).norm() << std::endl;
					if((cauchy_samples_.front() - cauchy_samples_.back()).norm() < this->GetTracker().TrackingTolerance())
					{
						std::cout << "loop is closed, after refinement\n";
						return true;
					}
					std::cout << "loop not closed\n";
					return false;	

				}//end CheckClosedLoop

				/*
				Input: No input, all data needed are class data members. 

				Output: Boolean declaring if our ratios are close enough or not. 

				Details: Finds minimum and maximum norms from tracking around the origin. Then we check against the user defined setting minimum_c_over_k_for_stabilization to see if 
				we have good ratios or not. 
					 
				*/

				/**
				\brief After we have used CircleTrack and have successfully closed the loop using CheckClosedLoop we need to check the maximum and minimum norms of the samples collected. 

				If the ratio of the maximum and minimum norm are within the threshold maximum_cauchy_ratio, and the difference is greater than final tolerance than we are successful. 
				*/
				bool RatioEGOperatingZoneTest()
				{	std::cout << "RatioEGOperatingZoneTest()\n";
					mpfr_float min(1e300);
					mpfr_float max(0);

					if(cauchy_times_.front().norm() < cauchy_settings_.ratio_cutoff_time)
					{
						return true;
					}
					else
					{
						mpfr_float norm;
						for(unsigned int ii=0; ii < this->EndgameSettings().num_sample_points; ++ii)
						{
							norm = cauchy_samples_[ii].norm();
							if(norm > max)
							{
								max = norm;
							}
							if(norm < min)
							{
								min = norm;
							}
						}

						if(min > this->Tolerances().final_tolerance && max > this->Tolerances().final_tolerance)
						{
							norm = min / max;
							if(norm < cauchy_settings_.maximum_cauchy_ratio && (max - min) > this->Tolerances().final_tolerance)
							{
								return false;
							}
						}
					}
					return true;

				}//end RatioEGOperatingZoneTest


				/*
				Input: All input needed is available as class data members. 

				Output: SuccessCode declaring if we are successful or not. 

				Details: Computes cauchy samples using Circle Track. Compares the ratios of the maximum and minimum norms of cauchy samples 
						 using RatioEGOperatingZoneTest. Then attempts to see if we have closed our loop, otherwise it continues the process outlined.
				*/

				/**
				\brief This function tracks into origin while computing loops around the origin. The function is checking to make sure we have reached a time the ratio of the 
				maximum and minimum norms are withing some tolerance. When this is the case we return success. If this does not happen we will return an error depending on the error 
				encountered. 
				*/		 
				SuccessCode InitialCauchyLoops()
				{	std::cout << "InitialCauchyLoops()\n";
					using std::max;
					bool continue_loop = true;

					auto fail_safe_max_cycle_number = max(cauchy_settings_.fail_safe_maximum_cycle_number,this->CycleNumber());

					auto return_value = SuccessCode::Success;

					while (continue_loop)
					{	
						cauchy_times_.clear();
						cauchy_samples_.clear();
						cauchy_samples_.push_back(pseg_samples_.back()); // cauchy samples and times should be empty before this point. 
						cauchy_times_.push_back(pseg_times_.back());

						auto next_time = pseg_times_.back();
						auto next_sample = pseg_samples_.back();

						// track around a circle once.  we'll use it to measure whether we believe we are in the eg operating zone, based on the ratio of ratios of norms of sample points around the circle
						auto tracking_success = CircleTrack(cauchy_times_.front(),cauchy_samples_.front());

						this->IncrementCycleNumber(1);

						if (tracking_success != SuccessCode::Success)
							return tracking_success;

						// find the ratio of the maximum and minimum coordinate wise for the loop. 
						if (RatioEGOperatingZoneTest())
						{ // then we believe we are in the EG operating zone, since the path is relatively flat.  i still disbelieve this is a good test (dab 20160310)
							while (true)
							{
								if (CheckClosedLoop())
								{//error is small enough, exit the loop with success. 
									return_value = SuccessCode::Success;
									continue_loop = false;
									break;
								}
								else if(this->CycleNumber() > fail_safe_max_cycle_number)
								{// too many iterations
									return_value = SuccessCode::CycleNumTooHigh;
									continue_loop = false;
									break;
								}

								//compute next loop, the last sample in times and samples is the sample our loop ended on. Either where we started or on another sheet at the same time value. 
								tracking_success = CircleTrack(cauchy_times_.back(),cauchy_samples_.back());

								this->IncrementCycleNumber(1);

								if(tracking_success != SuccessCode::Success)
									return tracking_success;
							}

							if(return_value == SuccessCode::CycleNumTooHigh)
							{//see if we should continue to the next sample point
								if(cauchy_times_.back().abs() < this->EndgameSettings().min_track_time)
									continue_loop = false;
								else
									continue_loop = true;
							}
						}//end if (RatioEGOperatingZoneTest())
						else 
						{
							std::cout << "not in EG zone yet, shrinking radius\n\n";
							//find the time for the next sample point
							next_time = pseg_times_.back() * this->EndgameSettings().sample_factor;

							SuccessCode tracking_success = this->GetTracker().TrackPath(next_sample,pseg_times_.back(),next_time,pseg_samples_.back());
							pseg_times_.pop_front();
							pseg_samples_.pop_front();

							pseg_times_.push_back(next_time);
							pseg_samples_.push_back(next_sample);

							if(tracking_success != SuccessCode::Success)
								return tracking_success;
						}
					} //end while(continue_loop)

					return return_value;
				}//end InitialCauchyLoops


				/*
				Input: Time and sample at the endgame boundary (usually 0.1). Along with address of time and sample for the approximation computed. 

				Output: SuccessCode deeming if we were successful. 

				Details: This function is in charge of finding the very first approximation of the origin. It does this by first computing some initial samples 
						 like what is done in the Power Series Endgame. We continue to track forward in this manner until we have stabilization of the cycle number being approximated. 
						 This prevents the unnecessary circle tracking if we are possibly not in the endgame operating zone. 

						 Once we have stabilization we then perform InitialCauchyLoops while get the accurate cycle number, and check the norms of the samples and make sure we are ready 
						 to approximate. When ready we call ComputePSEGApproximationOfXAtT0. This function will use a hermtie interpolater to get an approximation of the value at the origin. 
				*/
				/**
				\brief The Cauchy endgame will first find an initial approximation using the notion of the power series endgame. This function computes this approximation and returns a
				SuccessCode to let us know if an error was encountered. 

				\param[out] start_time The time value at which we start the endgame. 
				\param start_point The current space point at start_time.
				\param approximation_time the time at which we want to compute an approximation, usually the origin. 
				\param approximation The space value at approximation_time that we are computing. 

				\tparam ComplexType The complex number type.
				*/

				template<typename ComplexType>
				SuccessCode InitialPowerSeriesApproximation(ComplexType const& start_time, Vec<ComplexType> const& start_point, ComplexType & approximation_time, Vec<ComplexType> & approximation)
				{	std::cout << "InitialPowerSeriesApproximation()\n\n";
					//initialize array holding c_over_k estimates
					std::deque<mpfr_float> c_over_k; 

					//Compute initial samples for pseg
					auto initial_sample_success = this->ComputeInitialSamples(start_time, start_point, pseg_times_, pseg_samples_);
					if (initial_sample_success!=SuccessCode::Success)
						return initial_sample_success;

					c_over_k.push_back(ComputeCOverK());

					Vec<ComplexType> next_sample;
					ComplexType next_time;

					unsigned ii = 0;
					//track until for more c_over_k estimates or until we reach a cutoff time. 
					while ( (ii < cauchy_settings_.num_needed_for_stabilization) )
					{	std::cout << "getting next c/k\n";
						next_time = pseg_times_.back() * this->EndgameSettings().sample_factor;

						auto tracking_success = this->GetTracker().TrackPath(next_sample,pseg_times_.back(),next_time,pseg_samples_.back());
						if (tracking_success!=SuccessCode::Success)
							return tracking_success;

						pseg_samples_.pop_front();
						pseg_times_.pop_front();
						pseg_samples_.push_back(next_sample);
						pseg_times_.push_back(next_time);
						c_over_k.push_back(ComputeCOverK());

						++ii;
					}//end while

					std::cout << "have initial c_over_k:\n";
					for (auto iter : c_over_k)
						std::cout << iter << std::endl;
					//check to see if we continue. 

					//have we stabilized yet? 

					while(!CheckForCOverKStabilization(c_over_k) && pseg_times_.back().abs() > cauchy_settings_.cycle_cutoff_time)
					{
						next_time = pseg_times_.back() * this->EndgameSettings().sample_factor;

						auto tracking_success = this->GetTracker().TrackPath(next_sample,pseg_times_.back(),next_time,pseg_samples_.back());
						if(tracking_success != SuccessCode::Success)
							return tracking_success;

						c_over_k.pop_front();
						pseg_samples_.pop_front();
						pseg_times_.pop_front();

						pseg_samples_.push_back(next_sample);	
						pseg_times_.push_back(next_time);
						c_over_k.push_back(ComputeCOverK());

					}//end while
					std::cout << "have c/k stabilization\n";
	
					auto cauchy_loop_success = InitialCauchyLoops();
					if(cauchy_loop_success != SuccessCode::Success)
						return cauchy_loop_success;
					std::cout << "have initial cauchy loops\n";
					approximation = ComputePSEGApproximationAtT0(approximation_time);
					std::cout << "have first approximation\n";
					return SuccessCode::Success;

				}//end InitialPowerSeriesApproximation

				/*
				Input: time_t0 is the time value that we wish to interpolate at.


				Output: A new sample point that was found by interpolation to time_t0.

				Details: This function handles computing an approximation at the origin. 
					We compute the derivatives at the different times and samples. We then make sure all samples are to the same precision before refining them to final tolerance. 
					By InitialCauchyLoops we know what the cycle number is sow we convert derivatives and times to the s-plane where s = t^(1/(cyle number).
					We use the converted times and derivatives along with the samples to do a Hermite interpolation which is found in base_endgame.hpp.
				*/

				/**
				\brief This function takes the pseg_samples_ and pseg_times that have been collected and uses them to compute a Hermite interpolation to the time value time_t0. 


				\param[out] time_t0 is the time value that the power series approximation will interpolate to. 

				\tparam ComplexType The complex number type.
				*/
				template<typename ComplexType>
				Vec<ComplexType> ComputePSEGApproximationAtT0(const ComplexType time_t0)
				{
					//Compute dx_dt for each sample.
					std::deque< Vec<ComplexType> > pseg_derivatives;
					//initialize the derivative otherwise the first computation in the loop will be wrong. 

					Vec<ComplexType> pseg_derivative = -(this->GetSystem().Jacobian(pseg_samples_[0],pseg_times_[0]).inverse())*this->GetSystem().TimeDerivative(pseg_samples_[0],pseg_times_[0]);
					for(unsigned ii = 0; ii < this->EndgameSettings().num_sample_points; ++ii)
					{	
						// uses LU look at Eigen documentation on inverse in Eigen/LU.
					 	Vec<ComplexType> pseg_derivative = -(this->GetSystem().Jacobian(pseg_samples_[ii],pseg_times_[ii]).inverse())*this->GetSystem().TimeDerivative(pseg_samples_[ii],pseg_times_[ii]);
						pseg_derivatives.push_back(pseg_derivative);
					}

					//Checking to make sure all samples are of the same precision.
					unsigned max_precision = 0; 

					for(unsigned ii = 0; ii < this->EndgameSettings().num_sample_points;++ii)
						if(Precision(pseg_samples_[ii](0)) > max_precision)
							max_precision = Precision(pseg_samples_[ii](0));


					for(unsigned ii = 0; ii < this->EndgameSettings().num_sample_points;++ii)
						if(Precision(pseg_samples_[ii](0)) < max_precision)
							for(unsigned jj = 0; jj < this->GetSystem().NumVariables();++jj)
								pseg_samples_[ii](jj).precision(max_precision);


					this->GetSystem().precision(max_precision);

			 		//Conversion to S-plane.
					std::deque<ComplexType> s_times;
					std::deque< Vec<ComplexType> > s_derivatives;



					for(unsigned ii = 0; ii < pseg_samples_.size(); ++ii)
					{
						s_derivatives.push_back(pseg_derivatives[ii]*(ComplexType(this->CycleNumber())*pow(pseg_times_[ii],(ComplexType(this->CycleNumber()) - ComplexType(1))/this->CycleNumber())));
						s_times.push_back(pow(pseg_times_[ii],ComplexType(1)/ComplexType(this->CycleNumber())));

						this->GetTracker().Refine(cauchy_samples_[ii],cauchy_samples_[ii],cauchy_times_[ii],
						                          this->Tolerances().final_tolerance/100,
						                          this->EndgameSettings().max_num_newton_iterations);
					}

					return bertini::tracking::endgame::HermiteInterpolateAndSolve(time_t0, this->EndgameSettings().num_sample_points, s_times, pseg_samples_, s_derivatives);
				}//end ComputePSEGApproximationOfXAtT0

				/*
				Input: A deque of the cauchy samples needed to compute an approximation using the Cauchy Integral Formula. 

				Output: The approximation computed by using the samples collected from CircleTrack to sum and divide by the total number used to track back to the starting position. 

				Details: We can compute the Cauchy Integral Formula in this particular instance by computing the mean of the samples we have collected around the origin. 
				*/

				/**
				\brief Function that computes the mean of the samples that were collected while tracking around the origin. This value is the approximation of the value at the origin. 

				\param[out] Deque of cauchy_samples, used to help template the function with ComplexType.

				\tparam ComplexType The complex number type.
				*/
				template<typename ComplexType>
				Vec<ComplexType> ComputeCauchyApproximationOfXAtT0()
				{
					if (cauchy_samples_.size() != this->CycleNumber() * this->EndgameSettings().num_sample_points+1)
					{
						std::stringstream err_msg;
						err_msg << "to compute cauchy approximation, cauchy_samples must be of size " << this->CycleNumber() * this->EndgameSettings().num_sample_points+1 << " but is of size " << cauchy_samples_.size() << '\n';
						throw std::runtime_error(err_msg.str());
					}

					Vec<ComplexType> approximate_root = Vec<ComplexType>::Zero(this->GetSystem().NumVariables());//= (cauchy_samples_[0]+cauchy_samples_.back())/2; 

					for(unsigned int ii = 0; ii < this->CycleNumber() * this->EndgameSettings().num_sample_points; ++ii)
					{
						this->GetTracker().Refine(cauchy_samples_[ii],cauchy_samples_[ii],cauchy_times_[ii],
						                          this->Tolerances().final_tolerance/100,
						                          this->EndgameSettings().max_num_newton_iterations);
						approximate_root += cauchy_samples_[ii];
					}
					approximate_root /= (this->CycleNumber() * this->EndgameSettings().num_sample_points);
					return approximate_root;

				}


				/*
				Input: A time and sample value to start at for circle tracking around the origin. 

				Output: A SuccessCode deeming if we were successful or not. 

				Details: This function populates the deque cauchy_samples and cauchy_times. These are data members of the class and are not passed in. This function will continue to 
				call CircleTrack until we have closed the loop. 
				*/

				/**
				\brief Function that will utilize CircleTrack and CheckClosedLoop to collect all samples while tracking around the origin till we close the loop. 


				\param[out] starting_time The time value at which we start finding cauchy samples 
				\param starting_sample The current space point at starting_time, this will also be the first cauchy sample. 

				\tparam ComplexType The complex number type.
				*/
				template<typename ComplexType>
				SuccessCode ComputeCauchySamples(ComplexType const& starting_time, Vec<ComplexType> const& starting_sample)
				{
					std::cout << "ComputeCauchySamples()\n";
					cauchy_times_.clear();
					cauchy_samples_.clear();
					cauchy_times_.push_back(starting_time);
					cauchy_samples_.push_back(starting_sample);
					this->CycleNumber(0);


					while( this->CycleNumber() < cauchy_settings_.fail_safe_maximum_cycle_number )
					{
						std::cout << "ComputeCauchySamples() iteration " << this->CycleNumber() << "\n";
						//track around the origin once.
						auto tracking_success = CircleTrack(cauchy_times_.back(),cauchy_samples_.back());
						this->IncrementCycleNumber(1);

						std::cout << "loop ended at time \n" << cauchy_times_.back() << " and point\n" << cauchy_samples_.back() << "\n\n";
						if(tracking_success != SuccessCode::Success)
						{
							std::cout << "Cauchy loop fail tracking\n\n";
							return tracking_success;
						}
						else if(CheckClosedLoop())
						{
							std::cout << "Cauchy loop is CLOSED\n\n";
							return SuccessCode::Success;
						}
					} 

					return SuccessCode::CycleNumTooHigh;
				}//end ComputeCauchySamples
				

				/*
				Input: A time and sample to start the cauchy endgame.

				Output: SuccessCode deeming if we were successful or not.

				Details: This function runs the entire Cauchy Endgame. We first take our endgame boundary time value and sample to find a first approximation of the origin. This is done by
				using the idea for the power series endgame. We check for stabilization of the cycle number, and check to see when the ratios of the maximum and minimum norm of samples collected
				by CircleTrack are withing a tolerance. When both of these conditions are met we do a Hermite interpolation. 

				At this point we can start tracking in to the origin while using CircleTrack to compute samples and calculating their mean to get an approximation of the origin using the Cauchy
				Integral Formula. 
				*/

				/**
				\brief Function that runs the Cauchy endgame.

				To begin, this function will compute a first approximation using the power series endgame notion. This approximation is made after a heuristic on the stabilization of 
				the cyle number is made, and after the maximum and minimum norms of tracked space values around the origin are withing a certain tolerance. 			

				\param[out] start_time The time value at which we start the endgame. 
				\param start_point The current space point at start_time.

				\tparam ComplexType The complex number type.
				*/
				template<typename ComplexType>
				SuccessCode CauchyEG(ComplexType start_time, Vec<ComplexType> start_point)
				{	std::cout << "CauchyEG()\n";
					ClearTimesAndSamples(); //clear times and samples before we begin.
					this->CycleNumber(0);
					ComplexType origin(0,0);

					ComplexType next_time;
					Vec<ComplexType> next_sample;

					Vec<ComplexType> prev_approx, latest_approx;
					mpfr_float approximate_error;

				

					//Compute the first approximation using the power series approximation technique. 
					auto initial_ps_success = InitialPowerSeriesApproximation(start_time, start_point, origin, prev_approx);  // last argument is output here
					if(initial_ps_success != SuccessCode::Success)
						return initial_ps_success;

					next_time = pseg_times_.back();
					mpfr_float norm_of_dehom_of_prev_approx, norm_of_dehom_of_latest_approx;

					if(this->SecuritySettings().level <= 0)
						norm_of_dehom_of_prev_approx = this->GetSystem().DehomogenizePoint(prev_approx).norm();
					std::cout << "entering do-while in CauchyEG()\n";
					do
					{
						//Compute a cauchy approximation.
						latest_approx = ComputeCauchyApproximationOfXAtT0<ComplexType>();
						std::cout << "latest_approx\n" << latest_approx << "\n\n";
						if(this->SecuritySettings().level <= 0)
							norm_of_dehom_of_latest_approx = this->GetSystem().DehomogenizePoint(latest_approx).norm();

						approximate_error = (latest_approx - prev_approx).norm(); //Calculate the error between approximations. 
						std::cout << "approximate_error " << approximate_error << "\n\n";
						// dehom of prev approx and last approx not used because they are not updated with the most current information. However, prev approx and last approx are 
						// the most current. 

						if(approximate_error < this->Tolerances().final_tolerance)
						{
							this->final_approximation_at_origin_ = latest_approx;
							return SuccessCode::Success;
						}
						else if(cauchy_times_.front().abs() < this->EndgameSettings().min_track_time)
						{//we are too close to t = 0 but we do not have the correct tolerance - so we exit
							this->final_approximation_at_origin_ = latest_approx;
							return SuccessCode::FailedToConverge;
						}
						else if(this->SecuritySettings().level <= 0 && 
						   norm_of_dehom_of_prev_approx   > this->SecuritySettings().max_norm &&  
						   norm_of_dehom_of_latest_approx > this->SecuritySettings().max_norm  )
						{//we are too large, break out of loop to return error.
							this->final_approximation_at_origin_ = latest_approx;
							return SuccessCode::SecurityMaxNormReached;
						}


						prev_approx = latest_approx;
						norm_of_dehom_of_prev_approx = norm_of_dehom_of_latest_approx;

						next_time *= this->EndgameSettings().sample_factor;

						auto tracking_success = this->GetTracker().TrackPath(next_sample,cauchy_times_.back(),next_time,cauchy_samples_.front());
						if (tracking_success != SuccessCode::Success)
							return tracking_success;

						pseg_times_.push_back(next_time);  pseg_times_.pop_front();
						pseg_samples_.push_back(next_sample); pseg_samples_.pop_front();

						auto cauchy_samples_success = ComputeCauchySamples(next_time,next_sample);

						//Added because of Griewank osborne test case where tracker returns GoingToInfinity.  
						//The cauchy endgame should stop instead of attempting to continue. 
						//This is because the tracker will not track to any meaningful space values. 
						if (cauchy_samples_success == SuccessCode::GoingToInfinity)
							return SuccessCode::GoingToInfinity;
						else if (cauchy_samples_success != SuccessCode::Success && cauchy_times_.front().abs() < this->EndgameSettings().min_track_time)
						{// we are too close to t = 0 but we do have the correct tolerance -so we exit.
							this->final_approximation_at_origin_ = latest_approx;
							return SuccessCode::MinTrackTimeReached;
						}
						else if(cauchy_samples_success != SuccessCode::Success)
						{
							this->final_approximation_at_origin_ = latest_approx;
							return cauchy_samples_success;
						}

					} while (approximate_error > this->Tolerances().final_tolerance);

					this->final_approximation_at_origin_ = latest_approx;
					return SuccessCode::Success;
				} //end CauchyEG
			};


			
		}//namespace endgame

	}//namespace tracking

} // namespace bertini
#endif
