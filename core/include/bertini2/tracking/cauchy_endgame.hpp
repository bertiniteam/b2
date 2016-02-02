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

\brief 

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
		
			template<typename TrackerType> 
			class CauchyEndgame : public Endgame
			{
				public:
				config::Cauchy cauchy_settings_; 

				std::deque<mpfr> pseg_times_; //times that are used to find the first approximation of the origin using the Power series approximation. 
				std::deque< Vec<mpfr> > pseg_samples_; //samples used for the first approximation. 
				std::deque<mpfr> cauchy_times_; //times of samples gathered while using CircleTrack.
				std::deque< Vec<mpfr> > cauchy_samples_; //samples gathered while using CircleTrack.

				mpfr_float min_closed_loop_tolerance_;
				mpfr_float max_closed_loop_tolerance_;

				const TrackerType & endgame_tracker_;


				void ClearTimesAndSamples(){pseg_times_.clear(); pseg_samples_.clear(); cauchy_times_.clear(); cauchy_samples_.clear();}

				void SetPSEGTimes(std::deque<mpfr> pseg_times_to_set) { pseg_times_ = pseg_times_to_set;}
				std::deque< mpfr > GetPSEGTimes() {return pseg_times_;}

				void SetPSEGSamples(std::deque< Vec<mpfr> > pseg_samples_to_set) { pseg_samples_ = pseg_samples_to_set;}
				std::deque< Vec<mpfr> > GetPSEGSamples() {return pseg_samples_;}

				void SetCauchyTimes(std::deque<mpfr> cauchy_times_to_set) { cauchy_times_ = cauchy_times_to_set;}
				std::deque< mpfr > GetCauchyTimes() {return cauchy_times_;}
	
				const TrackerType & GetEndgameTracker(){return endgame_tracker_;}

				void SetEndgameSettings(config::EndGame new_endgame_settings){endgame_settings_ = new_endgame_settings;}
				config::EndGame GetEndgameStruct(){ return endgame_settings_;}

				void SetCauchySettings(config::PowerSeries new_cauchy_settings){cauchy_settings_ = new_cauchy_settings;}
				config::PowerSeries GetCauchySettings(){return cauchy_settings_;}

				void SetSecuritySettings(config::Security new_endgame_security_settings){ endgame_security_ = new_endgame_security_settings;}
				config::Security GetSecuritySettings(){return endgame_security_;}

				void SetToleranceSettings(config::Tolerances new_tolerances_settings){endgame_tolerances_ = new_tolerances_settings;}
				config::Tolerances GetTolerancesSettings(){return endgame_tolerances_;}

				CauchyEndgame(TrackerType const& tracker) : endgame_tracker_(tracker)
				{
					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				
				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings ,config::EndGame new_endgame_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances

					SetCauchySettings(new_cauchy_settings);
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings ,config::EndGame new_endgame_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings ,config::EndGame new_endgame_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetEndgameSettings(new_endgame_settings);
					SetToleranceSettings(new_tolerances_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker,config::Cauchy new_cauchy_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings ,config::EndGame new_endgame_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetEndgameSettings(new_endgame_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetToleranceSettings(new_tolerances_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetSecuritySettings(new_security_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetToleranceSettings(new_tolerances_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker,config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetEndgameSettings(new_endgame_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetSecuritySettings(new_security_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetToleranceSettings(new_tolerances_settings);

					min_closed_loop_tolerance_ = std::max(mpfr_float(1e-10),endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
					max_closed_loop_tolerance_ = std::max(endgame_tolerances_.newton_during_endgame,endgame_tolerances_.final_tolerance_times_final_tolerance_multiplier);
				} 

				~CauchyEndgame() {};


				/*
				Input: 
					A starting time and starting sample space to start tracking. 

				Output: The sample space above starting time after we have tracked around the origin. 

				Details:
					This function uses the number of sample points to track in a polygonal path around the origin. 
					This function should be called the number of times it takes to loop around the origin and get back the original 
					point we started with. In essence this function will help determine the cycle number. 
				*/

				template<typename ComplexType> 
				Vec<ComplexType> CircleTrack(ComplexType starting_time, Vec<ComplexType> starting_sample)
				{
					bool first_step_of_path = true;
					SuccessCode tracking_success = SuccessCode::Success;
					auto current_angle = ComplexType(-2 * M_PI);
					auto number_of_steps = 0;
					ComplexType i("0.0","1.0"); 

					Vec<ComplexType> current_sample = starting_sample;
					Vec<ComplexType> next_sample;
					ComplexType current_time = starting_time;
					ComplexType next_time;

					cauchy_times_.push_back(current_time);
					cauchy_samples_.push_back(current_sample);

					if(endgame_settings_.num_sample_points < 3) // need to make sure we won't track right through the origin.
					{
						std::cout << "ERROR: The number of sample points " << endgame_settings_.num_sample_points << " for circle tracking must be >= 3!" << '\n';
					}	

					if(starting_time.norm() <= 0) // error checking on radius of circle to track around. 
					{
						std::cout << "ERROR: The radius of the circle " << starting_time << " needs to be positive! " << '\n';
					}

					//if this is the first step of the path, initialize the current step size as the maximum step size. 
					if(first_step_of_path)
					{
						//endgame_tracker_.current_step_size = max_step_size;			
						first_step_of_path = false;			
					}

					auto consecutive_success = 0; 
					auto next_angle_to_stop_at = mpfr_float(2 * M_PI * (1.0 / (endgame_settings_.num_sample_points) - 1));
					auto absolute_value_of_distance_left = (next_angle_to_stop_at - current_angle) * starting_time;

					for(unsigned ii = 0; ii < endgame_settings_.num_sample_points && tracking_success == bertini::tracking::SuccessCode::Success; ++ii)
					{
						// std::cout << "ii is " << ii << '\n';

					// 	//calculate the next angle to stop at -2 * PI + 2 * PI * i / M

						next_angle_to_stop_at = 2 * M_PI * ((ii + 1.0) / (endgame_settings_.num_sample_points) - 1);
						absolute_value_of_distance_left = (next_angle_to_stop_at - current_angle) * starting_time;

						// std::cout << "type of real part is " << typeid(next_angle_to_stop_at).name() << '\n';	
						// std::cout << "cos(angle) is " << mpfr_float( starting_time.abs() * cos(next_angle_to_stop_at)) << '\n';
						// std::cout << "type of imag part is " << typeid(next_angle_to_stop_at).name() << '\n';
						// std::cout << "sin(angle) is " << mpfr_float( starting_time.abs() * sin(next_angle_to_stop_at)) << '\n';	

						next_time = ComplexType(mpfr_float(starting_time.abs() * cos(next_angle_to_stop_at)),mpfr_float( starting_time.abs() * sin(next_angle_to_stop_at)));	

						// std::cout << "next time is " << next_time << '\n';

						while((next_angle_to_stop_at - current_angle).norm() > endgame_tolerances_.track_tolerance_during_endgame)
						{	
							if ((next_angle_to_stop_at - current_angle).norm() > endgame_tolerances_.track_tolerance_during_endgame)
							{
								// auto next_time = ComplexType(mpfr_float( starting_time.abs() * sin(next_angle_to_stop_at),mpfr_float( starting_time.abs() * cos(next_angle_to_stop_at));	

								SuccessCode tracking_success = endgame_tracker_.TrackPath(next_sample,current_time,next_time,current_sample);	

								if (tracking_success == bertini::tracking::SuccessCode::Success)
								{//successful track
									// consecutive_success++;
									current_angle = next_angle_to_stop_at;	
									current_sample = next_sample;
									current_time = next_time;
									cauchy_times_.push_back(current_time);
									cauchy_samples_.push_back(current_sample);

								}
								else if(tracking_success == SuccessCode::HigherPrecisionNecessary)
								{ 
									std::cout << "Higher precision necessary. " << '\n';
									return next_sample;
								}
								else if(tracking_success == SuccessCode::GoingToInfinity)
								{
									std::cout << "Going to infinity. " << '\n';
									return next_sample;
								}
								else if(tracking_success == SuccessCode::MatrixSolveFailure)
								{
									std::cout << "Matrix solve failure. " << '\n';
									return next_sample;
								}
								else if(tracking_success == SuccessCode::MaxNumStepsTaken)
								{
									std::cout << "Max num steps taken. " << '\n';
									return next_sample;
								}
								else if(tracking_success == SuccessCode::MaxPrecisionReached)
								{
									std::cout << "Max precision reached. " << '\n';
									return next_sample;
								}
								else if(tracking_success == SuccessCode::MinStepSizeReached)
								{
									std::cout << "Min step size reached." << '\n';
									return next_sample;
								}
								else if(tracking_success == SuccessCode::Failure)
								{
									std::cout << "Failure. " << '\n';
									return next_sample;
								}
								else if(tracking_success == SuccessCode::SingularStartPoint)
								{
									std::cout << "Singular start point. " << '\n';
									return next_sample;
								}				
								else
								{	

								}		
							}		
						}
					}

					//final leg of tracking around the circle
					if ((next_time - starting_time).norm() < endgame_tolerances_.track_tolerance_during_endgame )
					{
						// std::cout << "in final leg" <<'\n';

						ComplexType next_time = ComplexType(mpfr_float(starting_time.abs() * cos(next_angle_to_stop_at)),mpfr_float( starting_time.abs() * sin(next_angle_to_stop_at)));	
						SuccessCode tracking_success = endgame_tracker_.TrackPath(next_sample,current_time,next_time,current_sample);
					}
					if ( tracking_success == bertini::tracking::SuccessCode::Success)
					{
						return next_sample;
					}
					else
					{
						std::cout << "ERROR: See tracking_success for error." << '\n';
						return next_sample;
					}
					return next_sample;
				}

				/*
				Input: 
					All input is available through class data members. 

				Output: An mpfr_float that is an estimate on the c/k ratio shown in the book under Cauchy Endgame.

				Details:
					Using the formula for the cycle test outline in the Bertini Book pg. 53, we can compute an estimate for c/k.
					This estimate is used to help find stabilization of the cycle number. 
				*/
				
				mpfr_float ComputeCOverK()
				{
					std::deque<mpfr_float> c_over_k; // belongs in preEG_d2..... should be ComplexType

					const Vec<mpfr> & sample0 = pseg_samples_[0];
					const Vec<mpfr> & sample1 = pseg_samples_[1];
					const Vec<mpfr> & sample2 = pseg_samples_[2];

					Vec<mpfr> rand_vector = Vec<mpfr>::Random(sample0.size()); //should be a row vector for ease in multiplying.


					// //DO NOT USE Eigen .dot() it will do conjugate transpose which is not what we want.
					// //Also, the .transpose*rand_vector returns an expression template that we do .norm of since abs is not available for that expression type. 
					mpfr_float estimate = abs(log(abs((((sample2 - sample1).transpose()*rand_vector).norm())/(((sample1 - sample0).transpose()*rand_vector).norm()))));
					estimate = abs(log(endgame_settings_.sample_factor))/estimate;
					// std::cout << "estimate is " << estimate << '\n';
					if (estimate < 1)
					{
					  	return mpfr_float("1");
					}
					else
					{	
						return estimate;
					}
				}

				/*
				Input: An array of c/k estimates. If we have stabilization than these estimates will be withing some treshold. 

				Output: Boolean of true or false. 

				Details: We consider the ratio of consecutive c/k estimates. If we are less than the minimum needed for stabilization we return false.
				Otherwise, we return true. 
	
				*/

				bool CheckForCOverKStabilization(std::deque<mpfr_float> c_over_k_array)
				{
					bool return_value = true;

					for(unsigned ii = 1; ii < cauchy_settings_.num_needed_for_stabilization ; ++ii)
					{
						// std::cout << " a is " << c_over_k_array[ii-1] << '\n';
						// std::cout << "abs(a) is " << abs(c_over_k_array[ii-1]) << '\n';


						auto a = abs(c_over_k_array[ii-1]);
						auto b = abs(c_over_k_array[ii]);

						auto divide = a;

						if(a < b)
						{
							divide = a/b;
						}
						else
						{
							divide = b/a;
						}

						// std::cout << "divide is " << divide << '\n';

						if(divide <  cauchy_settings_.minimum_for_c_over_k_stabilization)
						{
							return_value = false;
						}
					}
					return return_value;
				}


				/*
				Input: A time value and the space value above that time.
					
				Output: An mpfr_float representing a tolerance threshold for declaring a loop to be closed. 


				Details: Why/how the heck does this work?!
					
				*/


				template<typename ComplexType>
				mpfr_float FindToleranceForClosedLoop(ComplexType x_time, Vec<ComplexType> x_sample)
				{
					auto degree_max = std::max(endgame_tracker_.AMP_config_.degree_bound,mpfr_float("2.0")); 
					auto K = endgame_tracker_.AMP_config_.coefficient_bound;
					mpfr_float N;
					mpfr_float M;
					mpfr_float L;

					if(max_closed_loop_tolerance_ < min_closed_loop_tolerance_)
					{
						max_closed_loop_tolerance_ = min_closed_loop_tolerance_;
					}

					// if(samples_.back().precision < 64) //use double
					// {
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
					auto jacobian_at_current_time = endgame_tracker_.GetSystem().Jacobian(x_sample,x_time);


					auto minimum_singular_value = Eigen::JacobiSVD< Mat<ComplexType> >(jacobian_at_current_time).singularValues()(endgame_tracker_.GetSystem().NumVariables() - 1 );

					auto norm_of_sample = x_sample.norm();
					L = pow(norm_of_sample,degree_max - 2);
					auto tol = K * L * M;
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
					// }
				}

				/*
				Input: An mpfr_float for the maximum error allowed for a loop to be consider closed. This tolerance is found above.
					
				Output: A boolean declaring if we have closed the loop or not. 

				Details: At this point a very simple function, look below. 
					
				*/

				//Need to check with Dan if the extra parts of this are necessary. 
				bool CheckClosedLoop(mpfr_float closed_loop_max_error)
				{
					//auto mininimum_digits = ceil(-1 * log(closed_loop_max_error) + mpfr_float(0.5));
					//is binary operations going to work with different precisions correctly
					if(abs(cauchy_samples_.front() - cauchy_samples_.back()) < closed_loop_max_error )
					{
						return true;
					}

					return false;	
				}

				/*
				Input: All input needed is available as class data members. 

				Output: Boolean declaring if our ratios are close enough or not. 

				Details: Finds minimum and maximum norms from tracking around the origin. Then we check a few heuristics to see if 
				we have good ratios or not. 
					 
				*/

				bool CompareCauchyRatios()
				{
					mpfr_float min = mpfr_float(1e300);
					mpfr_float max = 0;
					mpfr_float norm;


					if(cauchy_times_.front().norm() < cauchy_settings_.ratio_cutoff_time)
					{
						return true;
					}
					else
					{
						for(unsigned int ii=0; ii <= endgame_settings_.num_sample_points; ++ii)
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

						// std::cout << "max is " << max <<  '\n';
						// std::cout << "min is " << min << '\n';
						// std::cout << "final tol is " << endgame_tolerances_.final_tolerance << '\n';
						// std::cout << "min / max is " << min / max << '\n';
						// std::cout << "max cauchy ratio is " << cauchy_settings_.maximum_cauchy_ratio << '\n';
						// std::cout << "max - min is " << max - min << '\n';
						// std::cout << "1st condition " << (norm < cauchy_settings_.maximum_cauchy_ratio) << '\n';
						// std::cout << "2nd condition " << ((max - min ) > endgame_tolerances_.final_tolerance) << '\n';


						if(min > endgame_tolerances_.final_tolerance && max > endgame_tolerances_.final_tolerance)
						{
							norm = min / max;
							if(norm < cauchy_settings_.maximum_cauchy_ratio && (max - min) > endgame_tolerances_.final_tolerance)
							{
								return false;
							}
						}
					}
					return true;
				}

				void preCauchyLoops()
				{
					bool continue_loop = true;

					while(continue_loop)
					{
						
					}

				}

				





			};



		}//namespace endgame

	}//namespace tracking

} // namespace bertini
#endif