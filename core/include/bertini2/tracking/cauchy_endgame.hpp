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

#include <deque>
#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include "tracking/base_endgame.hpp"
#include <cstdio>
#include <typeinfo> 


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

				std::deque<mpfr> times_;
				std::deque< Vec<mpfr> > samples_;

				const TrackerType & endgame_tracker_;


				void ClearTimesAndSamples(){times_.clear(); samples_.clear();}

				void SetTimes(std::deque<mpfr> times_to_set) { times_ = times_to_set;}
				std::deque< mpfr > GetTimes() {return times_;}

				void SetSamples(std::deque< Vec<mpfr> > samples_to_set) { samples_ = samples_to_set;}
				std::deque< Vec<mpfr> > GetSamples() {return samples_;}
	
				const TrackerType & GetEndgameTracker(){return endgame_tracker_;}

				void SetEndgameSettings(config::EndGame new_endgame_settings){endgame_settings_ = new_endgame_settings;}
				config::EndGame GetEndgameStruct(){ return endgame_settings_;}

				void SetCauchySettings(config::PowerSeries new_cauchy_settings){cauchy_settings_ = new_cauchy_settings;}
				config::PowerSeries GetCauchySettings(){return cauchy_settings_;}

				void SetSecuritySettings(config::Security new_endgame_security_settings){ endgame_security_ = new_endgame_security_settings;}
				config::Security GetSecuritySettings(){return endgame_security_;}

				void SetToleranceSettings(config::Tolerances new_tolerances_settings){endgame_tolerances_ = new_tolerances_settings;}
				config::Tolerances GetTolerancesSettings(){return endgame_tolerances_;}

				CauchyEndgame(TrackerType const& tracker) : endgame_tracker_(tracker){} //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				
				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings ,config::EndGame new_endgame_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances

					SetCauchySettings(new_cauchy_settings);
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings ,config::EndGame new_endgame_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings ,config::EndGame new_endgame_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetEndgameSettings(new_endgame_settings);
					SetToleranceSettings(new_tolerances_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);
				} 

				CauchyEndgame(TrackerType const& tracker,config::Cauchy new_cauchy_settings, config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings ,config::EndGame new_endgame_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetEndgameSettings(new_endgame_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetSecuritySettings(new_security_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetEndgameSettings(new_endgame_settings);
					SetToleranceSettings(new_tolerances_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetSecuritySettings(new_security_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
					SetToleranceSettings(new_tolerances_settings);
				} 

				CauchyEndgame(TrackerType const& tracker,config::Security new_security_settings, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetSecuritySettings(new_security_settings);
					SetToleranceSettings(new_tolerances_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::EndGame new_endgame_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetEndgameSettings(new_endgame_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Cauchy new_cauchy_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetCauchySettings(new_cauchy_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Security new_security_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetSecuritySettings(new_security_settings);
				} 

				CauchyEndgame(TrackerType const& tracker, config::Tolerances new_tolerances_settings) : endgame_tracker_(tracker)  //constructor specifying the system. Member initialization list used to initialize tracker of TrackerType
				{// Order of settings is in alphebetical order Cauchy Endgame Security Tolerances
					SetToleranceSettings(new_tolerances_settings);
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

			};



		}//namespace endgame

	}//namespace tracking

} // namespace bertini
#endif