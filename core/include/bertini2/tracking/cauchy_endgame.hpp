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


namespace bertini{ 

	namespace tracking {

		namespace endgame {
		
			template<typename TrackerType> 
			class CauchyEndgame : public Endgame{
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



		}//namespace endgame

	}//namespace tracking

} // namespace bertini
#endif