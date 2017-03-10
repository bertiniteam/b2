//This file is part of Bertini 2.
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
//along with base_endgame.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Tim Hodges, Colorado State University



#ifndef BERTINI_TRACKING_BASE_ENDGAME_HPP
#define BERTINI_TRACKING_BASE_ENDGAME_HPP

/**
\file base_endgame.hpp

\brief Contains base class, EndgameBase.
*/

#include <iostream>
#include <typeinfo>

#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include "bertini2/mpfr_complex.hpp"
#include "bertini2/limbo.hpp"

#include "bertini2/system/system.hpp"

#include "bertini2/detail/enable_permuted_arguments.hpp"

#include "bertini2/trackers/config.hpp"
#include "bertini2/endgames/config.hpp"
#include "bertini2/endgames/interpolation.hpp"

#include "bertini2/logging.hpp"

#include "bertini2/detail/configured.hpp"
namespace bertini{ 

	namespace tracking {

		namespace endgame {

			
			/**
			\class EndgameBase

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
			template<class TrackerType, class FinalEGT, typename... UsedNumTs>
			class EndgameBase : public detail::Configured<
				config::Endgame<typename TrackerTraits<TrackerType>::BaseRealType>,
				config::Security<typename TrackerTraits<TrackerType>::BaseRealType>>
			{
			protected:

				/**
				\brief convert the base endgame into the derived type.

				This enables the CRPT as used by the endgames
				*/
				const FinalEGT& AsDerived() const
				{
					return static_cast<const FinalEGT&>(*this);
				}

				using BaseComplexType = typename TrackerTraits<TrackerType>::BaseComplexType;
				using BaseRealType = typename TrackerTraits<TrackerType>::BaseRealType;

				using BCT = BaseComplexType;
				using BRT = BaseRealType;

				using Config = detail::Configured<
					config::Endgame<typename TrackerTraits<TrackerType>::BaseRealType>,
					config::Security<typename TrackerTraits<TrackerType>::BaseRealType>>;

				// state variables
				mutable std::tuple<Vec<UsedNumTs>... > final_approximation_; 
				mutable std::tuple<Vec<UsedNumTs>... > previous_approximation_; 
				mutable unsigned int cycle_number_ = 0; 
				mutable std::tuple<BaseRealType> approximate_error_;

				/**
				\brief A tracker that must be passed into the endgame through a constructor. This tracker is what will be used to track to all time values during the endgame. 
				*/
				std::reference_wrapper<const TrackerType> tracker_;

			public:

				inline
				const auto & EndgameSettings() const
				{
					return Config::template Get<config::Endgame<BRT>>();
				}
				
				inline
				const auto & SecuritySettings() const
				{
					return this->template Get<config::Security<BRT>>();
				}

				explicit EndgameBase(TrackerType const& tr, const std::tuple< const config::Endgame<BRT>&, const config::Security<BRT>&>& settings )
			      : tracker_(std::ref(tr)),
			      	Config( std::get<0>(settings),  std::get<1>(settings) )
			   	{}

			    template< typename... Ts >
   				EndgameBase(TrackerType const& tr, const Ts&... ts ) : EndgameBase(tr, Unpermute< config::Endgame<BRT>, config::Security<BRT> >( ts... ) ) 
   				{}


   				inline unsigned CycleNumber() const { return cycle_number_;}
				inline void CycleNumber(unsigned c) { cycle_number_ = c;}
				inline void IncrementCycleNumber(unsigned inc) { cycle_number_ += inc;}

				

				// const auto& EndgameSettings() const
				// {
				// 	return endgame_settings_;
				// }

				inline
				const auto& FinalTolerance() const
				{
					return this->template Get<config::Endgame<BRT>>().final_tolerance;
				}

				// const auto& SecuritySettings() const
				// {
				// 	return security_;
				// }



				// /**
				// \brief Setter for the general settings.
				// */
				// void SetEndgameSettings(config::Endgame<BRT> new_endgame_settings){endgame_settings_ = new_endgame_settings;}


				

				/**
				\brief Setter for the security settings.
				*/
				// void SetSecuritySettings(config::Security<BRT> new_endgame_security_settings){ security_ = new_endgame_security_settings;}


				/**
				\brief Setter for the final tolerance.
				*/
				inline
				void SetFinalTolerance(BRT const& ft){this->template Get<config::Endgame<BRT>>().final_tolerance = ft;}


				/**
				Tell the endgame to use the given tracker.  Takes a reference.  

				\note Ensure the tracker you are using doesn not go out of scope!
				*/
				inline
				void SetTracker(TrackerType const& new_tracker)
				{
					tracker_ = std::ref(new_tracker); // rebind the reference
				}


				/**
				\brief Getter for the tracker used inside an instance of the endgame. 
				*/
				inline
				const TrackerType & GetTracker() const
				{
					return tracker_.get();
				}

				/**
				\brief Get the most-recent approximation
				*/
				template<typename CT>
				inline
				const Vec<CT>& FinalApproximation() const 
				{
					return std::get<Vec<CT> >(final_approximation_);
				}

				/**
				\brief Get the second-most-recent approximation
				*/
				template<typename CT>
				inline
				const Vec<CT>& PreviousApproximation() const 
				{
					return std::get<Vec<CT> >(final_approximation_);
				}

				/**
				\brief Get the most recent accuracy estimate
				*/
				template<typename RT>
				inline
				const RT& ApproximateError() const
				{
					return std::get<RT>(approximate_error_);
				}

				/**
				\brief Get the system being tracked on, which is referred to by the tracker.
				*/
				inline
				const System& GetSystem() const 
				{ return GetTracker().GetSystem();}


				/**
				\brief Populates time and space samples so that we are ready to start the endgame. 

				## Input

				  		start_time: is the time when we start the endgame process usually this is .1
						x_endgame: is the space value at start_time
						times: a deque of time values. These values will be templated to be CT 
						samples: a deque of sample values that are in correspondence with the values in times. These values will be vectors with entries of CT. 

				## Output

					SuccessCode indicating whether tracking to all samples was successful.


				## Details

					The first sample will be (x_endgame) and the first time is start_time.
					From there we do a geometric progression using the sample factor (which by default is 1/2).
					Hence, next_time = start_time * sample_factor.
					We track then to the next_time and construct the next_sample.

				\param start_time The time value at which we start the endgame. 
				\param target_time The time value that we are trying to find a solution to. 
				\param x_endgame The current space point at start_time.
				\param times A deque that will hold all the time values of the samples we are going to use to start the endgame. 
				\param samples a deque that will hold all the samples corresponding to the time values in times. 

				\tparam CT The complex number type.
				*/	
				template<typename CT>
				SuccessCode ComputeInitialSamples(const CT & start_time,const CT & target_time, const Vec<CT> & x_endgame, TimeCont<CT> & times, SampCont<CT> & samples) // passed by reference to allow times to be filled as well.
				{	
					using RT = typename Eigen::NumTraits<CT>::Real;
					assert(this->template Get<config::Endgame<BRT>>().num_sample_points>0 && "number of sample points must be positive");

					if (TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						assert(Precision(start_time)==Precision(x_endgame) && "Computing initial samples requires input time and space with uniform precision");
					}

					samples.clear();
					times.clear();

					samples.push_back(x_endgame);
					times.push_back(start_time);

					auto num_vars = GetSystem().NumVariables();
					//start at 1, because the input point is the 0th element.
					for(int ii=1; ii < this->template Get<config::Endgame<BRT>>().num_sample_points; ++ii)
					{ 
						times.emplace_back((times[ii-1] + target_time) * RT(this->template Get<config::Endgame<BRT>>().sample_factor)); // next time is a point between the previous time and target time.
						samples.emplace_back(Vec<CT>(num_vars));											   // sample_factor gives us some point between the two, usually the midpoint.		

						auto tracking_success = GetTracker().TrackPath(samples[ii],times[ii-1],times[ii],samples[ii-1]);
						AsDerived().EnsureAtPrecision(times[ii],Precision(samples[ii]));

						if (tracking_success!=SuccessCode::Success)
							return tracking_success;
					}

					return SuccessCode::Success;
				}

			};
			
		}// end namespace endgame
	}// end namespace tracking 
}// end namespace bertini
				

#endif
