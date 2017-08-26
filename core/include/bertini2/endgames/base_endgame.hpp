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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University



#ifndef BERTINI_TRACKING_BASE_ENDGAME_HPP
#define BERTINI_TRACKING_BASE_ENDGAME_HPP

#pragma once
/**
\file base_endgame.hpp

\brief Contains base class, Endgame.
*/

#include <iostream>
#include <typeinfo>


#include "bertini2/mpfr_complex.hpp"
#include "bertini2/limbo.hpp"

#include "bertini2/system/system.hpp"

#include "bertini2/detail/enable_permuted_arguments.hpp"

#include "bertini2/trackers/config.hpp"
#include "bertini2/endgames/config.hpp"
#include "bertini2/endgames/interpolation.hpp"

#include "bertini2/logging.hpp"

#include "bertini2/detail/configured.hpp"


namespace bertini{ namespace endgame {

			
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
template<class FlavorT, class PrecT>
class EndgameBase : 
	public detail::Configured< typename AlgoTraits<FlavorT>::NeededConfigs >,
	public PrecT
{
public:
	using TrackerType = typename PrecT::TrackerType;

	using BaseComplexType = typename tracking::TrackerTraits<TrackerType>::BaseComplexType;
	using BaseRealType = typename tracking::TrackerTraits<TrackerType>::BaseRealType;


protected:

	using BCT = BaseComplexType;
	using BRT = BaseRealType;


	using Configured = detail::Configured< typename AlgoTraits<FlavorT>::NeededConfigs >;
	using Configs = typename AlgoTraits<FlavorT>::NeededConfigs;
	using ConfigsAsTuple = typename Configs::ToTuple;

	// a list of all the needed arithemtic types (complex for complex trackers)
	using NeededTypes = detail::TypeList<BCT>;
	using TupOfVec = typename NeededTypes::ToTupleOfVec;
	using TupOfReal = typename NeededTypes::ToTupleOfReal;
	using TupleOfTimes = typename NeededTypes::template ToTupleOfCont<TimeCont>;
	using TupleOfSamps = typename NeededTypes::template ToTupleOfCont<SampCont>;



	// universal endgame state variables
	mutable TupOfVec final_approximation_; 
	mutable TupOfVec previous_approximation_; 
	mutable unsigned int cycle_number_ = 0; 
	mutable NumErrorT approximate_error_;





	/**
	\brief convert the base endgame into the derived type.

	This enables the CRPT as used by the endgames
	*/
	const FlavorT& AsFlavor() const
	{
		return static_cast<const FlavorT&>(*this);
	}

	/**
	\brief Non-const version of AsFlavor
	*/
	FlavorT& AsFlavor()
	{
		return static_cast<FlavorT&>(*this);
	}

public:

	/**
	\brief The main function for running an endgame, from time to time, from a given point to a possibly singular solution.
	*/
	SuccessCode Run(const BCT & start_time, const Vec<BCT> & start_point, BCT const& target_time)
	{
		return this->AsFlavor().RunImpl(start_time, start_point, target_time);
	}

	/**
	\brief Run the endgame, shooting for default time of t=0.

	\see Run
	*/
	SuccessCode Run(BCT const& start_time, Vec<BCT> const& start_point)
	{
		return Run(start_time, start_point, static_cast<BCT>(0));
	}


	template<typename CT>
	SuccessCode RefineAllSamples(SampCont<CT> & samples, TimeCont<CT> & times)
	{
		for (size_t ii=0; ii<samples.size(); ++ii)
		{
			auto refine_success = this->RefineSample(samples[ii], samples[ii],  times[ii], 
										this->FinalTolerance() * this->EndgameSettings().sample_point_refinement_factor,
										this->EndgameSettings().max_num_newton_iterations);
			if (refine_success != SuccessCode::Success)
			{
				// BOOST_LOG_TRIVIAL(severity_level::trace) << "refining failed, code " << int(refine_success);
				return refine_success;
			}
		}

		if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec) // known at compile time
		{
			auto max_precision = this->EnsureAtUniformPrecision(times, samples);
			this->GetSystem().precision(max_precision);
		}

		return SuccessCode::Success;
	}


	/**
	A function passed off to the precision-specific endgame part
	*/
	SuccessCode RefineSample(Vec<BCT> & result, Vec<BCT> const& current_sample, BCT const& current_time, NumErrorT tol, unsigned max_iterations) const
	{
		return this->RefineSampleImpl(result, current_sample, current_time, tol, max_iterations);
	}

	void ChangePrecision(unsigned p)
	{
		AsFlavor().ChangePrecision(p);
		PrecT::ChangePrecision();
	}


	/**
	\brief This function is inferior to the templated Get<ConfigT> function provided
	*/
	inline
	const auto & EndgameSettings() const
	{
		return Configured::template Get<EndgameConfig>();
	}
	
	inline
	const auto & SecuritySettings() const
	{
		return this->template Get<SecurityConfig>();
	}

	explicit EndgameBase(TrackerType const& tr, const ConfigsAsTuple& settings ) :
      	Configured( settings ), PrecT(tr)
   	{}

    template< typename... Ts >
		EndgameBase(TrackerType const& tr, const Ts&... ts ) : EndgameBase(tr, Configs::Unpermute( ts... ) ) 
		{}


		inline unsigned CycleNumber() const { return cycle_number_;}
	inline void CycleNumber(unsigned c) { cycle_number_ = c;}
	inline void IncrementCycleNumber(unsigned inc) { cycle_number_ += inc;}

	
	/**
	\brief Get the final tolerance to which we are tracking the solution.
	*/
	inline
	const auto& FinalTolerance() const
	{
		return this->template Get<EndgameConfig>().final_tolerance;
	}



	/**
	\brief Setter for the final tolerance.
	*/
	inline
	void SetFinalTolerance(BRT const& ft){this->template Get<EndgameConfig>().final_tolerance = ft;}




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
	inline
	NumErrorT ApproximateError() const
	{
		return approximate_error_;
	}


	/**
	Get the latest time at which a point on the path was computed
	*/
	inline 
	const BCT& LatestTime() const
	{
		return AsFlavor().LatestTimeImpl();
	}

	// /**
	// \brief Get the system being tracked on, which is referred to by the tracker.
	// */
	// inline
	// const System& GetSystem() const 
	// { 
	// 	return GetTracker().GetSystem();
	// }


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
		assert(this->template Get<EndgameConfig>().num_sample_points>0 && "number of sample points must be positive");

		if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			assert(Precision(start_time)==Precision(x_endgame) && "Computing initial samples requires input time and space with uniform precision");
		}

		samples.clear();
		times.clear();

		samples.push_back(x_endgame);
		times.push_back(start_time);

		auto num_vars = this->GetSystem().NumVariables();
		//start at 1, because the input point is the 0th element.
		for(int ii=1; ii < this->template Get<EndgameConfig>().num_sample_points; ++ii)
		{ 
			times.emplace_back((times[ii-1] + target_time) * RT(this->template Get<EndgameConfig>().sample_factor)); // next time is a point between the previous time and target time.
			samples.emplace_back(Vec<CT>(num_vars));											   // sample_factor gives us some point between the two, usually the midpoint.		

			auto tracking_success = this->GetTracker().TrackPath(samples[ii],times[ii-1],times[ii],samples[ii-1]);
			this->EnsureAtPrecision(times[ii],Precision(samples[ii]));

			if (tracking_success!=SuccessCode::Success)
				return tracking_success;
		}

		return SuccessCode::Success;
	}

};
			
} }// end namespaces bertini
				

#endif
