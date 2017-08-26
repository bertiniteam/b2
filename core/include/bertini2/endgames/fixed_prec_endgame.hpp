//This file is part of Bertini 2.
//
//fixed_prec_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fixed_prec_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fixed_prec_endgame.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University




#pragma once

/**
\file fixed_prec_endgame.hpp

\brief Contains the policy for fixed precision endgame types.
*/

#include "bertini2/trackers/fixed_precision_tracker.hpp"
#include "bertini2/trackers/fixed_precision_utilities.hpp"

#include "bertini2/endgames/config.hpp"
#include "bertini2/endgames/prec_base.hpp"


namespace bertini{ namespace endgame {

template<typename TrackerT>
class FixedPrecEndgame : public EndgamePrecPolicyBase<TrackerT>
{
public:
	using TrackerType = TrackerT;
	using BaseComplexType = typename tracking::TrackerTraits<TrackerType>::BaseComplexType;
	using BaseRealType = typename tracking::TrackerTraits<TrackerType>::BaseRealType;

	using BCT = BaseComplexType;
	using BRT = BaseRealType;

	template<typename... T>
	static
	unsigned EnsureAtUniformPrecision(T& ...args)
	{
		return bertini::tracking::fixed::EnsureAtUniformPrecision(args...);
	}

	template<typename T>
	static
	void EnsureAtPrecision(T const & obj, unsigned prec)
	{
		using bertini::Precision;
		if (Precision(obj)!=prec)
		{
			std::stringstream err_msg;
			err_msg << "ensuring precision of object failed; precision is " << Precision(obj) << " and required precision is " << prec;
			throw std::runtime_error(err_msg.str());
		}
	}

	SuccessCode RefineSampleImpl(Vec<BCT> & result, Vec<BCT> const& current_sample, BCT const& current_time, double tol, unsigned max_iterations) const
	{
		using RT = mpfr_float;
		using std::max;
		auto& TR = this->GetTracker();

		auto refinement_success = this->GetTracker().Refine(result,current_sample,current_time,
		                          	tol,
		                          	max_iterations);

		return SuccessCode::Success;
	}

	FixedPrecEndgame(TrackerT const& new_tracker) : EndgamePrecPolicyBase<TrackerT>(new_tracker)
	{}
}; // re: fixed prec endgame policy

template<>
struct EGPrecSelector<tracking::DoublePrecisionTracker>
{
	using type = FixedPrecEndgame<tracking::DoublePrecisionTracker>;
};

template<>
struct EGPrecSelector<tracking::MultiplePrecisionTracker>
{
	using type = FixedPrecEndgame<tracking::MultiplePrecisionTracker>;
};

template<class D>
struct EGPrecSelector<tracking::FixedPrecisionTracker<D>>
{
	using type = FixedPrecEndgame<tracking::FixedPrecisionTracker<D>>;
};

}} //re: namespaces




