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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
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

/// \brief Fixed-precision endgame precision policy (refines and ensures precision at a single working precision).
template<typename TrackerT>
class FixedPrecEndgame : public virtual EndgamePrecPolicyBase<TrackerT>
{
public:
	using TrackerType = TrackerT;  ///< The path-tracker type.
	using BaseComplexT = typename tracking::TrackerTraits<TrackerType>::BaseComplexT;  ///< The complex number type.
	using BaseRealT = typename tracking::TrackerTraits<TrackerType>::BaseRealT;  ///< The real number type.

	using BCT = BaseComplexT;  ///< The complex number type.
	using BRT = BaseRealT;  ///< The real number type.

	/// \brief Bring a set of values to a common precision (a no-op of effect under fixed precision).
	template<typename... T>
	static
	unsigned EnsureAtUniformPrecision(T& ...args)
	{
		return bertini::tracking::fixed::EnsureAtUniformPrecision(args...);
	}

	/// \brief Assert that an object is already at the given precision (throws otherwise).
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

	/// \brief Refine a single sample point at the (fixed) working precision.
	SuccessCode RefineSampleImpl(Vec<BCT> & result, Vec<BCT> const& current_sample, BCT const& current_time, double tol, unsigned max_iterations) const
	{
		auto refinement_success = this->GetTracker().Refine(result,current_sample,current_time,
		                          	tol,
		                          	max_iterations);

		return refinement_success;
	}

	/// \brief Construct the fixed-precision endgame policy for a tracker.
	explicit
	FixedPrecEndgame(TrackerT const& new_tracker) : EndgamePrecPolicyBase<TrackerT>(new_tracker)
	{}

	virtual ~FixedPrecEndgame() = default;
}; // re: fixed prec endgame policy

/// \brief Selects the fixed-precision endgame policy for the double-precision tracker.
template<>
struct EGPrecSelector<tracking::DoublePrecisionTracker>
{
	using type = FixedPrecEndgame<tracking::DoublePrecisionTracker>;  ///< The endgame precision-policy type.
};

/// \brief Selects the fixed-precision endgame policy for the multiple-precision tracker.
template<>
struct EGPrecSelector<tracking::MultiplePrecisionTracker>
{
	using type = FixedPrecEndgame<tracking::MultiplePrecisionTracker>;  ///< The endgame precision-policy type.
};

/// \brief Selects the fixed-precision endgame policy for a fixed-precision tracker.
template<class D>
struct EGPrecSelector<tracking::FixedPrecisionTracker<D>>
{
	using type = FixedPrecEndgame<tracking::FixedPrecisionTracker<D>>;  ///< The endgame precision-policy type.
};

}} //re: namespaces




