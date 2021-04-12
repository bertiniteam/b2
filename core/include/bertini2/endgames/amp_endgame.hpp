//This file is part of Bertini 2.
//
//amp_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_endgame.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
// Tim Hodges, Colorado State University



#pragma once

/**
\file base_endgame.hpp

\brief Contains parent class, Endgame, the parent class for all endgames.
*/

#include "bertini2/trackers/amp_tracker.hpp"
#include "bertini2/trackers/adaptive_precision_utilities.hpp"

#include "bertini2/endgames/config.hpp"
#include "bertini2/endgames/prec_base.hpp"

#include "bertini2/endgames/events.hpp"
#include "bertini2/detail/observable.hpp"
namespace bertini{ namespace endgame {



/**
\brief Specifies some necessaries for AMP style endgame implementations, which differ from the fixed precision ones.
*/
class AMPEndgame : public virtual EndgamePrecPolicyBase<tracking::AMPTracker>
{
public:
	using TrackerT = tracking::AMPTracker;
	using EmitterType = AMPEndgame;

	template<typename... T>
	static
	unsigned EnsureAtUniformPrecision(T& ...args)
	{
		return tracking::adaptive::EnsureAtUniformPrecision(args...);
	}

	static
	void EnsureAtPrecision(double & obj, unsigned prec)
	{
		if (prec!=DoublePrecision())
			throw std::runtime_error("attempting to adjust precision of double to non-double precision");
	}

	static
	void EnsureAtPrecision(std::complex<double> & obj, unsigned prec)
	{
		if (prec!=DoublePrecision())
			throw std::runtime_error("attempting to adjust precision of std::complex<double> to non-double precision");
	}

	static
	void EnsureAtPrecision(mpfr_float & obj, unsigned prec)
	{
		using bertini::Precision;
		Precision(obj,prec);
	}

	static
	void EnsureAtPrecision(mpfr_complex & obj, unsigned prec)
	{
		using bertini::Precision;
		Precision(obj,prec);
	}



	SuccessCode RefineSampleImpl(Vec<mpfr_complex> & result, Vec<mpfr_complex> const& current_sample, mpfr_complex const& current_time, NumErrorT tol, unsigned max_iterations) const
	{

		using bertini::Precision;
		assert(Precision(current_time)==Precision(current_sample) && "precision of sample and time to be refined in AMP endgame must match");

		using RT = mpfr_float;
		using std::max;
		auto& TR = this->GetTracker();
		TR.ChangePrecision(Precision(current_time));

		auto refinement_success = this->GetTracker().Refine(result,current_sample,current_time,
		                          	tol,
		                          	max_iterations);

		
		if (refinement_success==SuccessCode::HigherPrecisionNecessary ||
		    refinement_success==SuccessCode::FailedToConverge)
		{
			

			using bertini::Precision;

			auto prev_precision = DefaultPrecision();
			auto higher_precision = max(prev_precision,LowestMultiplePrecision())+ PrecisionIncrement();
			DefaultPrecision(higher_precision);
			this->GetTracker().ChangePrecision(higher_precision);

			NotifyObservers(PrecisionChanged<EmitterType>(*this, prev_precision, higher_precision));

			auto next_sample_higher_prec = current_sample;
			Precision(next_sample_higher_prec, higher_precision);

			auto result_higher_prec = Vec<mpfr_complex>(current_sample.size());

			auto time_higher_precision = current_time;
			Precision(time_higher_precision,higher_precision);

			assert(time_higher_precision.precision()==DefaultPrecision());
			refinement_success = this->GetTracker().Refine(result_higher_prec,
			                                               next_sample_higher_prec,
			                                               time_higher_precision,
		                          							tol,
		                          							max_iterations);

			Precision(result, higher_precision);
			result = result_higher_prec;
			
			assert(Precision(result)==DefaultPrecision());
		}
		return refinement_success;
	}




	explicit
	AMPEndgame(TrackerT const& new_tracker) : EndgamePrecPolicyBase<TrackerT>(new_tracker)
	{}

	virtual ~AMPEndgame() = default;

}; // re: class AMPEndgame
		

template<>
struct EGPrecSelector<tracking::AMPTracker>
{
	using type = AMPEndgame;
};


} } // end namespaces 
				


