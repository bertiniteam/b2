//This file is part of Bertini 2.
//
//amp_powerseries_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_powerseries_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_powerseries_endgame.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Tim Hodges, Colorado State University




#pragma once

/**
\file amp_powerseries_endgame.hpp

\brief Contains the adaptive precision power series endgame type.
*/

#include "bertini2/tracking/powerseries_endgame.hpp"
#include "bertini2/tracking/amp_endgame.hpp"

namespace bertini{ namespace tracking { namespace endgame {


/**
\brief The Adaptive Precision Power Series Endgame, in a class.
*/
class AMPPowerSeriesEndgame : public PowerSeriesEndgame<AMPTracker,AMPPowerSeriesEndgame, dbl,mpfr>, 
							  public AMPEndgamePolicyBase
{
	
protected:

	void ChangeTemporariesPrecisionImpl(unsigned new_precision) const override
	{ }

	void MultipleToMultipleImpl(unsigned new_precision) const override
	{
		// first, change precision on the final approximation
		const auto& source_point = std::get<Vec<mpfr> >(this->final_approximation_at_origin_);
		auto& target_point = std::get<Vec<mpfr> >(this->final_approximation_at_origin_);
		target_point.resize(source_point.size());

		for (unsigned ii=0; ii<source_point.size(); ii++)
		{
			target_point(ii).precision(new_precision);
			target_point(ii) = mpfr(source_point(ii));
		}

		// then change precision of the permanent temporaries
		auto& times = std::get<TimeCont<mpfr> >(times_);
		for (auto& t : times)
			t.precision(new_precision);

		auto& samples = std::get< SampCont<mpfr> >(samples_);
		for (auto& s : samples)
			for (unsigned ii=0; ii<s.size(); ++ii)
				s(ii).precision(new_precision);
	}

	void DoubleToMultipleImpl(unsigned new_precision) const override
	{
		const auto& source_point = std::get<Vec<dbl> >(this->final_approximation_at_origin_);
		auto& target_point = std::get<Vec<mpfr> >(this->final_approximation_at_origin_);

		target_point.resize(source_point.size());
		for (unsigned ii=0; ii<source_point.size(); ii++)
		{
			target_point(ii).precision(new_precision);
			target_point(ii) = mpfr(source_point(ii));
		}


		auto& times_m = std::get<TimeCont<mpfr> >(times_);
		auto& times_d = std::get<TimeCont<dbl> >(times_);

		for (unsigned ii=0; ii<times_m.size(); ++ii)
		{
			times_m[ii].precision(new_precision);
			times_m[ii] = mpfr(times_d[ii]);
		}


		auto& samples_m = std::get< SampCont<mpfr> >(samples_);
		auto& samples_d = std::get< SampCont<dbl > >(samples_);

		for (unsigned ii=0; ii<times_m.size(); ++ii)
		{
			auto& s = samples_m[ii];
			auto& sd = samples_d[ii];
			for (unsigned jj=0; jj<s.size(); ++jj)
			{
				s(jj).precision(new_precision);
				s(jj) = mpfr(sd(jj));
			}
		}
	}

	void MultipleToDoubleImpl() const override
	{
		const auto& source_point = std::get<Vec<mpfr> >(this->final_approximation_at_origin_);
		auto& target_point = std::get<Vec<dbl> >(this->final_approximation_at_origin_);
		target_point.resize(source_point.size());

		for (unsigned ii=0; ii<source_point.size(); ii++)
			target_point(ii) = dbl(source_point(ii));


		auto& times_m = std::get<TimeCont<mpfr> >(times_);
		auto& times_d = std::get<TimeCont<dbl> >(times_);

		for (unsigned ii=0; ii<times_m.size(); ++ii)
		{
			times_d[ii] = dbl(times_m[ii]);
		}


		auto& samples_m = std::get< SampCont<mpfr> >(samples_);
		auto& samples_d = std::get< SampCont<dbl > >(samples_);

		for (unsigned ii=0; ii<times_m.size(); ++ii)
		{
			auto& s = samples_m[ii];
			auto& sd = samples_d[ii];
			for (unsigned jj=0; jj<s.size(); ++jj)
			{
				sd(jj) = dbl(s(jj));
			}
		}
	}
public:
	using TrackerType = AMPTracker;
	using EGType = PowerSeriesEndgame<TrackerType, AMPPowerSeriesEndgame, dbl, mpfr>;
	
	SuccessCode RefineSample(Vec<mpfr> & result, Vec<mpfr> const& current_sample, mpfr const& current_time) const
	{
		using RT = mpfr_float;
		using std::max;
		auto& TR = this->GetTracker();
		RT refinement_tolerance = this->Tolerances().final_tolerance/100;
		auto refinement_success = this->GetTracker().Refine(result,current_sample,current_time,
		                          	refinement_tolerance,
		                          	this->EndgameSettings().max_num_newton_iterations);

		
		if (refinement_success==SuccessCode::HigherPrecisionNecessary ||
		    refinement_success==SuccessCode::FailedToConverge)
		{
			using bertini::Precision;

			auto prev_precision = this->Precision();
			auto temp_higher_prec = max(prev_precision,LowestMultiplePrecision())+ PrecisionIncrement();
			mpfr_float::default_precision(temp_higher_prec);
			this->GetTracker().ChangePrecision(temp_higher_prec);


			auto next_sample_higher_prec = current_sample;  Precision(next_sample_higher_prec, temp_higher_prec);
			auto result_higher_prec = Vec<mpfr>(current_sample.size());
			auto time_higher_precision = current_time; Precision(time_higher_precision, temp_higher_prec);

			assert(time_higher_precision.precision()==mpfr_float::default_precision());
			RT refinement_tolerance = this->Tolerances().final_tolerance/100;
			refinement_success = this->GetTracker().Refine(result_higher_prec,
			                                               next_sample_higher_prec,
			                                               time_higher_precision,
		                          							refinement_tolerance,
		                          							this->EndgameSettings().max_num_newton_iterations);

			mpfr_float::default_precision(prev_precision);
			this->GetTracker().ChangePrecision(prev_precision);
			result = result_higher_prec;
			Precision(result,prev_precision);
			assert(result(0).precision()==mpfr_float::default_precision());
		}

		return refinement_success;
	}

	SuccessCode RefineSample(Vec<dbl> & result, Vec<dbl> const& current_sample, dbl const& current_time) const
	{
		using RT = double;

		auto refinement_success = this->GetTracker().Refine(result,current_sample,current_time,
		                          	RT(this->Tolerances().final_tolerance)/100,
		                          	this->EndgameSettings().max_num_newton_iterations);

		
		if (refinement_success==SuccessCode::HigherPrecisionNecessary ||
		    refinement_success==SuccessCode::FailedToConverge)
		{
			auto prev_precision = this->Precision();
			auto temp_higher_prec = LowestMultiplePrecision();
			mpfr_float::default_precision(temp_higher_prec);
			this->GetTracker().ChangePrecision(temp_higher_prec);


			auto next_sample_higher_prec = Vec<mpfr>(current_sample.size());
			for (int ii=0; ii<current_sample.size(); ++ii)
				next_sample_higher_prec(ii) = mpfr(current_sample(ii));

			auto result_higher_prec = Vec<mpfr>(current_sample.size());
			mpfr time_higher_precision(current_time);
			mpfr_float refinement_tolerance = this->Tolerances().final_tolerance/100;
			refinement_success = this->GetTracker().Refine(result_higher_prec,
			                                               next_sample_higher_prec,
			                                               time_higher_precision,
		                          							refinement_tolerance,
		                          							this->EndgameSettings().max_num_newton_iterations);


			mpfr_float::default_precision(prev_precision);
			this->GetTracker().ChangePrecision(prev_precision);
			for (unsigned ii(0); ii<current_sample.size(); ++ii)
				result(ii) = dbl(result_higher_prec(ii));
		}

		return refinement_success;
	}
	
	explicit AMPPowerSeriesEndgame(TrackerType const& tr, 
	                               const std::tuple< const config::PowerSeries &,
	                               					 const config::Endgame<BRT>&, 
	                               				     const config::Security<BRT>&, 
	                               				     const config::Tolerances<BRT>& 
	                               				    > & settings )
      : EGType(tr, settings)
   	{}

    template< typename... Ts >
		AMPPowerSeriesEndgame(TrackerType const& tr, const Ts&... ts ) : AMPPowerSeriesEndgame(tr, Unpermute<config::PowerSeries, config::Endgame<BRT>, config::Security<BRT>, config::Tolerances<BRT> >( ts... ) ) 
		{}


}; // re: AMPPowerSeriesEndgame

} } } //namespaces


