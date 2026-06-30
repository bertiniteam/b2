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
	using TrackerT = tracking::AMPTracker;  ///< The path-tracker type.
	using EmitterType = AMPEndgame;  ///< The event-emitter type.

	/// \brief Bring a set of values to a common (highest) precision.
	template<typename... T>
	static
	unsigned EnsureAtUniformPrecision(T& ...args)
	{
		return tracking::adaptive::EnsureAtUniformPrecision(args...);
	}

	/// \brief No-op precision check for a double (throws if asked for non-double precision).
	static
	void EnsureAtPrecision(double & /*obj*/, unsigned prec)
	{
		if (prec!=DoublePrecision())
			throw std::runtime_error("attempting to adjust precision of double to non-double precision");
	}

	/// \brief No-op precision check for a std::complex<double> (throws if asked for non-double precision).
	static
	void EnsureAtPrecision(std::complex<double> & /*obj*/, unsigned prec)
	{
		if (prec!=DoublePrecision())
			throw std::runtime_error("attempting to adjust precision of std::complex<double> to non-double precision");
	}

	/// \brief Set a multiprecision real to the given precision.
	static
	void EnsureAtPrecision(real_mp & obj, unsigned prec)
	{
		using bertini::Precision;
		Precision(obj,prec);
	}

	/// \brief Set a multiprecision complex to the given precision.
	static
	void EnsureAtPrecision(complex_mp & obj, unsigned prec)
	{
		using bertini::Precision;
		Precision(obj,prec);
	}



	/// \brief Refine a single sample point, escalating precision if Newton refinement needs it.
	SuccessCode RefineSampleImpl(Vec<complex_mp> & result, Vec<complex_mp> const& current_sample, complex_mp const& current_time, NumErrorT tol, unsigned max_iterations) const
	{

		using bertini::Precision;
		assert(Precision(current_time)==Precision(current_sample) && "precision of sample and time to be refined in AMP endgame must match");

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

			// thread-local precision: endgames may run on std::thread workers, so
			// neither read nor write the global default precision.
			auto prev_precision = ThreadPrecision();
			auto higher_precision = max(prev_precision,LowestMultiplePrecision())+ PrecisionIncrement();
			SetThreadPrecision(higher_precision);
			this->GetTracker().ChangePrecision(higher_precision);

			NotifyObservers(PrecisionChanged<EmitterType>(*this, prev_precision, higher_precision));

			auto next_sample_higher_prec = current_sample;
			Precision(next_sample_higher_prec, higher_precision);

			auto result_higher_prec = Vec<complex_mp>(current_sample.size());

			auto time_higher_precision = current_time;
			Precision(time_higher_precision,higher_precision);

			assert(time_higher_precision.precision()==ThreadPrecision());
			refinement_success = this->GetTracker().Refine(result_higher_prec,
			                                               next_sample_higher_prec,
			                                               time_higher_precision,
		                          							tol,
		                          							max_iterations);

			Precision(result, higher_precision);
			result = result_higher_prec;
			
			assert(Precision(result)==ThreadPrecision());
		}
		return refinement_success;
	}


	/**
	Hardware-double fast-lane refine for the adaptive-numeric-type (double-first) endgame.

	Pure-(i) escalation: this does NOT raise precision in place.  The endgame is computing in the
	complex_dbl slot, with the tracker and system already at double precision; we refine there.  If
	double cannot reach the tolerance the tracker returns HigherPrecisionNecessary / FailedToConverge,
	and we hand that code straight back to RunImplAMP, which migrates every container to mpfr and
	retries in the complex_mp slot (where the in-place mp->higher-mp escalation above takes over).
	*/
	SuccessCode RefineSampleImpl(Vec<complex_dbl> & result, Vec<complex_dbl> const& current_sample, complex_dbl const& current_time, NumErrorT tol, unsigned max_iterations) const
	{
		return this->GetTracker().Refine(result, current_sample, current_time, tol, max_iterations);
	}



	// ================================================================================================
	//   Adaptive-numeric-type endgame state + the precision/numeric-type helpers shared by both AMP
	//   endgame flavors (Cauchy and PowerSeries).  These live here, in the AMP precision policy, because
	//   they are entirely about the complex_dbl <-> complex_mp crossing and the AMP tracker's authority --
	//   not generic endgame concerns.  The flavors inherit them (AMPEndgame is the precision-policy base
	//   of EndgameBase) and call them through this->.
	// ================================================================================================

	// current_endgame_precision_ is the precision the endgame is *currently computing in*: DoublePrecision()
	// is the hardware-complex_dbl fast lane, higher is the mpfr slot.  adaptive_numeric_type_active_ gates
	// the escalation hooks in the flavors' shared phase methods so they fire only while the double-first
	// driver is orchestrating (never for the AMP-PowerSeries forwarding path, were one to exist).
	mutable unsigned current_endgame_precision_ = DoublePrecision();
	mutable bool     adaptive_numeric_type_active_ = false;

	// Choose the next (higher) working precision: at least LowestMultiplePrecision(), following the
	// tracker's authority if it climbed higher, and strictly above the current precision.
	unsigned NextEscalatedPrecision() const
	{
		using std::max;
		unsigned from_tracker = this->GetTracker().GetCurrentPrecision();
		unsigned newprec = max(static_cast<unsigned>(LowestMultiplePrecision()), from_tracker);
		if (newprec <= current_endgame_precision_)
			newprec = current_endgame_precision_ + PrecisionIncrement();
		return newprec;
	}

	// Cross a deque of times / samples, or a single vector, from the complex_dbl slot to the complex_mp
	// slot, then set the mpfr precision via the existing SetPrecision helper.  Templated on the tuple type
	// so this header need not name the endgames' TupleOfTimes / TupleOfSamps / TupOfVec aliases.  The
	// element crossing is the same one the AMP tracker performs; no library Vec-cast exists for these types.
	template<typename TimesTuple>
	void CrossTimesUp(TimesTuple& times, unsigned newprec) const
	{
		auto& from = std::get<TimeCont<complex_dbl> >(times);
		auto& to   = std::get<TimeCont<complex_mp> >(times);
		to.clear();
		for (auto const& t : from) to.push_back(complex_mp(t));
		from.clear();
		tracking::adaptive::SetPrecision(to, newprec);
	}

	template<typename SampsTuple>
	void CrossSampsUp(SampsTuple& samps, unsigned newprec) const
	{
		auto& from = std::get<SampCont<complex_dbl> >(samps);
		auto& to   = std::get<SampCont<complex_mp> >(samps);
		to.clear();
		for (auto const& v : from)
		{
			Vec<complex_mp> w(v.size());
			for (Eigen::Index i = 0; i < v.size(); ++i) w(i) = complex_mp(v(i));
			to.push_back(std::move(w));
		}
		from.clear();
		tracking::adaptive::SetPrecision(to, newprec);
	}

	template<typename VecTuple>
	void CrossVecUp(VecTuple& vec, unsigned newprec) const
	{
		using bertini::Precision;
		auto& from = std::get<Vec<complex_dbl> >(vec);
		auto& to   = std::get<Vec<complex_mp> >(vec);
		to.resize(from.size());
		for (Eigen::Index i = 0; i < from.size(); ++i) to(i) = complex_mp(from(i));
		if (to.size() > 0) Precision(to, newprec);
		from.resize(0);
	}

	// Downcast a complex_mp vector to the hardware complex_dbl fast lane.
	Vec<complex_dbl> DowncastToDouble(Vec<complex_mp> const& v) const
	{
		Vec<complex_dbl> out(v.size());
		for (Eigen::Index i = 0; i < v.size(); ++i) out(i) = complex_dbl(v(i));
		return out;
	}

	// Copy a complex_mp scalar / vector at the endgame's current working precision.
	complex_mp AtActivePrecisionScalar(complex_mp const& x) const
	{ using bertini::Precision; complex_mp r = x; Precision(r, current_endgame_precision_); return r; }

	Vec<complex_mp> AtActivePrecisionVec(Vec<complex_mp> const& v) const
	{ using bertini::Precision; Vec<complex_mp> r = v; if (r.size() > 0) Precision(r, current_endgame_precision_); return r; }

	// Widen an active-type vector (a fast-lane extrapolation result) to complex_mp at the given precision:
	// the single boundary conversion of an approximation back to the ambient type.
	template<typename ComplexT>
	Vec<complex_mp> ToBCT(Vec<ComplexT> const& v, unsigned prec) const
	{
		using bertini::Precision;
		Vec<complex_mp> out(v.size());
		for (Eigen::Index i = 0; i < v.size(); ++i) out(i) = complex_mp(v(i));
		if (out.size() > 0) Precision(out, prec);
		return out;
	}



	/// \brief Construct the AMP endgame precision policy for a tracker.
	explicit
	AMPEndgame(TrackerT const& new_tracker) : EndgamePrecPolicyBase<TrackerT>(new_tracker)
	{}

	virtual ~AMPEndgame() = default;

}; // re: class AMPEndgame


/// \brief Selects the AMP endgame precision policy for the adaptive-precision tracker.
template<>
struct EGPrecSelector<tracking::AMPTracker>
{
	using type = AMPEndgame;  ///< The endgame precision-policy type for this tracker.
};


} } // end namespaces 
				


