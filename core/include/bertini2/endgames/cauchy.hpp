//This file is part of Bertini 2.
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
//along with cauchy_endgame.hpp.  If not, see <http://www.gnu.org/licenses/>.
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

#include "bertini2/endgames/base_endgame.hpp"
#include "bertini2/trackers/adaptive_precision_utilities.hpp"  // tracking::adaptive::SetPrecision, for container migration


/**
\file bertini2/endgames/cauchy.hpp

\brief The concrete Cauchy endgame type

*/

namespace bertini{ namespace endgame{


/**
\class CauchyEndgame
\brief Class used to finish tracking paths during Homotopy Continuation.
\ingroup endgame

## Explanation
The bertini::CauchyEngame class enables us to finish tracking on possibly singular paths on an arbitrary square homotopy.
The intended usage is to:
1. Create a system, tracker, and instantiate some settings.
2. Using the tracker created track to the engame boundary.
3. Create a CauchyEndgame, associating it to the system you are going to solve or track on.
4. For each path being tracked send the CauchyEndgame the time value and other variable values at that time.
5. The CauchyEndgame, if successful, will store the target systems solution at $t = 0$.



## Example Usage
Below we demonstrate a basic usage of the CauchyEndgame class to find the singularity at $t = 0$.
The pattern is as described above: create an instance of the class, feeding it the system to be used, and the endgame boundary time and other variable values at the endgame boundary.
\code{.cpp}
using namespace bertini::tracking;
using RealT = tracking::TrackerTraits<TrackerType>::BaseRealT; // Real types
using ComplexT = tracking::TrackerTraits<TrackerType>::BaseComplexT; Complex types

// 1. Define the polynomial system that we wish to solve.
System target_sys;
Var x = Variable::Make("x"), t = Variable::Make("t"), y = Variable::Make("y");

VariableGroup vars{x,y};
target_sys.AddVariableGroup(vars);

target_sys.AddFunction((pow(x-1,3));
target_sys.AddFunction((pow(y-1,2));

// 1b. Homogenize and patch the polynomial system to work over projective space.
sys.Homogenize();
sys.AutoPatch();

// 2. Create a start system, for us we will use a total degree start system.
auto TD_start_sys = bertini::start_system::TotalDegreeLinearProduct(target_sys);

// 2b. Creating homotopy between the start system and system we wish to solve.
auto my_homotopy = (1-t)*target_sys + t*TD_start_sys*Rational::Rand(); //the random number is our gamma for a random path between t = 1 and t = 0.
my_homotopy.AddPathVariable(t);

//Sets up configuration settings for our particular system.
auto precision_config = PrecisionConfig(my_homotopy);


// 3. Creating a tracker. For us this is an AMPTracker.
AMPTracker tracker(my_homotopy);

//Tracker setup of settings.
SteppingConfig<RealT> stepping_preferences;
stepping_preferences.initial_step_size = RealT(1)/RealT(5);// change a stepping preference
NewtonConfig newton_preferences;
tracker.Setup(TestedPredictor,
            RealFromString("1e-6"),
            RealFromString("1e5"),
    stepping_preferences,
    newton_preferences);
tracker.PrecisionSetup(precision_config);

//We start at t = 1, and will stop at t = 0.1 before starting the endgames.
ComplexT t_start(1), t_endgame_boundary(0.1);

//This will hold our solutions at t = 0.1
std::vector<Vec<ComplexT> > my_homotopy_solutions_at_endgame_boundary;

// result holds the value we track to at 0.1, and tracking success will report if we are unsucessful.
Vec<ComplexT> result;

//4. Track all points to 0.1
for (unsigned ii = 0; ii < TD_start_sys.NumStartPoints(); ++ii)
{
    DefaultPrecision(ambient_precision);
    my_homotopy.precision(ambient_precision); // making sure our precision is all set up
    auto start_point = TD_start_sys.StartPoint<ComplexT>(ii);

    tracker.TrackPath(result,t_start,t_endgame_boundary,start_point);

    my_homotopy_solutions_at_endgame_boundary.push_back(result);
}


//Settings for the endgames.

config::Tolerances<RealT> tolerances;

CauchyConfig<RealT> cauchy_settings;
cauchy_settings.fail_safe_maximum_cycle_number = 6;


// 5. Create a cauchy endgame, and use them to get the soutions at t = 0.
EndgameSelector<TrackerType>::Cauchy my_cauchy_endgame(tracker,cauchy_settings,tolerances);


std::vector<Vec<ComplexT> > my_homotopy_solutions;

std::vector<Vec<ComplexT> > my_homotopy_divergent_paths;

for(auto s : my_homotopy_solutions_at_endgame_boundary)
{
    SuccessCode endgame_success = my_cauchy_endgame.Run(t_endgame_boundary,s);

    if(endgame_success == SuccessCode::Success)
    {
        my_homotopy_solutions.push_back(my_homotopy.DehomogenizePoint(my_endgame.FinalApproximation<ComplexT>()));
    }
    else
    {
        my_homotopy_divergent_paths.push_back(my_homotopy.DehomogenizePoint(my_endgame.FinalApproximation<ComplexT>()));
    }
}
\endcode


If this documentation is insufficient, please contact the authors with suggestions, or get involved!  Pull requests welcomed.

## Testing
Test suite driving this class: endgames_test.
File: test/endgames/generic_cauchy_test.hpp
File: test/endgames/amp_cauchy_test.cpp
File: test/endgames/fixed_double_cauchy_test.cpp
FIle: test/endgames/fixed_multiple_cauchy_test.cpp
*/
template<typename PrecT>
class CauchyEndgame :
	public virtual EndgameBase<CauchyEndgame<PrecT>, PrecT>
{
public:
	using BaseEGT = EndgameBase<CauchyEndgame<PrecT>, PrecT>;  ///< The base endgame type.
	using FinalEGT = CauchyEndgame<PrecT>;  ///< The final (derived) endgame type.
	using TrackerType = typename PrecT::TrackerType;  ///< The path-tracker type.

	using BaseComplexT = typename tracking::TrackerTraits<TrackerType>::BaseComplexT;  ///< The complex number type.
	using BaseRealT = typename tracking::TrackerTraits<TrackerType>::BaseRealT;  ///< The real number type.

	using EmitterType = CauchyEndgame<PrecT>;  ///< The event-emitter type.

protected:

	using EndgameBase<CauchyEndgame<PrecT>, PrecT>::NotifyObservers;




	using TupleOfTimes = typename BaseEGT::TupleOfTimes;  ///< A tuple of time containers, one per precision.
	using TupleOfSamps = typename BaseEGT::TupleOfSamps;  ///< A tuple of sample containers, one per precision.
	using TupOfVec = typename BaseEGT::TupOfVec;  ///< A tuple of vector containers, one per precision.

	using BCT = BaseComplexT;  ///< The complex number type.
	using BRT = BaseRealT;  ///< The real number type.

	using Configs = typename AlgoTraits<FinalEGT>::NeededConfigs;  ///< The configuration bundle (Configured base).
	using ConfigsAsTuple = typename Configs::ToTuple;  ///< The configuration structs as a tuple.

	/**
	\brief A deque of times that are specifically used to compute the power series approximation for the Cauchy endgame.
	*/
	mutable TupleOfTimes pseg_times_;
	/**
	\brief A deque of samples that are in correspondence with the pseg_times_. These samples are also used to compute the first power series approximation for the Cauchy endgame.
	*/
	mutable TupleOfSamps pseg_samples_; //samples used for the first approximation.
	/**
	\brief A deque of times that are collected by CircleTrack. These samples are used to compute all approximations of the origin after the first.
	*/
	mutable TupleOfTimes cauchy_times_;
	/**
	\brief A deque of samples collected by CircleTrack. Computed a mean of the values of this deque, after a loop has been closed, will give an approximation of the origin.
	*/
	mutable TupleOfSamps cauchy_samples_;

	/**
	\brief A fixed random probe vector used by ComputeCOverK to project sample differences to scalars.
	Generated ONCE (per precision) and reused across the whole endgame, so the c/k estimate is
	deterministic and consecutive estimates differ only because the samples differ -- not because the
	probe changed.  A fresh random probe every call made CheckForCOverKStabilization noisy (it could
	certify the operating zone spuriously) and churned mpfr allocations.  See z_notes/20260629.
	*/
	mutable TupOfVec c_over_k_probe_;

	// Scratch for LatestTimeImpl to return a BCT reference when the endgame is computing in the
	// complex_dbl fast lane (the latest time then lives in the complex_dbl slot, not the BCT slot).
	mutable BCT latest_time_cache_; ///< Scratch so LatestTimeImpl can return a BCT reference while in the double fast lane.





public:




	/**
	\brief Function that clears all samples and times from data members for the Cauchy endgame
	*/
	template<typename ComplexT>
	void ClearTimesAndSamples()
	{
		std::get<TimeCont<ComplexT> >(pseg_times_).clear();
		std::get<TimeCont<ComplexT> >(cauchy_times_).clear();
		std::get<SampCont<ComplexT> >(pseg_samples_).clear();
		std::get<SampCont<ComplexT> >(cauchy_samples_).clear();}
	/**
	\brief Setter for the time values for the power series approximation of the Cauchy endgame.
	*/
	template<typename ComplexT>
	void SetPSEGTimes(TimeCont<ComplexT> pseg_times_to_set)
	{ std::get<TimeCont<ComplexT> >(pseg_times_) = pseg_times_to_set;}

	/**
	\brief Getter for the time values for the power series approximation of the Cauchy endgame.
	*/
	template<typename ComplexT>
	TimeCont<ComplexT>& GetPSEGTimes() {return std::get<TimeCont<ComplexT> >(pseg_times_);}
	/// \brief Const overload returning the power-series time values.
	template<typename ComplexT>
	const TimeCont<ComplexT>& GetPSEGTimes() const {return std::get<TimeCont<ComplexT> >(pseg_times_);}

	/**
	\brief Setter for the space values for the power series approximation of the Cauchy endgame.
	*/
	template<typename ComplexT>
	void SetPSEGSamples(SampCont<ComplexT> const& pseg_samples_to_set) { std::get<SampCont<ComplexT> >(pseg_samples_) = pseg_samples_to_set;}

	/**
	\brief Getter for the space values for the power series approximation of the Cauchy endgame.

	Available in const and non-const flavors
	*/
	template<typename ComplexT>
	SampCont<ComplexT>& GetPSEGSamples() {return std::get<SampCont<ComplexT> >(pseg_samples_);}
	/// \brief Const overload returning the power-series sample values.
	template<typename ComplexT>
	const SampCont<ComplexT>& GetPSEGSamples() const {return std::get<SampCont<ComplexT> >(pseg_samples_);}
	/**
	\brief Setter for the space values for the Cauchy endgame.
	*/
	template<typename ComplexT>
	void SetCauchySamples(SampCont<ComplexT> const& cauchy_samples_to_set)
	{
		std::get<SampCont<ComplexT> >(cauchy_samples_) = cauchy_samples_to_set;
	}

	/**
	\brief Getter for the space values for the Cauchy endgame.

	Available in const and non-const flavors
	*/
	template<typename ComplexT>
	SampCont<ComplexT>& GetCauchySamples()
	{
		return std::get<SampCont<ComplexT> >(cauchy_samples_);
	}
	/// \brief Const overload returning the Cauchy sample values.
	template<typename ComplexT>
	const SampCont<ComplexT>& GetCauchySamples() const { return std::get<SampCont<ComplexT> >(cauchy_samples_); }


	/**
	\brief Setter for the time values for the Cauchy endgame.
	*/
	template<typename ComplexT>
	void SetCauchyTimes(TimeCont<ComplexT> const& cauchy_times_to_set)
	{
		std::get<TimeCont<ComplexT> >(cauchy_times_) = cauchy_times_to_set;
	}

	/**
	\brief Getter for the time values for the Cauchy endgame.
	*/
	template<typename ComplexT>
	TimeCont<ComplexT>& GetCauchyTimes()
	{
		return std::get<TimeCont<ComplexT> >(cauchy_times_);
	}
	/// \brief Const overload returning the Cauchy time values.
	template<typename ComplexT>
	const TimeCont<ComplexT>& GetCauchyTimes() const
	{
		return std::get<TimeCont<ComplexT> >(cauchy_times_);
	}


	/// \return The most recent power-series time value.
	const BCT& LatestTimeImpl() const
	{
		// In the adaptive-numeric-type endgame the latest time may live in the complex_dbl slot (the
		// fast lane), with the BCT slot empty.  Dispatch on which slot actually holds data, so this is
		// correct both during the run (observer events) and after it (solution metadata), regardless of
		// the adaptive_numeric_type_active_ flag.  Fixed precision compiles to the original BCT read.
		if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			if (GetPSEGTimes<BCT>().empty())
			{
				auto const& dbl_times = GetPSEGTimes<complex_dbl>();
				latest_time_cache_ = dbl_times.empty() ? BCT(0) : BCT(dbl_times.back());
				return latest_time_cache_;
			}
		}
		return GetPSEGTimes<BCT>().back();
	}


	/**
	\brief Setter for the specific settings in tracking_conifg.hpp under Cauchy.
	*/
	void SetCauchySettings(CauchyConfig const& new_cauchy_settings)
	{
		this->template Set<CauchyConfig>(new_cauchy_settings);
	}

	/**
	\brief Getter for the specific settings in tracking_conifg.hpp under Cauchy.
	*/
	const auto& GetCauchySettings() const
	{
		return this->template Get<CauchyConfig>();
	}


	/// \brief Construct the Cauchy endgame for a tracker, with its configuration as a tuple.
	explicit CauchyEndgame(TrackerType const& tr,
                            const ConfigsAsTuple& settings )
      : EndgamePrecPolicyBase<TrackerType>(tr), BaseEGT(tr, settings)
   	{ }

	/// \brief Construct the Cauchy endgame for a tracker, with configs given in any order.
    template< typename... Ts >
		CauchyEndgame(TrackerType const& tr, const Ts&... ts ) : CauchyEndgame(tr, Configs::Unpermute( ts... ) )
		{}


	virtual ~CauchyEndgame() = default;

	/// \brief Validate the endgame configuration (e.g. require >= 3 sample points for circle tracking).
	void ValidateConfigs()
	{
		if (this->EndgameSettings().num_sample_points < 3) // need to make sure we won't track right through the origin.
		{
			std::stringstream err_msg;
			err_msg << "ERROR: The number of sample points " << this->EndgameSettings().num_sample_points << " for circle tracking must be >= 3";
			throw std::runtime_error(err_msg.str());
		}
	}

	/**
		\brief Function to track around the origin

		## Input:
				starting_time: time value that we start from to track around the origin
				target_time: the time value that we are centering out loops around, default is t = 0
				starting_sample: an approximate solution of the homotopy at t = starting_time

		## Output:
				SuccessCode: This reports back if we were successful in advancing time.


		##Details:
	\tparam ComplexT The complex number type.
				Depeding on the number of samples points, we make a polgon around the origin with that many vertices. This function should be called the same number of times
				as paths converging to the solution we are approximating.
	*/
	template<typename ComplexT>
	SuccessCode CircleTrack(ComplexT const& target_time)
	{
		using bertini::Precision;
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;
		using std::acos;

		ValidateConfigs();

		auto& circle_times = std::get<TimeCont<ComplexT> >(cauchy_times_);
		auto& circle_samples = std::get<SampCont<ComplexT> >(cauchy_samples_);

		ComplexT starting_time = circle_times.back();  // take a COPY here, so won't invalidate it later

		// the initial sample has already been added to the sample repo... so don't do that here, please

		const auto num_vars = this->GetSystem().NumVariables();

		for (unsigned ii = 0; ii < this->EndgameSettings().num_sample_points; ++ii)
		{
			const Vec<ComplexT>& current_sample = circle_samples.back();
			const ComplexT& current_time = circle_times.back();

#ifndef BERTINI_DISABLE_PRECISION_CHECKS
			if (Precision(current_time)!=Precision(current_sample)){
				std::stringstream err_msg;
				err_msg << "current time and sample for circle track must be of same precision.  respective precisions: " << Precision(current_time) << " " << Precision(current_sample) << std::endl;
				throw std::runtime_error(err_msg.str());
			}
#endif

			//set up the time value for the next sample.
			using std::polar;

#ifndef USE_BMP_COMPLEX
			using bertini::polar;
#endif

			//Generalized since we could have a nonzero target time.
			using std::arg;
			RealT radius = abs(starting_time - target_time), angle = arg(starting_time - target_time); // generalized for nonzero target_time.

			auto next_sample = Vec<ComplexT>(num_vars);
			ComplexT next_time = (ii==this->EndgameSettings().num_sample_points-1)
								?
							  starting_time
								:
							  polar(radius, (ii+1)*2*acos(static_cast<RealT>(-1)) / (this->EndgameSettings().num_sample_points) + angle) + target_time;
			// If we are tracking to a nonzero target time we need to shift our values to track to. This is a step that may not be needed if target_time = 0
							  ;


			auto tracking_success = this->EndgameTrackPath(next_sample, current_time, next_time, current_sample);
			if (tracking_success != SuccessCode::Success)
			{
				return tracking_success;
			}

			// Pure-(i) numeric-type escalation: if the tracker's authoritative precision climbed above
			// the precision this endgame is computing in, double no longer suffices for this circle.
			// Bail to the migrate-and-retry driver (RunImplAMP), which discards this partial circle,
			// migrates the durable state to mpfr, and re-tracks the circle in mpfr.  Compile-time elided
			// for fixed precision (and shields GetCurrentPrecision(), which fixed trackers lack).
			if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				if (this->adaptive_numeric_type_active_ &&
				    this->GetTracker().GetCurrentPrecision() > this->current_endgame_precision_)
					return SuccessCode::HigherPrecisionNecessary;

			// CircleAdvanced carries the new point/time at BaseComplexT.  The fixed/mpfr lane emits exactly
			// as before.  The complex_dbl fast lane would have to convert the point to mpfr for the event,
			// so we only pay that when something is actually observing (temporaries live through the
			// synchronous NotifyObservers).
			if constexpr (std::is_same<ComplexT, BCT>::value)
				NotifyObservers(CircleAdvanced<EmitterType>(*this, next_sample, next_time));
			else if (this->HasObservers())
			{
				Vec<BCT> ev_pt(next_sample.size());
				for (Eigen::Index i = 0; i < next_sample.size(); ++i) ev_pt(i) = BCT(next_sample(i));
				BCT ev_t(next_time);
				NotifyObservers(CircleAdvanced<EmitterType>(*this, ev_pt, ev_t));
			}


			this->EnsureAtPrecision(next_time,Precision(next_sample)); assert(Precision(next_time)==Precision(next_sample));

			// auto refinement_success = this->RefineSample(next_sample, next_sample, next_time,
			// 							this->FinalTolerance() * this->EndgameSettings().sample_point_refinement_factor,
			// 							this->EndgameSettings().max_num_refinements);
			// if (refinement_success != SuccessCode::Success)
			// {
			// 	return refinement_success;
			// }

			// this->EnsureAtPrecision(next_time,Precision(next_sample));
			// assert(Precision(next_time)==Precision(next_sample));

			AddToCauchyData(next_time, next_sample);

			// NotifyObservers(SampleRefined<EmitterType>(*this));
		}

		return SuccessCode::Success;

	}//end CircleTrack

	/// \brief Append a (time, sample) pair to the Cauchy data containers.
	template<typename ComplexT>
	void AddToCauchyData(ComplexT const& time, Vec<ComplexT> const& sample)
	{
		std::get<TimeCont<ComplexT>>(cauchy_times_).push_back(time);
		std::get<SampCont<ComplexT>>(cauchy_samples_).push_back(sample);
	}

	/// \brief Append a (time, sample) pair to the power-series data containers.
	template<typename ComplexT>
	void AddToPSData(ComplexT const& time, Vec<ComplexT> const& sample)
	{
		std::get<TimeCont<ComplexT>>(pseg_times_).push_back(time);
		std::get<SampCont<ComplexT>>(pseg_samples_).push_back(sample);
	}

	/**
		\brief Lazily generate (once) and return the fixed random probe vector used by ComputeCOverK.
		The probe is generated the first time it is needed at a given size, then reused for the life of
		the endgame so the c/k estimate is deterministic.  For adaptive precision it is re-precisioned
		in place to match the working samples (the random direction is preserved).
	*/
	template<typename ComplexT>
	Vec<ComplexT> const& GetCOverKProbe(unsigned size, unsigned prec) const
	{
		using bertini::Precision;
		auto& probe = std::get<Vec<ComplexT> >(c_over_k_probe_);
		if (static_cast<unsigned>(probe.size()) != size)
		{
			probe.resize(size);
			for (unsigned ii = 0; ii < size; ++ii)
				probe(ii) = RandomUnit<ComplexT>();
		}
		if (Precision(probe) != prec)
			Precision(probe, prec);
		return probe;
	}

	/**
		\brief A function that uses the assumption of being in the endgame operating zone to compute an approximation of the ratio c over k.
			When the cycle number stabilizes we will see that the different approximations of c over k will stabilize.
			Returns the computed value of c over k.


		## Input:
				None: all data needed are class data members.

		## Output:
				estimate: The approximation of the ratio of the two numbers C and K heuristically signifying we are in the cauchy endgame operating zone.


		##Details:
				\tparam ComplexT The complex number type.
				Consult page 53 of \cite bertinibook, for the reasoning behind this heuristic.
	*/
	template<typename ComplexT>
	auto ComputeCOverK() const -> typename Eigen::NumTraits<ComplexT>::Real
	{//Obtain samples for computing C over K.
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;
		using std::abs;
		using std::log;

		const auto& pseg_samples = std::get<SampCont<ComplexT> >(pseg_samples_);

		assert(pseg_samples.size()>=3);
		const Vec<ComplexT> & sample0 = pseg_samples[0];
		const Vec<ComplexT> & sample1 = pseg_samples[1];
		const Vec<ComplexT> & sample2 = pseg_samples[2];

		// Use a fixed random probe vector, generated once and reused across the whole endgame, so this
		// estimate is deterministic.  A fresh random vector per call made consecutive c/k estimates
		// disagree by probe noise alone, which could trip (or stall) CheckForCOverKStabilization.
		const Vec<ComplexT> & rand_vector = GetCOverKProbe<ComplexT>(static_cast<unsigned>(sample0.size()), Precision(sample0));

		// //DO NOT USE Eigen .dot() it will do conjugate transpose which is not what we want.
		// //Also, the .transpose*rand_vector returns an expression template that we do .norm of since abs is not available for that expression type.
		RealT estimate = abs(log(abs((((sample2 - sample1).transpose()*rand_vector).template lpNorm<Eigen::Infinity>())/(((sample1 - sample0).transpose()*rand_vector).template lpNorm<Eigen::Infinity>()))));
		estimate = abs(log(RealT(this->EndgameSettings().sample_factor)))/estimate;
		if (estimate < 1)
		  	return RealT(1);
		else
			return estimate;

	}//end ComputeCOverK


	/**
		\brief Function to determine if ratios of c/k estimates are withing a user defined threshold.

		## Input:
				c_over_k_array: A container holding all previous computed C over K ratios. The stabilization of these ratios is key to the convergence of the cauchy endgame.

		## Output:
				true: if we have stabilized and can proceed with the endgame.
				false: if our ratios are not withing tolerances set by the user or by default.

		##Details:
				\tparam ComplexT The complex number type.

	*/
	template<typename ComplexT>
	bool CheckForCOverKStabilization(TimeCont<ComplexT> const& c_over_k_array) const
	{
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;
		using std::abs;

		assert(c_over_k_array.size()>=GetCauchySettings().num_needed_for_stabilization);
		for(unsigned ii = 1; ii < GetCauchySettings().num_needed_for_stabilization ; ++ii)
		{
			RealT a = abs(c_over_k_array[ii-1]);
			RealT b = abs(c_over_k_array[ii]);

			typename Eigen::NumTraits<ComplexT>::Real divide = a;

			if(a < b)
				divide = a/b;
			else
				divide = b/a;

			if(divide <  GetCauchySettings().minimum_for_c_over_k_stabilization)
				return false;
		}
		return true;

	}//end CheckForCOverKStabilization


	/*
	Input: A time value and the space value above that time.

	Output: An real_mp representing a tolerance threshold for declaring a loop to be closed.
	Details: Used in Bertini 1 as a heuristic for computing separatedness of roots. Decided to not be used since assumptions for this tolerance are not usually met.
	template<typename ComplexT>
	real_mp FindToleranceForClosedLoop(ComplexT x_time, Vec<ComplexT> x_sample)
	{
		auto degree_max = std::max(this->GetTracker().AMP_config_.degree_bound,real_mp("2.0"));
		auto K = this->GetTracker().AMP_config_.coefficient_bound;
		real_mp N;
		real_mp M;
		real_mp L;
		if(max_closed_loop_tolerance_ < min_closed_loop_tolerance_)
		{
			max_closed_loop_tolerance_ = min_closed_loop_tolerance_;
		}
		auto error_tolerance = real_mp("1e-13");
		if(x_sample.size() <= 1)
		{
			N = degree_max;
		}
		else
		{
			N = ComputeCombination(degree_max + x_sample[0].precision() - 1, x_sample[0].precision() - 1);
		}
		M = degree_max * (degree_max - 1) * N;
		auto jacobian_at_current_time = this->GetSystem().Jacobian(x_sample,x_time);
		auto minimum_singular_value = Eigen::JacobiSVD< Mat<ComplexT> >(jacobian_at_current_time).singularValues()(this->GetSystem().NumVariables() - 1 );
		auto norm_of_sample = x_sample.norm();
		L = pow(norm_of_sample,degree_max - 2);
		auto tol = K * L * M;
		if (tol == 0) // fail-safe
					tol = minimum_singular_value;
			else
			{
			tol = real_mp("2.0") / tol;
		tol = tol * minimum_singular_value;
			}
			// make sure that tol is between min_tol & max_tol
		if (tol > max_closed_loop_tolerance_) // tol needs to be <= max_tol
		    tol = max_closed_loop_tolerance_;
		if (tol < min_closed_loop_tolerance_) // tol needs to be >= min_tol
			tol = min_closed_loop_tolerance_;
			return tol;
	}// end FindToleranceForClosedLoop
	*/



	/**
		\brief Function that determines if we have closed a loop after calling CircleTrack().


		## Input:
				None: all data needed are class data members

		## Output:
				true: if we have closed the loop
				false: if we have not closed the loop

		##Details:
				\tparam ComplexT The complex number type
	*/
	template<typename ComplexT>
	bool CheckClosedLoop()
	{
		auto& times = std::get<TimeCont<ComplexT> >(cauchy_times_);
		auto& samples = std::get<SampCont<ComplexT> >(cauchy_samples_);

		if((samples.front() - samples.back()).template lpNorm<Eigen::Infinity>() < this->GetTracker().TrackingTolerance())
		{
			return true;
		}

		if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			//Ensure all samples are of the same precision.
			auto new_precision = this->EnsureAtUniformPrecision(times, samples);
			this->GetSystem().precision(new_precision);
		}

		this->GetTracker().Refine(samples.front(),samples.front(),times.front(),this->FinalTolerance(),this->EndgameSettings().max_num_refinements);
		this->GetTracker().Refine(samples.back(),samples.back(),times.back(),this->FinalTolerance(),this->EndgameSettings().max_num_refinements);

		if((samples.front() - samples.back()).template lpNorm<Eigen::Infinity>() < this->GetTracker().TrackingTolerance())
		{
			return true;
		}
		return false;

	}//end CheckClosedLoop




	/**
		\brief 	After we have used CircleTrack and have successfully closed the loop using CheckClosedLoop we need to check the maximum and minimum norms of the samples collected.
				If the ratio of the maximum and minimum norm are within the threshold maximum_cauchy_ratio, and the difference is greater than final tolerance than we are successful.


		## Input:
			target_time: Used to stay correct if we are using the endgame for a non-zero target time.

		## Output:
			true: If we are within the ratio cutoff time, or have a ratio within the thresholds.
			false: otherwise


		##Details:
				\tparam ComplexT The complex number type.
				It is important to know if we are within the endgame operating zone. This function allows us to have a check that
				heuristcially will tell us if we are.
	*/
	template<typename ComplexT>
	bool RatioEGOperatingZoneTest(ComplexT const& target_time) const
	{
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;
		RealT min(1e300);
		RealT max(0);
		auto& times = std::get<TimeCont<ComplexT> >(cauchy_times_);
		auto& samples = std::get<SampCont<ComplexT> >(cauchy_samples_);
		if(norm(times.front() - target_time) < GetCauchySettings().ratio_cutoff_time)
		{
			return true;
		}
		else
		{
			RealT norm;
			for(unsigned int ii=0; ii < this->EndgameSettings().num_sample_points; ++ii)
			{
				norm = samples[ii].template lpNorm<Eigen::Infinity>();
				if(norm > max)
				{
					max = norm;
				}
				if(norm < min)
				{
					min = norm;
				}
			}

			if(min > this->FinalTolerance() && max > this->FinalTolerance())
			{
				norm = min / max;
				if(norm < GetCauchySettings().maximum_cauchy_ratio && (max - min) > this->FinalTolerance())
				{
					return false; // bad ratio and too far apart.
				}
			}
		}
		return true;

	}//end RatioEGOperatingZoneTest


	/*
	Input: All input needed is available as class data members.
	Output: SuccessCode declaring if we are successful or not.
	Details: Computes cauchy samples using Circle Track. Compares the ratios of the maximum and minimum norms of cauchy samples
			 using RatioEGOperatingZoneTest. Then attempts to see if we have closed our loop, otherwise it continues the process outlined.
	*/



	/**
		\brief This function tracks into origin while computing loops around the origin. The function is checking to make sure we have reached a time when the ratio of the
	maximum and minimum norms are withing some tolerance. When this is the case we return success. If this does not happen we will return an error depending on the error
	encountered.

		## Input:
			target_time: the time value that we are creating loops around, default set to t= 0

		## Output:
			initial_cauchy_loop_success: This variable returns any error code we have encountered or success if we have sucessfully
			tracked to an appropriate time.

		##Details:
				\tparam ComplexT The complex number type.
	*/
	template<typename ComplexT>
	SuccessCode InitialCauchyLoops(ComplexT const& target_time)
	{
		using std::max;
		// auto fail_safe_max_cycle_number = max(GetCauchySettings().fail_safe_maximum_cycle_number,this->CycleNumber());
		auto fail_safe_max_cycle_number = GetCauchySettings().fail_safe_maximum_cycle_number;

		auto initial_cauchy_loop_success = SuccessCode::Success;

		bool loop_hasnt_closed = true;
		while (loop_hasnt_closed)
		{
			this->CycleNumber(0);
			ClearAndSeedCauchyData<ComplexT>();

			// track around a circle once.  we'll use it to measure whether we believe we are in the eg operating zone, based on the ratio of norms of sample points around the circle
			auto tracking_success = CircleTrack(target_time);
			this->IncrementCycleNumber(1);
			if (tracking_success != SuccessCode::Success)
				return tracking_success;

			// find the ratio of the maximum and minimum coordinate wise for the loop.
			if (RatioEGOperatingZoneTest<ComplexT>(target_time))
			{ // then we believe we are in the EG operating zone, since the path is relatively flat.  i still disbelieve this is a good test (dab 20160310)
				while (true)
				{
					if (CheckClosedLoop<ComplexT>())
					{//error is small enough, exit the loop with success.
						NotifyObservers(ClosedLoop<EmitterType>(*this));
						initial_cauchy_loop_success = SuccessCode::Success;
						loop_hasnt_closed = false;
						break;
					}
					else if(this->CycleNumber() > fail_safe_max_cycle_number)
					{// too many iterations
						initial_cauchy_loop_success = SuccessCode::CycleNumTooHigh;
						loop_hasnt_closed = false;
						break;
					}

					//compute next loop, the last sample in times and samples is the sample our loop ended on. Either where we started or on another sheet at the same time value.
					tracking_success = CircleTrack(target_time);
					this->IncrementCycleNumber(1);
					if(tracking_success != SuccessCode::Success)
						return tracking_success;
				}
			}//end if (RatioEGOperatingZoneTest())
			else
			{
				auto advance_success = AdvanceTime<ComplexT>(target_time);
				if (advance_success!=SuccessCode::Success)
					return advance_success;
			}
		} //end while(loop_hasnt_closed)

		return initial_cauchy_loop_success;
	}//end InitialCauchyLoops



	/// \brief Shift the power-series sample window forward by one (time, sample) pair.
	template <typename ComplexT>
	void RotateOntoPS(ComplexT const& next_time, Vec<ComplexT> const& next_sample)
	{
		auto& ps_times = std::get<TimeCont<ComplexT> >(pseg_times_);
		auto& ps_samples = std::get<SampCont<ComplexT> >(pseg_samples_);

		ps_times.pop_front();
		ps_samples.pop_front();

		ps_times.push_back(next_time);
		ps_samples.push_back(next_sample);
	}

	/// \brief Clear the Cauchy data and re-seed it from the current power-series samples.
	template <typename ComplexT>
	void ClearAndSeedCauchyData()
	{
		auto& cau_times = std::get<TimeCont<ComplexT> >(cauchy_times_);
		auto& cau_samples = std::get<SampCont<ComplexT> >(cauchy_samples_);
		auto& ps_times = std::get<TimeCont<ComplexT> >(pseg_times_);
		auto& ps_samples = std::get<SampCont<ComplexT> >(pseg_samples_);

		cau_times.clear();
		cau_samples.clear();
		cau_samples.push_back(ps_samples.back());
		cau_times.push_back(ps_times.back());
	}


	/**
		\brief 	Tracks til we believe we are in the Endgame Operating Zone, and then does a cauchy approximation

		## Input:
				start_time: time value for which we start to make a power series approximation
				start_point: approximate solution to our homotopy H at the start_time
				approximation_time: time at which we are trying to find the solution, usually t = 0
				approximation approximate solution to our homotopy H at the approxmation_time


		## Output:
			SuccessCode reporting if any errors had occurred. All data collected is stored in class data members.


		##Details:
	\tparam ComplexT The complex number type.

	This function is in charge of finding the very first approximation of the origin. It does this by first computing some initial samples
	like what is done in the Power Series Endgame. We continue to track forward in this manner until we have stabilization of the cycle number being approximated.
	This prevents the unnecessary circle tracking if we are possibly not in the endgame operating zone.
	Once we have stabilization we then perform InitialCauchyLoops while getting the accurate cycle number, and check the norms of the samples and make sure we are ready
	to approximate.
	*/
	template<typename ComplexT>
	SuccessCode InitialApproximation(ComplexT const& start_time, Vec<ComplexT> const& start_point,
	                                            ComplexT const& target_time, Vec<ComplexT> & approximation)
	{
		auto init_success = GetIntoEGZone(start_time, start_point, target_time);
		if (init_success!= SuccessCode::Success)
			return init_success;

		auto cauchy_loop_success = InitialCauchyLoops<ComplexT>(target_time);
		if (cauchy_loop_success != SuccessCode::Success)
			return cauchy_loop_success;

		return ComputeCauchyApproximationOfXAtT0(approximation);

	}//end InitialApproximation



	/// \brief Track from the boundary inward until the path enters the endgame (Cauchy) operating zone.
	template<typename ComplexT>
	SuccessCode GetIntoEGZone(ComplexT const& start_time, Vec<ComplexT> const& start_point, ComplexT const& target_time)
	{
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;

		//initialize array holding c_over_k estimates
		std::deque<RealT> c_over_k;

		auto& ps_times = std::get<TimeCont<ComplexT> >(pseg_times_);
		auto& ps_samples = std::get<SampCont<ComplexT> >(pseg_samples_);

		//Compute initial samples for pseg
		auto initial_sample_success = this->ComputeInitialSamples(start_time, target_time, start_point, ps_times, ps_samples);
		if (initial_sample_success!=SuccessCode::Success)
			return initial_sample_success;

		c_over_k.push_back(ComputeCOverK<ComplexT>());


		//track until for more c_over_k estimates or until we reach a cutoff time.
		for (unsigned ii = 0; ii < GetCauchySettings().num_needed_for_stabilization; ++ii)
		{
			auto advance_success = AdvanceTime<ComplexT>(target_time);
			if (advance_success!=SuccessCode::Success)
				return advance_success;
			c_over_k.push_back(ComputeCOverK<ComplexT>());
		}//end while


		//have we stabilized yet?
		while(!CheckForCOverKStabilization(c_over_k) && abs(ps_times.back()-target_time) > GetCauchySettings().cycle_cutoff_time)
		{
			auto advance_success = AdvanceTime<ComplexT>(target_time);
			if (advance_success!=SuccessCode::Success)
				return advance_success;

			c_over_k.pop_front();
			c_over_k.push_back(ComputeCOverK<ComplexT>());

		}//end while

		NotifyObservers(InEGOperatingZone<EmitterType>(*this));

		return SuccessCode::Success;
	}

	/**
	\brief Function that computes the mean of the samples that were collected while tracking around the origin. This value is the approximation of the value at the origin.

		## Input:
			result: This vector, passed by reference, holds the approximation that we calculate.

		## Output:
			SuccessCode deeming if we were suceessful, or if we encountered an error.

		##Details:
	\tparam ComplexT The complex number type.
				We can compute the Cauchy Integral Formula in this particular instance by computing the mean of the samples we have collected around the origin.

				/todo i believe this function works incorrectly when the target time is not 0.  hence, the target time needs to be passed in.
	*/
	template<typename ComplexT>
	SuccessCode ComputeCauchyApproximationOfXAtT0(Vec<ComplexT>& result)
	{
		auto& cau_times = std::get<TimeCont<ComplexT> >(cauchy_times_);
		auto& cau_samples = std::get<SampCont<ComplexT> >(cauchy_samples_);

		if (cau_samples.size() != this->CycleNumber() * this->EndgameSettings().num_sample_points+1)
		{
			std::stringstream err_msg;
			err_msg << "to compute cauchy approximation, cau_samples must be of size " << this->CycleNumber() * this->EndgameSettings().num_sample_points+1 << " but is of size " << cau_samples.size() << '\n';
			throw std::runtime_error(err_msg.str());
		}


		//Ensure all samples are of the same precision.
		if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			auto new_precision = this->EnsureAtUniformPrecision(cau_times, cau_samples);
			this->GetSystem().precision(new_precision);
		}


		auto total_num_pts = this->CycleNumber() * this->EndgameSettings().num_sample_points;
		auto refine_code = this->template RefineAllSamples<ComplexT>(cau_samples, cau_times);
		// Pure-(i): when the adaptive-numeric-type driver is orchestrating, a refine that double cannot
		// satisfy is a request to cross to mpfr -- propagate it so RunImplAMP migrates and retries.
		// Fixed precision (and the AMP-PowerSeries path) keep ignoring the code, exactly as before.
		if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
			if (this->adaptive_numeric_type_active_ &&
			    (refine_code == SuccessCode::HigherPrecisionNecessary || refine_code == SuccessCode::FailedToConverge))
				return refine_code;

		Precision(result, Precision(cau_samples.back()));

		result = Vec<ComplexT>::Zero(static_cast<Eigen::Index>(this->GetSystem().NumVariables()));
		for(unsigned int ii = 0; ii < total_num_pts; ++ii)
			result += cau_samples[ii];
		result /= this->CycleNumber() * this->EndgameSettings().num_sample_points;

		return SuccessCode::Success;

	}

	/**
	\brief Collects samples while tracking around the target time, until we close the loop, or exceed the limit on # of loops.

		## Input:
			the target time, at which we are computing roots of a system


		## Output:
			SuccessCode deeming if we were able to collect all samples around the origin, or if we encounted an error at some point.

		##Details:

			the starting time and point for this routine are the most recent power series samples.
	\tparam ComplexT The complex number type.
				This function populates the deque cauchy_samples and cauchy_times. These are data members of the class and are not passed in. This function will continue to
				call CircleTrack until we have closed the loop.

	*/
	template<typename ComplexT>
	SuccessCode ComputeCauchySamples(ComplexT const& target_time)
	{
		using bertini::Precision;

		ClearAndSeedCauchyData<ComplexT>();
		this->CycleNumber(0);


		while( this->CycleNumber() < GetCauchySettings().fail_safe_maximum_cycle_number )
		{
			//track around the origin once.
			auto tracking_success = CircleTrack(target_time);
			this->IncrementCycleNumber(1);

			if(tracking_success != SuccessCode::Success)
			{
				return tracking_success;
			}
			else if(CheckClosedLoop<ComplexT>())
			{
				return SuccessCode::Success;
			}
		}
		NotifyObservers(CycleNumTooHigh<EmitterType>(*this));
		return SuccessCode::CycleNumTooHigh;
	}//end ComputeCauchySamples


	/**
	\brief Advances time, marching toward the target time.

	Works from the most recent time-sample pair stored in the power series data.

	If the distance between next and target is too small, dies (returns not success).
	*/
	template<typename ComplexT>
	SuccessCode AdvanceTime(ComplexT const& target_time)
	{
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;

		auto& ps_times = std::get<TimeCont<ComplexT> >(pseg_times_);
		auto& ps_samples = std::get<SampCont<ComplexT> >(pseg_samples_);

		auto& current_time = ps_times.back();
		auto& current_sample = ps_samples.back();

		//Generalized next_time in case if we are not trying to converge to the t = 0.
		ComplexT next_time = (target_time-current_time) * static_cast<RealT>(this->EndgameSettings().sample_factor)+current_time;

		if (abs(next_time - target_time) < this->EndgameSettings().min_track_time)//we are too close to t = 0 but we do not have the correct tolerance - so we exit
			return SuccessCode::MinTrackTimeReached;

		// advance in time
		Vec<ComplexT> next_sample;
		auto time_advance_success = this->EndgameTrackPath(next_sample,current_time, next_time, current_sample);
		if (time_advance_success != SuccessCode::Success)
		{
			NotifyObservers(EndgameFailure<EmitterType>(*this));
			return time_advance_success;
		}

		// Pure-(i) escalation: return BEFORE RotateOntoPS so the PSEG window stays an untouched
		// checkpoint -- RunImplAMP migrates it to mpfr and retries this advance.  Elided for fixed prec.
		if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
			if (this->adaptive_numeric_type_active_ &&
			    this->GetTracker().GetCurrentPrecision() > this->current_endgame_precision_)
				return SuccessCode::HigherPrecisionNecessary;

		this->EnsureAtPrecision(next_time,Precision(next_sample));
		RotateOntoPS(next_time, next_sample);

		NotifyObservers(TimeAdvanced<EmitterType>(*this));
		return SuccessCode::Success;
	}

	/**
	\brief The pole-component mass of the current Cauchy loop: the endgame
	operating-zone measurement.

	The loop samples are a discrete Fourier series of the Puiseux/Laurent expansion of
	the path around the target time.  Every CLEAN Puiseux term t^{k/c} with k >= 0
	contributes exactly zero to the twisted means below (roots-of-unity
	orthogonality), so any surviving mass in the negative modes k = -1..-c measures
	structure the operating zone forbids:

	- a pole AT the target time (an affine path to genuine infinity -- there is no
	  patch to give it a finite place to go), or
	- any other singularity of the ramified cover INSIDE the loop (a branch point
	  t* != 0 with |t*| < r), where the Cauchy integral assumption fails and the
	  plain mean is garbage.

	The two cases separate as the radius shrinks: an inside singularity's mass dies
	once r < |t*| (then the endgame may proceed -- the zone finally reached); a pole
	at the target time has mass GROWING like 1/r, so it never reaches the operating
	zone and never accepts convergence (the acceptance gate in RunImpl).

	Uses exactly the samples the mean uses: the closed c-circuit loop, uniform in
	angle, the duplicate closing sample excluded.  Coordinates only -- no function
	values, so system scaling cannot affect the verdict.

	\tparam ComplexT The complex number type of the samples.
	\return max over m = 1..c of the infinity norm of the e^{i m theta / c}-twisted
	sample mean (the estimated magnitude of the t^{-m/c} coefficient at the current
	radius).
	*/
	template<typename ComplexT>
	NumErrorT PoleComponentMass() const
	{
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;
		auto const& samples = std::get<SampCont<ComplexT>>(cauchy_samples_);
		auto const& times   = std::get<TimeCont<ComplexT>>(cauchy_times_);
		auto const c = this->CycleNumber();
		auto const N = this->EndgameSettings().num_sample_points;
		auto const M = c * N;
		if (c == 0 || samples.size() < M + 1 || times.size() < 2)
			return NumErrorT(0);

		// orientation of the loop (which way theta advances); the frequency sign follows it
		RealT const orientation_test = imag(times[1] * conj(times[0]));
		RealT const orientation = (orientation_test > RealT(0)) ? RealT(1) : RealT(-1);

		RealT const two_pi = RealT(2) * acos(RealT(-1));
		NumErrorT mass(0);
		for (unsigned m = 1; m <= c; ++m)
		{
			Vec<ComplexT> twisted = Vec<ComplexT>::Zero(samples[0].size());
			for (unsigned j = 0; j < M; ++j)
			{
				// weight e^{+i m theta_j / c} with theta_j = orientation * 2 pi j / N:
				// picks out exactly the k = -m Puiseux mode
				RealT const angle = orientation * two_pi * RealT(j * m) / RealT(N * c);
				twisted += samples[j] * ComplexT(cos(angle), sin(angle));
			}
			twisted /= RealT(M);
			auto const twisted_norm = static_cast<NumErrorT>(twisted.template lpNorm<Eigen::Infinity>());
			if (twisted_norm > mass)
				mass = twisted_norm;
		}
		return mass;
	}

	/// \brief PoleComponentMass read from whichever numeric lane the adaptive
	/// (double-first) endgame is currently computing in.
	template<typename Dummy = void>
	NumErrorT PoleComponentMassAMP() const
	{
		if (this->current_endgame_precision_ == DoublePrecision())
			return PoleComponentMass<complex_dbl>();
		return PoleComponentMass<complex_mp>();
	}

	/**
	\brief Is a pole-component mass SIGNIFICANT, i.e. above the noise floor of the
	refined samples, relative to the approximation's own scale?

	\param mass The pole-component mass (PoleComponentMass).
	\param approx_scale The infinity norm of the current approximation.
	\return true when the mass indicates the operating zone has not been reached.
	*/
	bool PoleMassSignificant(NumErrorT mass, NumErrorT approx_scale) const
	{
		return mass > NumErrorT(1e3) * this->FinalTolerance() * (NumErrorT(1) + approx_scale);
	}

	/**
	\brief Primary function that runs the Cauchy endgame.
	To begin, this function will compute a first approximation using the power series endgame notion. This approximation is made after a heuristic on the stabilization of
	the cyle number is made, and after the maximum and minimum norms of tracked space values around the origin are withing a certain tolerance.

		## Input:
			start_time: the time value at which we start the endgame
			start_point: an approximate solution to our homotopy H at start_time
			target_time: the time value that we are using the endgame to interpolate to.

		## Output:
			SuccessCode: reporting if we were successful in the endgame or if we encountered an error

		##Details:
	\tparam ComplexT The complex number type.
				This function runs the entire Cauchy Endgame. We first take our endgame boundary time value and sample to find a first approximation of the origin. This is done by
					using the idea for the power series endgame. We check for stabilization of the cycle number, and check to see when the ratios of the maximum and minimum norm of samples collected
					by CircleTrack are withing a tolerance. When both of these conditions are met we do a Hermite interpolation.
					At this point we can start tracking in to the origin while using CircleTrack to compute samples and calculating their mean to get an approximation of the origin using the Cauchy
					Integral Formula.
	*/
	template<typename ComplexT>
	SuccessCode RunImpl(ComplexT const& start_time, Vec<ComplexT> const& start_point, ComplexT const& target_time)
	{
		if (start_point.size()!=static_cast<Eigen::Index>(this->GetSystem().NumVariables()))
		{
			std::stringstream err_msg;
			err_msg << "number of variables in start point for CauchyEG, " << start_point.size() << ", must match the number of variables in the system, " << this->GetSystem().NumVariables();
			throw std::runtime_error(err_msg.str());
		}

		if (Precision(start_time)!=Precision(start_point))
		{
			std::stringstream ss;
			ss << "CauchyEG Run time and point must be of matching precision. (" << Precision(start_time) << "!=" << Precision(start_point) << ")";
			throw std::runtime_error(ss.str());
		}

		using RealT = typename Eigen::NumTraits<ComplexT>::Real;

		Vec<ComplexT>& latest_approx = this->final_approximation_;
		Vec<ComplexT>& prev_approx = this->previous_approximation_;
		NumErrorT& approx_error = this->approximate_error_;


		ClearTimesAndSamples<ComplexT>(); //clear times and samples before we begin.
		this->CycleNumber(0);
		prev_approx = start_point;

		auto init_success = GetIntoEGZone(start_time, start_point, target_time);
		if (init_success!= SuccessCode::Success)
			return init_success;

		auto cauchy_loop_success = InitialCauchyLoops<ComplexT>(target_time);
		if (cauchy_loop_success != SuccessCode::Success)
			return cauchy_loop_success;


		// The security check watches the extrapolated ENDPOINT for divergence to
		// infinity.  A genuinely-infinite endpoint -- a patched projective path leaving
		// the affine chart -- has an approximation whose dehomogenized norm grows to
		// infinity, and it grows FASTER than any finite loop sample along the way (the
		// endpoint IS the infinity; the samples are all finite), so the endpoint is the
		// earliest, strongest divergence signal.  Watching the loop samples instead
		// over-truncates finite paths whose samples transiently spike above max_norm
		// before settling (found via cyclic-6: 156 -> 153 lost, all nonsingular).  The
		// pole-annihilation case -- a Laurent pole whose Cauchy mean is a STATIONARY
		// FINITE number, which no norm check can catch -- is handled separately by the
		// pole-component operating-zone check below, so this check need only see honest
		// infinity.  Initialized to 0 so the check never reads indeterminate values.
		RealT norm_of_dehom_prev(0), norm_of_dehom_latest(0);

		if(this->SecuritySettings().level <= 0)
			norm_of_dehom_prev = this->GetSystem().InfinityNormOfDehomogenized(prev_approx);

		// Cycle-number consistency: refuse to accept a converged approximation until the cycle number
		// has reported the SAME value for num_consecutive_same_cycle_number consecutive approximations.
		// When the working precision is too low to close the loop accurately (NOT monodromy -- a
		// nonsingular endpoint is cycle 1 at every radius in high precision), the cycle number thrashes
		// (41, 14, 36, ...) and a coincidental approx_error dip must not be mistaken for convergence.
		// prev_cycle == 0 means "no prior measurement".
		unsigned prev_cycle = 0, same_cycle_count = 0;

		do
		{
			//Compute a cauchy approximation.  Uses the previously computed samples,
			//either from InitialCauchyLoops, or ComputeCauchySamples
			auto extrapolation_success = ComputeCauchyApproximationOfXAtT0<ComplexT>(latest_approx);
			if (extrapolation_success!=SuccessCode::Success)
				return extrapolation_success;

			unsigned cur_cycle = this->CycleNumber();
			if (cur_cycle == prev_cycle)
				++same_cycle_count;
			else
				same_cycle_count = 1;
			prev_cycle = cur_cycle;

			approx_error = static_cast<NumErrorT>((latest_approx - prev_approx).template lpNorm<Eigen::Infinity>());
			NotifyObservers(ApproximatedRoot<EmitterType>(*this));

			// the operating-zone ACCEPTANCE GATE: negative-mode (pole) mass in the loop.
			// Nonzero mass means the disk between here and the target time is not clean
			// -- a pole at the target, or another branch point inside the loop -- and
			// the mean is not to be trusted, however stationary it looks.  So convergence
			// is refused until the zone is reached (pole mass insignificant).  This alone
			// cures the junk-success bug: an affine Laurent pole has a stationary FINITE
			// Cauchy mean (roots-of-unity annihilation), so no norm check can catch it,
			// but its pole mass never vanishes, so this gate never accepts it -- the path
			// instead runs to the endgame's natural terminal condition (min track time /
			// fail-safe max cycle) and returns non-Success.  There is deliberately NO
			// active pole-growth truncation: distinguishing a genuine pole (mass grows
			// unboundedly) from a branch point merely inside the loop (mass grows a round
			// or two, then dies as the radius shrinks past it) by counting growing rounds
			// produced false positives that truncated CLEAN convergent paths -- cyclic-6
			// lost genuine nonsingular solutions (156 -> 153).  The gate suffices.
			auto const pole_mass = PoleComponentMass<ComplexT>();
			bool const in_operating_zone = !PoleMassSignificant(pole_mass,
				static_cast<NumErrorT>(latest_approx.template lpNorm<Eigen::Infinity>()));

			if (in_operating_zone
			    && approx_error < this->FinalTolerance()
			    && same_cycle_count >= GetCauchySettings().num_consecutive_same_cycle_number)
			{
				NotifyObservers(Converged<EmitterType>(*this));
				return SuccessCode::Success;
			}

			if (this->SecuritySettings().level <= 0)
			{//the endpoint out of bounds twice running: the path is diverging; truncate.
				norm_of_dehom_latest = this->GetSystem().InfinityNormOfDehomogenized(latest_approx);

				if (norm_of_dehom_prev   > this->SecuritySettings().max_norm &&
					norm_of_dehom_latest > this->SecuritySettings().max_norm  )
				{
					NotifyObservers(SecurityMaxNormReached<EmitterType>(*this));
					return SuccessCode::SecurityMaxNormReached;
				}
			}

			prev_approx = latest_approx;
			norm_of_dehom_prev = norm_of_dehom_latest;

			auto advance_success = AdvanceTime<ComplexT>(target_time);
			if (advance_success != SuccessCode::Success)
				return advance_success;

			// then compute the next set of cauchy samples used for extrapolating the point at target time
			auto cauchy_samples_success = ComputeCauchySamples(target_time);
			if (cauchy_samples_success != SuccessCode::Success)
				return cauchy_samples_success;

		} while (true);

		return SuccessCode::Success;
	} //end main CauchyEG function


	// ================================================================================================
	//   Adaptive-numeric-type (double-first) Cauchy endgame.
	//
	//   Computes in the hardware-complex_dbl fast lane while the AMP tracker's authoritative precision
	//   stays double, and crosses to complex_mp only when the tracker escalates a TrackPath/Refine
	//   (pure-(i)).  Fixed precision never enters here -- base Run() sends it to RunImpl<BCT>.  Every
	//   method below is a member template so the explicit fixed-precision class instantiations do not
	//   force-compile them (which would std::get a complex_mp/complex_dbl slot a single-type endgame
	//   does not have).  The complex_dbl->complex_mp container migration is the one piece with no
	//   pre-existing analog; mpfr->higher-mpfr co-vary already happens via EnsureAtUniformPrecision.
	// ================================================================================================

	/// \brief Run the Cauchy endgame via the double-first adaptive-numeric-type driver, migrating up to mpfr only on tracker escalation.
	template<typename Dummy = void>
	SuccessCode RunImplAMP(BCT const& start_time, Vec<BCT> const& start_point, BCT const& target_time)
	{
		using bertini::Precision;
		using RealT = typename Eigen::NumTraits<BCT>::Real;

		if (start_point.size()!=static_cast<Eigen::Index>(this->GetSystem().NumVariables()))
		{
			std::stringstream err_msg;
			err_msg << "number of variables in start point for CauchyEG, " << start_point.size() << ", must match the number of variables in the system, " << this->GetSystem().NumVariables();
			throw std::runtime_error(err_msg.str());
		}

		// Arm the escalation hooks in the shared phase methods (CircleTrack / AdvanceTime / refine);
		// disarm on every exit path.
		this->adaptive_numeric_type_active_ = true;
		struct Disarmer { bool& flag; ~Disarmer(){ flag = false; } } disarm{this->adaptive_numeric_type_active_};

		// Start in the precision the tracker handed us at the endgame boundary: double for the easy
		// majority, already-mpfr for the few paths that escalated before the endgame.  Migrate UP only.
		this->current_endgame_precision_ = std::max(DoublePrecision(), Precision(start_point));
		this->CycleNumber(0);

		// ---- INITIALIZATION (GetIntoEGZone + InitialCauchyLoops).  On escalation, restart init from the
		//      boundary at the higher precision: init lives at the well-conditioned large-|t| end where
		//      escalation is rare, and a fresh init clears its own containers so it is self-consistent. ----
		while (true)
		{
			SuccessCode init_code;
			if (this->current_endgame_precision_ == DoublePrecision())
				init_code = RunInitSegmentT<complex_dbl>(complex_dbl(start_time), this->DowncastToDouble(start_point), complex_dbl(target_time));
			else
				init_code = RunInitSegmentT<complex_mp>(this->AtActivePrecisionScalar(start_time), this->AtActivePrecisionVec(start_point), this->AtActivePrecisionScalar(target_time));

			if (init_code == SuccessCode::HigherPrecisionNecessary)
			{
				this->current_endgame_precision_ = this->NextEscalatedPrecision();
				SetThreadPrecision(this->current_endgame_precision_);
				this->GetSystem().precision(this->current_endgame_precision_);
				continue;
			}
			if (init_code != SuccessCode::Success)
				return init_code;
			break;
		}

		// ---- MAIN CONVERGENCE LOOP.  Each phase self-heals: on escalation it migrates the durable state
		//      up to mpfr (widening the retained PSEG window, never re-tracking it) and retries in mpfr.
		//      The approximations live in BCT, so the per-iteration bookkeeping arithmetic is mpfr -- one
		//      small vector op next to the tracking, which itself stays in the fast lane. ----
		// the security check watches the extrapolated ENDPOINT for honest divergence to
		// infinity (which grows faster than any finite loop sample); acceptance is gated
		// on the pole-component operating-zone measurement, which handles the
		// finite-mean pole case -- see RunImpl for the full reasoning
		RealT norm_of_dehom_prev(0), norm_of_dehom_latest(0);
		if (this->SecuritySettings().level <= 0)
			norm_of_dehom_prev = this->GetSystem().InfinityNormOfDehomogenized(this->previous_approximation_);

		unsigned prev_cycle = 0, same_cycle_count = 0;

		while (true)
		{
			auto extrap_code = ComputeCauchyApproxAMP();
			if (extrap_code != SuccessCode::Success)
				return extrap_code;

			unsigned cur_cycle = this->CycleNumber();
			if (cur_cycle == prev_cycle) ++same_cycle_count; else same_cycle_count = 1;
			prev_cycle = cur_cycle;

			Precision(this->previous_approximation_, Precision(this->final_approximation_));
			this->approximate_error_ = static_cast<NumErrorT>((this->final_approximation_ - this->previous_approximation_).template lpNorm<Eigen::Infinity>());
			NotifyObservers(ApproximatedRoot<EmitterType>(*this));

			// the operating-zone ACCEPTANCE GATE (see RunImpl for the full reasoning):
			// convergence is refused while pole mass is significant, which alone cures
			// junk-success; there is deliberately no active pole-growth truncation (it
			// false-positived on clean convergent paths -- cyclic-6, 156 -> 153).
			auto const pole_mass = PoleComponentMassAMP();
			bool const in_operating_zone = !PoleMassSignificant(pole_mass,
				static_cast<NumErrorT>(this->final_approximation_.template lpNorm<Eigen::Infinity>()));

			if (in_operating_zone
			    && this->approximate_error_ < this->FinalTolerance()
			    && same_cycle_count >= GetCauchySettings().num_consecutive_same_cycle_number)
			{
				NotifyObservers(Converged<EmitterType>(*this));
				return SuccessCode::Success;
			}

			if (this->SecuritySettings().level <= 0)
			{
				norm_of_dehom_latest = this->GetSystem().InfinityNormOfDehomogenized(this->final_approximation_);
				if (norm_of_dehom_prev   > this->SecuritySettings().max_norm &&
				    norm_of_dehom_latest > this->SecuritySettings().max_norm)
				{
					NotifyObservers(SecurityMaxNormReached<EmitterType>(*this));
					return SuccessCode::SecurityMaxNormReached;
				}
			}

			this->previous_approximation_ = this->final_approximation_;
			norm_of_dehom_prev = norm_of_dehom_latest;

			auto advance_code = AdvanceTimeAMP(target_time);
			if (advance_code != SuccessCode::Success)
				return advance_code;

			auto samples_code = ComputeCauchySamplesAMP(target_time);
			if (samples_code != SuccessCode::Success)
				return samples_code;
		}

		return SuccessCode::Success;
	}


	// Initialization at one numeric type: seed the PSEG window and reach the EG operating zone, exactly
	// as the head of RunImpl does, but reporting HigherPrecisionNecessary up to the driver on escalation.
	/// \brief Run the initialization segment (seed the PSEG window and reach the EG zone) in a given numeric type.
	template<typename ComplexT>
	SuccessCode RunInitSegmentT(ComplexT const& start_time, Vec<ComplexT> const& start_point, ComplexT const& target_time)
	{
		ClearTimesAndSamples<ComplexT>();
		this->CycleNumber(0);
		this->previous_approximation_ = this->ToBCT(start_point, this->current_endgame_precision_);

		auto init_success = GetIntoEGZone(start_time, start_point, target_time);
		if (init_success != SuccessCode::Success)
			return init_success;

		return InitialCauchyLoops<ComplexT>(target_time);
	}


	// Cauchy mean (extrapolation) at the active numeric type, written into final_approximation_ (BCT).
	/// \brief Compute the Cauchy mean (extrapolation) at the active numeric type, migrating-and-retrying in mpfr on escalation.
	template<typename Dummy = void>
	SuccessCode ComputeCauchyApproxAMP()
	{
		unsigned guard = 0;
		while (true)
		{
			SuccessCode code;
			if (this->current_endgame_precision_ == DoublePrecision())
			{
				Vec<complex_dbl> r;
				code = ComputeCauchyApproximationOfXAtT0<complex_dbl>(r);
				if (code == SuccessCode::Success)
					this->final_approximation_ = this->ToBCT(r, this->current_endgame_precision_);
			}
			else
			{
				code = ComputeCauchyApproximationOfXAtT0<complex_mp>(this->final_approximation_);
			}

			if (code == SuccessCode::HigherPrecisionNecessary || code == SuccessCode::FailedToConverge)
			{
				if (this->template EscalateAndMigrate<>(++guard) != SuccessCode::Success)
					return SuccessCode::HigherPrecisionNecessary;
				continue;
			}
			return code;
		}
	}


	// Advance time at the active numeric type; migrate-and-retry in mpfr on escalation.  AdvanceTime
	// returns HigherPrecisionNecessary BEFORE it rotates the PSEG window, so the window stays a clean
	// checkpoint and the retry continues from the widened window with no double-advance.
	/// \brief Advance time at the active numeric type, migrating-and-retrying in mpfr on escalation.
	template<typename Dummy = void>
	SuccessCode AdvanceTimeAMP(BCT const& target_time)
	{
		unsigned guard = 0;
		while (true)
		{
			SuccessCode code = (this->current_endgame_precision_ == DoublePrecision())
				? AdvanceTime<complex_dbl>(complex_dbl(target_time))
				: AdvanceTime<complex_mp>(this->AtActivePrecisionScalar(target_time));

			if (code == SuccessCode::HigherPrecisionNecessary)
			{
				if (this->template EscalateAndMigrate<>(++guard) != SuccessCode::Success)
					return SuccessCode::HigherPrecisionNecessary;
				continue;
			}
			return code;
		}
	}


	// Build a closed Cauchy loop's samples at the active numeric type; migrate-and-retry in mpfr on
	// escalation.  ComputeCauchySamples clears and re-seeds its cauchy data from the PSEG window, so the
	// mpfr retry simply re-tracks the circle from the (migrated) window -- it never reuses the lossy,
	// double-tracked partial circle that triggered the escalation.
	/// \brief Build a closed Cauchy loop's samples at the active numeric type, migrating-and-retrying in mpfr on escalation.
	template<typename Dummy = void>
	SuccessCode ComputeCauchySamplesAMP(BCT const& target_time)
	{
		unsigned guard = 0;
		while (true)
		{
			SuccessCode code = (this->current_endgame_precision_ == DoublePrecision())
				? ComputeCauchySamples<complex_dbl>(complex_dbl(target_time))
				: ComputeCauchySamples<complex_mp>(this->AtActivePrecisionScalar(target_time));

			if (code == SuccessCode::HigherPrecisionNecessary)
			{
				if (this->template EscalateAndMigrate<>(++guard) != SuccessCode::Success)
					return SuccessCode::HigherPrecisionNecessary;
				continue;
			}
			return code;
		}
	}


	// Flavor-specific: cross every Cauchy container from the complex_dbl slot to complex_mp (or, if
	// already mpfr, raise its precision uniformly), via the shared base Cross* / SetPrecision helpers.
	// Widen-only by default -- retained samples were already tracked/refined to final_tolerance
	// (pure-(i)/B).  Called by the base EscalateAndMigrate.
	/// \brief Widen this flavor's Cauchy containers from the complex_dbl slot up to complex_mp at the new precision.
	template<typename Dummy = void>
	void MigrateContainersToPrecision(unsigned newprec)
	{
		using bertini::Precision;
		if (this->current_endgame_precision_ == DoublePrecision())
		{
			this->template CrossTimesUp<>(pseg_times_,     newprec);
			this->template CrossSampsUp<>(pseg_samples_,   newprec);
			this->template CrossTimesUp<>(cauchy_times_,   newprec);
			this->template CrossSampsUp<>(cauchy_samples_, newprec);
			this->template CrossVecUp<>  (c_over_k_probe_, newprec);
		}
		else
		{
			tracking::adaptive::SetPrecision(std::get<TimeCont<complex_mp>>(pseg_times_),     newprec);
			tracking::adaptive::SetPrecision(std::get<SampCont<complex_mp>>(pseg_samples_),   newprec);
			tracking::adaptive::SetPrecision(std::get<TimeCont<complex_mp>>(cauchy_times_),   newprec);
			tracking::adaptive::SetPrecision(std::get<SampCont<complex_mp>>(cauchy_samples_), newprec);
			auto& pm = std::get<Vec<complex_mp>>(c_over_k_probe_);
			if (pm.size() > 0) Precision(pm, newprec);
		}

		if (this->final_approximation_.size()    > 0) Precision(this->final_approximation_,    newprec);
		if (this->previous_approximation_.size() > 0) Precision(this->previous_approximation_, newprec);
		this->GetSystem().precision(newprec);

		if (this->EndgameSettings().refine_when_increasing_precision)
		{
			auto& cau_t = std::get<TimeCont<complex_mp>>(cauchy_times_);
			auto& cau_s = std::get<SampCont<complex_mp>>(cauchy_samples_);
			if (!cau_s.empty())
				this->template RefineAllSamples<complex_mp>(cau_s, cau_t);
		}
	}

};


}} // namespaces
