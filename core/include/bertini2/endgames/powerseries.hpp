//This file is part of Bertini 2.
//
//powerseries_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//powerseries_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with powerseries_endgame.hpp.  If not, see <http://www.gnu.org/licenses/>.
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


namespace bertini{ namespace endgame{


/** 
\class PowerSeriesEndgame

\brief class used to finish tracking paths during Homotopy Continuation

\ingroup endgame

## Explanation

The bertini::PowerSeriesEndgame class enables us to finish tracking on possibly singular paths on an arbitrary square homotopy.  

The intended usage is to:

1. Create a system, tracker, and instantiate some settings.
2. Using the tracker created track to the engame boundary, by default this is t = 0.1. 
3. Create a PowerSeriesEndgame, associating it to the tracker you wish to use. The tracker knows the system being solved.
4. For each path being tracked send the PowerSeriesEndgame the time value and other variable values that it should use to start the endgame. 
5. The PowerSeriesEndgame, if successful, will store the homotopy solutions at t = 0.

## Example Usage

Below we demonstrate a basic usage of the PowerSeriesEndgame class to find the singularity at t = 0. 

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

PowerSeriesConfig power_series_settings;
power_series_settings.max_cycle_number = 4;


// 5. Create a power series endgame, and use them to get the soutions at t = 0. 
EndgameSelector<TrackerType>::PSEG my_pseg_endgame(tracker,power_series_settings,tolerances);


std::vector<Vec<ComplexT> > my_homotopy_solutions; 

std::vector<Vec<ComplexT> > my_homotopy_divergent_paths; 

for(auto s : my_homotopy_solutions_at_endgame_boundary) 
{
    SuccessCode endgame_success = my_pseg_endgame.Run(t_endgame_boundary,s);

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

File: test/endgames/generic_pseg_test.hpp
File: test/endgames/amp_powerseries_test.cpp
File: test/endgames/fixed_double_powerseries_test.cpp
FIle: test/endgames/fixed_multiple_powerseries_test.cpp
*/

template<typename PrecT> 
class PowerSeriesEndgame : 
	public virtual EndgameBase<PowerSeriesEndgame<PrecT>, PrecT>
{
public:
	using BaseEGT = EndgameBase<PowerSeriesEndgame<PrecT>, PrecT>;  ///< The base endgame type.
	using FinalEGT = PowerSeriesEndgame<PrecT>;  ///< The final (derived) endgame type.
	using TrackerType = typename BaseEGT::TrackerType;  ///< The path-tracker type.

	using BaseComplexT = typename BaseEGT::BaseComplexT;  ///< The complex number type.
	using BaseRealT = typename BaseEGT::BaseRealT;  ///< The real number type.

	using EmitterType = PowerSeriesEndgame<PrecT>;  ///< The event-emitter type.

protected:
	
	using EndgameBase<PowerSeriesEndgame<PrecT>, PrecT>::NotifyObservers;

	using TupleOfTimes = typename BaseEGT::TupleOfTimes;  ///< A tuple of time containers, one per numeric type.
	using TupleOfSamps = typename BaseEGT::TupleOfSamps;  ///< A tuple of sample containers, one per numeric type.
	using TupOfVec = typename BaseEGT::TupOfVec;  ///< A tuple of vector containers, one per numeric type.

	using BCT = BaseComplexT;  ///< The complex number type.
	using BRT = BaseRealT;  ///< The real number type.

	using Configs = typename AlgoTraits<FinalEGT>::NeededConfigs;  ///< The configuration bundle (Configured base).
	using ConfigsAsTuple = typename Configs::ToTuple;  ///< The configuration structs as a tuple.

	/**
	\brief State variable representing a computed upper bound on the cycle number.
	*/
	mutable unsigned upper_bound_on_cycle_number_;

	/**
	\brief Holds the time values for different space values used in the Power series endgame. 
	*/	
	mutable TupleOfTimes times_;

	/**
	\brief Holds the space values used in the Power series endgame. 
	*/			
	mutable TupleOfSamps samples_;

	/**
	\brief Holds the derivatives at each space point. 
	*/			
	mutable TupleOfSamps derivatives_;

	/**
	\brief Random vector used in computing an upper bound on the cycle number.

	Dual-slot (TupOfVec) so the adaptive-numeric-type endgame can form the cycle-number dot products in
	the active type (it is multiplied against Vec<ComplexT> sample differences, which must match scalar
	type).  Single-slot for fixed precision, as before.
	*/
	mutable TupOfVec rand_vector_;

	/// \brief Debug-assert that the stored time and sample containers have consistent, sufficient sizes.
	template<typename ComplexT>
	void AssertSizesTimeSpace() const
	{
#ifndef NDEBUG
		const auto num_sample_points = this->EndgameSettings().num_sample_points;
		assert(std::get<SampCont<ComplexT> >(samples_).size()==std::get<TimeCont<ComplexT> >(times_).size() && "must have same number of samples in times and spaces");
		assert(std::get<SampCont<ComplexT> >(samples_).size()>=num_sample_points && "must have sufficient number of samples");
#endif
	}

	/// \brief Debug-assert that the time, sample, and derivative containers have consistent, sufficient sizes.
	template<typename ComplexT>
	void AssertSizesTimeSpaceDeriv() const
	{
#ifndef NDEBUG
		const auto num_sample_points = this->EndgameSettings().num_sample_points;
		assert(std::get<SampCont<ComplexT> >(samples_).size()==std::get<TimeCont<ComplexT> >(times_).size() && "must have same number of samples in times and spaces");
		assert(std::get<SampCont<ComplexT> >(samples_).size()==std::get<SampCont<ComplexT> >(derivatives_).size() && "must have same number of samples in derivatives and spaces");
		assert(std::get<SampCont<ComplexT> >(samples_).size()>=num_sample_points && "must have sufficient number of samples");
#endif
	}

public:

	/// \return The computed upper bound on the cycle number.
	auto UpperBoundOnCycleNumber() const { return upper_bound_on_cycle_number_;}


	/**
	\brief Function that clears all samples and times from data members for the Power Series endgame
	*/	
	template<typename ComplexT>
	void ClearTimesAndSamples()
	{
		std::get<TimeCont<ComplexT> >(times_).clear(); 
		std::get<SampCont<ComplexT> >(samples_).clear();
	}

	/**
	\brief Function to set the times used for the Power Series endgame.
	*/	
	template<typename ComplexT>
	void SetTimes(TimeCont<ComplexT> const& times_to_set) { std::get<TimeCont<ComplexT> >(times_) = times_to_set;}

	/**
	\brief Function to get the times used for the Power Series endgame.
	*/	
	template<typename ComplexT>
	const auto& GetTimes() const {return std::get<TimeCont<ComplexT> >(times_);}


	mutable BCT latest_time_cache_;  ///< Scratch so LatestTimeImpl can return a BCT reference in the complex_dbl fast lane.

	/// \return The most recent time value in the sample sequence.
	const BCT& LatestTimeImpl() const
	{
		// In the adaptive-numeric-type endgame the latest time may live in the complex_dbl slot, with
		// the BCT slot empty.  Dispatch on which slot holds data so this is correct in either lane and
		// both during and after the run.  Fixed precision compiles to the original BCT read.
		if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			if (GetTimes<BCT>().empty())
			{
				auto const& dbl_times = GetTimes<complex_dbl>();
				latest_time_cache_ = dbl_times.empty() ? BCT(0) : BCT(dbl_times.back());
				return latest_time_cache_;
			}
		}
		return GetTimes<BCT>().back();
	}

	/**
	\brief Function to set the space values used for the Power Series endgame.
	*/	
	template<typename ComplexT>
	void SetSamples(SampCont<ComplexT> const& samples_to_set) { std::get<SampCont<ComplexT> >(samples_) = samples_to_set;}

	/**
	\brief Function to get the space values used for the Power Series endgame.
	*/	
	template<typename ComplexT>
	const auto& GetSamples() const {return std::get<SampCont<ComplexT> >(samples_);}

	/**
	\brief Function to set the times used for the Power Series endgame.
	*/	
	template<typename ComplexT>
	void SetRandVec(int size)
	{
		auto& rv = std::get<Vec<ComplexT> >(rand_vector_);
		rv.resize(size);
		for (int ii = 0; ii < size; ++ii)
			rv(ii) = RandomUnit<ComplexT>();
	}





	/// \brief Construct the power-series endgame for a tracker, with its configuration as a tuple.
	explicit PowerSeriesEndgame(TrackerType const& tr,
	                            const ConfigsAsTuple& settings )
      : EndgamePrecPolicyBase<TrackerType>(tr), BaseEGT(tr, settings)
   	{}

	/// \brief Construct the power-series endgame for a tracker, with configs given in any order.
    template< typename... Ts >
	explicit
	PowerSeriesEndgame(TrackerType const& tr, const Ts&... ts ) : PowerSeriesEndgame(tr, Configs::Unpermute( ts... ) )
		{}



	/**
	\brief Computes an upper bound on the cycle number. Consult page 53 of \cite bertinibook.

	## Input: 
			None: all data needed are class data members.

	## Output:
			upper_bound_on_cycle_number_: Used for an exhaustive search for the best cycle number for approimating the path to t = 0.

	##Details:
			\tparam ComplexT The complex number type.
	*/
	template<typename ComplexT>
	unsigned ComputeBoundOnCycleNumber()
	{ 
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;
		using std::log; using std::abs;

		const auto& samples = std::get<SampCont<ComplexT> >(samples_);
		
		AssertSizesTimeSpace<ComplexT>();

		auto num_samples = samples.size();
		const Vec<ComplexT> & sample0 = samples[num_samples-3];
		const Vec<ComplexT> & sample1 = samples[num_samples-2];
		const Vec<ComplexT> & sample2 = samples[num_samples-1]; // most recent sample.  oldest samples at front of the container


// should this only be if the system is homogenized?
		const auto& rand_vector = std::get<Vec<ComplexT> >(rand_vector_);
		ComplexT rand_sum1 = ((sample1 - sample0).transpose()*rand_vector).sum();
		ComplexT rand_sum2 = ((sample2 - sample1).transpose()*rand_vector).sum();

		if ( abs(rand_sum1)==0 || abs(rand_sum2)==0) // avoid division by 0
		{
			upper_bound_on_cycle_number_ = 1;
			return upper_bound_on_cycle_number_;
		}

		RealT estimate = log(static_cast<RealT>(this->EndgameSettings().sample_factor))/log(abs(rand_sum2/rand_sum1));

		const auto& ps_config = this->template Get<PowerSeriesConfig>();
		if (estimate < 1)
		  	upper_bound_on_cycle_number_ = 1;
		else
		{
			// max_cycle_number is a CEILING on the candidate search (each candidate costs a full
			// Hermite solve).  Clamp before any conversion to unsigned: near-unity sample ratios
			// (slow convergence -- high multiplicity, or a slow diverger) drive the estimate toward
			// +inf, and unsigned(inf) is undefined behavior.  The !(a < b) form also routes NaN to
			// the ceiling.  This was `max` -- the ceiling was a floor, and the search was unbounded.
			using std::round;
			RealT amplified = round(estimate) * static_cast<RealT>(ps_config.cycle_number_amplification);
			if (!(amplified < static_cast<RealT>(ps_config.max_cycle_number)))
				upper_bound_on_cycle_number_ = ps_config.max_cycle_number;
			else
				upper_bound_on_cycle_number_ = std::max(1u, static_cast<unsigned>(amplified));
		}

		return upper_bound_on_cycle_number_;
	}//end ComputeBoundOnCycleNumber




	/**
	\brief This function computes the cycle number using an exhaustive search up the upper bound computed by the above function BoundOnCyleNumber. 

		## Input: 
				None: all data needed are class data members.

		## Output:
				cycle_number_: Used to create a hermite interpolation to t = 0. 

		##Details:
				\tparam ComplexT The complex number type.
			This is done by an exhaustive search from 1 to upper_bound_on_cycle_number. There is a conversion to the s-space from t-space in this function. 
	As a by-product the derivatives at each of the samples is returned for further use. 
	*/

	template<typename ComplexT>
	unsigned ComputeCycleNumber(ComplexT const& t0)
	{
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;

		const auto& samples = std::get<SampCont<ComplexT> >(samples_);
		const auto& times   = std::get<TimeCont<ComplexT> >(times_);

		AssertSizesTimeSpaceDeriv<ComplexT>();
		
		const Vec<ComplexT> &most_recent_sample = samples.back();  
		const ComplexT& most_recent_time = times.back();

		//Compute upper bound for cycle number.
		ComputeBoundOnCycleNumber<ComplexT>();


		unsigned num_pts;
		if (samples.size() > this->EndgameSettings().num_sample_points)
			num_pts = this->EndgameSettings().num_sample_points;
		else 
			num_pts = this->EndgameSettings().num_sample_points-1;


		auto min_found_difference = Eigen::NumTraits<RealT>::highest();

		// Default the selection before the candidate loop: if every candidate's difference is
		// NaN (poisoned samples), the comparisons below are all false and nothing is assigned --
		// a fresh endgame would otherwise carry cycle number 0 into TransformToSPlane and throw.
		this->cycle_number_ = 1;

		TimeCont<ComplexT> s_times(num_pts);
		SampCont<ComplexT> s_derivatives(num_pts);

		for(unsigned int candidate = 1; candidate <= upper_bound_on_cycle_number_; ++candidate)
		{			
			using std::pow;

			std::tie(s_times, s_derivatives) = TransformToSPlane(static_cast<int>(candidate), t0, num_pts, ContStart::Front);
			RealT cand_power{1/static_cast<RealT>(candidate)};
			RealT curr_diff = (HermiteInterpolateAndSolve<ComplexT>(
								  pow((most_recent_time-t0)/(times[0]-t0),cand_power), // the target time
			                      num_pts,s_times,samples,s_derivatives, ContStart::Front) // the input data
			                 - 
			                 most_recent_sample).template lpNorm<Eigen::Infinity>();

			if (curr_diff < min_found_difference)
			{
				min_found_difference = curr_diff;
				this->cycle_number_ = candidate;
			}

		}// end cc loop over cycle number possibilities

		return this->cycle_number_;
	}//end ComputeCycleNumber



	/**
		\brief Compute a set of derivatives using internal data to the endgame.

		## Input: 
				None: all data needed are class data members.

		## Output:
				None: Derivatives are members of this class.

		##Details:
				\tparam ComplexT The complex number type.
	*/
	template<typename ComplexT>
	void ComputeAllDerivatives()
	{
		auto& samples = std::get<SampCont<ComplexT> >(samples_);
		auto& times   = std::get<TimeCont<ComplexT> >(times_);
		auto& derivatives = std::get<SampCont<ComplexT> >(derivatives_);

		assert((samples.size() == times.size()) && "must have same number of times and samples");

		if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec) // known at compile time
		{
			auto max_precision = this->EnsureAtUniformPrecision(times, samples);
			this->GetSystem().precision(max_precision);
		}

		//Compute dx_dt for each sample.
		derivatives.clear(); derivatives.resize(samples.size());
		for(unsigned ii = 0; ii < samples.size(); ++ii)
		{	
			derivatives[ii] = -this->GetSystem().Jacobian(samples[ii],times[ii]).lu().solve(this->GetSystem().TimeDerivative(samples[ii],times[ii]));
		}
	}


	/**
	\brief Transform the stored times and derivatives into the S-plane, scaled by the cycle number.

	This function transforms the times and derivatives into the S-plane, scaled by the cycle number, and maps them into the interval [0, 1].

	\param cycle_num The cycle number to use.
	\param t0 The base time value of the transform.
	\param num_pts The number of points to transform.
	\param shift_from Which end of the sample containers to transform from.
	*/
	template <typename ComplexT>
	std::tuple<TimeCont<ComplexT>, SampCont<ComplexT>> TransformToSPlane(int cycle_num, ComplexT const& t0, unsigned num_pts, ContStart shift_from)
	{
		if (cycle_num==0)
			throw std::runtime_error("cannot transform to s plane with cycle number 0");
		AssertSizesTimeSpaceDeriv<ComplexT>();


		using RealT = typename Eigen::NumTraits<ComplexT>::Real;

		const auto& times   = std::get<TimeCont<ComplexT> >(times_);
		const auto& derivatives  = std::get<SampCont<ComplexT> >(derivatives_);

		RealT c = static_cast<RealT>(cycle_num);
		RealT one_over_c = 1/c;

		size_t offset_t, offset_d;
		if (shift_from == ContStart::Back)
		{
			offset_t = times.size()-num_pts;
			offset_d = derivatives.size()-num_pts;
		}
		else
			offset_t = offset_d = 0;


		TimeCont<ComplexT> s_times(num_pts);
		SampCont<ComplexT> s_derivatives(num_pts);

		ComplexT time_shift = times[offset_t] - t0;

		for(unsigned ii = 0; ii < num_pts; ++ii){
			s_times[ii] = pow((times[ii+offset_t]-t0)/time_shift, one_over_c); 
			s_derivatives[ii] = derivatives[ii+offset_d]*( c*pow(s_times[ii],cycle_num-1))*time_shift;
		}

		return std::make_tuple(s_times, s_derivatives);
	}



	/**
	\brief This function computes an approximation of the space value at the time time_t0. 

		## Input: 
				result: Passed by reference this holds the value of the approximation we compute
				t0: This is the time value for which we wish to compute an approximation at. 

		## Output:
				SuccessCode: This reports back if we were successful in making an approximation.

		##Details:
	\tparam ComplexT The complex number type.
				This function handles computing an approximation at the origin. 
				We compute the cycle number best for the approximation, and convert derivatives and times to the s-plane where s = t^(1/c).
				We use the converted times and derivatives along with the samples to do a Hermite interpolation.
	*/
	template<typename ComplexT>
	SuccessCode ComputeApproximationOfXAtT0(Vec<ComplexT>& result, const ComplexT & t0)
	{	
		const auto c = ComputeCycleNumber<ComplexT>(t0);

		auto num_pts = this->EndgameSettings().num_sample_points;

		TimeCont<ComplexT> s_times;
		SampCont<ComplexT> s_derivatives;

		std::tie(s_times, s_derivatives) = TransformToSPlane(static_cast<int>(c), t0, num_pts, ContStart::Back);
		// the data was transformed to be on the interval [0 1] so we can hard-code the time-to-solve as 0 here.

		Precision(result, Precision(s_derivatives.back()));
		result = HermiteInterpolateAndSolve(ComplexT(0), num_pts, s_times, std::get<SampCont<ComplexT> >(samples_), s_derivatives, ContStart::Back);
		// A NaN extrapolation must be a FAILURE code.  The run loop converges on
		// `approx_error > FinalTolerance()` becoming false, and every IEEE comparison against
		// NaN is false -- so a NaN approximation would exit the loop down the SUCCESS path,
		// reporting Converged with a poisoned answer.
		if (bertini::ContainsNaN(result))
			return SuccessCode::FailedToConverge;
		return SuccessCode::Success;
	}//end ComputeApproximationOfXAtT0



	/**
		\brief The samples used in the power series endgame are collected by advancing time to t = 0, by multiplying the current time by the sample factor. 

		## Input: 
				target_time: This is the time we are trying to approximate, default is t = 0.

		## Output:
				SuccessCode: This reports back if we were successful in advancing time. 

		##Details:
				\tparam ComplexT The complex number type.
				This function computes the next time value for the power series endgame. After computing this time value, 
				it will track to it and compute the derivative at this time value for further appoximations to be made during the
				endgame.
	*/
	template<typename ComplexT>
	SuccessCode AdvanceTime(const ComplexT & target_time)
	{
		using RealT = typename Eigen::NumTraits<ComplexT>::Real;
		
		auto& samples = std::get<SampCont<ComplexT> >(samples_);
		auto& times   = std::get<TimeCont<ComplexT> >(times_);

		AssertSizesTimeSpaceDeriv<ComplexT>();

		Vec<ComplexT> next_sample;
		ComplexT next_time = (times.back() + target_time) * static_cast<RealT>(this->EndgameSettings().sample_factor); //setting up next time value using the midpoint formula, sample_factor will give us some 

  		if (abs(next_time - target_time) < this->EndgameSettings().min_track_time) // generalized for target_time not equal to 0.
  		{
  			NotifyObservers(MinTrackTimeReached<EmitterType>(*this));
  			return SuccessCode::MinTrackTimeReached;
  		}

		SuccessCode tracking_success = this->EndgameTrackPath(next_sample,times.back(),next_time,samples.back());
			if (tracking_success != SuccessCode::Success)
				return tracking_success;

		// Pure-(i) escalation: return BEFORE pushing the new sample so the window stays an untouched
		// checkpoint; the driver migrates to mpfr and retries this advance.  Elided for fixed precision.
		if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
			if (this->adaptive_numeric_type_active_ &&
			    this->GetTracker().GetCurrentPrecision() > this->current_endgame_precision_)
				return SuccessCode::HigherPrecisionNecessary;

		NotifyObservers(InEGOperatingZone<EmitterType>(*this));

		this->EnsureAtPrecision(next_time,Precision(next_sample));


		times.push_back(next_time);
		samples.push_back(next_sample);




		auto refine_success = this->RefineSample(samples.back(), next_sample,  times.back(),
										this->FinalTolerance() * this->EndgameSettings().sample_point_refinement_factor,
										this->EndgameSettings().max_num_refinements);
		if (refine_success != SuccessCode::Success)
		{
			// Adaptive-numeric-type: a refine that double cannot satisfy is a request to cross to mpfr.
			// Roll back the just-pushed sample so the window is a clean checkpoint, then signal the driver.
			if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				if (this->adaptive_numeric_type_active_ &&
				    (refine_success == SuccessCode::HigherPrecisionNecessary || refine_success == SuccessCode::FailedToConverge))
				{
					times.pop_back();
					samples.pop_back();
					return SuccessCode::HigherPrecisionNecessary;
				}
			NotifyObservers(RefiningFailed<EmitterType>(*this));
			return refine_success;
		}
		
		this->EnsureAtPrecision(times.back(),Precision(samples.back()));

		NotifyObservers(SampleRefined<EmitterType>(*this));

		// we keep one more samplepoint than needed around, for estimating the cycle number
		if (times.size() > this->EndgameSettings().num_sample_points+1)
		{
			times.pop_front();
			samples.pop_front();
		}

 		return SuccessCode::Success;
	}


	/**
	\brief Primary function running the Power Series endgame. 

		## Input: 
				start_time: This is the time value for which the endgame begins, by default this is t = 0.1
				start_point: An approximate solution of the homotopy at t = start_time
				target_time: The time value that we are wishing to approximate to. This is default set to t = 0. 

		## Output:
				SuccessCode: This reports back if we were successful in advancing time. 

		##Details:
	\tparam ComplexT The complex number type.
				Tracking forward with the number of sample points, this function will make approximations using Hermite interpolation. This process will continue until two consecutive
				approximations are withing final tolerance of each other. 
	*/		
	template<typename ComplexT>
	SuccessCode RunImpl(const ComplexT & start_time, const Vec<ComplexT> & start_point, ComplexT const& target_time)
	{
		if (start_point.size()!=static_cast<Eigen::Index>(this->GetSystem().NumVariables()))
		{
			std::stringstream err_msg;
			err_msg << "number of variables in start point for PSEG, " << start_point.size() << ", must match the number of variables in the system, " << this->GetSystem().NumVariables();
			throw std::runtime_error(err_msg.str());
		}

		// thread-local only: the endgame may run on a std::thread worker.
		SetThreadPrecision(Precision(start_point));

		using RealT = typename Eigen::NumTraits<ComplexT>::Real;
		//Set up for the endgame.
		ClearTimesAndSamples<ComplexT>();


		// unpack some references for easy use
		auto& samples = std::get<SampCont<ComplexT> >(samples_);
		auto& times   = std::get<TimeCont<ComplexT> >(times_);
		Vec<ComplexT>& latest_approx = this->final_approximation_;
		Vec<ComplexT>& prev_approx = this->previous_approximation_;

		// this is for estimating a ... norm?
		SetRandVec<ComplexT>(static_cast<int>(start_point.size()));
		
		
		auto initial_sample_success = this->ComputeInitialSamples(start_time, target_time, start_point, times, samples);

		if (initial_sample_success!=SuccessCode::Success)
		{
			NotifyObservers(EndgameFailure<EmitterType>(*this));
			return initial_sample_success;
		}

		this->template RefineAllSamples<ComplexT>(samples, times);
		ComputeAllDerivatives<ComplexT>();



	 	auto extrapolation_code = ComputeApproximationOfXAtT0(prev_approx, target_time);
	 	latest_approx = prev_approx;

	 	if (extrapolation_code != SuccessCode::Success)
	 		return extrapolation_code;

	 	RealT norm_of_dehom_of_latest_approx(0); // initialized to 0 so the security check never reads indeterminate values
	 	RealT norm_of_dehom_of_prev_approx(0);
	 	if (this->SecuritySettings().level <= 0)
	 	 	norm_of_dehom_of_prev_approx = this->GetSystem().InfinityNormOfDehomogenized(prev_approx);


	 	NumErrorT& approx_error = this->approximate_error_;
		approx_error = 1;

		while (approx_error > this->FinalTolerance())
		{
	  		auto advance_code = AdvanceTime<ComplexT>(target_time);
	  		if (advance_code!=SuccessCode::Success)
	 		{
	 			NotifyObservers(EndgameFailure<EmitterType>(*this));
	 			return advance_code;
	 		}

	 		// this code is what bertini1 does... it refines all samples, like, all the time.
	 		this->template RefineAllSamples<ComplexT>(samples, times);
	 		ComputeAllDerivatives<ComplexT>();

	 		extrapolation_code = ComputeApproximationOfXAtT0(latest_approx, target_time);
	 		if (extrapolation_code!=SuccessCode::Success)
	 		{
	 			NotifyObservers(EndgameFailure<EmitterType>(*this));
	 			return extrapolation_code;
	 		}

	 		approx_error = static_cast<NumErrorT>((latest_approx - prev_approx).template lpNorm<Eigen::Infinity>());
	 		NotifyObservers(ApproximatedRoot<EmitterType>(*this));


	 		if(this->SecuritySettings().level <= 0)
	 		{
	 			norm_of_dehom_of_latest_approx = this->GetSystem().InfinityNormOfDehomogenized(latest_approx);
		 		if(this->BeyondSecurityMaxNorm(norm_of_dehom_of_latest_approx) && this->BeyondSecurityMaxNorm(norm_of_dehom_of_prev_approx))
		 		{
		 			NotifyObservers(SecurityMaxNormReached<EmitterType>(*this));
	 				return SuccessCode::SecurityMaxNormReached;
		 		}
	 			norm_of_dehom_of_prev_approx = norm_of_dehom_of_latest_approx;
	 		}

	 		


	 		Precision(prev_approx, Precision(latest_approx));
	 		prev_approx = latest_approx;
		} //end while	

		NotifyObservers(Converged<EmitterType>(*this));
		return SuccessCode::Success;

	} //end PSEG


	// ================================================================================================
	//   Adaptive-numeric-type (double-first) PowerSeries endgame.  Mirrors the Cauchy flavor: compute in
	//   the hardware complex_dbl fast lane while the tracker's authoritative precision stays double, and
	//   cross to complex_mp only when the tracker escalates (pure-(i)).  Member templates so the explicit
	//   fixed-precision class instantiations never force-compile them.  Fixed precision uses RunImpl<BCT>.
	// ================================================================================================

	/// \brief Run the power-series endgame via the double-first adaptive-numeric-type driver, migrating up to mpfr only on tracker escalation.
	template<typename Dummy = void>
	SuccessCode RunImplAMP(BCT const& start_time, Vec<BCT> const& start_point, BCT const& target_time)
	{
		using bertini::Precision;
		using RealT = typename Eigen::NumTraits<BCT>::Real;

		if (start_point.size()!=static_cast<Eigen::Index>(this->GetSystem().NumVariables()))
		{
			std::stringstream err_msg;
			err_msg << "number of variables in start point for PSEG, " << start_point.size() << ", must match the number of variables in the system, " << this->GetSystem().NumVariables();
			throw std::runtime_error(err_msg.str());
		}

		this->adaptive_numeric_type_active_ = true;
		struct Disarmer { bool& flag; ~Disarmer(){ flag = false; } } disarm{this->adaptive_numeric_type_active_};

		this->current_endgame_precision_ = std::max(DoublePrecision(), Precision(start_point));

		// ---- SETUP (restart-in-mpfr on escalation; bounded) ----
		while (true)
		{
			SuccessCode code = (this->current_endgame_precision_==DoublePrecision())
				? SetupSegmentT<complex_dbl>(complex_dbl(start_time), this->DowncastToDouble(start_point), complex_dbl(target_time))
				: SetupSegmentT<complex_mp>(this->AtActivePrecisionScalar(start_time), this->AtActivePrecisionVec(start_point), this->AtActivePrecisionScalar(target_time));
			if (code == SuccessCode::HigherPrecisionNecessary)
			{
				this->current_endgame_precision_ = this->NextEscalatedPrecision();
				SetThreadPrecision(this->current_endgame_precision_);
				this->GetSystem().precision(this->current_endgame_precision_);
				continue;
			}
			if (code != SuccessCode::Success)
				return code;
			break;
		}

		// previous_approximation_ (BCT) holds the first extrapolation, set by SetupSegmentT.
		RealT norm_prev(0), norm_latest(0);
		if (this->SecuritySettings().level <= 0)
			norm_prev = this->GetSystem().InfinityNormOfDehomogenized(this->previous_approximation_);

		this->approximate_error_ = 1;

		// ---- MAIN LOOP: each tracker-touching phase self-heals (migrate-and-retry in mpfr) ----
		while (this->approximate_error_ > this->FinalTolerance())
		{
			auto adv = AdvanceTimeAMP(target_time);
			if (adv != SuccessCode::Success) { NotifyObservers(EndgameFailure<EmitterType>(*this)); return adv; }

			auto ref = RefineAllSamplesAMP();
			if (ref != SuccessCode::Success) return ref;

			ComputeAllDerivativesAMP();

			auto ext = ComputeApproxAMP(target_time);   // -> final_approximation_ (BCT)
			if (ext != SuccessCode::Success) { NotifyObservers(EndgameFailure<EmitterType>(*this)); return ext; }

			Precision(this->previous_approximation_, Precision(this->final_approximation_));
			this->approximate_error_ = static_cast<NumErrorT>((this->final_approximation_ - this->previous_approximation_).template lpNorm<Eigen::Infinity>());
			NotifyObservers(ApproximatedRoot<EmitterType>(*this));

			if (this->SecuritySettings().level <= 0)
			{
				norm_latest = this->GetSystem().InfinityNormOfDehomogenized(this->final_approximation_);
				if (this->BeyondSecurityMaxNorm(norm_latest) && this->BeyondSecurityMaxNorm(norm_prev))
				{
					NotifyObservers(SecurityMaxNormReached<EmitterType>(*this));
					return SuccessCode::SecurityMaxNormReached;
				}
				norm_prev = norm_latest;
			}

			this->previous_approximation_ = this->final_approximation_;
		}

		NotifyObservers(Converged<EmitterType>(*this));
		return SuccessCode::Success;
	}


	// Pre-loop work + first extrapolation, at one numeric type.  Reports HigherPrecisionNecessary up to
	// the driver (which restarts setup in mpfr) on escalation.
	/// \brief Run the pre-loop setup and first extrapolation in a given numeric type, reporting escalation up to the driver.
	template<typename ComplexT>
	SuccessCode SetupSegmentT(ComplexT const& start_time, Vec<ComplexT> const& start_point, ComplexT const& target_time)
	{
		SetThreadPrecision(Precision(start_point));
		ClearTimesAndSamples<ComplexT>();

		auto& samples = std::get<SampCont<ComplexT> >(samples_);
		auto& times   = std::get<TimeCont<ComplexT> >(times_);

		SetRandVec<ComplexT>(static_cast<int>(start_point.size()));

		auto init = this->ComputeInitialSamples(start_time, target_time, start_point, times, samples);
		if (init != SuccessCode::Success) { NotifyObservers(EndgameFailure<EmitterType>(*this)); return init; }
		if (this->GetTracker().GetCurrentPrecision() > this->current_endgame_precision_)
			return SuccessCode::HigherPrecisionNecessary;

		auto ref = this->template RefineAllSamples<ComplexT>(samples, times);
		if (ref == SuccessCode::HigherPrecisionNecessary || ref == SuccessCode::FailedToConverge)
			return SuccessCode::HigherPrecisionNecessary;
		if (ref != SuccessCode::Success) return ref;

		ComputeAllDerivatives<ComplexT>();

		Vec<ComplexT> first_approx;
		auto ext = ComputeApproximationOfXAtT0<ComplexT>(first_approx, target_time);
		if (ext != SuccessCode::Success) return ext;
		this->previous_approximation_ = this->ToBCT(first_approx, this->current_endgame_precision_);
		return SuccessCode::Success;
	}


	/// \brief Advance time at the active numeric type, migrating-and-retrying in mpfr on escalation.
	template<typename Dummy = void>
	SuccessCode AdvanceTimeAMP(BCT const& target_time)
	{
		unsigned guard = 0;
		while (true)
		{
			SuccessCode code = (this->current_endgame_precision_==DoublePrecision())
				? AdvanceTime<complex_dbl>(complex_dbl(target_time))
				: AdvanceTime<complex_mp>(this->AtActivePrecisionScalar(target_time));
			if (code == SuccessCode::HigherPrecisionNecessary)
			{
				if (this->template EscalateAndMigrate<>(++guard) != SuccessCode::Success) return SuccessCode::HigherPrecisionNecessary;
				continue;
			}
			return code;
		}
	}

	/// \brief Refine all retained samples at the active numeric type, migrating-and-retrying in mpfr on escalation.
	template<typename Dummy = void>
	SuccessCode RefineAllSamplesAMP()
	{
		unsigned guard = 0;
		while (true)
		{
			SuccessCode code = (this->current_endgame_precision_==DoublePrecision())
				? this->template RefineAllSamples<complex_dbl>(std::get<SampCont<complex_dbl> >(samples_), std::get<TimeCont<complex_dbl> >(times_))
				: this->template RefineAllSamples<complex_mp>(std::get<SampCont<complex_mp> >(samples_), std::get<TimeCont<complex_mp> >(times_));
			if (code == SuccessCode::HigherPrecisionNecessary || code == SuccessCode::FailedToConverge)
			{
				if (this->template EscalateAndMigrate<>(++guard) != SuccessCode::Success) return SuccessCode::HigherPrecisionNecessary;
				continue;
			}
			return code;
		}
	}

	/// \brief Compute all sample derivatives at the active numeric type.
	template<typename Dummy = void>
	void ComputeAllDerivativesAMP()
	{
		if (this->current_endgame_precision_==DoublePrecision())
			ComputeAllDerivatives<complex_dbl>();
		else
			ComputeAllDerivatives<complex_mp>();
	}

	/// \brief Compute the power-series extrapolation at the active numeric type, writing into final_approximation_ (BCT).
	template<typename Dummy = void>
	SuccessCode ComputeApproxAMP(BCT const& target_time)
	{
		if (this->current_endgame_precision_==DoublePrecision())
		{
			Vec<complex_dbl> r;
			auto code = ComputeApproximationOfXAtT0<complex_dbl>(r, complex_dbl(target_time));
			if (code == SuccessCode::Success)
				this->final_approximation_ = this->ToBCT(r, this->current_endgame_precision_);
			return code;
		}
		else
		{
			return ComputeApproximationOfXAtT0<complex_mp>(this->final_approximation_, this->AtActivePrecisionScalar(target_time));
		}
	}


	// Flavor-specific: cross every PowerSeries container from the complex_dbl slot to complex_mp (or, if
	// already mpfr, raise its precision uniformly), via the shared AMP-policy Cross* / SetPrecision
	// helpers.  Called by the base EscalateAndMigrate.
	/// \brief Widen this flavor's PowerSeries containers from the complex_dbl slot up to complex_mp at the new precision.
	template<typename Dummy = void>
	void MigrateContainersToPrecision(unsigned newprec)
	{
		using bertini::Precision;
		if (this->current_endgame_precision_ == DoublePrecision())
		{
			this->CrossTimesUp(times_,       newprec);
			this->CrossSampsUp(samples_,     newprec);
			this->CrossSampsUp(derivatives_, newprec);
			this->CrossVecUp  (rand_vector_, newprec);
		}
		else
		{
			tracking::adaptive::SetPrecision(std::get<TimeCont<complex_mp> >(times_),       newprec);
			tracking::adaptive::SetPrecision(std::get<SampCont<complex_mp> >(samples_),     newprec);
			tracking::adaptive::SetPrecision(std::get<SampCont<complex_mp> >(derivatives_), newprec);
			auto& rv = std::get<Vec<complex_mp> >(rand_vector_);
			if (rv.size() > 0) Precision(rv, newprec);
		}
		if (this->final_approximation_.size()    > 0) Precision(this->final_approximation_,    newprec);
		if (this->previous_approximation_.size() > 0) Precision(this->previous_approximation_, newprec);
		this->GetSystem().precision(newprec);
	}


	virtual ~PowerSeriesEndgame() = default;
}; // end powerseries class




}} // re: namespaces
