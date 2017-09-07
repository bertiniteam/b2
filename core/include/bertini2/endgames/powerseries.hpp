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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University


#pragma once

#include "bertini2/endgames/base_endgame.hpp"


namespace bertini{ namespace endgame{


/** 
\class PowerSeriesEndgame

\brief class used to finish tracking paths during Homotopy Continuation


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
using RealT = tracking::TrackerTraits<TrackerType>::BaseRealType; // Real types
using ComplexT = tracking::TrackerTraits<TrackerType>::BaseComplexType; Complex types

// 1. Define the polynomial system that we wish to solve. 
System target_sys;
Var x = MakeVariable("x"), t = MakeVariable("t"), y = MakeVariable("y");

VariableGroup vars{x,y};
target_sys.AddVariableGroup(vars); 

target_sys.AddFunction((pow(x-1,3));
target_sys.AddFunction((pow(y-1,2));

// 1b. Homogenize and patch the polynomial system to work over projective space. 
sys.Homogenize();
sys.AutoPatch();

// 2. Create a start system, for us we will use a total degree start system.
auto TD_start_sys = bertini::start_system::TotalDegree(target_sys);

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
    mpfr_float::default_precision(ambient_precision);
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
	public EndgameBase<PowerSeriesEndgame<PrecT>, PrecT>
{
public:
	using BaseEGT = EndgameBase<PowerSeriesEndgame<PrecT>, PrecT>;
	using FinalEGT = PowerSeriesEndgame<PrecT>;
	using TrackerType = typename BaseEGT::TrackerType;

	using BaseComplexType = typename BaseEGT::BaseComplexType;
	using BaseRealType = typename BaseEGT::BaseRealType;


protected:
	

	using TupleOfTimes = typename BaseEGT::TupleOfTimes;
	using TupleOfSamps = typename BaseEGT::TupleOfSamps;

	using BCT = BaseComplexType;
	using BRT = BaseRealType;

	using Configs = typename AlgoTraits<FinalEGT>::NeededConfigs;
	using ConfigsAsTuple = typename Configs::ToTuple;

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
	*/
	mutable Vec<BCT> rand_vector_;

	template<typename CT>
	void AssertSizesTimeSpace() const
	{
		const auto num_sample_points = this->EndgameSettings().num_sample_points;
		assert(std::get<SampCont<CT> >(samples_).size()==std::get<TimeCont<CT> >(times_).size() && "must have same number of samples in times and spaces");
		assert(std::get<SampCont<CT> >(samples_).size()>=num_sample_points && "must have sufficient number of samples");
	}

	template<typename CT>
	void AssertSizesTimeSpaceDeriv() const
	{
		const auto num_sample_points = this->EndgameSettings().num_sample_points;
		assert(std::get<SampCont<CT> >(samples_).size()==std::get<TimeCont<CT> >(times_).size() && "must have same number of samples in times and spaces");
		assert(std::get<SampCont<CT> >(samples_).size()==std::get<SampCont<CT> >(derivatives_).size() && "must have same number of samples in derivatives and spaces");
		assert(std::get<SampCont<CT> >(samples_).size()>=num_sample_points && "must have sufficient number of samples");
	}

public:

	auto UpperBoundOnCycleNumber() const { return upper_bound_on_cycle_number_;}


	/**
	\brief Function that clears all samples and times from data members for the Power Series endgame
	*/	
	template<typename CT>
	void ClearTimesAndSamples()
	{
		std::get<TimeCont<CT> >(times_).clear(); 
		std::get<SampCont<CT> >(samples_).clear();
	}

	/**
	\brief Function to set the times used for the Power Series endgame.
	*/	
	template<typename CT>
	void SetTimes(TimeCont<CT> times_to_set) { std::get<TimeCont<CT> >(times_) = times_to_set;}

	/**
	\brief Function to get the times used for the Power Series endgame.
	*/	
	template<typename CT>
	const auto& GetTimes() const {return std::get<TimeCont<CT> >(times_);}


	const BCT& LatestTimeImpl() const
	{
		return GetTimes<BCT>().back();
	}

	/**
	\brief Function to set the space values used for the Power Series endgame.
	*/	
	template<typename CT>
	void SetSamples(SampCont<CT> samples_to_set) { std::get<SampCont<CT> >(samples_) = samples_to_set;}

	/**
	\brief Function to get the space values used for the Power Series endgame.
	*/	
	template<typename CT>
	const auto& GetSamples() const {return std::get<SampCont<CT> >(samples_);}

	/**
	\brief Function to set the times used for the Power Series endgame.
	*/	
	template<typename CT>
	void SetRandVec(int size) {rand_vector_ = Vec<CT>::Random(size);}





	explicit PowerSeriesEndgame(TrackerType const& tr, 
	                            const ConfigsAsTuple& settings )
      : BaseEGT(tr, settings)
   	{}

    template< typename... Ts >
		PowerSeriesEndgame(TrackerType const& tr, const Ts&... ts ) : PowerSeriesEndgame(tr, Configs::Unpermute( ts... ) ) 
		{}


	~PowerSeriesEndgame() {};


	/**
	\brief Computes an upper bound on the cycle number. Consult page 53 of \cite bertinibook.

	## Input: 
			None: all data needed are class data members.

	## Output:
			upper_bound_on_cycle_number_: Used for an exhaustive search for the best cycle number for approimating the path to t = 0.

	##Details:
			\tparam CT The complex number type.
	*/
	template<typename CT>
	unsigned ComputeBoundOnCycleNumber()
	{ 
		using RT = typename Eigen::NumTraits<CT>::Real;
		using std::log; using std::abs;

		const auto& samples = std::get<SampCont<CT> >(samples_);
		
		AssertSizesTimeSpace<CT>();

		auto num_samples = samples.size();
		const Vec<CT> & sample0 = samples[num_samples-3];
		const Vec<CT> & sample1 = samples[num_samples-2];
		const Vec<CT> & sample2 = samples[num_samples-1]; // most recent sample.  oldest samples at front of the container


		CT rand_sum1 = ((sample1 - sample0).transpose()*rand_vector_).sum();
		CT rand_sum2 = ((sample2 - sample1).transpose()*rand_vector_).sum();

		if ( abs(rand_sum1)==0 || abs(rand_sum2)==0) // avoid division by 0
		{
			upper_bound_on_cycle_number_ = 1;
			return upper_bound_on_cycle_number_;
		}

		RT estimate = log(static_cast<RT>(this->EndgameSettings().sample_factor))/log(abs(rand_sum2/rand_sum1));

		if (estimate < 1) // would be nan if sample points are same as each other
		  	upper_bound_on_cycle_number_ = 1;
		else
		{
			using std::max;
			auto upper_bound = unsigned(round(estimate)*this->template Get<PowerSeriesConfig>().cycle_number_amplification);
			upper_bound_on_cycle_number_ = max(upper_bound,this->template Get<PowerSeriesConfig>().max_cycle_number);
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
				\tparam CT The complex number type.
			This is done by an exhaustive search from 1 to upper_bound_on_cycle_number. There is a conversion to the s-space from t-space in this function. 
	As a by-product the derivatives at each of the samples is returned for further use. 
	*/

	template<typename CT>
	unsigned ComputeCycleNumber(CT const& t0)
	{
		using RT = typename Eigen::NumTraits<CT>::Real;

		const auto& samples = std::get<SampCont<CT> >(samples_);
		const auto& times   = std::get<TimeCont<CT> >(times_);
		const auto& derivatives = std::get<SampCont<CT> >(derivatives_);

		AssertSizesTimeSpaceDeriv<CT>();
		
		const Vec<CT> &most_recent_sample = samples.back();  
		const CT& most_recent_time = times.back();

		//Compute upper bound for cycle number.
		ComputeBoundOnCycleNumber<CT>();


		unsigned num_pts;
		if (samples.size() > this->EndgameSettings().num_sample_points)
			num_pts = this->EndgameSettings().num_sample_points;
		else 
			num_pts = this->EndgameSettings().num_sample_points-1;


		auto min_found_difference = Eigen::NumTraits<RT>::highest();

		TimeCont<CT> s_times(num_pts);
		SampCont<CT> s_derivatives(num_pts);

		auto offset = samples.size() - num_pts - 1; // -1 here to shift away from the back of the container
		for(unsigned int candidate = 1; candidate <= upper_bound_on_cycle_number_; ++candidate)
		{			
			// BOOST_LOG_TRIVIAL(severity_level::trace) << "testing cycle candidate " << candidate;

			std::tie(s_times, s_derivatives) = TransformToSPlane(candidate, t0, num_pts, ContStart::Front);

//TODO move these to an observer on the endgame.  of course, you have to implement the event types first, ha.	
// for (int ii=0; ii<times.size(); ++ii)
// {
// 	BOOST_LOG_TRIVIAL(severity_level::trace) << "t-plane data " << ii << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(times[ii])) << " time " << times[ii] << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(samples[ii])) << " sample " << samples[ii] << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(derivatives[ii])) << " derivative " << derivatives[ii] << "\n\n";
// }

// for (int ii=0; ii<s_times.size(); ++ii)
// {
// 	BOOST_LOG_TRIVIAL(severity_level::trace) << "s-plane data " << ii << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(s_times[ii])) << " s_time " << s_times[ii] << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(s_derivatives[ii])) << " s_derivative" << s_derivatives[ii] << "\n\n\n";
// }



			RT curr_diff = (HermiteInterpolateAndSolve(
								  pow((most_recent_time-t0)/(times[0]-t0), 1/static_cast<RT>(candidate)), // the target time
			                      num_pts,s_times,samples,s_derivatives, ContStart::Front) // the input data
			                 - 
			                 most_recent_sample).template lpNorm<Eigen::Infinity>();
			// BOOST_LOG_TRIVIAL(severity_level::trace) << "residual " << curr_diff;
			if (curr_diff < min_found_difference)
			{
				min_found_difference = curr_diff;
				this->cycle_number_ = candidate;
			}

		}// end cc loop over cycle number possibilities
		// BOOST_LOG_TRIVIAL(severity_level::trace) << "cycle number computed to be " << this->CycleNumber();
		return this->cycle_number_;
	}//end ComputeCycleNumber



	/**
		\brief Compute a set of derivatives using internal data to the endgame.

		## Input: 
				None: all data needed are class data members.

		## Output:
				None: Derivatives are members of this class.

		##Details:
				\tparam CT The complex number type.
	*/
	template<typename CT>
	void ComputeAllDerivatives()
	{
		auto& samples = std::get<SampCont<CT> >(samples_);
		auto& times   = std::get<TimeCont<CT> >(times_);
		auto& derivatives = std::get<SampCont<CT> >(derivatives_);

		assert((samples.size() == times.size()) && "must have same number of times and samples");

		if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec) // known at compile time
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
	\param c The cycle number you want to use

	This function transforms the times and derivatives into the S-plane, scaled by the cycle number.

	this function also transforms them into the interval [0 1]
	*/
	template <typename CT>
	std::tuple<TimeCont<CT>, SampCont<CT>> TransformToSPlane(int cycle_num, CT const& t0, unsigned num_pts, ContStart shift_from)
	{
		if (cycle_num==0)
			throw std::runtime_error("cannot transform to s plane with cycle number 0");
		AssertSizesTimeSpaceDeriv<CT>();


		using RT = typename Eigen::NumTraits<CT>::Real;

		const auto& times   = std::get<TimeCont<CT> >(times_);
		const auto& derivatives  = std::get<SampCont<CT> >(derivatives_);

		RT c = static_cast<RT>(cycle_num);
		RT one_over_c = 1/c;

		unsigned offset_t, offset_d;
		if (shift_from == ContStart::Back)
		{
			offset_t = times.size()-num_pts;
			offset_d = derivatives.size()-num_pts;
		}
		else
			offset_t = offset_d = 0;


		TimeCont<CT> s_times(num_pts);
		SampCont<CT> s_derivatives(num_pts);

		CT time_shift = times[offset_t] - t0;

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
	\tparam CT The complex number type.
				This function handles computing an approximation at the origin. 
				We compute the cycle number best for the approximation, and convert derivatives and times to the s-plane where s = t^(1/c).
				We use the converted times and derivatives along with the samples to do a Hermite interpolation.
	*/
	template<typename CT>
	SuccessCode ComputeApproximationOfXAtT0(Vec<CT>& result, const CT & t0)
	{	
		const auto c = ComputeCycleNumber<CT>(t0);

		auto num_pts = this->EndgameSettings().num_sample_points;

		TimeCont<CT> s_times;
		SampCont<CT> s_derivatives;

		std::tie(s_times, s_derivatives) = TransformToSPlane(c, t0, num_pts, ContStart::Back);

//TODO move these to an observer on the endgame.  of course, you have to implement the event types first, ha.	
// const auto& times   = std::get<TimeCont<CT> >(times_);
// const auto& samples  = std::get<SampCont<CT> >(samples_);
// const auto& derivatives  = std::get<SampCont<CT> >(derivatives_);
// for (int ii=0; ii<times.size(); ++ii)
// {
// 	BOOST_LOG_TRIVIAL(severity_level::trace) << "t-plane data " << ii << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(times[ii])) << " time " << times[ii] << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(samples[ii])) << " sample " << samples[ii] << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(derivatives[ii])) << " derivative " << derivatives[ii] << "\n\n";
// }

// for (int ii=0; ii<s_times.size(); ++ii)
// {
// 	BOOST_LOG_TRIVIAL(severity_level::trace) << "s-plane data " << ii << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(s_times[ii])) << " s_time " << s_times[ii] << '\n';
// 		BOOST_LOG_TRIVIAL(severity_level::trace) << std::setprecision(Precision(s_derivatives[ii])) << " s_derivative" << s_derivatives[ii] << "\n\n\n";
// }

		// the data was transformed to be on the interval [0 1] so we can hard-code the time-to-solve as 0 here.

		Precision(result, Precision(s_derivatives.back()));
		result = HermiteInterpolateAndSolve(CT(0), num_pts, s_times, std::get<SampCont<CT> >(samples_), s_derivatives, ContStart::Back);
		return SuccessCode::Success;
	}//end ComputeApproximationOfXAtT0



	/**
		\brief The samples used in the power series endgame are collected by advancing time to t = 0, by multiplying the current time by the sample factor. 

		## Input: 
				target_time: This is the time we are trying to approximate, default is t = 0.

		## Output:
				SuccessCode: This reports back if we were successful in advancing time. 

		##Details:
				\tparam CT The complex number type.
				This function computes the next time value for the power series endgame. After computing this time value, 
				it will track to it and compute the derivative at this time value for further appoximations to be made during the
				endgame.
	*/
	template<typename CT>
	SuccessCode AdvanceTime(const CT & target_time)
	{
		using RT = typename Eigen::NumTraits<CT>::Real;
		
		auto& samples = std::get<SampCont<CT> >(samples_);
		auto& times   = std::get<TimeCont<CT> >(times_);
		auto& derivatives  = std::get<SampCont<CT> >(derivatives_);

		AssertSizesTimeSpaceDeriv<CT>();

		Vec<CT> next_sample;
		CT next_time = (times.back() + target_time) * static_cast<RT>(this->EndgameSettings().sample_factor); //setting up next time value using the midpoint formula, sample_factor will give us some 

  		if (abs(next_time - target_time) < this->EndgameSettings().min_track_time) // generalized for target_time not equal to 0.
  		{
//TODO move these to an observer on the endgame.  of course, you have to implement the event types first, ha.	
  			// BOOST_LOG_TRIVIAL(severity_level::trace) << "Current time norm is less than min track time." << '\n';
  			return SuccessCode::MinTrackTimeReached;
  		}

//TODO move these to an observer on the endgame.  of course, you have to implement the event types first, ha.	
  		// BOOST_LOG_TRIVIAL(severity_level::trace) << "tracking to t = " << next_time << ", default precision: " << DefaultPrecision() << "\n";
		SuccessCode tracking_success = this->GetTracker().TrackPath(next_sample,times.back(),next_time,samples.back());
			if (tracking_success != SuccessCode::Success)
				return tracking_success;

		this->EnsureAtPrecision(next_time,Precision(next_sample));
	

		times.push_back(next_time);
		samples.push_back(next_sample);

	


		auto refine_success = this->RefineSample(samples.back(), next_sample,  times.back(), 
										this->FinalTolerance() * this->EndgameSettings().sample_point_refinement_factor,
										this->EndgameSettings().max_num_newton_iterations);
			if (refine_success != SuccessCode::Success)
			{
				// BOOST_LOG_TRIVIAL(severity_level::trace) << "refining failed, code " << int(refine_success);
				return refine_success;
			}
		
		this->EnsureAtPrecision(times.back(),Precision(samples.back()));

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
	\tparam CT The complex number type.
				Tracking forward with the number of sample points, this function will make approximations using Hermite interpolation. This process will continue until two consecutive
				approximations are withing final tolerance of each other. 
	*/		
	template<typename CT>
	SuccessCode RunImpl(const CT & start_time, const Vec<CT> & start_point, CT const& target_time)
	{
		if (start_point.size()!=this->GetSystem().NumVariables())
		{
			std::stringstream err_msg;
			err_msg << "number of variables in start point for PSEG, " << start_point.size() << ", must match the number of variables in the system, " << this->GetSystem().NumVariables();
			throw std::runtime_error(err_msg.str());
		}

		// BOOST_LOG_TRIVIAL(severity_level::trace) << "\n\nPSEG(), default precision: " << DefaultPrecision() << "\n\n";
		// BOOST_LOG_TRIVIAL(severity_level::trace) << "start point precision: " << Precision(start_point) << "\n\n";

		DefaultPrecision(Precision(start_point));

		using RT = typename Eigen::NumTraits<CT>::Real;
		//Set up for the endgame.
		ClearTimesAndSamples<CT>();

		auto& samples = std::get<SampCont<CT> >(samples_);
		auto& times   = std::get<TimeCont<CT> >(times_);
		auto& derivatives  = std::get<SampCont<CT> >(derivatives_);
		Vec<CT>& latest_approx = std::get<Vec<CT> >(this->final_approximation_);
		Vec<CT>& prev_approx = std::get<Vec<CT> >(this->previous_approximation_);


		SetRandVec<CT>(start_point.size());	 	
		
		
		auto initial_sample_success = this->ComputeInitialSamples(start_time, target_time, start_point, times, samples);

		if (initial_sample_success!=SuccessCode::Success)
		{
			// BOOST_LOG_TRIVIAL(severity_level::trace) << "initial sample gathering failed, code " << int(initial_sample_success) << std::endl;
			return initial_sample_success;
		}

		this->template RefineAllSamples<CT>(samples, times);
		ComputeAllDerivatives<CT>();


	   


	 	auto extrapolation_code = ComputeApproximationOfXAtT0(prev_approx, target_time);
	 	latest_approx = prev_approx;

	 	if (extrapolation_code != SuccessCode::Success)
	 		return extrapolation_code;

	 	RT norm_of_dehom_of_latest_approx;
	 	RT norm_of_dehom_of_prev_approx;
	 	if (this->SecuritySettings().level <= 0)
	 	 	norm_of_dehom_of_prev_approx = this->GetSystem().DehomogenizePoint(prev_approx).template lpNorm<Eigen::Infinity>();

	 	

	  	

	 	NumErrorT& approx_error = this->approximate_error_;
		approx_error = 1;

		while (approx_error > this->FinalTolerance())
		{
	  		auto advance_code = AdvanceTime<CT>(target_time);
	  		if (advance_code!=SuccessCode::Success)
	 		{
	 			// BOOST_LOG_TRIVIAL(severity_level::trace) << "unable to advance time, code " << int(advance_code);
	 			return advance_code;
	 		}

	 		// this code is what bertini1 does... it refines all samples, like, all the time.
	 		this->template RefineAllSamples<CT>(samples, times);
	 		ComputeAllDerivatives<CT>();

	 		extrapolation_code = ComputeApproximationOfXAtT0(latest_approx, target_time);
	 		if (extrapolation_code!=SuccessCode::Success)
	 		{
	 			// BOOST_LOG_TRIVIAL(severity_level::trace) << "failed to compute the approximation at " << target_time << "\n\n";
	 			return extrapolation_code;
	 		}
	 		// BOOST_LOG_TRIVIAL(severity_level::trace)  << std::setprecision(Precision(latest_approx)) << "latest approximation:\n" << latest_approx << '\n';

	 		if(this->SecuritySettings().level <= 0)
	 		{
	 			norm_of_dehom_of_latest_approx = this->GetSystem().DehomogenizePoint(latest_approx).template lpNorm<Eigen::Infinity>();
		 		if(norm_of_dehom_of_latest_approx > this->SecuritySettings().max_norm && norm_of_dehom_of_prev_approx > this->SecuritySettings().max_norm)
	 				return SuccessCode::SecurityMaxNormReached;

	 			//update
	 			norm_of_dehom_of_prev_approx = norm_of_dehom_of_latest_approx;
	 		}

	 		//update
	 		approx_error = static_cast<NumErrorT>((latest_approx - prev_approx).template lpNorm<Eigen::Infinity>());
	 		// BOOST_LOG_TRIVIAL(severity_level::trace) << "consecutive approximation error:\n" << approx_error << '\n';

	 		Precision(prev_approx, Precision(latest_approx));
	 		prev_approx = latest_approx;

		} //end while	

		return SuccessCode::Success;

	} //end PSEG

}; // end powerseries class




}} // re: namespaces
