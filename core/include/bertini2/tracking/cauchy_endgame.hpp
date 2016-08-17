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
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Tim Hodges, Colorado State University



#pragma once

#include "bertini2/tracking/base_endgame.hpp"


namespace bertini{ namespace tracking { namespace endgame{


/** 
\class CauchyEndgame
\brief Class used to finish tracking paths during Homotopy Continuation.

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
using RealT = TrackerTraits<TrackerType>::BaseRealType; // Real types
using ComplexT = TrackerTraits<TrackerType>::BaseComplexType; Complex types

// 1. Define the polynomial system that we wish to solve. 
System target_sys;
Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t"), y = std::make_shared<Variable>("y");

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
config::Stepping<RealT> stepping_preferences;
stepping_preferences.initial_step_size = RealT(1)/RealT(5);// change a stepping preference
config::Newton newton_preferences;
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

config::Tolerances<RealT> tolerances;
tolerances.final_tolerance_multiplier = RealT(100);

config::Cauchy<RealT> cauchy_settings;
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
template<typename TrackerType, typename FinalEGT, typename... UsedNumTs> 
class CauchyEndgame : public EndgameBase<TrackerType, FinalEGT>
{
private:
	// convert the base endgame into the derived type.
	const FinalEGT& AsDerived() const
	{
		return static_cast<const FinalEGT&>(*this);
	}


protected:

	using BaseComplexType = typename TrackerTraits<TrackerType>::BaseComplexType;
	using BaseRealType = typename TrackerTraits<TrackerType>::BaseRealType;

	using BCT = BaseComplexType;
	using BRT = BaseRealType;

	/**
	\brief Settings that are specific to the Cauchy endgame. 
	*/	
	config::Cauchy<BRT> cauchy_settings_; 

	/**
	\brief A deque of times that are specifically used to compute the power series approximation for the Cauchy endgame. 
	*/
	mutable std::tuple<TimeCont<UsedNumTs>...> pseg_times_;
	/**
	\brief A deque of samples that are in correspondence with the pseg_times_. These samples are also used to compute the first power series approximation for the Cauchy endgame.
	*/
	mutable std::tuple<SampCont<UsedNumTs>...> pseg_samples_; //samples used for the first approximation. 
	/**
	\brief A deque of times that are collected by CircleTrack. These samples are used to compute all approximations of the origin after the first. 
	*/
	mutable std::tuple<TimeCont<UsedNumTs>...> cauchy_times_;
	/**
	\brief A deque of samples collected by CircleTrack. Computed a mean of the values of this deque, after a loop has been closed, will give an approximation of the origin.
	*/
	mutable std::tuple<SampCont<UsedNumTs>...> cauchy_samples_;



	
	

public:
	/**
	\brief Function that clears all samples and times from data members for the Cauchy endgame
	*/
	template<typename CT>
	void ClearTimesAndSamples()
	{
		std::get<TimeCont<CT> >(pseg_times_).clear(); 
		std::get<TimeCont<CT> >(cauchy_times_).clear(); 
		std::get<SampCont<CT> >(pseg_samples_).clear(); 
		std::get<SampCont<CT> >(cauchy_samples_).clear();}
	/**
	\brief Setter for the deque holding time values for the power series approximation of the Cauchy endgame. 
	*/
	template<typename CT>
	void SetPSEGTimes(TimeCont<CT> pseg_times_to_set) 
	{ std::get<TimeCont<CT> >(pseg_times_) = pseg_times_to_set;}

	/**
	\brief Getter for the deque holding time values for the power series approximation of the Cauchy endgame.
	*/
	template<typename CT>
	TimeCont<CT>& GetPSEGTimes() {return std::get<TimeCont<CT> >(pseg_times_);}

	/**
	\brief Setter for the deque holding space values for the power series approximation of the Cauchy endgame. 
	*/
	template<typename CT>
	void SetPSEGSamples(SampCont<CT> pseg_samples_to_set) { std::get<SampCont<CT> >(pseg_samples_) = pseg_samples_to_set;}

	/**
	\brief Getter for the deque holding space values for the power series approximation of the Cauchy endgame. 
	*/
	template<typename CT>
	SampCont<CT>& GetPSEGSamples() {return std::get<SampCont<CT> >(pseg_samples_);}

	/**
	\brief Setter for the deque holding space values for the Cauchy endgame. 
	*/
	template<typename CT>
	void SetCauchySamples(SampCont<CT> cauchy_samples_to_set) { std::get<SampCont<CT> >(cauchy_samples_) = cauchy_samples_to_set;}

	/**
	\brief Getter for the deque holding space values for the Cauchy endgame. 
	*/
	template<typename CT>
	SampCont<CT>& GetCauchySamples() {return std::get<SampCont<CT> >(cauchy_samples_);}

	/**
	\brief Setter for the deque holding time values for the Cauchy endgame. 
	*/
	template<typename CT>
	void SetCauchyTimes(TimeCont<CT> cauchy_times_to_set) { std::get<TimeCont<CT> >(cauchy_times_) = cauchy_times_to_set;}

	/**
	\brief Getter for the deque holding time values for the Cauchy endgame. 
	*/
	template<typename CT>
	TimeCont<CT>& GetCauchyTimes() {return std::get<TimeCont<CT> >(cauchy_times_);}

	
	/**
	\brief Setter for the specific settings in tracking_conifg.hpp under Cauchy.
	*/
	void SetCauchySettings(config::Cauchy<BRT> new_cauchy_settings){cauchy_settings_ = new_cauchy_settings;}

	/**
	\brief Getter for the specific settings in tracking_conifg.hpp under Cauchy.
	*/
	auto& GetCauchySettings(){return cauchy_settings_;}
	

	explicit CauchyEndgame(TrackerType const& tr, 
	                            const std::tuple< config::Cauchy<BRT> const&, 
	                            					const config::Endgame<BRT>&, 
	                            					const config::Security<BRT>&, 
	                            					const config::Tolerances<BRT>& 
	                            				>& settings )
      : EndgameBase<TrackerType, FinalEGT>(tr, std::get<1>(settings), std::get<2>(settings), std::get<3>(settings) ), 
          cauchy_settings_( std::get<0>(settings) )
   	{ }

    template< typename... Ts >
		CauchyEndgame(TrackerType const& tr, const Ts&... ts ) : CauchyEndgame(tr, Unpermute<config::Cauchy<BRT>, config::Endgame<BRT>, config::Security<BRT>, config::Tolerances<BRT> >( ts... ) ) 
		{}


	~CauchyEndgame() {};



	/**
		\brief Function to track around the origin 

		## Input: 
				starting_time: time value that we start from to track around the origin
				target_time: the time value that we are centering out loops around, default is t = 0
				starting_sample: an approximate solution of the homotopy at t = starting_time

		## Output:
				SuccessCode: This reports back if we were successful in advancing time. 


		##Details:
	\tparam CT The complex number type.
				Depeding on the number of samples points, we make a polgon around the origin with that many vertices. This function should be called the same number of times  
				as paths converging to the solution we are approximating. 
	*/
	template<typename CT> 
	SuccessCode CircleTrack(CT const& starting_time, CT const& target_time, Vec<CT> const& starting_sample)
	{	
		using bertini::Precision;
		assert(Precision(starting_time)==Precision(starting_sample) && "starting time and sample for circle track must be of same precision");
		DefaultPrecision(Precision(starting_time));

		using RT = typename Eigen::NumTraits<CT>::Real;
		using std::acos;

		if (this->EndgameSettings().num_sample_points < 3) // need to make sure we won't track right through the origin.
		{
			std::stringstream err_msg;
			err_msg << "ERROR: The number of sample points " << this->EndgameSettings().num_sample_points << " for circle tracking must be >= 3";
			throw std::runtime_error(err_msg.str());
		}	

		auto& circle_times = std::get<TimeCont<CT> >(cauchy_times_);
		auto& circle_samples = std::get<SampCont<CT> >(cauchy_samples_);

		// the initial sample has already been added to the sample repo... so don't do that here, please
		
		const auto num_vars = this->GetSystem().NumVariables();

		for (unsigned ii = 0; ii < this->EndgameSettings().num_sample_points; ++ii)
		{
			const Vec<CT>& current_sample = circle_samples.back();
			const CT& current_time = circle_times.back();
			assert(Precision(current_time)==Precision(current_sample) && "current time and sample for circle track must be of same precision");

			//set up the time value for the next sample. 
			using std::polar;
			using bertini::polar;

			//Generalized since we could have a nonzero target time. 
			RT radius = abs(starting_time - target_time), angle = arg(starting_time - target_time); // generalized for nonzero target_time.

			auto next_sample = Vec<CT>(num_vars);
			auto next_time = (ii==this->EndgameSettings().num_sample_points-1) 
								?
							  starting_time
								:
							  polar(radius, (ii+1)*2*acos(static_cast<RT>(-1)) / (this->EndgameSettings().num_sample_points) + angle) + target_time;
			// If we are tracking to a nonzero target time we need to shift our values to track to. This is a step that may not be needed if target_time = 0
							  ; 


			auto tracking_success = this->GetTracker().TrackPath(next_sample, current_time, next_time, current_sample);	
			if (tracking_success != SuccessCode::Success)
			{
				std::cout << "tracker fail in circle track, radius " << radius << ", type " << int(tracking_success) << std::endl;
				return tracking_success;
			}

			AsDerived().EnsureAtPrecision(next_time,Precision(next_sample)); assert(Precision(next_time)==Precision(next_sample));

			auto refinement_success = AsDerived().RefineSample(next_sample, next_sample, next_time);
			if (refinement_success != SuccessCode::Success)
			{
				std::cout << "refinement fail in circle track, type " << int(refinement_success) << std::endl;
				return refinement_success;
			}

			AsDerived().EnsureAtPrecision(next_time,Precision(next_sample)); assert(Precision(next_time)==Precision(next_sample));

			circle_times.push_back(next_time);
			circle_samples.push_back(next_sample);

			// down here next_sample and next_time should have the same precision.
		}

		return SuccessCode::Success;

	}//end CircleTrack


	/**
		\brief A function that uses the assumption of being in the endgame operating zone to compute an approximation of the ratio c over k. 
			When the cycle number stabilizes we will see that the different approximations of c over k will stabilize. 
			Returns the computed value of c over k. 
	

		## Input: 
				None: all data needed are class data members.

		## Output:
				estimate: The approximation of the ratio of the two numbers C and K heuristically signifying we are in the cauchy endgame operating zone. 


		##Details:
				\tparam CT The complex number type.
				Consult page 53 of \cite bertinibook, for the reasoning behind this heuristic.
	*/
	template<typename CT>
	auto ComputeCOverK() -> typename Eigen::NumTraits<CT>::Real
	{//Obtain samples for computing C over K.
		using RT = typename Eigen::NumTraits<CT>::Real;
		using std::abs;
		using std::log;

		const auto& pseg_samples = std::get<SampCont<CT> >(pseg_samples_);

		assert(pseg_samples.size()>=3);
		const Vec<CT> & sample0 = pseg_samples[0];
		const Vec<CT> & sample1 = pseg_samples[1];
		const Vec<CT> & sample2 = pseg_samples[2];

		Vec<CT> rand_vector = Vec<CT>::Random(sample0.size()); //should be a row vector for ease in multiplying.


		// //DO NOT USE Eigen .dot() it will do conjugate transpose which is not what we want.
		// //Also, the .transpose*rand_vector returns an expression template that we do .norm of since abs is not available for that expression type. 
		RT estimate = abs(log(abs((((sample2 - sample1).transpose()*rand_vector).norm())/(((sample1 - sample0).transpose()*rand_vector).norm()))));
		estimate = abs(log(RT(this->EndgameSettings().sample_factor)))/estimate;
		if (estimate < 1)
		  	return RT(1);
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
				\tparam CT The complex number type.

	*/
	template<typename CT>
	bool CheckForCOverKStabilization(TimeCont<CT> const& c_over_k_array)
	{	
		using std::abs;

		assert(c_over_k_array.size()>=cauchy_settings_.num_needed_for_stabilization);
		for(unsigned ii = 1; ii < cauchy_settings_.num_needed_for_stabilization ; ++ii)
		{
			auto a = abs(c_over_k_array[ii-1]);
			auto b = abs(c_over_k_array[ii]);

			typename Eigen::NumTraits<CT>::Real divide = a;

			if(a < b)
				divide = a/b;
			else
				divide = b/a;

			if(divide <  cauchy_settings_.minimum_for_c_over_k_stabilization)
				return false;
		}
		return true;

	}//end CheckForCOverKStabilization


	/*
	Input: A time value and the space value above that time.
		
	Output: An mpfr_float representing a tolerance threshold for declaring a loop to be closed. 
	Details: Used in Bertini 1 as a heuristic for computing separatedness of roots. Decided to not be used since assumptions for this tolerance are not usually met. 
	template<typename CT>
	mpfr_float FindToleranceForClosedLoop(CT x_time, Vec<CT> x_sample)
	{
		auto degree_max = std::max(this->GetTracker().AMP_config_.degree_bound,mpfr_float("2.0")); 
		auto K = this->GetTracker().AMP_config_.coefficient_bound;
		mpfr_float N;
		mpfr_float M;
		mpfr_float L;
		if(max_closed_loop_tolerance_ < min_closed_loop_tolerance_)
		{
			max_closed_loop_tolerance_ = min_closed_loop_tolerance_;
		}
		auto error_tolerance = mpfr_float("1e-13");
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
		auto minimum_singular_value = Eigen::JacobiSVD< Mat<CT> >(jacobian_at_current_time).singularValues()(this->GetSystem().NumVariables() - 1 );
		auto norm_of_sample = x_sample.norm();
		L = pow(norm_of_sample,degree_max - 2);
		auto tol = K * L * M;
		// std::cout << "TOL IS " << tol << '\n';
		if (tol == 0) // fail-safe
					tol = minimum_singular_value;
			else
			{
			tol = mpfr_float("2.0") / tol;
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
				\tparam CT The complex number type
	*/
	template<typename CT>
	bool CheckClosedLoop()
	{	
		using RT = typename Eigen::NumTraits<CT>::Real;
		auto& times = std::get<TimeCont<CT> >(cauchy_times_);
		auto& samples = std::get<SampCont<CT> >(cauchy_samples_);

		if((samples.front() - samples.back()).norm() < this->GetTracker().TrackingTolerance())
		{
			return true;
		}

		if (TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			//Ensure all samples are of the same precision.
			auto new_precision = AsDerived().EnsureAtUniformPrecision(times, samples);
		}

		this->GetTracker().Refine(samples.front(),samples.front(),times.front(),RT(this->Tolerances().final_tolerance),this->EndgameSettings().max_num_newton_iterations);
		this->GetTracker().Refine(samples.back(),samples.back(),times.back(),RT(this->Tolerances().final_tolerance),this->EndgameSettings().max_num_newton_iterations);

		if((samples.front() - samples.back()).norm() < this->GetTracker().TrackingTolerance())
		{
			return true;
		}
		return false;	

	}//end CheckClosedLoop
		 



	/**
		\brief 	After we have used CircleTrack and have successfully closed the loop using CheckClosedLoop we need to check the maximum and minimum norms of the samples collected. 
				If the ratio of the maximum and minimum norm are within the threshold maximum_cauchy_ratio, and the difference is greater than final tolerance than we are successful. 


		## Input: 
			None: all data needed are class data members
		 
		## Output:
			true: If we are within the ratio cutoff time, or have a ratio within the thresholds.
			false: otherwise


		##Details:
				\tparam CT The complex number type.
				It is important to know if we are within the endgame operating zone. This function allows us to have a check that 
				heuristcially will tell us if we are. 
	*/
	template<typename CT>
	bool RatioEGOperatingZoneTest()
	{	
		using RT = typename Eigen::NumTraits<CT>::Real;
		RT min(1e300);
		RT max(0);
		auto& times = std::get<TimeCont<CT> >(cauchy_times_);
		auto& samples = std::get<SampCont<CT> >(cauchy_samples_);
		if(norm(times.front()) < cauchy_settings_.ratio_cutoff_time)
		{
			return true;
		}
		else
		{
			RT norm;
			for(unsigned int ii=0; ii < this->EndgameSettings().num_sample_points; ++ii)
			{
				norm = samples[ii].norm();
				if(norm > max)
				{
					max = norm;
				}
				if(norm < min)
				{
					min = norm;
				}
			}
			//compute the ratio if they are large enough.
			if(min > this->Tolerances().final_tolerance && max > this->Tolerances().final_tolerance)
			{
				norm = min / max;
				if(norm < cauchy_settings_.maximum_cauchy_ratio && (max - min) > this->Tolerances().final_tolerance)
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
				\tparam CT The complex number type.
	*/		 
	template<typename CT>
	SuccessCode InitialCauchyLoops(CT const& target_time)
	{	
		using RT = typename Eigen::NumTraits<CT>::Real;
		auto& cau_times = std::get<TimeCont<CT> >(cauchy_times_);
		auto& cau_samples = std::get<SampCont<CT> >(cauchy_samples_);
		auto& ps_times = std::get<TimeCont<CT> >(pseg_times_);
		auto& ps_samples = std::get<SampCont<CT> >(pseg_samples_);

		using std::max;
		bool continue_loop = true;

		auto fail_safe_max_cycle_number = max(cauchy_settings_.fail_safe_maximum_cycle_number,this->CycleNumber());

		auto initial_cauchy_loop_success = SuccessCode::Success;

		while (continue_loop)
		{	
			this->CycleNumber(0);
			cau_times.clear();
			cau_samples.clear();
			cau_samples.push_back(ps_samples.back()); // cauchy samples and times should be empty before this point. 
			cau_times.push_back(ps_times.back());

			auto next_time = ps_times.back();
			auto next_sample = ps_samples.back();

			// track around a circle once.  we'll use it to measure whether we believe we are in the eg operating zone, based on the ratio of norms of sample points around the circle
			auto tracking_success = CircleTrack(cau_times.front(),target_time,cau_samples.front());

			this->IncrementCycleNumber(1);

			if (tracking_success != SuccessCode::Success)
				return tracking_success;

			// find the ratio of the maximum and minimum coordinate wise for the loop. 
			if (RatioEGOperatingZoneTest<CT>())
			{ // then we believe we are in the EG operating zone, since the path is relatively flat.  i still disbelieve this is a good test (dab 20160310)
				while (true)
				{
					if (CheckClosedLoop<CT>())
					{//error is small enough, exit the loop with success. 
						initial_cauchy_loop_success = SuccessCode::Success;
						continue_loop = false;
						break;
					}
					else if(this->CycleNumber() > fail_safe_max_cycle_number)
					{// too many iterations
						initial_cauchy_loop_success = SuccessCode::CycleNumTooHigh;
						continue_loop = false;
						break;
					}

					//compute next loop, the last sample in times and samples is the sample our loop ended on. Either where we started or on another sheet at the same time value. 
					tracking_success = CircleTrack(cau_times.back(),target_time,cau_samples.back());

					this->IncrementCycleNumber(1);

					if(tracking_success != SuccessCode::Success)
						return tracking_success;
				}

				if(initial_cauchy_loop_success == SuccessCode::CycleNumTooHigh)
				{//see if we should continue to the next sample point
					if(abs(cau_times.back() - target_time) < this->EndgameSettings().min_track_time)
						continue_loop = false;
					else
						continue_loop = true;
				}
			}//end if (RatioEGOperatingZoneTest())
			else
			{
				//compute the time for the next sample point
				AsDerived().EnsureAtPrecision(next_time,Precision(ps_samples.back()));
				next_time = (ps_times.back() + target_time) * static_cast<RT>(this->EndgameSettings().sample_factor);

				SuccessCode tracking_success = this->GetTracker().TrackPath(next_sample,ps_times.back(),next_time,ps_samples.back());
				AsDerived().EnsureAtPrecision(next_time,Precision(next_sample));

				ps_times.pop_front();
				ps_samples.pop_front();

				ps_times.push_back(next_time);
				ps_samples.push_back(next_sample);

				if(tracking_success != SuccessCode::Success)
					return tracking_success;
			}
		} //end while(continue_loop)

		return initial_cauchy_loop_success;
	}//end InitialCauchyLoops


	/**
		\brief 	The Cauchy endgame will first find an initial approximation using the notion of the power series endgame. This function computes this approximation and returns a
	SuccessCode to let us know if an error was encountered. 

		## Input: 
				start_time: time value for which we start to make a power series approximation
				start_point: approximate solution to our homotopy H at the start_time
				approximation_time: time at which we are trying to find the solution, usually t = 0
				approximation approximate solution to our homotopy H at the approxmation_time


		## Output:
			SuccessCode reporting if any errors had occurred. All data collected is stored in class data members. 


		##Details:
	\tparam CT The complex number type.
				This function is in charge of finding the very first approximation of the origin. It does this by first computing some initial samples 
							 like what is done in the Power Series Endgame. We continue to track forward in this manner until we have stabilization of the cycle number being approximated. 
							 This prevents the unnecessary circle tracking if we are possibly not in the endgame operating zone. 
							 Once we have stabilization we then perform InitialCauchyLoops while getting the accurate cycle number, and check the norms of the samples and make sure we are ready 
							 to approximate. When ready we call ComputePSEGApproximationOfXAtT0. This function will use a hermtie interpolater to get an approximation of the value at the origin. 
	*/
	template<typename CT>
	SuccessCode InitialPowerSeriesApproximation(CT const& start_time, Vec<CT> const& start_point, 
	                                            CT const& approximation_time, Vec<CT> & approximation)
	{	
		using RT = typename Eigen::NumTraits<CT>::Real;

		//initialize array holding c_over_k estimates
		std::deque<RT> c_over_k; 

		auto& ps_times = std::get<TimeCont<CT> >(pseg_times_);
		auto& ps_samples = std::get<SampCont<CT> >(pseg_samples_);

		//Compute initial samples for pseg
		auto initial_sample_success = this->ComputeInitialSamples(start_time, approximation_time, start_point, ps_times, ps_samples);
		if (initial_sample_success!=SuccessCode::Success)
			return initial_sample_success;

		c_over_k.push_back(ComputeCOverK<CT>());

		Vec<CT> next_sample;
		CT next_time;

		unsigned ii = 0;
		//track until for more c_over_k estimates or until we reach a cutoff time. 
		while ( (ii < cauchy_settings_.num_needed_for_stabilization) )
		{	
			next_time = (ps_times.back() + approximation_time) * static_cast<RT>(this->EndgameSettings().sample_factor); // using general midpoint formula with sample_factor to give us a time
																												  // between ps_times.back() and the approximation_time.

			auto tracking_success = this->GetTracker().TrackPath(next_sample,ps_times.back(),next_time,ps_samples.back());
			if (tracking_success!=SuccessCode::Success)
				return tracking_success;

			AsDerived().EnsureAtPrecision(next_time, Precision(next_sample));

			ps_samples.pop_front();
			ps_times.pop_front();
			ps_samples.push_back(next_sample);
			ps_times.push_back(next_time);
			c_over_k.push_back(ComputeCOverK<CT>());

			++ii;
		}//end while


		//have we stabilized yet? 
		while(!CheckForCOverKStabilization(c_over_k) && abs(ps_times.back()) > cauchy_settings_.cycle_cutoff_time)
		{
			next_time = (ps_times.back() + approximation_time) * static_cast<RT>(this->EndgameSettings().sample_factor); // using general midpoint formula with sample_factor to give us a time
																												  // between ps_times.back() and the approximation_time.

			auto tracking_success = this->GetTracker().TrackPath(next_sample,ps_times.back(),next_time,ps_samples.back());

			if(tracking_success != SuccessCode::Success)
				return tracking_success;

			AsDerived().EnsureAtPrecision(next_time, Precision(next_sample));

			c_over_k.pop_front();
			ps_samples.pop_front();
			ps_times.pop_front();

			ps_samples.push_back(next_sample);	
			ps_times.push_back(next_time);
			c_over_k.push_back(ComputeCOverK<CT>());

		}//end while

		auto cauchy_loop_success = InitialCauchyLoops<CT>(approximation_time);
		if (cauchy_loop_success != SuccessCode::Success)
			return cauchy_loop_success;

		return ComputePSEGApproximationAtT0(approximation, approximation_time);

	}//end InitialPowerSeriesApproximation

	/**
		\brief 	This function takes the pseg_samples_ and pseg_times that have been collected and uses them to compute a Hermite interpolation to the time value time_t0. 

		## Input: 
				result: This vector holds the approxmation that we end up calculating
				time_t0: The time value at which we are trying to approximate, usually t = 0

		## Output:
				SuccessCode deeming if we were successful or if we encountered an error. 

		##Details:
				\tparam CT The complex number type. 
				This function handles computing an approximation at the origin. 
		We compute the derivatives at the different times and samples. We then make sure all samples are to the same precision before refining them to final tolerance. 
						By InitialCauchyLoops we know what the cycle number is so we convert derivatives and times to the s-plane where s = t^(1/(cyle number).
		We use the converted times and derivatives along with the samples to do a Hermite interpolation which is found in base_endgame.hpp.
	*/
	template<typename CT>
	SuccessCode ComputePSEGApproximationAtT0(Vec<CT>& result, const CT & time_t0)
	{
		using RT = typename Eigen::NumTraits<CT>::Real;

		auto& ps_times = std::get<TimeCont<CT> >(pseg_times_);
		auto& ps_samples = std::get<SampCont<CT> >(pseg_samples_);

		if (TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			//Ensure all samples are of the same precision.
			auto new_precision = AsDerived().EnsureAtUniformPrecision(ps_times, ps_samples);
		}


		for (unsigned ii=0; ii<ps_samples.size(); ++ii)
		{
			auto refine_code = AsDerived().RefineSample(ps_samples[ii],ps_samples[ii],ps_times[ii]);
			if (refine_code != SuccessCode::Success)
				return refine_code;
		}
		
		if (TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			//Ensure all samples are of the same precision.
			auto new_precision = AsDerived().EnsureAtUniformPrecision(ps_times, ps_samples);
		}
		

		auto num_sample_points = this->EndgameSettings().num_sample_points;
		//Compute dx_dt for each sample.
		SampCont<CT> pseg_derivatives;
		for(unsigned ii = 0; ii < num_sample_points; ++ii)
		{	
			// the inverse() call uses LU look at Eigen documentation on inverse in Eigen/LU.
			pseg_derivatives.push_back(
			                           -(this->GetSystem().Jacobian(ps_samples[ii],ps_times[ii]).inverse())*this->GetSystem().TimeDerivative(ps_samples[ii],ps_times[ii])
			                           );
		}

 		//Conversion to S-plane.
		TimeCont<CT> s_times(num_sample_points);
		SampCont<CT> s_derivatives(num_sample_points);
		auto c = static_cast<RT>(this->CycleNumber());
		for(unsigned ii = 0; ii < num_sample_points; ++ii)
		{
			s_derivatives[ii] = pseg_derivatives[ii]*
			                        (c*pow(ps_times[ii],(c-1)/c));
			s_times[ii] =  pow(ps_times[ii], 1/c);
		}

		result = HermiteInterpolateAndSolve(time_t0, num_sample_points, 
                                            s_times, ps_samples, s_derivatives);
		return SuccessCode::Success;
	}//end ComputePSEGApproximationOfXAtT0


	/**
	\brief Function that computes the mean of the samples that were collected while tracking around the origin. This value is the approximation of the value at the origin. 

		## Input: 
			result: This vector, passed by reference, holds the approximation that we calculate. 

		## Output:
			SuccessCode deeming if we were suceessful, or if we encountered an error. 

		##Details:
	\tparam CT The complex number type.
				We can compute the Cauchy Integral Formula in this particular instance by computing the mean of the samples we have collected around the origin. 
	*/
	template<typename CT>
	SuccessCode ComputeCauchyApproximationOfXAtT0(Vec<CT>& result)
	{	
		using RT = typename Eigen::NumTraits<CT>::Real;
		auto& cau_times = std::get<TimeCont<CT> >(cauchy_times_);
		auto& cau_samples = std::get<SampCont<CT> >(cauchy_samples_);

		if (cau_samples.size() != this->CycleNumber() * this->EndgameSettings().num_sample_points+1)
		{
			std::stringstream err_msg;
			err_msg << "to compute cauchy approximation, cauchy_samples must be of size " << this->CycleNumber() * this->EndgameSettings().num_sample_points+1 << " but is of size " << cau_samples.size() << '\n';
			throw std::runtime_error(err_msg.str());
		}

		result = Vec<CT>::Zero(this->GetSystem().NumVariables());//= (cau_samples[0]+cau_samples.back())/2; 

		if (TrackerTraits<TrackerType>::IsAdaptivePrec)
		{
			//Ensure all samples are of the same precision.
			auto new_precision = AsDerived().EnsureAtUniformPrecision(cau_times, cau_samples);
		}

		for(unsigned int ii = 0; ii < this->CycleNumber() * this->EndgameSettings().num_sample_points; ++ii)
		{
			auto refine_code = AsDerived().RefineSample(cau_samples[ii],cau_samples[ii],cau_times[ii]);
			if (refine_code!=SuccessCode::Success)
				return refine_code;
			result += cau_samples[ii];
		}
		result /= static_cast<RT>(this->CycleNumber() * this->EndgameSettings().num_sample_points);
		return SuccessCode::Success;

	}

	/**
	\brief Function that will utilize CircleTrack and CheckClosedLoop to collect all samples while tracking around the origin till we close the loop. 

		## Input: 
			starting_time: the time value at which we start finding cauchy samples 
			target_time: the time value we are computing cauchy samples around
			starting_sample: the current space point at starting_time, this will also be the first cauchy sample


		## Output: 
			SuccessCode deeming if we were able to collect all samples around the origin, or if we encounted an error at some point. 

		##Details:
	\tparam CT The complex number type.
				This function populates the deque cauchy_samples and cauchy_times. These are data members of the class and are not passed in. This function will continue to 
				call CircleTrack until we have closed the loop. 

	*/
	template<typename CT>
	SuccessCode ComputeCauchySamples(CT const& starting_time,CT const& target_time, Vec<CT> const& starting_sample)
	{
		using bertini::Precision;
		assert(Precision(starting_time)==Precision(starting_sample));

		auto& cau_times = std::get<TimeCont<CT> >(cauchy_times_);
		auto& cau_samples = std::get<SampCont<CT> >(cauchy_samples_);

		cau_times.clear();
		cau_samples.clear();
		cau_times.push_back(starting_time);
		cau_samples.push_back(starting_sample);
		this->CycleNumber(0);


		while( this->CycleNumber() < cauchy_settings_.fail_safe_maximum_cycle_number )
		{
			//track around the origin once.
			auto tracking_success = CircleTrack(cau_times.back(),target_time,cau_samples.back());
			this->IncrementCycleNumber(1);

			if(tracking_success != SuccessCode::Success)
			{
				std::cout << "Cauchy loop fail tracking "<< int(tracking_success) <<"\n\n";
				return tracking_success;
			}
			else if(CheckClosedLoop<CT>())
			{
				return SuccessCode::Success;
			}
		} 

		return SuccessCode::CycleNumTooHigh;
	}//end ComputeCauchySamples



	/**
	\brief Primary function that runs the Cauchy endgame.
	To begin, this function will compute a first approximation using the power series endgame notion. This approximation is made after a heuristic on the stabilization of 
	the cyle number is made, and after the maximum and minimum norms of tracked space values around the origin are withing a certain tolerance. 			

		## Input: 
			start_time: the time value at which we start the endgame
			start_point: an approximate solution to our homotopy H at start_time
			target_time: the time value that we are using the endgame to interpolate to, default set to t = 0.

		## Output: 
			SuccessCode: reporting if we were successful in the endgame or if we encountered an error

		##Details:
	\tparam CT The complex number type.
				This function runs the entire Cauchy Endgame. We first take our endgame boundary time value and sample to find a first approximation of the origin. This is done by
					using the idea for the power series endgame. We check for stabilization of the cycle number, and check to see when the ratios of the maximum and minimum norm of samples collected
					by CircleTrack are withing a tolerance. When both of these conditions are met we do a Hermite interpolation. 
					At this point we can start tracking in to the origin while using CircleTrack to compute samples and calculating their mean to get an approximation of the origin using the Cauchy
					Integral Formula. 
	*/
	template<typename CT>
	SuccessCode Run(CT const& start_time, Vec<CT> const& start_point, CT const& target_time = static_cast<CT>(0))
	{	
		if (start_point.size()!=this->GetSystem().NumVariables())
		{
			std::stringstream err_msg;
			err_msg << "number of variables in start point for CauchyEG, " << start_point.size() << ", must match the number of variables in the system, " << this->GetSystem().NumVariables();
			throw std::runtime_error(err_msg.str());
		}

		assert(Precision(start_time)==Precision(start_time) && ("CauchyEG Run time and point must be of matching precision"));

		using RT = typename Eigen::NumTraits<CT>::Real;

		auto& cau_times = std::get<TimeCont<CT> >(cauchy_times_);
		auto& cau_samples = std::get<SampCont<CT> >(cauchy_samples_);
		auto& ps_times = std::get<TimeCont<CT> >(pseg_times_);
		auto& ps_samples = std::get<SampCont<CT> >(pseg_samples_);


		ClearTimesAndSamples<CT>(); //clear times and samples before we begin.
		this->CycleNumber(0);

		Vec<CT> prev_approx, latest_approx;
		RT approximate_error;
		Vec<CT>& final_approx = std::get<Vec<CT> >(this->final_approximation_at_origin_);
	

		//Compute the first approximation using the power series approximation technique. 
		auto initial_ps_success = InitialPowerSeriesApproximation(start_time, start_point, target_time, prev_approx);  // last argument is output here
		if(initial_ps_success != SuccessCode::Success)
			return initial_ps_success;

		//Generalized next_time in case if we are not trying to converge to the t = 0.
		CT next_time = (ps_times.back() + target_time) * static_cast<RT>(this->EndgameSettings().sample_factor);
		RT norm_of_dehom_of_prev_approx, norm_of_dehom_of_latest_approx;

		if(this->SecuritySettings().level <= 0)
			norm_of_dehom_of_prev_approx = this->GetSystem().DehomogenizePoint(prev_approx).norm();
		do
		{
			//Compute a cauchy approximation.  Uses the previously computed samples, 
			//either from InitialCauchyLoops, or ComputeCauchySamples
			auto extrapolation_success = ComputeCauchyApproximationOfXAtT0<CT>(latest_approx);
			if (extrapolation_success!=SuccessCode::Success)
				return extrapolation_success;

			if (this->SecuritySettings().level <= 0)
				norm_of_dehom_of_latest_approx = this->GetSystem().DehomogenizePoint(latest_approx).norm();

			approximate_error = (latest_approx - prev_approx).norm(); //Calculate the error between approximations. 

			// dehom of prev approx and last approx not used because they are not updated with the most current information. However, prev approx and last approx are 
			// the most current. 

			if (approximate_error < this->Tolerances().final_tolerance)
			{
				final_approx = latest_approx;
				return SuccessCode::Success;
			}
			else if (abs(cau_times.front() - target_time) < this->EndgameSettings().min_track_time) // generalized for target_time not equal to 0. 
			{//we are too close to t = 0 but we do not have the correct tolerance - so we exit
				final_approx = latest_approx;
				return SuccessCode::FailedToConverge;
			}
			else if (this->SecuritySettings().level <= 0 && 
			   norm_of_dehom_of_prev_approx   > this->SecuritySettings().max_norm &&  
			   norm_of_dehom_of_latest_approx > this->SecuritySettings().max_norm  )
			{//we are too large, break out of loop to return error.
				final_approx = latest_approx;
				return SuccessCode::SecurityMaxNormReached;
			}


			prev_approx = latest_approx;
			norm_of_dehom_of_prev_approx = norm_of_dehom_of_latest_approx;

			//Generalized next_time in case if we are not trying to converge to the t = 0.
			next_time = (next_time + target_time) * static_cast<RT>(this->EndgameSettings().sample_factor);

			

			Vec<CT> next_sample;
			auto tracking_success = this->GetTracker().TrackPath(next_sample,cau_times.back(),next_time,cau_samples.front());
			if (tracking_success != SuccessCode::Success)
				return tracking_success;

			AsDerived().EnsureAtPrecision(next_time,Precision(next_sample));

			ps_times.push_back(next_time);  ps_times.pop_front();
			ps_samples.push_back(next_sample); ps_samples.pop_front();

			auto cauchy_samples_success = ComputeCauchySamples(next_time,target_time,next_sample);

			//Added because of Griewank osborne test case where tracker returns GoingToInfinity.  
			//The cauchy endgame should stop instead of attempting to continue. 
			//This is because the tracker will not track to any meaningful space values. 
			if (cauchy_samples_success == SuccessCode::GoingToInfinity)
				return SuccessCode::GoingToInfinity;
			else if (cauchy_samples_success != SuccessCode::Success && abs(cau_times.front()-target_time) < this->EndgameSettings().min_track_time)
			{// we are too close to t = 0 but we do have the correct tolerance -so we exit.
				final_approx = latest_approx;
				return SuccessCode::MinTrackTimeReached;
			}
			else if(cauchy_samples_success != SuccessCode::Success)
			{
				final_approx = latest_approx;
				return cauchy_samples_success;
			}
		} while (approximate_error > this->Tolerances().final_tolerance);

		final_approx = latest_approx;
		return SuccessCode::Success;
	} //end main CauchyEG function
};


}}} // namespaces
